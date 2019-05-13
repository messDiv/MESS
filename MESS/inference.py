from __future__ import print_function

import MESS.stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from boruta import BorutaPy
from MESS.util import progressbar
from skgarden import RandomForestQuantileRegressor
from sklearn import metrics
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier,\
                            GradientBoostingClassifier, GradientBoostingRegressor,\
                            AdaBoostRegressor, AdaBoostClassifier
from sklearn.model_selection import train_test_split, cross_val_score, RandomizedSearchCV
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline

import logging
LOGGER = logging.getLogger(__name__)

## The base Ensemble class takes care of reading in the empirical dataframe,
## calculating summary stats, reading the simulated date, and reshaping the sim 
## sumstats to match the stats of the real data.
class Ensemble(object):
    def __init__(self, empirical_df, simfile, algorithm="rf", verbose=False):
        self.algorithm = algorithm
        ## Fetch empirical data and calculate the sumstats for the datatypes passed in
        self.empirical_df = empirical_df
        try:
            self.empirical_sumstats = MESS.stats.calculate_sumstats(empirical_df)
            if verbose: print("Got empirical summary statistics: {}".format(self.empirical_sumstats))
        except Exception as inst:
            print("  Malformed input dataframe: {}".format(inst))

        try:
            ## Read in simulated data, split it into features and targets,
            ## and prune the features to correspond to the loaded empirical data
            self.sim_df = pd.read_csv(simfile, sep="\t", header=0)
            ## Species richness is the first summary statistic, so split to
            ## features and targets here. We're keeping all features and
            ## targets in the _X and _y variables, and can selectively prune
            ## them later
            S_idx = list(self.sim_df.columns).index("S")
            self._X = self.sim_df.iloc[:, S_idx:]
            self._y = self.sim_df.iloc[:, :S_idx]
        except Exception as inst:
            print("  Failed loading simulations file: {}".format(inst))
            raise

        ## Take the default targets intelligently depending on whether this
        ## is a classifier or a regressor.
        self.set_targets()
        ## Set features to correspond to real data. Trims self.X to only
        ## the features necessary, leaves self._X intact.
        self.set_features(self.empirical_sumstats.columns)

        ## If you want to estimate parameters independently then
        ## we keep track of which features are most relevant per target.
        ## By default we'll just assume you want to use all features,
        ## if you mod this variable or call the feature_selection method
        ## then this will get over-written.
        self.model_by_target = {x:{"features":list(self.X.columns)} for x in self.y.columns}


    def __repr__(self):
        return "{}: nsims - {}\n\tFeatures - {}\n\tTargets - {}".format(type(self), len(self._X), list(self.features), self.targets)


    def set_features(self, feature_list=''):
        if not len(feature_list):
            self.features = self._X.columns
        else:
            self.features = feature_list
        try:
            self.X = self._X[self.features]
        except exception as inst:
            print("  Failed setting features: {}".format(inst))
            raise


    def set_targets(self, target_list=''):
        if isinstance(target_list, str) and target_list:
            target_list = [target_list]
        if not len(target_list):
            self.targets = self._default_targets
        else:
            self.targets = target_list
        try:
            self.y = self._y[self.targets]
        except Exception as inst:
            print("  Failed setting targets: {}".format(inst))
            raise


    ## The feature selection method only takes one target vector
    ## at a time, so we'll do each target seperately and then
    ## take the union of the results.
    ## Uses BorutaPy, an all-relevant feature selection method
    ## https://github.com/scikit-learn-contrib/boruta_py
    ## http://danielhomola.com/2015/05/08/borutapy-an-all-relevant-feature-selection-method/
    def feature_selection(self, quick=False, verbose=False):
        mask = np.zeros(len(self.X.columns), dtype=bool)
        if verbose: print("Selecting features:")
        for t in self.y:
            if verbose: print("  {}\t".format(t), end='')
            model = self._base_model(n_estimators=600, max_depth=5)

            if quick:
                ## Random subset the data and run fewer iterations. If you don't have
                ## at least 1000 samples then don't bother
                nsamps = min(len(self.y), 1000)
                idxs = np.random.choice(len(self.y), nsamps, replace=False)
                tmpX = self.X.iloc[idxs].values
                tmpy = self.y[t].iloc[idxs].values
                max_iter=10
            else:
                tmpX = self.X.values
                tmpy = self.y[t].values 
                max_iter=100

            # define Boruta feature selection method
            feat_selector = BorutaPy(model, max_iter=max_iter, n_estimators='auto', verbose=0, random_state=1)
            ## Turn off np warnings
            np.warnings.filterwarnings('ignore')
            feat_selector.fit(tmpX, tmpy)

            # check ranking of features
            features = list(self.features[feat_selector.support_])
            if verbose: print("{}".format(features))
            if not features:
                print("  NB: No features found relevant for target {}\n      Fall back to using all parameters.".format(t))
                feat_selector.support_ = feat_selector.support_ | True
                features = self.features
            ## Remember the relevant features for this target, for later prediction
            ## Use feat_selector.support_weak_ to include more variables
            self.model_by_target[t]["features"] = features

            mask += feat_selector.support_

        selected_features = self.features[mask]
        if verbose: print("All selected features: {}".format(" ".join(selected_features)))
        ## Filter on the selected features
        self.X = self.X[selected_features]   
        ## Also filter the empirical sumstats
        self.empirical_sumstats = self.empirical_sumstats[selected_features]

        self.features = selected_features


    def param_search_cv(self, by_target=False, quick=False, verbose=False):

        if verbose: print("Finding best model parameters.")
        if quick:
            n_iter=5
            nsamps = min(len(self.y), 1000)
            idxs = np.random.choice(len(self.y), nsamps, replace=False)
            tmpX = self.X.iloc[idxs]
            tmpy = self.y.iloc[idxs]
        else:
            n_iter=100
            tmpX = self.X
            tmpy = self.y
        ## Randomly search 100 different parameter combinations and take the
        ## one that reduces CV error
        cvsearch = RandomizedSearchCV(estimator=self._base_model(),\
                                       param_distributions=self._param_grid,
                                       n_iter=n_iter, cv=4, verbose=verbose, n_jobs=-1)
        if by_target:
            for t in self.targets:
                ## Munge the tmpX to fit the features selected for each target.
                ## This is annoying housekeeping.
                tmpX_filt = tmpX[self.model_by_target[t]["features"]]
                if verbose: print("\t{} - Finding best params using features: {}".format(t, list(tmpX_filt.columns)))
                cvsearch.fit(tmpX_filt, tmpy[t])
                if verbose: print("Best params for {}: {}".format(t, cvsearch.best_params_))
                self.model_by_target[t]["cvsearch"] = cvsearch
                self.model_by_target[t]["model"] = cvsearch.best_estimator_
                self.model_by_target[t]["feature_importances"] = self.model_by_target[t]["model"].feature_importances_
            if len(self.targets) == 1:
                ## For convenience. If there's only one target then set self.best_model
                self.best_model = self.model_by_target[self.targets[0]]["model"]
        else:
            cvsearch.fit(tmpX, tmpy)
            if verbose: print(cvsearch.best_params_)
            self._cvsearch = cvsearch
            self.best_model = cvsearch.best_estimator_


    ## The magic method to just do-it-all
    def predict(self, select_features=True, param_search=True, by_target=False, quick=False, verbose=False):
        ## Select relevant features
        if select_features:
            self.feature_selection(quick=quick, verbose=verbose)

        ## If using either the RFQuantileRegressor or GradientBoosting, then
        ## we force by_targets to be true, since they'll only do one target at at time
        ## RandomForest can handle multi-target regression, but not rfq or gb, so we
        ## allow plain rf to optionally do by_target
        self._by_target = by_target
        if self.algorithm in ["rfq", "gb", "ab"]:
            self._by_target = True

        if param_search:
            self.param_search_cv(by_target=self._by_target, quick=quick, verbose=verbose)
        else:
            ## If you don't want to do the param search, then just take the default param
            ## Also, assuming if you don't care about setting model parameters, then you
            ## won't care about feature selection, so just use all features for the
            ## the by_target code
            if self._by_target:
                for t in self.targets:
                    self.model_by_target[t]["model"] = self._base_model()
                    self.model_by_target[t]["model"].fit(self.X[self.model_by_target[t]["features"]], self.y[t])
                    self.model_by_target[t]["feature_importances"] = self.model_by_target[t]["model"].feature_importances_
            else:
                ## TODO: Make default base_model params smarter
                self.best_model = self._base_model(n_jobs=-1)
                self.best_model.fit(self.X, self.y)


    ## Return a summary of feature importances, either summed over all
    ## parameters, or individually per parameter
    def feature_importances(self):
        importances = []
        if self._by_target:
            for t in self.targets:
                importances.append(pd.DataFrame(self.model_by_target[t]["feature_importances"],\
                                            columns=[t], index=self.model_by_target[t]["features"]).T)
            importances = pd.concat(importances)
        else:
            ## If model_by_target dictionary isn't populated then get all
            ## feature importances for the joint inference.
            importances = pd.DataFrame(self.best_model.feature_importances_, columns=["Feature importance"], index=self.features).T
        return importances


class Classifier(Ensemble):
    _default_targets = ["community_assembly_model"]

    def __init__(self, empirical_df, simfile, algorithm="rf", verbose=False):
        super(Classifier, self).__init__(empirical_df, simfile, algorithm=algorithm, verbose=verbose)

        if algorithm == "rf":
            self._base_model = RandomForestClassifier
        elif algorithm == "gb":
            self._base_model = GradientBoostingClassifier
        elif algorithm == "ab":
            self._base_model = AdaBoostClassifier
        else:
            raise Exception(" Unsupported classification algorithm: {}".format(algorithm))
        self.algorithm = algorithm
        self._param_grid = _get_param_grid(algorithm)


    def predict(self, select_features=True, param_search=True, by_target=False, quick=False, verbose=False):

        super(Classifier, self).predict(select_features=select_features, param_search=param_search,\
                                        by_target=by_target, quick=quick, verbose=verbose)

        ## Force by_target to be true for GradientBoosting and AdaBoost
        self._by_target = by_target
        if self.algorithm in ["gb", "ab"]:
            self._by_target = True

        if self._by_target:
            ## Predict each target independently using it's own trained RF
            preds = [self.model_by_target[t]["model"].predict(self.empirical_sumstats[self.model_by_target[t]["features"]]) for t in self.targets]
            self.empirical_pred = pd.DataFrame(np.array(preds).T, columns=self.targets, index=["estimate"])

            proba = [self.model_by_target[t]["model"].predict_proba(self.empirical_sumstats[self.model_by_target[t]["features"]])[0] for t in self.targets]
            ## Somewhat annoying, but we grab the classes vector from the model of the first target
            self.empirical_proba = pd.DataFrame(proba, columns=self.model_by_target.values()[0]["model"].classes_, index=self.targets)
        else:
            ## Do all targets at once. Also, you don't get prediction intervls
            ## if you don't do by_target
            self.empirical_pred = pd.DataFrame(self.best_model.predict(self.empirical_sumstats),\
                                                columns=self.targets, index=["estimate"])
            self.empirical_proba = pd.DataFrame(self.best_model.predict_proba(self.empirical_sumstats),\
                                                columns=self.best_model.classes_, index=self.targets)

        return self.empirical_pred, self.empirical_proba


class Regressor(Ensemble):
    _default_targets = ["alpha", "J_m", "ecological_strength", "m", "speciation_prob", "_lambda"]

    def __init__(self, empirical_df, simfile, algorithm="rfq", verbose=False):
        super(Regressor, self).__init__(empirical_df, simfile, algorithm="rf", verbose=False)

        if algorithm == "rf":
            self._base_model = RandomForestRegressor
        elif algorithm == "rfq":
            self._base_model = RandomForestQuantileRegressor
        elif algorithm == "gb":
            self._base_model = GradientBoostingRegressor
        elif self.algorithm == "ab":
            self._base_model = AdaBoostRegressor
        else:
            raise Exception(" Unsupported regression algorithm: {}".format(algorithm))
        self.algorithm = algorithm
        self._param_grid = _get_param_grid(algorithm)


    ## Add upper and lower prediction interval for algorithms that support
    ## quantile regression (rfq, gq)
    ## The quick parameter doesn't do anything for RFQuantileRegressor because
    ## it's already really fast (you don't have to refit the model).
    def prediction_interval(self, interval=0.95, quick=False, verbose=False):
        if verbose: print("Calculating prediction interval(s)")
        upper = 1.0 - ((1.0 - interval)/2.)
        lower = 1.0 - upper
        if self.algorithm == "rfq":
            y_lower = [self.model_by_target[t]["model"].predict(self.empirical_sumstats[self.model_by_target[t]["features"]], lower*100) for t in self.targets]
            y_upper = [self.model_by_target[t]["model"].predict(self.empirical_sumstats[self.model_by_target[t]["features"]], upper*100) for t in self.targets]
        elif self.algorithm == "gb":
            if quick:
                nsamps = min(len(self.y), 1000)
                idxs = np.random.choice(len(self.y), nsamps, replace=False)
                tmpX = self.X.iloc[idxs]
                tmpy = self.y.iloc[idxs]
            else:
                tmpX = self.X
                tmpy = self.y
            y_lower = []
            y_upper = []
            for t in self.targets:
                if verbose: print("\t{}".format(t))
                tmp_gb = self.model_by_target[t]["model"].set_params(loss="quantile").set_params(alpha=lower)
                tmp_gb.fit(tmpX, tmpy[t])
                y_lower.append(tmp_gb.predict(self.empirical_sumstats))
                tmp_gb = self.model_by_target[t]["model"].set_params(loss="quantile").set_params(alpha=upper)
                tmp_gb.fit(tmpX, tmpy[t])
                y_upper.append(tmp_gb.predict(self.empirical_sumstats))
        else:
            print("Unsupported algorithm for prediction intervals - {}".format(self.algorithm))
            return self.empirical_pred

        ## Concatenate lower and upper quartiles onto the prediction df and name the rows nicely
        self.y_lower = pd.DataFrame(y_lower, columns=["lower {}".format(lower)], index=self.targets).T
        self.y_upper = pd.DataFrame(y_upper, columns=["upper {}".format(upper)], index=self.targets).T
        self.empirical_pred = pd.concat([self.empirical_pred, self.y_lower, self.y_upper])

        return self.empirical_pred


    def predict(self, select_features=True, param_search=True, by_target=False, quick=False, verbose=False):
        super(Regressor, self).predict(select_features=select_features, param_search=param_search,\
                                        by_target=by_target, quick=quick, verbose=verbose)

        self._by_target = by_target
        if self.algorithm in ["rfq", "gb", "ab"]:
            self._by_target = True

        if self._by_target:
            ## Predict each target independently using it's own trained RF
            preds = [self.model_by_target[t]["model"].predict(self.empirical_sumstats[self.model_by_target[t]["features"]]) for t in self.targets]
            self.empirical_pred = pd.DataFrame(np.array(preds).T, columns=self.targets, index=["estimate"])

            ## If using one of the algorithms that supports quantile regression then
            ## return the prediction intervals as well
            if self.algorithm in ["rfq", "gb"]:
                self.empirical_pred = self.prediction_interval(interval=0.95, quick=quick, verbose=verbose)
        else:
            ## Do all targets at once. Also, you don't get prediction intervls
            ## if you don't do by_target
            self.empirical_pred = pd.DataFrame(self.best_model.predict(self.empirical_sumstats),\
                                                columns=self.targets, index=["estimate"])

        return self.empirical_pred


def posterior_predictive_check(empirical_df,\
                                parameter_estimates,\
                                ax='',\
                                est_only=False,\
                                nsims=100,\
                                outfile='',\
                                verbose=False):

    ## If est_only, drop the lower and upper prediction interval and just use
    ## the mean estimated parameters for generating posterior predictive simes
    if est_only:
        parameter_estimates = pd.DataFrame(parameter_estimates.loc["estimate", :]).T

    ## Pre-fetch nsims worth of  sampled parameters from the prediction interval
    param_df = pd.DataFrame(data=np.zeros(shape=(nsims,len(parameter_estimates.columns))),\
                            columns=parameter_estimates.columns)
    ## You only have the point estimate that's all you can use
    if len(parameter_estimates) == 1:
        for param in parameter_estimates:
            param_df[param] = [parameter_estimates[param].values[0]] * nsims
    ## If you have prediction intervals sample uniform between them
    elif len(parameter_estimates) == 3:
        for param in parameter_estimates:
            param_df[param] = np.random.uniform(low=parameter_estimates[param]["lower 0.025"],
                                                high=parameter_estimates[param]["upper 0.975"],
                                                size=nsims)
    else:
        raise Exception("Shape of parameter estimate df not understood: {}".format(parameter_estimates.shape))

    r = MESS.Region("ppc")
    r.set_param("project_dir", "./ppc")

    if verbose: progressbar(nsims, 0, "Performing simulations")
    for i in range(nsims):
        for param in parameter_estimates:
            if param == "_lambda":
                #r.set_param("generations", est[param]["estimate"])
                r.set_param("generations", param_df[param].iloc[i])
            else:
                #r.set_param(param, est[param]["estimate"])
                r.set_param(param, param_df[param].iloc[i])
        r.run(sims=1, quiet=True)

        if verbose: progressbar(nsims, i, "Performing simulations")
    if verbose: progressbar(nsims, nsims, "Performing simulations")
    if verbose: print("\nCalculating PCs and plotting")

    simfile = "./ppc/SIMOUT.txt"
    sim_df = pd.read_csv(simfile, sep="\t", header=0)

    ## Chop the sims down to only include ss contained in the observed data
    obs_ss = MESS.util.calculate_sumstats(empirical_df)
    sim_df = sim_df[obs_ss.columns]

    if not ax:
        fig, ax = plt.subplots(figsize=(5, 5))

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(pd.concat([obs_ss, sim_df]))

    ax.scatter(pcs[:, 0], pcs[:, 1])
    ## Plot the observed ss in red
    ax.scatter(pcs[:, 0][0], pcs[:, 1][0], c='r')

    if outfile:
        try:
            ax.savefig(outfile)
        except Exception as inst:
            raise Exception("Failed saving figure: {}".format(inst))

    return ax


####################################################################
## Convenience functions for wrapping ML parameter grid construction
####################################################################
def _get_param_grid(algorithm):
    if algorithm in ["rfq", "rf"]:
        return _rf_param_grid()
    elif algorithm == "gb":
        return _gb_param_grid()
    elif algorithm == "ab":
        return _ab_param_grid()


def _ab_param_grid():
    n_estimators = [int(x) for x in np.linspace(10, 1000, 10)]
    learning_rate = np.linspace(0.01, 2, 6)
    random_grid = {'n_estimators':n_estimators,
                    'learning_rate':learning_rate}


def _rf_param_grid():
    n_estimators = [int(x) for x in np.linspace(200, 2000, 10)]
    max_depth = [int(x) for x in np.linspace(10, 110, 11)]
    max_depth.append(None)
    min_samples_split = [2, 5, 10]
    min_samples_leaf = [1, 2, 4]
    bootstrap = [True, False]
    random_grid = {'n_estimators': n_estimators,
                   'max_depth': max_depth,
                   'min_samples_split': min_samples_split,
                   'min_samples_leaf': min_samples_leaf,
                   'bootstrap': bootstrap}
    return random_grid


def _gb_param_grid():
    learning_rate = np.logspace(-4, -0.5, 6)
    n_estimators = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)]
    min_samples_split = [2, 5, 10]
    min_samples_leaf = [1, 2, 4]
    max_features = ['auto', 'sqrt']
    max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
    max_depth.append(None)
    random_grid = {'learning_rate': learning_rate,
                   'n_estimators': n_estimators,
                   'min_samples_split': min_samples_split,
                   'min_samples_leaf': min_samples_leaf,
                   'max_features': max_features,
                   'max_depth': max_depth}
    return random_grid
