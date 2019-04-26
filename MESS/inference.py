from __future__ import print_function

import MESS.stats
import numpy as np
import pandas as pd

from boruta import BorutaPy
from skgarden import RandomForestQuantileRegressor
from sklearn import metrics
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier, GradientBoostingClassifier, GradientBoostingRegressor
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
            if verbose: print("{}".format(list(self.features[feat_selector.support_])))

            ## Remember the relevant features for this target, for later prediction
            ## Use feat_selector.support_weak_ to include more variables
            self.model_by_target[t]["features"] = list(self.features[feat_selector.support_])

            mask += feat_selector.support_
            #X_filtered = pd.DataFrame(X_filtered, columns = X.columns[feat_selector.support_])

        selected_features = self.features[mask]
        if verbose: print("All selected features: {}".format(" ".join(selected_features)))
        ## Filter on the selected features
        self.X = self.X[selected_features]   
        ## Also filter the empirical sumstats
        self.empirical_sumstats = self.empirical_sumstats[selected_features]

        self.features = selected_features


class Classifier(Ensemble):
    _default_targets = ["community_assembly_model"]

    def __init__(self, empirical_df, simfile, algorithm="rf", verbose=False):
        super(Classifier, self).__init__(empirical_df, simfile, algorithm=algorithm, verbose=verbose)

        if algorithm == "rf":
            self._base_model = RandomForestClassifier
        elif algorithm == "gb":
            self._base_model = GradientBoostingClassifier
        self._param_grid = _get_param_grid(algorithm)


class Regressor(Ensemble):
    _default_targets = ["alpha", "J_m", "ecological_strength", "m", "speciation_prob", "_lambda"]

    def __init__(self, empirical_df, simfile, algorithm="rfq", verbose=False):
        super(Regressor, self).__init__(empirical_df, simfile, algorithm="rf", verbose=False)

        if algorithm == "rf":
            self.algorithm = "rf"
            self._base_model = RandomForestRegressor
        elif algorithm == "rfq":
            self.algorithm = "rfq"
            self._base_model = RandomForestQuantileRegressor
        elif algorithm == "gb":
            self.algorithm = "gb"
            self._base_model = GradientBoostingRegressor
        else:
            raise Exception(" Unsupported regression algorithm: {}".format(algorithm))
        self._param_grid = _get_param_grid(algorithm)


    def param_search_cv(self, by_target=False, quick=False, verbose=False):

        if verbose: print("Finding best model parameters.")
        if quick:
            n_iter=5
            tmpX = self.X[:1000]
            tmpy = self.y[:1000]
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
                ## TODO: Munge the tmpX to fit the features selected for each target
                tmpX_filt = tmpX[self.model_by_target[t]["features"]]
                if verbose: print("\t{} - Finding best params using features: {}".format(t, list(tmpX_filt.columns)))
                cvsearch.fit(tmpX, tmpy[t])
                if verbose: print("Best params for {}: {}".format(t, cvsearch.best_params_))
                self.model_by_target[t]["cvsearch"] = cvsearch
                self.model_by_target[t]["model"] = cvsearch.best_estimator_
        else:
            cvsearch.fit(tmpX, tmpy)
            if verbose: print(cvsearch.best_params_)
            self._cvsearch = cvsearch
            self.best_model = cvsearch.best_estimator_


    ## Add upper and lower prediction interval for algorithms that support
    ## quantile regression (rfq, gq)
    def prediction_interval(self, interval=0.95):
        upper = 1.0 - ((1.0 - interval)/2.)
        lower = 1.0 - upper
        if self.algorithm == "rfq":
            y_lower = [self.model_by_target[t]["model"].predict(self.empirical_sumstats, lower*100) for t in self.targets]
            y_upper = [self.model_by_target[t]["model"].predict(self.empirical_sumstats, upper*100) for t in self.targets]
        elif self.algorithm == "gb":
            y_lower = []
            y_upper = []
            for t in self.targets:
                tmp_gb = self.model_by_target[t]["model"].set_params(loss="quantile").set_params(alpha=lower)
                tmp_gb.fit(self.X, self.y[t])
                y_lower.append(tmp_gb.predict(self.empirical_sumstats))
                tmp_gb = self.model_by_target[t]["model"].set_params(loss="quantile").set_params(alpha=upper)
                tmp_gb.fit(self.X, self.y[t])
                y_upper.append(tmp_gb.predict(self.empirical_sumstats))
        else:
            print("Unsupported algorithm for prediction intervals - {}".format(self.algorithm))
            return self.empirical_pred

        ## Concatenate lower and upper quartiles onto the prediction df and name the rows nicely
        self.y_lower = pd.DataFrame(y_lower, columns=["lower {}".format(lower)], index=self.targets).T
        self.y_upper = pd.DataFrame(y_upper, columns=["upper {}".format(upper)], index=self.targets).T
        self.empirical_pred = pd.concat([self.empirical_pred, self.y_lower, self.y_upper])

        return self.empirical_pred


    ## The magic method to just do-it-all
    def predict(self, select_features=True, param_search=True, by_target=False, quick=False, verbose=False):

        ## Select relevant features
        if select_features:
            self.feature_selection(quick=quick, verbose=verbose)
       
        ## If using either the RFQuantileRegressor or GradientBoosting, then
        ## we force by_targets to be true, since they'll only do one target at at time
        ## RandomForest can handle multi-target regression, but not rfq or gb, so we
        ## allow plain rf to optionally do by_target
        if self.algorithm in ["rfq", "gb"]:
            by_target = True

        if param_search:
            self.param_search_cv(by_target=by_target, quick=quick, verbose=verbose)
        else:
            ## If you don't want to do the param search, then just take the default param
            ## Also, assuming if you don't care about setting model parameters, then you
            ## won't care about feature selection, so just use all features for the
            ## the by_target code
            if by_target:
                for t in self.targets:
                    self.model_by_target[t]["model"] = self._base_model()
                    self.model_by_target[t]["model"].fit(self.X, self.y[t])
            else:
                self.best_model = self._base_model(n_jobs=-1) 
                self.best_model.fit(self.X, self.y)

        if by_target:
            ## Predict each target independently using it's own trained RF
            preds = [self.model_by_target[t]["model"].predict(self.empirical_sumstats) for t in self.targets]
            self.empirical_pred = pd.DataFrame(np.array(preds).T, columns=self.targets, index=["estimate"])

            ## If using one of the algorithms that supports quantile regression then
            ## return the prediction intervals as well
            if self.algorithm in ["rfq", "gb"]:
                self.empirical_pred = self.prediction_interval(interval=0.95)
        else:
            ## Do all targets at once. Also, you don't get prediction intervls
            ## if you don't do by_target
            self.empirical_pred = pd.DataFrame(self.best_model.predict(self.empirical_sumstats),\
                                                columns=self.targets, index=["estimate"])

        return self.empirical_pred


####################################################################
## Convenience functions for wrapping ML parameter grid construction
####################################################################
def _get_param_grid(algorithm):
    if algorithm in ["rfq", "rf"]:
        return _rf_param_grid()
    elif algorithm == "gb":
        return _gb_param_grid()


def _rf_param_grid():
    n_estimators = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)]
    max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
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
