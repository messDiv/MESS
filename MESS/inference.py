from __future__ import print_function

import MESS.stats
import numpy as np
import pandas as pd

from boruta import BorutaPy
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

        ## Take the default targets
        self.set_targets()
        ## Set features to correspond to real data
        self.set_features(self.empirical_sumstats.columns)


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

    def __init__(self, empirical_df, simfile, algorithm="rf", verbose=False):
        super(Regressor, self).__init__(empirical_df, simfile, algorithm="rf", verbose=False)

        ## If you want to estimate parameters independently then
        ## we keep track of which features are most relevant per target
        self.relevant_features_per_target = {x:[] for x in self.y.columns}

        if algorithm == "rf":
            self._base_model = RandomForestRegressor
        elif algorithm == "gb":
            self._base_model = GradientBoostingRegressor
        self._param_grid = _get_param_grid(algorithm)


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
            tmpX = self.X.values
            tmpy = self.y[t].values
            model = self._base_model(n_jobs=-1, n_estimators=600, max_depth=5)

            if quick:
                ## Subset the data and run fewer iterations
                tmpX = tmpX[:1000]
                tmpy = tmpy[:1000]
                max_iter=10
            else:
                max_iter=100

            # define Boruta feature selection method
            feat_selector = BorutaPy(model, max_iter=max_iter, n_estimators='auto', verbose=0, random_state=1)
            ## Turn off np warnings
            np.warnings.filterwarnings('ignore')
            feat_selector.fit(tmpX, tmpy)

            # check ranking of features
            if verbose: print("{}".format(list(self.features[feat_selector.support_])))
            ## Use feat_selector.support_weak_ to include more variables
            self.relevant_features_per_target[t] = list(self.features[feat_selector.support_])
            mask += feat_selector.support_
            #X_filtered = pd.DataFrame(X_filtered, columns = X.columns[feat_selector.support_])

        selected_features = self.features[mask]
        if verbose: print("All selected features: {}".format(" ".join(selected_features)))
        ## Filter on the selected features
        self.X = self.X[selected_features]
        ## Also filter the empirical sumstats
        self.empirical_sumstats = self.empirical_sumstats[selected_features]

        self.features = selected_features


    def param_search_cv(self, quick=False, verbose=False):

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
        cvsearch.fit(tmpX, tmpy)
        if verbose: print(cvsearch.best_params_)
        self._cvsearch = cvsearch
        self.best_model = cvsearch.best_estimator_


    ## The magic method to just do-it-all
    def predict(self, select_features=True, param_search=True, quick=False, verbose=False):

        ## Select relevant features
        if select_features:
            self.feature_selection(quick=quick, verbose=verbose)
       
        if param_search:
            self.param_search_cv(quick=quick, verbose=verbose)
        else:
            ## If you don't want to do the param search, then just take the default params
            self.best_model = model(n_jobs=-1) 

        self.empirical_pred = pd.DataFrame(self.best_model.predict(self.empirical_sumstats), columns=self.targets)
        print(self.empirical_pred)


####################################################################
## Convenience functions for wrapping ML parameter grid construction
####################################################################
def _get_param_grid(algorithm):
    if algorithm == "rf":
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
