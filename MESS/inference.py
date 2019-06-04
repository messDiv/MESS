"""
.. module:: inference
   :synopsis: Machine learning model selection and parameter estimation.

"""

from __future__ import print_function

import MESS.stats
import matplotlib.pyplot as plt
import joblib
import numpy as np
import pandas as pd

from boruta import BorutaPy
from MESS.util import progressbar, MESSError
from skgarden import RandomForestQuantileRegressor
from sklearn import metrics
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier,\
                            GradientBoostingClassifier, GradientBoostingRegressor,\
                            AdaBoostRegressor, AdaBoostClassifier
from sklearn.model_selection import train_test_split, cross_val_predict, cross_val_score, RandomizedSearchCV
from sklearn.multioutput import MultiOutputRegressor
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline

import logging
LOGGER = logging.getLogger(__name__)

## The base Ensemble class takes care of reading in the empirical dataframe,
## calculating summary stats, reading the simulated date, and reshaping the sim
## sumstats to match the stats of the real data.
class Ensemble(object):
    """
    The Ensemble class is a parent class from which Classifiers and Regressors
    inherit shared methods. You normally will not want to create an Ensemble
    class directly, but the methods documented here are inherited by both
    Classifier() and Regressor() so may be called on either of them.

    :attention: Ensemble objects should never be created directly. It is a base class that provides functionality to Classifier() and Regressor().
    """
    def __init__(self, empirical_df, simfile, target_model=None, algorithm="rf", metacommunity_traits=None, verbose=False):
        self.empirical_df = empirical_df
        self.simfile = simfile
        self.algorithm = algorithm
        self.metacommunity_traits = metacommunity_traits

        self.set_data(empirical_df,\
                        metacommunity_traits=metacommunity_traits,\
                        verbose=verbose)

        try:
            ## Read in simulated data, split it into features and targets,
            ## and prune the features to correspond to the loaded empirical data
            ## If target_model is passed in then prune out all other sims. This
            ## is only used by the Regressor code, but it's here because this
            ## is where the data splitting happens.
            if target_model:
                self.sim_df = pd.read_csv(simfile, sep="\t", header=0)
                self.sim_df = self.sim_df[self.sim_df["community_assembly_model"] == target_model]
            else:
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


    def set_data(self, empirical_df, metacommunity_traits=None, verbose=False):
        """
        A convenience function to allow using pre-trained models to make
        predictions on new datasets without retraining the model. This will
        calculate summary statistics on input data (recycling metacommunity
        traits if these were previously input), and reshape the statistics
        to match the features selected during initial model construction.

        This is only sensible if the data from the input community consists
        of identical axes as the data used to build the model. This will be
        useful if you have community data from mutiple islands in the same
        archipelago, different communities that share a common features,
        and share a metacommunity.

        :param pandas.DataFrame empirical_df: A DataFrame containing the empirical
            data. This df has a very specific format which is documented here.
        :param array-like metacommunity_traits: A list or np.array of the trait values
            from the metacommunity. Used for calculating some of the trait based
            summary statistics.
        :param bool verbose: Print progress information.
        """
        ## Fetch empirical data and calculate the sumstats for the datatypes passed in
        self.empirical_df = empirical_df

        if metacommunity_traits is None:
            try:
                metacommunity_traits = self.metacommunity_traits
            except AttributeError:
                ## If you don't pass it in, and it isn't already set on the
                ## on the model, then just keep on truckin'.
                pass

        try:
            self.empirical_sumstats = MESS.stats.calculate_sumstats(empirical_df, metacommunity_traits=metacommunity_traits)
            if verbose: print("Got empirical summary statistics:\n{}".format(self.empirical_sumstats))
        except Exception as inst:
            raise MESSError("  Malformed input dataframe: {}".format(inst))

        try:
            ## Reshape the input sumstats to fit the previously trained model
            ## if there is one. If self.features is unset (on Ensemble.__init__)
            ## this will raise, which is fine.
            self.empirical_sumstats = self.empirical_sumstats[self.features]
        except AttributeError:
            if verbose: print("No features previously selected, using all.")


    def dump(self, outfile):
        """
        Save the model to a file on disk. Useful for saving trained models
        to prevent having to retrain them.

        :param str outfile: The file to save the model to.
        """
        joblib.dump(self, outfile)


    @staticmethod
    def load(infile):
        """
        Load a MESS.inference model from disk. This is complementary to the
        MESS.inference.Ensemble.dump() method.

        :param str infile: The file to load a trained model from.
        """
        return joblib.load(infile)


    def set_features(self, feature_list=''):
        """
        Specify the feature list to use for classification/regression. By
        default the methods use all features, but if you want to specify exact
        feature sets to use you may call this method.

        :param list feature_list: The list of features (summary statistics)
            to retain for downstream analysis. Items in this list should
            correspond exactly to summary statistics in the simulations or
            else it will complain.
        """
        if not len(feature_list):
            self.features = self._X.columns
        else:
            self.features = feature_list
        try:
            self.X = self._X[self.features]
        except Exception as inst:
            print("  Failed setting features: {}".format(inst))
            raise


    def set_targets(self, target_list=''):

        """
        Specify the target (parameter) list to use for classification/regression. By
        default the classifier will only consider :ref:`community_assembly_model`
        and the regressor will use all targets, but if you want to specify
        exact target sets to use you may call this method.

        :param list target_list: The list of targets (model parameters)
            to retain for downstream analysis. Items in this list should
            correspond exactly to parameters in the simulations or
            else it will complain.
        """

        ## If just one target then convert it to a list for compatibility
        if isinstance(target_list, str) and target_list:
            target_list = [target_list]
        ## If no targets passed in then just use the defaults (all),
        ## otherwise use the list of targets passed in.
        if not len(target_list):
            self.targets = self._default_targets
        else:
            self.targets = target_list
        ## Filter (and do not modify) the internal full targets dataset (_y).
        ## This makes it so you can easily change self.y, without reloading
        ## the data.
        try:
            self.y = self._y[self.targets]
        except Exception as inst:
            print("  Failed setting targets: {}".format(inst))
            raise


    ## The feature selection method only takes one target vector
    ## at a time, so we'll do each target seperately and then
    ## take the union of the results.
    def feature_selection(self, quick=False, verbose=False):
        """
        Access to the feature selection routine. Uses BorutaPy,
        an all-relevant feature selection method:
        https://github.com/scikit-learn-contrib/boruta_py
        http://danielhomola.com/2015/05/08/borutapy-an-all-relevant-feature-selection-method/

        :hint: Normally you will not run this on your own, but will use it indirectly through the predict() methods.

        :param bool quick: Run fast but do a bad job.
        :param bool verbose: Print lots of quasi-informative messages.
        """
        mask = np.zeros(len(self.X.columns), dtype=bool)
        if verbose: print("Selecting features:")
        for t in self.y:
            if verbose: print("  {}\t".format(t), end='')
            try:
                model = self._base_model(n_estimators=600, max_depth=5)
            except TypeError as inst:
                msg = "Ensemble model must support max_depth parameter to enable feature selection. AdaBoost can't be run with feature_selection."
                raise MESSError(msg)

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
    def predict(self, select_features=True, param_search=True, by_target=False, quick=False, force=False, verbose=False):

        try:
            _ = self.best_model
            if not force:
                if verbose: "Model already trained. Use 'force' to retrain."
                return
        except AttributeError:
            pass

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
                ## Set the best_model variable just using the model for the first target
                self.best_model = self.model_by_target[self.targets[0]]["model"]
            else:
                ## TODO: Make default base_model params smarter
                self.best_model = self._base_model(n_jobs=-1)
                self.best_model.fit(self.X, self.y)


    def feature_importances(self):
        """
        Assuming predict() has already been called, this method will return
        the feature importances of all features used for prediction.

        :return: A pandas.DataFrame of feature importances.
        """
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


    def _cv_check(self, quick=False, verbose=False):
        """
        A convenience method for checking whether Ensemble.best_model exists
        as this is what's used by both Ensemble.cross_val_* methods. If it
        doesn't exist this method will call self.predict() to generate it.
        """
        try:
            _ = self.best_model
        except:
            msg = """
Estimator has not been optimized prior to CV. Performing estimator CV skipping
feature selection and hyperparameter optimization. For better performance call
predict() on the estimator prior to calling the cv_predict/cv_score methods.
"""
            if verbose: print(msg)
            self.predict(select_features=(not quick), param_search=(not quick),\
                        quick=quick, verbose=verbose) 


    def cross_val_predict(self, cv=5, features='', quick=False, verbose=False):
        """
        Perform K-fold cross-validation prediction. For each of the `cv` folds,
        simulations will be split into sets of `K - (1/K)` training simulations
        and 1/K test simulations.

        .. note:: CV predictions are not appropriate for evaluating model
            generalizability, these should only be used for visualization and
            exploration.

        :param int cv: The number of K-fold cross-validation splits to perform.
        :param bool quick: If `True` skip feature selection and
            hyper-parameter tuning, and subset simulations. Runs fast but does
            a bad job. For testing.
        :param bool verbose: Report on progress. Depending on the number of CV
            folds this will be more or less chatty (mostly useless except for
            debugging).

        :return: The array of predicted targets for each set of features when
            it was a member of the held-out testing set. Also saves the results
            in the Estimator.cv_preds variable.
        """
        self._cv_check(quick=quick, verbose=verbose)

        ## If using rf then you don't need to use the multioutputregressor
        ## structure and in fact the RF regressor will take advantage of
        ## correlations between the targets.
        if self.algorithm != "rf":
            best_model = MultiOutputRegressor(self.best_model)
        else:
            best_model = self.best_model

        self.cv_preds = cross_val_predict(best_model, self.X, self.y, cv=cv, n_jobs=-1)
        self.cv_preds = pd.DataFrame(self.cv_preds, columns=self.targets)

        return self.cv_preds


    def cross_val_score(self, cv=5, quick=False, verbose=False):
        """
        Perform K-fold cross-validation scoring. For each of `cv` folds
        simulations will be split into sets of `K - (1/K)` training simulations
        and 1/K test simulations.

        :param int cv: The number of K-fold cross-validation splits to perform.
        :param bool quick: If `True` skip feature selection and
            hyper-parameter tuning, and subset simulations. Runs fast but does
            a bad job. For testing.
        :param bool verbose: Report on progress. Depending on the number of CV
            folds this will be more or less chatty (mostly useless except for
            debugging).

        :return: The array of scores of the estimator for each K-fold. Also
            saves the results in the Estimator.cv_scores variable.
        """
        self._cv_check(quick=quick, verbose=verbose)
        self.cv_scores = cross_val_score(self.best_model, self.X, self.y, cv=cv, n_jobs=-1)

        return self.cv_scores


    ####################################################################
    ## Plotting functions
    ####################################################################

    def plot_feature_importance(self,\
                                cutoff=0.05,\
                                figsize=(10, 12),\
                                layout=None,\
                                subplots=True,\
                                legend=False):
        """
        Construct a somewhat crude plot of feature importances, useful for a
        quick and dirty view of these values. If more than one feature present
        in the model then a grid-layout is constructed and each individual
        feature is displayed within a subplot. This function is a thin wrapper
        around pandas.DataFrame.plot.barh().

        :param float cutoff: Remove any features that do not have greater importance
            than this value across all plotted features. Just remove uninteresting
            features to reduce the amount of visual noise in the figures.
        :param tuple figsize: A tuple specifying figure width, height in inches.
        :param tuple layout: A tuple specifying the row, column layout of the
            sub-panels. By default we do our best, and it's normally okay.
        :param bool subplots: Whether to plot each feature individually, or just
            cram them all into one huge plot. Unless you have only a few features,
            setting this option to `False` will look insane.
        :param bool legend: Whether to plot the legend.

        :return: Returns all the matplotlib axes
        """
        imps = self.feature_importances()
        ## Drop any features that are relatively unimportant
        imps = imps[imps.columns[imps.max() > cutoff]]

        ## Just do our best here to figure out a good layout. These should
        ## work in most normal cases, but you an fall
        if not layout is None:
            if len(imps) == 1:
                layout = (1, 1)
            else:
                layout = (2, len(imps)/2)
        axs = imps.T.plot.barh(figsize=figsize, subplots=subplots, legend=legend, width=0.9, layout=layout, sharey=True)
        return axs


class Classifier(Ensemble):
    """
    This class wraps all the model selection machinery.

    :param pandas.DataFrame empirical_df: A DataFrame containing the empirical
        data. This df has a very specific format which is documented here.
    :param string simfile: The path to the file containing all the simulations.
    :param string algorithm: One of the :ref:`ensemble_methods` to use for
        parameter estimation.
    :param array-like metacommunity_traits: A list or np.array of the trait values
        from the metacommunity. Used for calculating some of the trait based
        summary statistics.
    :param bool verbose: Print detailed progress information.
    """


    _default_targets = ["community_assembly_model"]

    def __init__(self, empirical_df, simfile, algorithm="rf", metacommunity_traits=None, verbose=False):
        super(Classifier, self).__init__(empirical_df, simfile, algorithm=algorithm, metacommunity_traits=metacommunity_traits, verbose=verbose)

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


    def predict(self, select_features=True, param_search=True, by_target=False, quick=False, force=False, verbose=False):
        """
        Predict the most likely community assembly model class.

        :param bool select_features: Whether to perform relevant feature selection.
            This will remove features with little information useful for model
            prediction. Should improve classification performance, but does take
            time.
        :param bool param_search: Whether to perform ML classifier hyperparameter
            tuning. If ``False`` then classification will be performed with default
            classifier options, which will almost certainly result in poor performance,
            but it will run really fast!.
        :param bool by_target: Whether to predict multiple target variables
            simultaneously, or each individually and sequentially.
        :param bool quick: Reduce the number of retained simulations and the number
            of feature selection and hyperparameter tuning iterations to make the
            prediction step run really fast! Useful for testing.
        :param bool force: Force re-running feature selection and hyper-parameter
            tuning. This is basically here to prevent you from shooting yourself
            in the foot inside a for loop with `select_features=True` when
            really what you want (most of the time) is to just run this once,
            and call predict() multiple times without redoing this. 
        :param bool verbose: Print detailed progress information.

        :return: A tuple including the predicted model and the probabilities per model class.
        """
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
            self.empirical_proba = pd.DataFrame(proba, columns=list(self.model_by_target.values())[0]["model"].classes_, index=self.targets)
        else:
            ## Do all targets at once. Also, you don't get prediction intervls
            ## if you don't do by_target
            self.empirical_pred = pd.DataFrame(self.best_model.predict(self.empirical_sumstats),\
                                                columns=self.targets, index=["estimate"])
            self.empirical_proba = pd.DataFrame(self.best_model.predict_proba(self.empirical_sumstats),\
                                                columns=self.best_model.classes_, index=self.targets)

        return self.empirical_pred, self.empirical_proba


    def plot_confusion_matrix(self,\
                                ax='',\
                                figsize=(8, 8),\
                                cmap=plt.cm.Greys,\
                                cbar=False,\
                                title="",\
                                normalize=False,
                                outfile=''):
        """
        Plot the confusion matrix for CV predictions. Assumes
        `Classifier.cross_val_predict()` has been called. If not it complains
        and tells you to do that first.

        :param matplotlib.pyploat.axis ax: The matplotlib axis to draw the plot
            on.
        :param tuple figsize: If not passing in an axis, specify the size of
            the figure to plot.
        :param matplotlib.pyplot.cm cmap: Specify the colormap to use.
        :param bool cbar: Whether to add a colorbar to the figure.
        :param str title: Add a title to the figure.
        :param bool normalize: Whether to normalize the bin values (scale to
            1/# simulations).
        :param str outfile: Where to save the figure. This parameter should
            include the desired output file format, e.g. `.png`, `.svg` or
            `.svg`.
        """
        try:
            _ = self.cv_preds
        except AttributeError:
            msg = "No CV predictions. You must call `Classifier.cross_val_predict()` first"
            raise MESSError(msg)

        if not ax:
            fig, ax = plt.subplots(figsize=figsize)

        ## By default the confusion_matrix() function lays out the matrix based
        ## on the alpha-sorted order of labels found in the data.
        labels = ["Competition", "Filtering", "Neutral"]

        conf_matrix = metrics.confusion_matrix(self.y, self.cv_preds)

        if normalize:
            conf_matrix = conf_matrix.astype(float)/conf_matrix.sum(axis=1)

        im = ax.imshow(conf_matrix, interpolation='nearest', cmap=cmap)
        if cbar:
            ax.figure.colorbar(im, ax=ax)

        ax.set_xticks(np.arange(conf_matrix.shape[1]))
        ax.set_yticks(np.arange(conf_matrix.shape[0]))
        ax.set_xticklabels(labels, fontsize=15)
        ax.set_yticklabels(labels, fontsize=15)
        ax.set_title(title, fontsize=25)
        ax.set_ylabel('True assembly model class', fontsize=20)
        ax.set_xlabel('Predicted assembly model class', fontsize=20)

        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                rotation_mode="anchor")

        fmt = '.2f' if normalize else 'd'
        thresh = conf_matrix.max() / 2.
        for i in range(conf_matrix.shape[0]):
            for j in range(conf_matrix.shape[1]):
                ax.text(j, i, format(conf_matrix[i, j], fmt),
                        ha="center", va="center", fontsize=15,
                        color="white" if conf_matrix[i, j] > thresh else "black")

        fig.tight_layout()

        if outfile:
            plt.savefig(outfile)

        return ax


class Regressor(Ensemble):
    """
    This class wraps all the parameter estimation machinery.

    :param pandas.DataFrame empirical_df: A DataFrame containing the empirical
        data. This df has a very specific format which is documented here.
    :param string simfile: The path to the file containing all the simulations.
    :param string algorithm: The ensemble method to use for parameter estimation.
    :param string target_model: The community assembly model to specifically use.
        If you include this then the simulations will be read and then filtered
        for only this `community_assembly_model`.
    :param array-like metacommunity_traits: A list or np.array of the trait values
        from the metacommunity. Used for calculating some of the trait based
        summary statistics.
    :param bool verbose: Print lots of status messages. Good for debugging,
        or if you're *really* curious about the process.
    """

    _default_targets = ["alpha", "S_m", "J_m", "speciation_rate", "death_proportion",\
                        "trait_rate_meta", "ecological_strength", "J", "m",\
                        "generation", "speciation_prob", "_lambda"]

    def __init__(self, empirical_df, simfile, target_model=None, algorithm="rfq", metacommunity_traits=None, verbose=False):
        super(Regressor, self).__init__(empirical_df, simfile, target_model=target_model, algorithm=algorithm, metacommunity_traits=metacommunity_traits, verbose=False)

        if algorithm == "rf":
            self._base_model = RandomForestRegressor
        elif algorithm == "rfq":
            self._base_model = RandomForestQuantileRegressor
        elif algorithm == "gb":
            self._base_model = GradientBoostingRegressor
        elif algorithm == "ab":
            self._base_model = AdaBoostRegressor
        else:
            raise Exception(" Unsupported regression algorithm: {}".format(algorithm))
        self.algorithm = algorithm
        self._param_grid = _get_param_grid(algorithm)

        ## Remove invariant targets (save time)
        self.y = self.y.loc[:, (self.y != self.y.iloc[0]).any()]
        self.targets = list(self.y.columns)
        if verbose: print("Removed invariant targets. Retained: {}".format(list(self.y.columns)))

    def prediction_interval(self, interval=0.95, quick=False, verbose=False):
        """
        Add upper and lower prediction interval for algorithms that support
        quantile regression (`rf`, `gb`).

        :hint: You normaly won't have to call this by hand, as it is incorporated automatically into the predict() methods. We allow access to in for experimental purposes.

        :param float interval: The prediction interval to generate.
        :param bool quick: Subsample the data to make it run fast, for testing.
            The `quick` parameter doesn't do anything for `rf` because it's
            already really fast (the model doesn't have to be refit).
        :param bool verbose: Print information about progress.

        :return: A pandas.DataFrame containing the model predictions and the
            prediction intervals.
        """
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


    def predict(self, select_features=True, param_search=True, by_target=False, quick=False, force=True, verbose=False):
        """
        Predict parameter estimates for selected targets.

        :param bool select_features: Whether to perform relevant feature selection.
            This will remove features with little information useful for parameter
            estimation. Should improve parameter estimation performance, but does
            take time.
        :param bool param_search: Whether to perform ML regressor hyperparamter
            tuning. If ``False`` then prediction will be performed with default
            options, which will almost certainly result in poor performance,
            but it will run really fast!.
        :param bool by_target: Whether to estimate all parameters simultaneously,
            or each individually and sequentially. Some ensemble methods are only
            capable of performing individual parameter estimation, in which case
            this parameter is forced to ``True``.
        :param bool quick: Reduce the number of retained simulations and the number
            of feature selection and hyperparameter tuning iterations to make the
            prediction step run really fast! Useful for testing.
        :param bool force: Force re-running feature selection and hyper-parameter
            tuning. This is basically here to prevent you from shooting yourself
            in the foot inside a for loop with `select_features=True` when
            really what you want (most of the time) is to just run this once,
            and call predict() multiple times without redoing this. 
        :param bool verbose: Print detailed progress information.

        :return: A pandas DataFrame including the predicted value per target
            parameter, and 95% prediction intervals if the ensemble method
            specified for this Regressor supports it.
        """
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


    def cross_val_predict(self, cv=5, quick=False, verbose=False):
        """
        A thin wrapper around Ensemble.cross_val_predict() that basically just
        calculates some Regressor specific statistics after the cross validation
        prodecure. This function will calculate and populate class variables:

        Regressor.MAE: Mean absolute error
        Regressor.RMSE: Root mean squared error
        Regressor.vscore: Explained variance score
        Regressor.r2: Coefficient of determination regression score

        As well as Regressor.cv_stats which is just a pandas.DataFrame of the
        above stats.

        :param int cv: The number of cross-validation folds to perform.
        :param bool quick: Whether to downsample to run fast but do a bad job.
        :param bool verbose: Whether to print progress messages.
        """
        super(Regressor, self).cross_val_predict(cv=cv, quick=quick, verbose=verbose)        

        ## Calculate mean absolute error and root mean squared error of estimates
        self.MAE = pd.Series(np.mean(np.abs((self.cv_preds - self.y))), name="MAE")
        self.RMSE = pd.Series(np.sqrt(np.mean((self.cv_preds - self.y)**2)/len(self.y)), name="RMSE")

        ## Here multioutput='raw_values' returns a score per target, rather
        ## aggregating all scores into one value
        self.vscore = metrics.explained_variance_score(self.y,\
                                                        self.cv_preds,\
                                                        multioutput='raw_values')
        self.vscore = pd.Series(self.vscore, index=self.targets, name="vscore")

        self.r2 = metrics.r2_score(self.y,\
                                    self.cv_preds,\
                                    multioutput='raw_values')
        self.r2 = pd.Series(self.r2, index=self.targets, name="R2")

        ## Write it all out as a DataFrame for convenience
        self.cv_stats = pd.concat([self.MAE, self.RMSE, self.vscore, self.r2], axis=1)

        if verbose:
            print("MAE:\n{}\nRMSE:\n{}\nvscores:\n{}\nR2:\n{}".format(self.MAE,\
                                                                        self.RMSE,\
                                                                        self.vscore,\
                                                                        self.r2))
        return self.cv_preds


    def plot_cv_predictions(self,\
                                ax='',\
                                figsize=(10, 5),\
                                figdims=(2, 3),\
                                n_cvs=1000,\
                                title="",\
                                outfile=''):
        """
        Plot the cross validation predictions for this Regressor. Assumes
        `Regressor.cross_val_predict()` has been called. If not it complains
        and tells you to do that first.

        :param matplotlib.pyploat.axis ax: The matplotlib axis to draw the plot
            on.
        :param tuple figsize: If not passing in an axis, specify the size of
            the figure to plot.
        :param tuple figdims: The number of rows and columns (specified in that
            order) of the output figure. There will be one plot per target 
            parameter, so there should be at least as many available cells in
            the specified grid.
        :param int n_cvs: The number of true/estimated points to plot on the
            figure.
        :param str title: Add a title to the figure.
        :param str outfile: Where to save the figure. This parameter should
            include the desired output file format, e.g. `.png`, `.svg` or
            `.svg`.

        :return: The flattened list of matplotlib axes on which the scatter
            plots were drawn, one per target.
        """
        try:
            _ = self.cv_preds
        except AttributeError:
            msg = "No CV predictions. You must call `Regressor.cross_val_predict()` first"
            raise MESSError(msg)

        fig, axs = plt.subplots(figdims[0], figdims[1], figsize=figsize)
        axs = axs.flatten()

        for t, ax in zip(self.targets, axs):
            self._plot_cv_prediction(t, ax, n_cvs=n_cvs)

        fig.tight_layout()
        plt.suptitle(title)

        if outfile:
            plt.savefig(outfile)

        return axs


    def _plot_cv_prediction(self,\
                            target,\
                            ax='',\
                            n_cvs=1000):
        """
        An internal function for plotting just one set of cv predictions
        for one target parameter. Normally you'll want to use the
        Regressor.plot_cv_predictions() method which wraps calls to this.

        :param str target: The target parameter of the set of CVs to plot.
        :param matplotlib.pyplot.axis ax: The axis to draw to, otherwise
            create a new one.
        :param int n_cvs: The number of true/estimated points to plot.
        """
        if not ax:
            fig, ax = plt.subplots()

        if target == "ecological_strength":
            y_true = np.log(self.y[target])
        else:
            y_true = self.y[target]

        ax.scatter(y_true[:n_cvs], self.cv_preds[target][:n_cvs], c='black', marker='.', s=2)

        ax.set_title(target)
        if target in ["m", "speciation_prob"]:
            ax.set_xlim(np.min(y_true), np.max(y_true))
            ax.set_ylim(np.min(y_true), np.max(y_true))


#############################
## Module methods
#############################

def parameter_estimation_cv(simfile, outdir='', target_model=None, data_axes='',
                            algorithm="rf", quick=True, verbose=False):

    ## Make outdir if it doesn't exist.
    if not os.path.exists(model_outdir):
        os.mkdir(model_outdir)
    importance_out_df = pd.DataFrame([], index=feature_sets["all"])

    if not data_axes:
        data_axes = ["abundance", "pi", "dxy", "trait"]

    ## Generate a synthetic community and filter only data axes specified.
    synthetic_community = MESS.util.synthetic_community()[data_axes]

    rgr = MESS.inference.Regressor(empirical_df=synthetic_community,\
                                    simfile=simfile,\
                                    target_model=target_model,\
                                    verbose=verbose)

    ## Feature selection is independent of the real data, so it doesn't matter
    ## that we passed in synthetic empirical data here. It chooses features
    ## that are only relevant for each summary statistic.
    rgr.feature_selection(quick=quick, verbose=verbose)
    ## Searching for the best model hyperparameters is the same, it is done
    ## independently of the observed data.
    rgr.param_search_cv(by_target=rgr._by_target, quick=quick, verbose=verbose)

    
    cv_preds = cross_val_predict(rgr.best_model, rgr.X, rgr.y, cv=5, n_jobs=-1)


def posterior_predictive_check(empirical_df,\
                                parameter_estimates,\
                                ax='',\
                                est_only=False,\
                                nsims=100,\
                                outfile='',\
                                use_lambda=True,\
                                verbose=False):
    """
    Perform posterior predictive simulations. This function will take
    parameter estimates and perform MESS simulations using these parameter
    values. It will then plot the resulting summary statistics in PC
    space, along with the summary statistics of the observed data. The
    logic of posterior predictive checks is that if the estimated parameters
    are a good fit to the data, then summary statistics generated using
    these parameters should resemble those of the real data.

    :param pandas.DataFrame empirical_df: A DataFrame containing the empirical
        data. This df has a very specific format which is documented here.
    :param pandas.DataFrame parameter_estimates: A DataFrame containing the
        the parameter estimates from a MESS.inference.Regressor.predict() call
        and optional prediction interval upper and lower bounds.
    :param bool ax: The matplotlib axis to use for plotting. If not specified
        then a new axis will be created.
    :param bool est_only: If True, drop the lower and upper prediction
        interval (PI) and just use the mean estimated parameters for generating
        posterior predictive simulations. If False, and PIs exist, then
        parameter values will be sampled uniformly between the lower and upper
        PI.
    :param bool nsims: The number of posterior predictive simulations to perform.
    :param bool outfile: A file path for saving the figure. If not specified
        the figure is simply not saved to the filesystem.
    :param bool use_lambda: Whether to generated simulations using time
        as measured in _lambda or in generations.
    :param bool verbose: Print detailed progress information.

    :return: A matplotlib axis containing the plot.
    """

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
    ## FIXME: Oh boy, this is not parallelized at all. Fix this.
    for i in range(nsims):
        for param in parameter_estimates:
            ## Only use one of either the estimated lambda values
            ## or the time in generations, not both.
            if param == "_lambda":
                if use_lambda == True:
                    r.set_param("generations", param_df[param].iloc[i])
                else: pass
            elif param == "generation":
                if use_lambda == False:
                    r.set_param("generations", param_df[param].iloc[i])
            else:
                r.set_param(param, param_df[param].iloc[i])
        r.run(sims=1, quiet=True)

        if verbose: progressbar(nsims, i, "Performing simulations")
    if verbose: progressbar(nsims, nsims, "Performing simulations")
    if verbose: print("\nCalculating PCs and plotting")

    simfile = "./ppc/SIMOUT.txt"
    sim_df = pd.read_csv(simfile, sep="\t", header=0)

    ## Chop the sims down to only include ss contained in the observed data
    obs_ss = MESS.stats.calculate_sumstats(empirical_df)
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
    return random_grid


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
