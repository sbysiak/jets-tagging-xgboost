from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
# import xgboost as xgb
#
# from sklearn.grid_search import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
# from sklearn.metrics import accuracy_score, recall_score, precision_score, roc_curve, roc_auc_score, auc
# from sklearn.metrics import make_scorer
#
# from sklearn.neural_network import MLPClassifier
# from sklearn.neighbors import KNeighborsClassifier
# from sklearn.svm import SVC
# from sklearn.tree import DecisionTreeClassifier
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.naive_bayes import GaussianNB

import keras as K
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Convolution2D, MaxPooling2D
from keras.layers import Conv1D, MaxPooling1D
from keras.utils import np_utils
from keras.datasets import mnist
from keras import optimizers
from keras.callbacks import Callback, TensorBoard
from keras.wrappers.scikit_learn import KerasClassifier



# for sake of google.colab
try:
    import ROOT
    import root_numpy
    import root_pandas
except:
    print 'ROOT, root_numpy, root_pandas has NOT been imported'

import os
import os.path
from time import time
from collections import OrderedDict
# os.getcwd()

# from IPython.core.interactiveshell import InteractiveShell
# InteractiveShell.ast_node_interactivity = "all"

# suppress ConvergenceWarnings
# import warnings
# from sklearn.exceptions import ConvergenceWarning
# warnings.filterwarnings('ignore',category=ConvergenceWarning)
# warnings.filterwarnings('ignore',category=DeprecationWarning)

pd.set_option('precision', 4)


#
#
#
#
#


# # # # # # # # # # # # # # # # # # # # # # # #
# Reading ROOT trees in various formats
#

def read_msc(fnames=['jetsTree_hadd_LHC17f8g_17_255540.root',],
             n_rows=1000,
             observables=None,
             add_variables=['E', 'Pt', 'Phi', 'Eta'],
             target='tagExp',
             verbose=True):
    """function reading ROOT trees
    (format used in msc-thesis: many tags True/Exp)

    Parameters
    ----------
    fnames : string or array_like of strings
        list of files' names
    nrows : int or None, default=1000
        no. of rows to be read from each file
        if None, then reads all
    observables : array_like of strings or None, default=None
        names of branch leafes to be read from jetsObservablesBranch,
        if None reads all
    add_variables : array_like of strings, default=['E', 'Pt', 'Phi', 'Eta']
        names of variables to be read from jetsBranch (jet as instance of TParticle)
    target : string, default='tagExp'
        name of variable to be used as the target with name 'target'
    verbose : bool, default=True

    Returns
    -------
    tree : pandas DataFrame
        tabular form of tree
    """
    start = time()

    if isinstance(fnames, str): # enabling single string as input
        fnames = [fnames, ]
    if not observables:
        observables = [v[0] for v in root_numpy.list_structures(fnames[0])['jetObservablesBranch']]
    add_variables = ['jetsBranch.'+var+'()' for var in add_variables]
    trees_lst = []
    for fname in fnames:
        nptree = root_numpy.root2array(fname, branches=observables+add_variables, stop=n_rows)
        pdtree = pd.DataFrame(nptree)
        trees_lst.append(pdtree)

    tree = pd.concat(trees_lst)
    tree.index = range(len(tree))

    tree.columns = [col.replace('jetsBranch.', '') for col in tree.columns]
    tree.columns = [col.replace('()', '') for col in tree.columns]
    tree['target'] = tree[target]
    tree = tree.drop([c for c in tree.columns if 'tag' in c], axis=1)

    if verbose:
        print 'read_msc: {} rows loaded, exec. time: {} sec'.format(len(tree), time()-start)
        #print pdtree.head()
    return tree



def read_inz(fnames=['INZYNIERKA/pythia/OUTPUTS/treeOfJets_q5_en200_nev50k.root',
                     'INZYNIERKA/pythia/OUTPUTS/treeOfJets_g_en200_nev100k.root'] ,
             n_rows=1000,
             observables=None,
             add_variables=['E', 'Pt', 'Phi', 'Eta'],
             verbose=True):
    """function reading ROOT trees
    (format used in inz-thesis: most basic,
    each flavour in separate file, included in file name)

    Parameters
    ----------
    fnames : string or array_like of strings
        list of files' names, target bases only on it
        must contain 'Jets_g/qN_en', e.g. 'treeOfJets_q5_en200.root'
    nrows : int or None, default=1000
        no. of rows to be read from each file
        if None, then reads all
    observables : array_like of strings or None, default=None
        names of branch leafes to be read from jetsObservablesBranch,
        if None reads all
    add_variables : array_like of strings, default=['E', 'Pt', 'Phi', 'Eta']
        names of variables to be read from jetsBranch (jet as instance of TParticle)
    verbose : bool, default=True

    Returns
    -------
    tree : pandas DataFrame
        tabular form of tree
    """
    start = time()

    if isinstance(fnames, str): # enabling single string as input
        fnames = [fnames, ]
    if not observables: observables = ['ptRel', 'radMom', 'angular', 'svR', 'svN', 'mult', 'eIn']
    add_variables = ['jetsBranch.'+var+'()' for var in add_variables]
    #n_rows = int(n_rows/2)
    trees_lst = []
    for fname in fnames:
        nptree = root_numpy.root2array(fname, branches=observables+add_variables, stop=n_rows)
        pdtree = pd.DataFrame(nptree)
        target = 0 if 'Jets_g_en' in fname else int(fname.split('_')[1][-1])
        print fname,target
        pdtree['target'] = np.ones_like(pdtree.index)*target
        trees_lst.append(pdtree)

    tree = pd.concat(trees_lst)
    tree.index = range(len(tree))

    tree.columns = [col.replace('jetsBranch.', '') for col in tree.columns]
    tree.columns = [col.replace('()', '') for col in tree.columns]

    if verbose:
        print 'read_inz: {} rows loaded, exec. time: {} sec'.format(len(tree), time()-start)
        #print pdtree.head()
    return tree



def read_jet_extractor_observables(fnames=['../alianalysis-tasks/JetExtractor_HFCJ/jetObserv.root', ],
                       n_rows=1000,
                       observables=None,
                       flavours=['light', 'b'],
                       verbose=True):
    """function reading ROOT trees (format of jetObserv generated from JetExtractor task)

    Parameters
    ----------
    fnames : string or array_like of strings
        list of files' names
    nrows : int or None, default=1000
        no. of rows to be read from each file
        if None, then reads all
    observables : array_like of strings or None, default=None
        names of branch leafes to be read from jetsObservablesBranch,
        if None reads all
    flavours : string or array_like of strings, default=['light', 'b']
        parton flavours to be read. Available: 'light', 'b', 'c'
    verbose : bool, default=True

    Returns
    -------
    tree : pandas DataFrame
        tabular form of tree
    """
    start = time()

    if isinstance(fnames, str): # enabling single string as input
        fnames = [fnames, ]
    if isinstance(flavours, str): # enabling single string as input
        flavours = [flavours, ]
    trees_flavour_lst = []
    for flavour in flavours:
        if not observables: observables = [v[0] for v in root_numpy.list_structures(fnames[0], 'light/jetObservTree')['jetObservablesBranch']]
        trees_lst = []
        for fname in fnames:
            nptree = root_numpy.root2array(fname, flavour+'/jetObservTree', branches=observables, stop=n_rows)
            pdtree = pd.DataFrame(nptree)
            map_flavour = {'light':0, 'b':5, 'c':4}
            pdtree['target'] = np.ones_like(pdtree[observables[0]])*map_flavour[flavour]
            trees_lst.append(pdtree)

        tree_flavour = pd.concat(trees_lst)
        tree_flavour.index = range(len(tree_flavour))

        tree_flavour = tree_flavour.drop([c for c in tree_flavour.columns if 'tag' in c], axis=1)
        trees_flavour_lst.append(tree_flavour)

    tree = pd.concat(trees_flavour_lst)
    tree.index = range(len(tree))

    if verbose:
        print 'read_jet_extractor: {} rows loaded, exec. time: {} sec'.format(len(tree), time()-start)
        #print pdtree.head()
    return tree



def read_branch_names(tree, branches_lst=[], verbose=False, indent=0):
    """ help function for reading tree with custom format / custom classes
    reads nested branches' names
    """
    for branch in tree.GetListOfBranches():
        if len(branch.GetListOfBranches()) == 0:
            if verbose: print '\t'*indent + branch.GetName()
            branches_lst.append(branch.GetName())
        else:
            if verbose: print'\t'*indent + branch.GetName(), ' ::'
            branches_lst = read_branch_names(branch, branches_lst=branches_lst, verbose=verbose, indent=indent+1)
    return branches_lst



def read_jet_extractor(fnames=['../benchmark-dataset/DATA/trains/train_LHC17f8g_11.root', ],
                       n_rows=1000,
                       observables=None,
                       flavours=['light', 'b'],
                       selection=None,
                       get_entries_only=False,
                       return_ttrees=False,
                       verbose=True):
    """function reading ROOT trees (format generated by JetExtractor task)

    Parameters
    ----------
    fnames : string or array_like of strings
        list of files' names
    nrows : int or None, default=1000
        no. of rows to be read from each file
        if None, then reads all
    observables : array_like of strings or None, default=None
        names of branch leafes to be read from jetsObservablesBranch,
        if None reads all extracted with read_branch_names()
    flavours : string or array_like of strings, default=['light', 'b']
        parton flavours to be read. Available: 'light', 'b', 'c'
    selection : string
        selection condition passed to root_numpy.tree2array()
    get_entries_only : bool
        will only print number of entries in each tree
    return_ttrees : bool
        will return ROOT.TTree objects
    verbose : bool, default=True

    Returns
    -------
    tree : pandas DataFrame
        tabular form of tree
    """
    start = time()

    if isinstance(fnames, str): # enabling single string as input
        fnames = [fnames, ]
    if isinstance(flavours, str): # enabling single string as input
        flavours = [flavours, ]

    trees_flavour_lst = []
    for flavour in flavours:
        trees_lst = []
        for fname in fnames:
            tfile = ROOT.TFile(fname)
            tdir_file = tfile.Get('ChargedJetsHadronCF')
            tlist = tdir_file.Get('AliAnalysisTaskJetExtractorHF_Jet_AKTChargedR040_tracks_pT0150_E_scheme_RhoR020KT_{}Jets_histos;1'.format(flavour))
            jettree = tlist.FindObject('ExtractedJets')
            if get_entries_only:
                print fname, flavour, jettree.GetEntries()
                continue
            if return_ttrees:
                trees_flavour_lst.append(jettree)
                continue
            if not observables: observables = read_branch_names(jettree, [])
            nptree = root_numpy.tree2array(jettree, branches=observables, stop=n_rows, selection=selection)
            pdtree = pd.DataFrame(nptree)
            map_flavour = {'light':0, 'b':5, 'c':4}
            pdtree['target'] = np.ones_like(pdtree[observables[0]])*map_flavour[flavour]
            trees_lst.append(pdtree)

        if get_entries_only or return_ttrees: continue

        tree_flavour = pd.concat(trees_lst)
        tree_flavour.index = range(len(tree_flavour))
        tree_flavour = tree_flavour.drop(['fMotherInitialCollision','fMotherHadronMatching'], axis=1)
        trees_flavour_lst.append(tree_flavour)

    if get_entries_only: return None
    if return_ttrees: return trees_flavour_lst

    tree = pd.concat(trees_flavour_lst)
    tree.index = range(len(tree))

    if verbose:
        print 'read_jet_extractor: {} rows loaded, exec. time: {} sec'.format(len(tree), time()-start)
    return tree



# # # # # # # # # # # # # # # # # # # # # # # #
# Data processing
#



def prepare_dataset(paramtree,
                 light_heavy_ratio=1,
                 light_heavy_targets=[0,5],
                 cut_var=None, cut_low=0, cut_high=999,
                 train_size=0.6,
                 test_size=None,
                 valid_size=None,
                 standarize_scales=True,
                 shuffle=True,
                 return_flavours=False,
                 verbose=True):
    """converts input DataFrame into (ready to be fed into algo) train/test datasets
    performing a couple of preprocessing steps

    Parameters:
    -----------
    paramtree : pandas DataFrame
        data in tabular format (e.g. read by read:_msc/_inz/_jet_extractor)
    light_heavy_ratio : numeric, default=1
        light to heavy jets ratio in output datasets,
        use value other than 0 to create unbalanced dataset
    light_heavy_targets : array_like of int, default=[0,5]
        targets values (labels) for light and heavy jets in input
        (will be mapped to 0/1)
    cut_var : string or array_like of strings or None, default=None
        name of variable/feature on which cuts should be performed
        if None, then no cuts will be performed
    cut_low : numeric or array_like of numerics, default=0
        lower limits for consecutive cuts, ignored if cut_var == None
    cut_high : numeric or array_like of numerics, default=999
        higher limits for consecutive cuts, ignored if cut_var == None
    train_size, test_size, valid_size : float, int, None
        values passed directly to sklearn's train_test_split
        if ``bool(valid_size)`` is True, then two splits are performed:
        first with train_size and test_size passed,
        second on train_subset with valid_size as test_size parameter.
        Otherwise only first split is performed.
        E.g. train/test/valid _size = 0.5,  None, 0.3
             fractions will be        0.35, 0.5,  0.15
        See: http://scikit-learn.org/stable/modules/generated/sklearn.model_selection.train_test_split.html
    standarize_scales : bool, default=True
        if scaling using StandardScaler should be performed
    shuffle : bool, default=True
        if shuffling of light flavour jets should be performed
        or simply first <demanded number> of jets should be selected
    return_flavours : bool, default=False
        if flvours should be returned as 3rd array for each datasubset: train/test/val
    verbose : bool, default=True

    Returns:
    --------
        (X_train, y_train, X_test, y_test, X_val, y_val) : numpy arrays
        self-explanatory datasets,
        X_val, y_val can be empty lists if only first split (to train/test) was performed

        (X_train, y_train, F_train, X_test, y_test, F_test, X_val, y_val, F_val) : numpy arrays
        if ``return_flavours`` is True then also 3rd arrays is returned (F_xxx),
        which contain true flavours or jets (from 'flavour' column in ``paramtree``)
    """



    pdtree = paramtree.copy()
    pdtree.index = range(len(pdtree))  # for cases when paramtree is result of concactination

    len0 = len(pdtree)
    pdtree = pdtree.dropna()
    len_noNA = len(pdtree)
    if cut_var:
        # multiple (try:) or single (except:) cut
        try:
            len_cuts = []
            for c_var, c_low, c_high in zip(cut_var, cut_low, cut_high):
                pdtree = pdtree.loc[(pdtree[c_var] >= c_low) & (pdtree[c_var] <= c_high)]
                len_cuts.append(len(pdtree))
        except TypeError:
            pdtree = pdtree.loc[(pdtree[cut_var] >= cut_low) & (pdtree[cut_var] <= cut_high)]
            len_cuts = len(pdtree)
        msg_suffix = 'after cuts: {}'.format(len_cuts)
    else: msg_suffix = 'no cuts applied'
    if verbose:
        print 'no. jets available: \tall: {}, after dropNA: {}, '.format(len0,len_noNA) + msg_suffix

    light_target, heavy_target = light_heavy_targets
    df_heavy = pdtree.loc[pdtree['target'] == heavy_target]
    df_light = pdtree.loc[pdtree['target'] == light_target]
    # shuffle
    idx = df_light.index.tolist()
    if shuffle:
        np.random.shuffle(idx)
    # light-to-heavy ratio
    idx = idx[:len(df_heavy)*light_heavy_ratio]
    df = pd.concat([df_light.loc[idx], df_heavy])

    # project into X,y
    y = df['target']
    X = df.drop([c for c in df.columns if 'tag' in c], axis=1)
    X = X.drop(['target'], axis=1)
    if return_flavours:
        if 'flavour' in X.columns:
            flavours = X['flavour']
            X = X.drop(['flavour'], axis=1)
        else:
            print '\nWARNING: \"flavour\" column absent in DataFrame passed\n'
            return_flavours = False

    # map light vs heavy
    y = y.map(lambda y: 0 if y == light_target else 1)
    y = np.array(y)

    # train-validation-test split + scaling
    if standarize_scales:
        scaler = StandardScaler()
        scaler.fit(X)
        X = scaler.transform(X)

    seed = np.random.randint(99999)

    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=train_size, test_size=test_size, stratify=y, random_state=seed)
    if valid_size: X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=valid_size, stratify=y_train, random_state=seed)
    else: X_val, y_val = [],[]

    if return_flavours:
        F_train, F_test, y_train, y_test = train_test_split(flavours, y, train_size=train_size, test_size=test_size, stratify=y, random_state=seed)
        if valid_size: F_train, F_val, y_train, y_val = train_test_split(F_train, y_train, test_size=valid_size, stratify=y_train, random_state=seed)
        else: F_val, y_val = [],[]

    if verbose:
        print 'no. jets selected (light/heavy) \t train: {}/{}, val: {}/{}, test: {}/{}'.format(
                                                len(y_train)-sum(y_train), sum(y_train),
                                                len(y_val)-sum(y_val), sum(y_val),
                                                len(y_test)-sum(y_test), sum(y_test))
    if return_flavours:
        return (X_train, y_train, F_train, X_test, y_test, F_test, X_val, y_val, F_val)
    else:
        return (X_train, y_train, X_test, y_test, X_val, y_val)





def unroll_df(paramtree,
            n_constit_sv = {'fConstituents':5, 'fSecondaryVertices':3},
            sorting_vars = {'fConstituents':'fpT_max', 'fSecondaryVertices':'fLxy_max'},
            relative_eta_phi=True,
            separate=True,
            verbose=True):
    """ Funtion rewriting columns filled with arrays into separate columns
    (those initial columns are 'fConstituents.*' and 'fSecondaryVertices.*').
    Only ``N`` first values from arrays + arrays' lengths are taken into account
    based on dict ```n_constit_sv``.
    If array's length is smaller than demanded number of new columns, then
    fields are filled with zeros.

    Parameters:
    -----------
    paramtree : pandas DataFrame
        data in tabular format (e.g. read by read:_msc/_inz/_jet_extractor)
    n_constit_sv : dict default={'fConstituents':5, 'fSecondaryVertices':3}
        dictionary with numbers of jet's fConstituents and fSecondaryVertices
        to be considered
    sorting_vars : dict default = {'fConstituents':'fpT_max', 'fSecondaryVertices':'fLxy_max'}
        dictionary with sorting variables for each type and information
        if increasing (min) or decreasing (max) sorting is to be performed and
        therefore N smallest or N largest values are to be selected respectively
    relative_eta_phi : bool, default=True
        if eta and phi in constituents should be relative or absolute
    separate : bool, default=True
        if should return only new columns + target column or
        previous dataframe + new columns - columns with arrays
    verbose : bool, default=True

    """

    start = time()
    # help function
    def fill_value_or_zero(row, variable, i):
        try:
            indices = row[variable.split('.')[0]+'_indices']
            return row[variable][indices[i]]
        except IndexError:
            return 0

    df = paramtree.copy()

    for constit_sv in ['fConstituents', 'fSecondaryVertices']:
        df[constit_sv+'_size'] = df.apply(lambda row: len(row[constit_sv+'.fVx']), axis=1)
        var, min_max = sorting_vars[constit_sv].split('_')
        to_be_sorted = df[constit_sv+'.'+var].tolist()
        if min_max == 'max':
            df[constit_sv+'_indices'] = [np.argsort(row)[::-1] for row in to_be_sorted]  ### reversed !!!
        elif min_max == 'min':
            df[constit_sv+'_indices'] = [np.argsort(row)       for row in to_be_sorted]

    for col in df.columns:
        for constit_sv in ['fConstituents', 'fSecondaryVertices']:
            if constit_sv+'.' not in col: continue
            if constit_sv not in n_constit_sv.keys(): n_constit_sv[constit_sv] = 0
            n = n_constit_sv[constit_sv]
            if verbose: print col,':\t'
            for i in range(n):
                if verbose: print i,
                df[constit_sv+'-{}.{}'.format(i, col.split('.')[1])] = df.apply(lambda row: fill_value_or_zero(row, col, i), axis=1)
            if verbose: print

    if relative_eta_phi:
        for i in range(n_constit_sv['fConstituents']):
            df['fConstituents-'+str(i)+'.fPhi'] = df['fConstituents-'+str(i)+'.fPhi'] - df['fPhi']
            df['fConstituents-'+str(i)+'.fEta'] = df['fConstituents-'+str(i)+'.fEta'] - df['fEta']

    # keep fConstituents-NN.XX and fConstituents_size but drop fConstituents.XX
    df = df.drop([col for col in df.columns
            if 'fConstituents.' in col
            or 'fSecondaryVertices.' in col
            or '_indices' in col],
            axis=1)

    if verbose:
        print 'unroll_df: exec. time: {} sec'.format(time()-start)

    if separate:
        sv_cols = [c for c in df.columns if 'fSecondaryVertices' in c or 'fVertex' in c]
        constit_cols = [c for c in df.columns if 'fConstituents' in c]
        if  n_constit_sv['fSecondaryVertices'] > 0 and n_constit_sv['fConstituents'] > 0:
            return df[sv_cols + constit_cols + ['target', 'flavour']]
        elif n_constit_sv['fSecondaryVertices'] > 0:
            return df[sv_cols +  ['target', 'flavour']]
        elif  n_constit_sv['fConstituents'] > 0:
            return df[constit_cols + ['target', 'flavour']]
    else:
        return df



# # # # # # # # # # # # # # # # # # # # # # # #
# KERAS callbacks
#

class StopHopeless(Callback):
    """ stops training if after ``epoch_min`` specified metric has not reached
    demanded level, either low or high, e.g.:
        StopHopeless(metric_name='acc', mertric_limit_low=0.55)
    but for loss it should be high limit:
        StopHopeless(metric_name='loss', metric_limit_high=0.3)

    """
    def __init__(self, metric_limit_low=None, metric_limit_high=None, metric_name='acc', epoch_min=10):
        self.epoch_min = epoch_min
        self.metric_limit_low = metric_limit_low
        self.metric_limit_high = metric_limit_high
        self.metric_name = metric_name

    def on_epoch_end(self, epoch, logs={}):
        if epoch >= self.epoch_min:
            if self.metric_limit_low:
                if logs.get(self.metric_name) < self.metric_limit_low:
                    print 'stopping HOPELESS training!\n\tAfter {} epochs {}={}, \
                           while low limit={}'.format(epoch, self.metric_name,
                           logs.get(self.metric_name), self.metric_limit_low)
                    self.model.stop_training = True
            if self.metric_limit_high:
                if logs.get(self.metric_name) > self.metric_limit_high:
                    print 'stopping HOPELESS training!\n\tAfter {} epochs {}={},\
                           while high limit={}'.format(epoch, self.metric_name,
                           logs.get(self.metric_name), self.metric_limit_high)
                    self.model.stop_training = True



class StopMaxETA(Callback):
    """stops training if ETA after ``epoch_min`` excesses ``max_ETA``,
    checked after each epoch
    """

    def __init__(self, max_ETA=3600, epoch_min=5):
        self.epoch_min = epoch_min
        self.max_ETA = max_ETA

    def on_epoch_end(self, epoch, logs={}):
        if epoch == 0:
            self.start_time = time()
        if epoch+1 > self.epoch_min:
            time_so_far = time() - self.start_time
            epochs_total = self.params['epochs']
            ETA = time_so_far / (epoch+1-1) * epochs_total
            if ETA > self.max_ETA:
                print 'stopping training due to MaxETA!\n\tAfter {} epochs finished in {} sec, ETA={}, while maxETA={}'.format(epoch, time_so_far, ETA, self.max_ETA)
                self.model.stop_training = True


class LogPerEpoch(Callback):
    """logs metrics per epoch instead of standard in comet.ml: per batch"""

    def __init__(self, experiment, metrics=None):
        self.experiment = experiment
        self.metrics = metrics

    def on_epoch_end(self, epoch, logs={}):
        if self.metrics:
            for m in self.metrics:
                self.experiment.log_metric(m+'_per_epoch', logs.get(m),  step=epoch+1)
        else:
            for m in logs.keys():
                self.experiment.log_metric(m+'_per_epoch', logs.get(m),  step=epoch+1)


# # # # # # # # # # # # # # # # # # # # # # # #
# KERAS models
#

def create_model_FC(activation='relu', dropout=0.0, neurons_layers='2x64', input_shape=1,
                 batch_size=64, lr=0.01, optimizer='Nadam',
                 return_descr=False, use_as_subnet=False):
    # for 1-layer nets dropout is not used !!!
    model = Sequential()
    if 'x' in neurons_layers:
        n_layers, n_hidden = [int(n) for n in neurons_layers.split('x')]
        model.add(Dense(n_hidden, activation=activation, input_shape=input_shape))
        for neurons in range(n_layers-1):
            # add dropout before, not after dense layer,
            # as there should be no dropout between two last layers
            model.add(Dropout(dropout))
            model.add(Dense(n_hidden, activation=activation))

        if use_as_subnet:
            return model
        model.add(Dense(2, activation='softmax'))
    else:
        hidden_lst = [int(h) for h in neurons_layers.split('-')]
        model.add(Dense(hidden_lst[0], activation=activation, input_shape=input_shape))
        for n_hidden in hidden_lst[1:]:
            # add dropout before, not after dense layer,
            # as there should be no dropout between two last layers
            model.add(Dropout(dropout))
            model.add(Dense(n_hidden, activation=activation))

        if use_as_subnet:
            return model
        model.add(Dense(2, activation='softmax'))

    if   optimizer == 'SGD': opt = optimizers.SGD(lr=lr)
    elif optimizer == 'Adam': opt = optimizers.Adam(lr=lr)
    elif optimizer == 'Nadam': opt = optimizers.Nadam(lr=lr)
    else: opt = optimizer  # accept also optimizer objects
    model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['accuracy'])

    model_descr = 'struct=FC:{struct}_lr={lr}_dropout={dropout}__opt={opt}_act={act}_batchsize={batch_size}'.format(
                                                                     struct=neurons_layers,
                                                                     opt=str(opt.__class__.__name__),
                                                                     lr=lr,
                                                                     act=activation,
                                                                     dropout=dropout,
                                                                     batch_size=batch_size)
    print model_descr
    if return_descr:
        return model, model_descr
    else:
        return model


def create_model_conv(n_conv_layers=2, n_filters_first=128, n_filters_change='down',
                      kernel_size=2, pool_size=2, strides=1,
                      n_fc_layers=4, n_fc_units=64,
                      activation='relu', dropout_fc=0.0, dropout_conv=0.0,
                      input_shape=1, batch_size=64, lr=0.0003, optimizer='Nadam',
                      return_descr=False, use_as_subnet=False):
    ''' n_filters_change : string 'down' or 'up' or 'flat'
            each next conv layer will have consecutively
            2x less or 2x more or same number of filters
        pool_size : int or None
            if None then no MaxPooling layers will be used,
            otherwise will be applied after each conv layer
        use_as_subnet : bool, default=False
            if True, then not compiled model without softmax layer is returned
    '''
    model = Sequential()
    # Conv Layers
    model.add( Conv1D(n_filters_first, kernel_size, strides=strides,  activation=activation, padding='same', input_shape=input_shape) )
    if pool_size:
        model.add( MaxPooling1D(pool_size=pool_size) )
    n_filters_prev = n_filters_first
    for i in range(n_conv_layers-1):

        if dropout_conv > 1e-5:
            model.add( Dropout(dropout_conv) )

        if n_filters_change == 'flat':
            n_filters = n_filters_prev
            n_filters_prev = n_filters
        elif n_filters_change == 'down':
            n_filters = int(n_filters_prev/2)
            n_filters_prev = n_filters
        elif n_filters_change == 'up':
            n_filters = int(n_filters_prev*2)
            n_filters_prev = n_filters
        model.add( Conv1D(n_filters, kernel_size, strides=strides,  activation=activation, padding='same') )
        if pool_size:
            model.add( MaxPooling1D(pool_size=pool_size) )


    model.add( Flatten() )

    # Fully-Connected Layers
    for i in range(n_fc_layers):
        if dropout_fc > 1e-5:
            model.add( Dropout(dropout_fc) )
        model.add( Dense(n_fc_units, activation=activation) )

    if use_as_subnet:
        return model


    model.add(Dense(2, activation='softmax'))

    if   optimizer == 'SGD': opt = optimizers.SGD(lr=lr)
    elif optimizer == 'Adam': opt = optimizers.Adam(lr=lr)
    elif optimizer == 'Nadam': opt = optimizers.Nadam(lr=lr)
    else: opt = optimizer  # accept also optimizer objects
    model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['accuracy'])

    index_flatten = ['Flatten' in str(l) for l in model.layers].index(True)
    strides_str = '' if strides == 1  else 's='+str(strides)
    kernel_str = '' if kernel_size == 2 else 'k='+str(kernel_size)
    kernel_strides_str = '' if not strides_str and not kernel_str else '({},{})'.format(kernel_str, strides_str)
    pool_str = '' if pool_size == 2 else '(p={})'.format(pool_size)
    structure = '{conv_maxpool}x{n_conv}[{shape_first}->{shape_last}]+Dense({n_fc_units})x{n_fc_layers}'.format(
                        conv_maxpool='(Conv1D{}+MaxPool{})'.format(kernel_strides_str, pool_str) if pool_size else 'Conv1D'+kernel_strides_str,
                        n_conv=n_conv_layers,
                        shape_first=model.layers[0].output_shape[1:],
                        shape_last=model.layers[index_flatten].input_shape[1:],
                        n_fc_units=n_fc_units,
                        n_fc_layers=n_fc_layers)
    model_descr = 'struct={struct}_lr={lr}_dropouts=({dropout_conv},{dropout_fc})'.format(
                                                                     struct=structure,
                                                                     lr=lr,
                                                                     dropout_conv=dropout_conv,
                                                                     dropout_fc=dropout_fc)
    print model_descr
    if return_descr:
        return model, model_descr
    else:
        return model



def create_model_conv_set(n_fc_layers=2, n_fc_units=64,
                 activation='relu', dropout_conv=0.0, dropout_fc=0.0, input_size=1,
                 batch_size=64, lr=0.0003, optimizer='Nadam',
                 return_descr=False, use_as_subnet=False):
    model = Sequential()
    model.add(Conv1D(128, 4, activation=activation, input_shape=input_size))
    model.add(MaxPooling1D(pool_size=4))
    model.add(Dropout(dropout_conv))
    model.add(Conv1D(64, 2, activation=activation))
    model.add(Dropout(dropout_conv))
    model.add(Conv1D(32, 2, activation=activation))
    model.add(MaxPooling1D(pool_size=2))

    model.add(Flatten())

    for _ in range(n_fc_layers):
        # add dropout before, not after dense layer,
        # as there should be no dropout between two last layers
        model.add(Dropout(dropout_fc))
        model.add(Dense(n_fc_units, activation=activation))

    if use_as_subnet:
        return model

    model.add(Dense(2, activation='softmax'))

    if   optimizer == 'SGD': opt = optimizers.SGD(lr=lr)
    elif optimizer == 'Adam': opt = optimizers.Adam(lr=lr)
    elif optimizer == 'Nadam': opt = optimizers.Nadam(lr=lr)
    else: opt = optimizer  # accept also optimizer objects
    model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['accuracy'])

    model_descr = 'struct={struct}_lr={lr}_dropouts=({dropout_conv},{dropout_fc})'.format(
                                                                     struct='~haake+FC:{}x{}'.format(n_fc_layers, n_fc_units),
                                                                     lr=lr,
                                                                     dropout_conv=dropout_conv,
                                                                     dropout_fc=dropout_fc)
    print model_descr
    if return_descr:
        return model, model_descr
    else:
        return model






# # # # # # # # # # # # # # # # # # # # # # # #
# Plotting
#

def log_roc_plots(y, y_pred_proba, datasubset, experiment=None, close=True):
    """ Function plotting ROC curve and optionally logging it to comet.ml

        Parameters:
        -----------
        y : array_like of int
            true labels
        y_pred_proba : array_like of floats
            score predicted by classifier, with sklearn: clf.predict_proba()
        dataset : string
            'train'/'valid'/'test' - used for naming plots
        experiment : comet_ml.Experiment object or None, default=None
            if None, plots will not be logged
        close : boolean, default=True
            if plt.close('all') should be executed at the end
    """
    fpr, tpr, _ = roc_curve(y, y_pred_proba)
    roc_auc = auc(fpr, tpr)
    plt.figure()
    plt.plot(fpr, tpr, color='darkorange',
             lw=2, label='ROC curve (area = {:.3f})'.format(roc_auc))
    plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
    fig_name = 'ROC cuve - ' + datasubset.upper() + ' set'
    plt.title(fig_name)
    if experiment:
        experiment.log_figure(figure_name=fig_name)

    plt.figure()
    plt.plot(tpr, fpr, color='darkorange',
             lw=2, label='ROC curve (area = {:.3f})'.format(roc_auc))
    plt.xlim([0.0, 1.0])
    plt.xlabel('b-jet efficiency')
    plt.ylabel('light mistagging efficiency')
    plt.ylim([1e-5,1.05])
    plt.semilogy()
    plt.legend(loc="lower right")
    fig_name = '(Mis)tagging efficiency - ' + datasubset.upper() + ' set'
    plt.title(fig_name)
    if experiment:
        experiment.log_figure(figure_name=fig_name)

    if close:
        plt.close('all')



def log_roc_plots_flavours(y, y_pred_proba, f, datasubset, experiment=None, close=True):
    """ Function plotting ROC curve and optionally logging it to comet.ml
    similar to log_roc_plots but plotting separate lines for each flavour

        Parameters:
        -----------
        y : array_like of int
            true labels
        y_pred_proba : array_like of floats
            score predicted by classifier, with sklearn: clf.predict_proba()
        f : array_like of int
            true flavours
        dataset : string
            'train'/'valid'/'test' - used for naming plots
        experiment : comet_ml.Experiment object or None, default=None
            if None, plots will not be logged
        close : boolean, default=True
            if plt.close('all') should be executed at the end
    """

    y_c, y_l, y_pred_proba_c, y_pred_proba_l = [],[],[],[]
    for y_i, y_pred_proba_i, f_i in zip(y, y_pred_proba, f):
        f_i = int(f_i)
        if f_i == 5 or f_i == 4:
            y_c.append(y_i)
            y_pred_proba_c.append(y_pred_proba_i)
        if f_i == 5 or f_i == 0:
            y_l.append(y_i)
            y_pred_proba_l.append(y_pred_proba_i)


    plt.figure()
    fpr_c, tpr_c, _ = roc_curve(y_c, y_pred_proba_c)
    fpr_l, tpr_l, _ = roc_curve(y_l, y_pred_proba_l)
    roc_auc_c = auc(fpr_c, tpr_c)
    roc_auc_l = auc(fpr_l, tpr_l)
    plt.plot(fpr_c, tpr_c, color='black',
             lw=2, label='b vs c ROC curve (area = {:.3f})'.format(roc_auc_c))
    plt.plot(fpr_l, tpr_l, color='blue',
             lw=2, label='b vs light ROC curve (area = {:.3f})'.format(roc_auc_l))
    plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
    fig_name = 'ROC cuve - ' + datasubset.upper() + ' set - flavours'
    plt.title(fig_name)
    if experiment:
        experiment.log_figure(figure_name=fig_name)

    plt.figure()
    plt.plot(tpr_c, fpr_c, color='black',
             lw=2, label='b vs c ROC curve (area = {:.3f})'.format(roc_auc_c))
    plt.plot(tpr_l, fpr_l, color='blue',
             lw=2, label='b vs light ROC curve (area = {:.3f})'.format(roc_auc_l))
    plt.xlim([0.0, 1.0])
    plt.xlabel('b-jet efficiency')
    plt.ylabel('c/light mistagging efficiency')
    plt.ylim([1e-5,1.05])
    plt.semilogy()
    plt.legend(loc="lower right")
    fig_name = '(Mis)tagging efficiency - ' + datasubset.upper() + ' set'
    plt.title(fig_name)
    if experiment:
        experiment.log_figure(figure_name=fig_name)

    if close:
        plt.close('all')



def log_score_plot(y_train, y_pred_proba_train,
                   y_test,  y_pred_proba_test,
                   experiment=None):
    """ Function plotting prediction score distributions,
        one with linear and second with log yscale

        Parameters:
        -----------
        y_train : array_like of int
            true labels for training set
        y_pred_proba_train : array_like of floats
            score predicted by classifier for training set, with sklearn: clf.predict_proba()
        y_test : array_like of int
            true labels for testing set
        y_pred_proba_test : array_like of floats
            score predicted by classifier for testing set, with sklearn: clf.predict_proba()
        experiment : comet_ml.Experiment object or None, default=None
            if None, plots will not be logged
    """

    for yscale in ['linear', 'log']:
        sns.distplot(np.array(y_pred_proba_train)[np.array(y_train) == 0],
                     hist=0,
                     kde_kws=dict(bw=0.05, color='r', alpha=1, linestyle=':', lw=2),
                     bins=50, norm_hist=1, label='udsg, train');
        sns.distplot(np.array(y_pred_proba_train)[np.array(y_train) == 1],
                     hist=0,
                     kde_kws=dict(bw=0.05, color='b', alpha=1, linestyle=':', lw=2),
                     label='b, train');

        sns.distplot(np.array(y_pred_proba_test)[np.array(y_test) == 0],
                     hist=0,
                     kde_kws=dict(bw=0.05, color='r', alpha=0.5, lw=2),
                     bins=50, norm_hist=1, label='udsg, test');
        sns.distplot(np.array(y_pred_proba_test)[np.array(y_test) == 1],
                     hist=0,
                     kde_kws=dict(bw=0.05, color='b', alpha=0.5, lw=2),
                     label='udsg, test');

        plt.xlim([-0.02,1.02]);
        plt.xlabel('score')
        plt.ylabel('probability density')

        if yscale == 'log':
            plt.ylim([1e-4,1e2])
        plt.yscale(yscale)
        fig_name = 'scores distr - ' + yscale
        plt.title(fig_name)
        if experiment:
            experiment.log_figure(figure_name=fig_name)
        plt.close('all')




def log_b_eff(y, y_pred_proba,
              experiment=None,
              mistag_thresholds=[1e-4, 1e-3, 1e-2, 1e-1]):
    """ Function logging b-jet efficiency levels (true positive rate) at specific
        mistagging efficiencies of light jets (false positive rate)

        Parameters:
        -----------
        y : array_like of int
            true labels
        y_pred_proba : array_like of floats
            score predicted by classifier with sklearn: clf.predict_proba()
        experiment : comet_ml.Experiment object or None, default=None
            if None, plots will not be logged
        mistag_thresholds : list of floats, default=[1e-4, 1e-3, 1e-2, 1e-1]
            list of mistagging efficiencies
            for which b-jet efficiencies should be logged
    """
    fpr, tpr, _ = roc_curve(y, y_pred_proba)
    prev_fpri = fpr[0]
    prev_tpri = tpr[0]
    for fpri, tpri in zip(fpr, tpr):
        for fpr_thresh in mistag_thresholds:
            if prev_fpri < fpr_thresh and fpri > fpr_thresh:
                print 'mistagging eff. = {} for b-jet eff. = {}'.format(fpr_thresh, tpri)
                if experiment:
                    experiment.log_parameter('bEff@mistagEff={}'.format(fpr_thresh), round(tpri,4))
        prev_fpri = fpri
        prev_tpri = tpri




def roc_uncert(y_pred_mat, y_test_mat, outname=None, make_plot=True, verbose=None):
    """ function plotting ROC curve with uncertainties for multiple trainings

    Parameters:
    -----------
    y_pred_mat : 2D array_like or array_like of ints
        y_pred for multiple (2D) or single (1D) training
    y_test_mat : 2D array_like or array_like of ints
        y_test for multiple (2D) or single (1D) training
    outname : string or None, default=None
        name of output file,
        if None, then plot is produced only in pop-up window
    make_plot : bool, default=True
        if False, then only calculations are performed
    verbose : bool or None, default=None
        if None, then it's negation of make_plot

    Returns:
    --------
    (mean_auc, std_auc) : tuple of floats
    """

    # enable 1D input
    if np.isscalar(y_pred_mat[0]): y_pred_mat = [y_pred_mat, ]
    if np.isscalar(y_test_mat[0]): y_pred_mat = [y_test_mat, ]

    if verbose == None: verbose = not make_plot

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    for i, (y_pred, y_test) in enumerate(zip(y_pred_mat, y_test_mat)):
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(y_test, y_pred)
        tprs.append(interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
        if make_plot:
            plt.plot(fpr, tpr,
                     lw=1, alpha=0.2, color='blue',
                     label='ROC fold {} (AUC = {:.3f})'.format(i, roc_auc))
        if verbose:
            print 'training no. {}: AUC = {:.4f}'.format(i, roc_auc)


    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)


    if make_plot:
        plt.plot([0, 1], [0, 1],
                 linestyle='--', lw=2,
                 color='r', alpha=.8,
                 label='Luck')
        plt.plot(mean_fpr, mean_tpr,
                 color='b', lw=2, alpha=.8,
                 label=r'Mean ROC (AUC = {:.3f} $\pm$ {:.3f})'.format(mean_auc, std_auc)
                 )
        plt.fill_between(mean_fpr, tprs_lower, tprs_upper,
                         color='grey', alpha=.2,
                         label=r'$\pm$ 1 std. dev.')
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('train_size = {}'.format(len(y_pred)))
        plt.legend(loc="lower right")

        if outname: plt.savefig(outname)
        else: plt.show()
        plt.clf()

    if verbose:
        print 'summary: AUC = {:.4f} +/- {:.4f}'.format(i, roc_auc)

    return (mean_auc, std_auc)
