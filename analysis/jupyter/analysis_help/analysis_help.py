from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# import seaborn as sns
# import xgboost as xgb
#
# from sklearn.grid_search import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
# from sklearn.metrics import accuracy_score, recall_score, precision_score, roc_curve, roc_auc_score, auc
# from sklearn.metrics import make_scorer
#
# from sklearn.neural_network import MLPClassifier
# from sklearn.neighbors import KNeighborsClassifier
# from sklearn.svm import SVC
# from sklearn.tree import DecisionTreeClassifier
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.naive_bayes import GaussianNB

# %matplotlib inline

import ROOT
import root_numpy
import root_pandas

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
            if not observables: observables = read_branch_names(jettree, [])

            nptree = root_numpy.tree2array(jettree, branches=observables, stop=n_rows)
            pdtree = pd.DataFrame(nptree)
            map_flavour = {'light':0, 'b':5, 'c':4}
            pdtree['target'] = np.ones_like(pdtree[observables[0]])*map_flavour[flavour]
            trees_lst.append(pdtree)

        tree_flavour = pd.concat(trees_lst)
        tree_flavour.index = range(len(tree_flavour))

        tree_flavour = tree_flavour.drop(['fMotherInitialCollision','fMotherHadronMatching'], axis=1)
        trees_flavour_lst.append(tree_flavour)

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
    verbose : bool, default=True

    Returns:
    --------
        (X_train, y_train, X_test, y_test, X_val, y_val) : numpy arrays
        self-explanatory datasets,
        X_val, y_val can be empty lists if only first split (to train/test) was performed
    """

    pdtree = paramtree.copy()

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
    np.random.shuffle(idx)
    # light-to-heavy ratio
    idx = idx[:len(df_heavy)*light_heavy_ratio]
    df = pd.concat([df_light.loc[idx], df_heavy])

    # project into X,y
    y = df['target']
    X = df.drop([c for c in df.columns if 'tag' in c], axis=1)
    X = df.drop(['target'], axis=1)
    # map light vs heavy
    y = y.map(lambda y: 0 if y == light_target else 1)
    y = np.array(y)

    # train-validation-test split + scaling
    if standarize_scales:
        scaler = StandardScaler()
        scaler.fit(X)
        X = scaler.transform(X)

    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=train_size, test_size=test_size, stratify=y)
    if valid_size: X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=valid_size, stratify=y_train)
    else: X_val, y_val = [],[]

    if verbose:
        print 'no. jets selected (light/heavy) \t train: {}/{}, val: {}/{}, test: {}/{}'.format(
                                                len(y_train)-sum(y_train), sum(y_train),
                                                len(y_val)-sum(y_val), sum(y_val),
                                                len(y_test)-sum(y_test), sum(y_test))
    return (X_train, y_train, X_test, y_test, X_val, y_val)



def unroll_df(paramtree,
            n_constit_sv = {'fConstituents':5, 'fSecondaryVertices':3},
            sorting_vars = {'fConstituents':'fpT_max', 'fSecondaryVertices':'fLxy_max'},
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
            for i in range(n):
                df[constit_sv+'-{}.{}'.format(i, col.split('.')[1])] = df.apply(lambda row: fill_value_or_zero(row, col, i), axis=1)

    # keep fConstituents-NN.XX and fConstituents_size but drop fConstituents.XX
    df = df.drop([col for col in df.columns
            if 'fConstituents.' in col
            or 'fSecondaryVertices.' in col
            or '_indices' in col],
            axis=1)

    if verbose:
        print 'unroll_df: exec. time: {} sec'.format(time()-start)

    if separate:
        return df[[col for col in df.columns if 'fConstituents' in col or 'fSecondaryVertices' in col or col == 'target']]
    else:
        return df


# # # # # # # # # # # # # # # # # # # # # # # #
# Plotting
#

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
