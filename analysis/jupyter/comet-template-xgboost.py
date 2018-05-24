from __future__ import division
from comet_ml import Experiment

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import xgboost as xgb
import scipy.stats as scistats

from sklearn.grid_search import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

import ROOT
import root_numpy

from time import time
import analysis_help as ahelp


pd.set_option('precision', 4)


def print_framed(msg):
    print '\n'+'='*20 + '\n{}\n'.format(msg) + '='*20+'\n'


np.random.seed(0)


start_0 = time()
df_orig = ahelp.read_jet_extractor(fnames='../../../benchmark-dataset/DATA/trains/train_LHC17f8g_11.root', n_rows=500)
df = df_orig.copy()
df_full = ahelp.unroll_df(df, n_constit_sv={'fConstituents':1, 'fSecondaryVertices':2}, separate=False)
df = df_full.copy()
# df = df[[c for c in df.columns if ('fSecondaryVertices-' in c and 'PID' not in c) or c == 'target']]
df = df[[c for c in df.columns if 'PID' not in c]]
df.head()
X_train, y_train, X_test, y_test, X_val, y_val = ahelp.prepare_dataset(df, light_heavy_ratio=1, train_size=0.8, standarize_scales=False)

experiment = Experiment(api_key="Vr3JAwshw0tYqpsQluIhQxdaL",
                        project_name='xgb-template-tests')

# PROVIDE DESCRIPTION OF EXPERIMENT
experiment.log_other('descr', 'if cols_used as a string is saved')
print 'cols used = ', df.columns.tolist()
experiment.log_parameter('N_col_total_sv_constit', [len(df.columns)-1,
                                                    len([c for c in df.columns if 'fSecondaryVertices' in c]),
                                                    len([c for c in df.columns if 'fConstituents' in c])])
experiment.log_parameter('input_shape', X_train.shape)
experiment.log_dataset_hash(X_train)

time_data_preparing = time()-start_0
print_framed('total data preparing time: {} sec'.format(time_data_preparing))
experiment.log_other('time_data_preparing [sec]', round(time_data_preparing,1))


param_dist =  {
                'max_depth':range(1,10)*3 + range(10,101,10), # enhance small numbers
                'n_estimators':[1,2,5,10,15,20,40,60,80,100],
                'learning_rate':[round(x,6) for x in np.logspace(-5,1,13)],
                'colsample_bytree':[round(x,2) for x in  np.arange(0.1, 1.01, 0.1)],
                'colsample_bylevel':[round(x,2) for x in np.arange(0.1, 1.01, 0.1)],
                'subsample':[round(x,2) for x in np.arange(0.1, 1.01, 0.1)],
                'gamma':[round(x,6) for x in np.logspace(-5,3,17)]+[0]*3
                }


start = time()
clf_rand = RandomizedSearchCV(xgb.XGBClassifier(), param_distributions=param_dist, n_iter=10,
                              scoring='roc_auc', cv=3, verbose=1, return_train_score=True)
clf_rand.fit(X_train, y_train)
experiment.log_parameter('n_iter', clf_rand.n_iter)
time_training = time()-start
print_framed('RandomSearch exec time: {} sec'.format(time_training))
experiment.log_other('time_model_training [sec]', round(time_training,1))


experiment.log_multiple_params(param_dist, prefix='grid')
experiment.log_multiple_params(clf_rand.best_params_, prefix='best')



# all tested paramsets
# print_framed('all tested paramsets:')
# gscores = clf_rand.grid_scores_
# for gs in sorted(gscores, key=lambda(x): x.mean_validation_score, reverse=True): print gs

print_framed('all tested paramsets:')
cv_res = clf_rand.cv_results_
print '    {}{}{}\t\t{}'.format('valid:', ' '*14,'train:', 'params:')
for te_m, te_s, tr_m, tr_s, par in sorted(zip(cv_res['mean_test_score'],
                                              cv_res['std_test_score'],
                                              cv_res['mean_train_score'],
                                              cv_res['std_train_score'],
                                              cv_res['params']),
                                          key=lambda x: -x[0]):
    print '{:.4f}+/-{:.4f}    {:.4f}+/-{:.4f}\t\t{}'.format(te_m, te_s, tr_m, tr_s, par)


# feature importance
print_framed('feature importances:')
col_names = [c for c in df.columns if c != 'target']
results = sorted(zip(clf_rand.best_estimator_.feature_importances_, col_names), key=lambda x: -x[0])
dict_vars = {}
dict_nums = {}
for imp, var in results:
    if '.' not in var: continue
    num,var_type = var.split('.')
    if num in dict_nums.keys(): dict_nums[num].append(imp)
    else: dict_nums[num] = [imp,]

    if var_type in dict_vars.keys(): dict_vars[var_type].append(imp)
    else: dict_vars[var_type] = [imp,]

res_lst = []
for k,v in dict_vars.items():
    res_lst.append([k, np.mean(v), np.std(v)])
for r in sorted(res_lst, key=lambda row: -row[1]):
    print '{}: {:.3f}+/-{:.3f}'.format(*r)

res_lst = []
for k,v in dict_nums.items():
    res_lst.append([k, np.mean(v), np.std(v)])
for r in sorted(res_lst, key=lambda row: -row[1]):
    print '{}: {:.3f}+/-{:.3f}'.format(*r)
print '\n'+' ='*15+'\n'
for r in results: print r


# metrics
print_framed('metrics:')
from sklearn.metrics import accuracy_score, auc, precision_score, recall_score, roc_auc_score

for mode, X, y in zip( [experiment.train, experiment.test], [X_train, X_test], [y_train, y_test] ):
    with mode():
        print mode.__name__
        y_pred = clf_rand.best_estimator_.predict(X)
        y_pred_proba = [i[1] for i in clf_rand.best_estimator_.predict_proba(X)]

        metrics = {}
        for metric in accuracy_score, precision_score, recall_score:
            print '{}: \t{}'.format( metric.__name__, metric(y, y_pred))
            metrics[metric.__name__.split('_')[0]] = metric(y, y_pred)
        print 'AUC: \t\t{}'.format(roc_auc_score(y, y_pred_proba))
        metrics['roc_auc'] = roc_auc_score(y, y_pred_proba)

        experiment.log_multiple_metrics(metrics)
