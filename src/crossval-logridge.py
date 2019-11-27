#### Exhaustive search for all scaling relations 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

import numpy as np
import pandas as pd
import sys 
import os
from uncertainties import ufloat, unumpy
from sklearn.linear_model import RidgeCV

np.random.seed(42) # for reproducibility

X        = 'age' if len(sys.argv) < 2 else sys.argv[1] 
combo_N  = 0     if len(sys.argv) < 3 else int(sys.argv[2]) # var select 
rem_idx  = -1    if len(sys.argv) < 4 else int(sys.argv[3]) # cv 
folds    = 10    if len(sys.argv) < 5 else int(sys.argv[4]) # cv 

nu_max_ref = 104.5
Dnu_ref    = 9.25
Teff_ref   = 4790
X_ref = {'age': 4.3, 'M':1.33, 'R':6.6}

# Directories and save files 
if rem_idx >= 0: # cross validate 
    out_dir  = "crossval_res"
    comp_dir = "crossval_comp"
    X_combo_fold = X + "_" + str(combo_N) + "_" + str(rem_idx) + ".dat"
    out_filename  = os.path.join(out_dir,  X_combo_fold)
    comp_filename = os.path.join(comp_dir, X_combo_fold)
    os.makedirs(comp_dir, exist_ok=True)
else: # run for all the data 
    out_dir = "scaling_res"
    out_filename  = os.path.join(out_dir, X + "_" + str(combo_N) + ".dat")

os.makedirs(out_dir, exist_ok=True)

DF = pd.read_table('rg.dat', sep='\t')
DF = DF[DF['ev_stage'] == 0]
#DF = DF.sort_values(by=X)
DF = DF.reset_index(drop=True)

def combos(x): 
    if len(x) == 0: 
        return [] 
    cs = combos(x[1:]) # recursion! 
    return [(x[0],)] + [(x[0],) + c for c in cs] + cs 

combinations = combos(range(4))
combinations = sorted(combinations, key=lambda x: -len(x))
combo = combinations[combo_N]

def get_stars(rem_idx=-1, folds=folds, train=True):
    stars = []
    np.random.seed(42) 
    for index, star in DF.iterrows():
        Fe_H_err = ufloat(np.random.normal(0, 0.05), 0.05)
        
        if ((index % folds == rem_idx and train) or
            (index % folds != rem_idx and not train)):
            continue 
        
        val    = ufloat(star[X],        star['e_'+X])     / X_ref[X] 
        nu_max = ufloat(star['nu_max'], star['e_nu_max']) / nu_max_ref 
        Dnu0   = ufloat(star['Dnu'],    star['e_Dnu'])    / Dnu_ref 
        Teff   = ufloat(star['Teff'],   star['e_Teff'])   / Teff_ref 
        Fe_H   = np.exp(1) ** (ufloat(star['FeH'], star['e_FeH']) + Fe_H_err)
        stars += [[nu_max, Dnu0, Teff, Fe_H, val]] 
    return np.array(stars) 

stars = get_stars(rem_idx)
stars = stars[stars[:,-1].argsort()]

train_set = unumpy.log(stars[:,combo])
train_set = np.array([[t.nominal_value for t in star]
    for star in train_set])
targets   = unumpy.log(stars[:,-1])
vals   = np.array([t.nominal_value for t in targets])
e_vals = np.array([t.std_dev for t in targets])

alphas = np.logspace(-10, 10, 10000)

cv = [[[x for i, x in enumerate(range(len(train_set))) if i%10 != k],
       [x for i, x in enumerate(range(len(train_set))) if i%10 == k]]
      for k in range(10)]

rlm = RidgeCV(alphas=alphas, fit_intercept=False, cv=cv)
rlm.fit(train_set, vals, sample_weight=1/e_vals**2)

def get_P_best(P, combo):
    P_best = P
    for ii in range(4):
        if ii not in combo: 
            P_best = list(P_best[:ii]) + [0] + list(P_best[ii:])
    return P_best 

P_best = get_P_best(rlm.coef_, combo)

# calculate R^2 
if rem_idx >= 0:
    test_set = get_stars(rem_idx, train=False) 
else:
    test_set = stars

alpha, beta, gamma, delta = P_best 
Xs      = []
scal_Xs = []
ratios  = []

for star in test_set: 
    nu_max, Dnu0, Teff, Fe_H, val = star 
    scaling_X = (nu_max**alpha * 
                   Dnu0**beta  * 
                   Teff**gamma * 
                   Fe_H**delta)
    
    val *= X_ref[X] 
    scaling_X *= X_ref[X] 
    
    Xs      += [val]
    scal_Xs += [scaling_X]
    ratios  += [scaling_X / val - 1]

Xs      = np.array(Xs)
scal_Xs = np.array(scal_Xs)

weights = [1/x.std_dev**2  for x in Xs]
SS_tot = ( (Xs - Xs.mean())**2 * weights ).sum()
SS_res = ( (Xs - scal_Xs)**2 * weights ).sum()
R2     = 1 - SS_res / SS_tot
R2adj  = 1 - (1-R2)*(len(test_set)-1) / (len(test_set)-len(combo)-1)

# calculate mean and stdev of relative differences 
ratios_ = [r.nominal_value for r in ratios]
weights = [1/r.std_dev**2  for r in ratios]
mean = np.average(ratios_, weights=weights) 
sd = np.sqrt( np.average((ratios_ - mean)**2, weights=weights) )

# calculate systematic errors 
diffs     = np.array([x-y for x,y in zip(Xs, scal_Xs)])
diffs_    = np.array([x.nominal_value for x in diffs])
weights   = np.array([1/x.std_dev**2  for x in diffs])
#w_mean    = np.average(diffs_, weights=weights)
sigma_sys = np.sqrt(np.average( (diffs_)**2 , weights=weights))

sigma_ran     = np.median([x.std_dev for x in scal_Xs])
sigma_ran_pct = np.median([x.std_dev / x.nominal_value * 100 
    for x in scal_Xs])
sigma_tot     = np.median([np.sqrt(x.std_dev**2 + sigma_sys**2) 
    for x in scal_Xs])
sigma_tot_pct = np.median([
        np.sqrt(x.std_dev**2 + sigma_sys**2) / x.nominal_value * 100
    for x in scal_Xs])

# save result 
out_str = '\t'.join(map(str, 
    [combo_N] + list(P_best) + \
    [mean, sd, R2.nominal_value, R2adj.nominal_value, 
     rlm.alpha_,
     sigma_sys,
     sigma_ran, sigma_ran_pct, 
     sigma_tot, sigma_tot_pct]
))
print(out_str)

with open(out_filename, 'w') as f:
    f.write('\t'.join(
        ['combo', 
         'alpha', 'beta', 'gamma', 'delta', 
         'mean', 'sd', 'r2', 'r2adj', 
         'lambda', 
         'sigma_sys', 
         'sigma_ran', 'sigma_ran_pct',
         'sigma_tot', 'sigma_tot_pct']) + '\n')
    f.write(out_str)
    f.write('\n')

if rem_idx >= 0: # save comparisons from cross validation 
    with open(comp_filename, 'w') as f:
        f.write('\t'.join(
            ['combo', 'rem_idx', 
             'X', 'e_X', 
             'scaling_X', 'scaling_e_X', 
             'ratio', 'e_ratio']) + '\n')
        for index in range(len(test_set)):
            f.write('\t'.join(map(str, 
               [combo_N, rem_idx, 
                     Xs[index].nominal_value,      Xs[index].std_dev,
                scal_Xs[index].nominal_value, scal_Xs[index].std_dev,
                 ratios[index].nominal_value,  ratios[index].std_dev])))
            f.write('\n')
