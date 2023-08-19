import os
import numpy as np
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from statsmodels.tools.tools import pinv_extended
import statsmodels.api as sm
import statsmodels
from statsmodels.stats.outliers_influence import variance_inflation_factor
import pandas as pd
import math
import sys
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import RandomForestRegressor
from random import sample
import tqdm
from matplotlib import pyplot

import torch
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def train(n,size1,size2,a):

    cwd = os.getcwd()

    xt = np.empty((size1,n)) 
    yt = np.empty(size1)

    x = np.empty((size2,n)) 
    y = np.empty(size2)

    '''
    xppi = np.empty((90,n))
    yppi = np.empty(90)
    '''

    # load data into arrays x and y

    f = open(cwd + "/Machine_Learning/output.txt", 'a')
    
    with open(cwd + "/Machine_Learning/data2.txt") as data:
        lines = data.readlines()
        count = 0
        for line in lines:
            line = line.split(' ')
            # print(line)
            for i in range(0,n):
                x[count][i] = float(line[i])
            y[count] = float(line[n])
            count += 1
    
    """ with open(cwd + "/Machine_Learning/data.txt") as data:
        lines = data.readlines()
        count = 0
        for line in lines:
            line = line.split(' ')
            for i in range(0,n):
                xt[count][i] = float(line[i])
            yt[count] = float(line[n])
            count += 1 """

    scaler = StandardScaler()
    # x = scaler.fit_transform(x)
    # print(x)
    # y = np.log(y*-1)
    
    rf = RandomForestRegressor()
    fi = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    for i in range(0,100):
        rf.fit(x,y)
        j = 0
        for fimp in rf.feature_importances_:
            fi[j] += fimp
            j += 1
        
        f.write("\n")
        for fimp in rf.feature_importances_:
            f.write(str(fimp) + " ")
        
    f.write("\n")
    for fimportance in fi:
        f.write(str(fimportance/100) + "\n")

    # print(pearsonr(np.reshape(x,141),np.reshape(y,141)))
    x = sm.add_constant(x)
    # xt = sm.add_constant(xt)
    # model = sm.Logit(y,x)
    model = sm.OLS(y,x)
    # results = model.fit_regularized(alpha=a,L1_wt=0)
    results = model.fit()
    

    # print(str(spearmanr(results.predict(x),y).correlation) + " " + str(spearmanr(results.predict(x),y).pvalue))

    """ for i in range(0,n+1):
        print(str(statsmodels.stats.outliers_influence.variance_inflation_factor(x, i))) """


    """ pinv_wexog,_ = pinv_extended(model.wexog)
    normalized_cov_params = np.dot(pinv_wexog, np.transpose(pinv_wexog))
    summary = sm.regression.linear_model.OLSResults(model,results.params,normalized_cov_params)
    f.write(str(summary.summary()) + "\n") """



    f.write(str(results.summary()) + "\n")

    
    
    """ sumtot = 0
    for i in range(0,10):
        kf = KFold(n_splits=4, shuffle=True)
        sum = 0
        for i, (train_index, test_index) in enumerate(kf.split(x)):
            trainX = [x[row] for row in train_index]
            trainY = [y[row] for row in train_index]
            testX = [x[row] for row in test_index]
            testY = [y[row] for row in test_index]
            m = sm.OLS(trainY,trainX)
            r = m.fit()
            # r = m.fit_regularized(alpha=a,L1_wt=0)
            sum += pearsonr(r.predict(testX),testY)[0]
        
        # print(str(lty) + " " + str(lowest) + "\n" + str(results.predict(ltx)) + "\n\n")
        sumtot += sum
        f.write(str(sum/4) + "\n")
    
    # f.write(str(sumtot/40) + "\n") 
    f.write("\n") """
    

    """ for pred in results.predict(x):
        f.write(str(pred) + "\n") """
    
    """ for pred in results.predict(xt):
        f.write(str(pred) + "\n")     """

    # f.write(str(results.aic) + "\t" + str(results.rsquared) + "\n")
    
    
    # f.write(str(results.rsquared) + "\n")

    # get best model

    """ x = sm.add_constant(x)

    model2 = sm.OLS(y,x)
    results = model2.fit()
    f.write(str(results.aic) + "\t" + str(results.rsquared) + "\n") """

# train(9,81,141,1.740)