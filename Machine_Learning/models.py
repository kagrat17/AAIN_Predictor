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

def mmd(x, y, sigma):
    # compare kernel MMD paper and code:
    # A. Gretton et al.: A kernel two-sample test, JMLR 13 (2012)
    # http://www.gatsby.ucl.ac.uk/~gretton/mmd/mmd.htm
    # x shape [n, d] y shape [m, d]
    # n_perm number of bootstrap permutations to get p-value, pass none to not get p-value
    n, d = x.shape
    m, d2 = y.shape
    assert d == d2
    xy = torch.cat([x.detach(), y.detach()], dim=0)
    dists = torch.cdist(xy, xy, p=2.0)
    # we are a bit sloppy here as we just keep the diagonal and everything twice
    # note that sigma should be squared in the RBF to match the Gretton et al heuristic
    k = torch.exp((-1/(2*sigma**2)) * dists**2) + torch.eye(n+m)*1e-5
    k_x = k[:n, :n]
    k_y = k[n:, n:]
    k_xy = k[:n, n:]
    # The diagonals are always 1 (up to numerical error, this is (3) in Gretton et al.)
    # note that their code uses the biased (and differently scaled mmd)
    mmd = k_x.sum() / (n * (n - 1)) + k_y.sum() / (m * (m - 1)) - 2 * k_xy.sum() / (n * m)
    return mmd

def MMD(x, y, kernel):
    """Emprical maximum mean discrepancy. The lower the result
       the more evidence that distributions are the same.

    Args:
        x: first sample, distribution P
        y: second sample, distribution Q
        kernel: kernel type such as "multiscale" or "rbf"
    """
    xx, yy, zz = torch.mm(x, x.t()), torch.mm(y, y.t()), torch.mm(x, y.t())
    rx = (xx.diag().unsqueeze(0).expand_as(xx))
    ry = (yy.diag().unsqueeze(0).expand_as(yy))

    dxx = rx.t() + rx - 2. * xx # Used for A in (1)
    dyy = ry.t() + ry - 2. * yy # Used for B in (1)
    dxy = rx.t() + ry - 2. * zz # Used for C in (1)

    XX, YY, XY = (torch.zeros(xx.shape).to(device),
                  torch.zeros(xx.shape).to(device),
                  torch.zeros(xx.shape).to(device))

    if kernel == "multiscale":

        bandwidth_range = [0.2, 0.5, 0.9, 1.3]
        for a in bandwidth_range:
            XX += a**2 * (a**2 + dxx)**-1
            YY += a**2 * (a**2 + dyy)**-1
            XY += a**2 * (a**2 + dxy)**-1

    if kernel == "rbf":

        bandwidth_range = [10, 15, 20, 50]
        for a in bandwidth_range:
            XX += torch.exp(-0.5*dxx/a)
            YY += torch.exp(-0.5*dyy/a)
            XY += torch.exp(-0.5*dxy/a)

    return torch.mean(XX + YY - 2. * XY)

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
    
    """ pca = PCA(n_components=13)
    pca.fit(x)

    PC = range(1, pca.n_components_+1)
    plt.bar(PC, pca.explained_variance_ratio_, color='gold')
    plt.xlabel('Principal Components')
    plt.ylabel('Variance %')
    plt.xticks(PC)
    plt.show()

    PCnames = ['PC'+str(i+1) for i in range(pca.n_components_)]
    Loadings = pd.DataFrame(pca.components_,columns=PCnames)

    print(Loadings.iloc[:,:3])
    print(pca.singular_values_)
    print(pca.n_features_in_)
    print(pca.n_samples_)

    pcr = make_pipeline(StandardScaler(), PCA(n_components=3), LinearRegression())
    pcr.fit(x, y)

    print(pcr.score(x,y)) """

    
    rf = RandomForestRegressor()
    fi = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    for i in range(0,100):
        rf.fit(x,y)
        j = 0
        for fimp in rf.feature_importances_:
            fi[j] += fimp
            j += 1
        '''
        f.write("\n")
        for fimp in rf.feature_importances_:
            f.write(str(fimp) + " ")
        '''
    for fimportance in fi:
        f.write(str(fimportance/100) + "\n")
    '''
    x = x[np.random.choice(x.shape[0],50,replace=False),:]
    xt = xt[np.random.choice(xt.shape[0],50,replace=False),:]

    xtr = torch.from_numpy(x)
    ytr = torch.from_numpy(xt)

    our_mmd = MMD(xtr,ytr,kernel="multiscale")
    f.write(str(our_mmd) + "\n")


    # bootstrapping
    xy = np.concatenate((x,xt),axis=0)

    mmds = []
    for i in tqdm.tqdm(range(1000)):
        xy1, xy2 = np.split(np.random.permutation(xy),2)
        mmds.append(MMD(torch.from_numpy(xy1), torch.from_numpy(xy2), kernel="multiscale").item())
    mmds = torch.tensor(mmds)

    pyplot.hist(mmds.numpy(), bins=20)

    print((our_mmd < mmds).float().mean())

    print(torch.quantile(mmds,0.95))

    '''

    # MMD
    """ x = x[np.random.choice(x.shape[0],50,replace=False),:]
    xt = xt[np.random.choice(xt.shape[0],50,replace=False),:]
    xy = np.concatenate((x,xt),axis=0)

    xtr = torch.from_numpy(x)
    ytr = torch.from_numpy(xt)
    xytr = torch.from_numpy(xy)

    sigma = xytr.median()/2
    our_mmd = mmd(xtr, ytr, sigma)
    print(our_mmd)

    N_X = len(xtr)
    N_Y = len(ytr)

    mmds = []
    for i in tqdm.tqdm(range(1000)):
        xy1, xy2 = np.split(np.random.permutation(xy),2)
        mmds.append(mmd(torch.from_numpy(xy1), torch.from_numpy(xy2), sigma).item())
    mmds = torch.tensor(mmds)
    print((our_mmd < mmds).float().mean()) """

    # f.write(str(pearsonr(x,y)[0]))
    
    
    x = sm.add_constant(x)
    xt = sm.add_constant(xt)
    model = sm.OLS(y,x)
    # results = model.fit_regularized(alpha=a,L1_wt=0)
    results = model.fit()
    

    # print(str(spearmanr(results.predict(x),y).correlation) + " " + str(spearmanr(results.predict(x),y).pvalue))

    """ for i in range(0,13):
        print(str(statsmodels.stats.outliers_influence.variance_inflation_factor(x, i)))

    pinv_wexog,_ = pinv_extended(model.wexog)
    normalized_cov_params = np.dot(pinv_wexog, np.transpose(pinv_wexog))
    summary = sm.regression.linear_model.OLSResults(model,results.params,normalized_cov_params)
    f.write(str(summary.summary()) + "\n") """



    """ f.write(str(results.summary()) + "\n")

    sumtot = 0
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
    f.write("\n")

    
   
    

    f.write("\n")
    

    for pred in results.predict(x):
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

train(14,81,141,1.740)