import os
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.model_selection import RepeatedKFold
from sklearn.model_selection import cross_val_score
from scipy.stats import pearsonr
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import pandas as pd
import math


def train(n,size1,size2):

    cwd = os.getcwd()

    x = np.empty((size2,n)) 
    y = np.empty(size2)

    xt = np.empty((size1,n)) 
    yt = np.empty(size1)

    '''
    xppi = np.empty((90,n))
    yppi = np.empty(90)
    '''
    

    # load data into arrays x and y
    '''
    with open(cwd + "\\Machine_Learning\\prodigy_data.txt") as data:
        lines = data.readlines()
        count = 0
        for line in lines:
            line = line.split(' ')
            for i in range(0,n):
                x[count][i] = float(line[i])
            y[count] = float(line[n])
            count += 1
    '''
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
    
    with open(cwd + "/Machine_Learning/data.txt") as data:
        lines = data.readlines()
        count = 0
        for line in lines:
            line = line.split(' ')
            for i in range(0,n):
                xt[count][i] = float(line[i])
            yt[count] = float(line[n])
            count += 1

    # scaler = MinMaxScaler()
    # x = scaler.fit_transform(x)
    # print(x)

    # model = LinearRegression()
    model = RandomForestRegressor()
    
    # random forest
    model.fit(x,y)
    for val in model.feature_importances_:
        f.write(str(val) + "\n")
    f.write(str(model.score(xt,yt)))

    # repeated cross validation
    """ cv = RepeatedKFold(n_splits=4,n_repeats=10)
    scores = cross_val_score(model,x,y,cv=cv)
    i = 0
    avg = 0
    tot = 0
    for score in scores:
        avg += score
        tot += score
        if i == 3:
            print(str(avg/4))
            avg = 0
            i = 0
        else:
            i += 1
    print("tot avg: " + str(tot/40))
    print(scores) """
    

    # Other model
    
    """ # Prodigy
    # coef = [-0.1830,0.0324,-0.0786,0.2035,-0.1525,-0.0530,0.1155,-0.0610,0.0136,-0.0196,4.2761,4.5334,4.4457,-446.5052]
    # coef = [-0.1595,0.2377,-0.1674,-0.0244,-0.1740,1.1235]
    coef = [0.3882,0.1459,-0.1269,0.0915,-0.0035,0.0840,-0.1280,0.0092,-0.0250,0.0051,-0.1402,-0.0375,-0.0449,-0.0022]

    pred = []
    for row in x:
        pred.append(coef[0]*row[0] + coef[1]*row[1] + coef[2]*row[2] + coef[3]*row[3] + coef[4]*row[4] + coef[5]*row[5] + coef[6]*row[6]+ coef[7]*row[7]+ coef[8]*row[8]+ coef[9]*row[9]+ coef[10]*row[10]+ coef[11]*row[11] + coef[12]*row[12] + coef[13])
        # pred.append(coef[0]*row[0] + coef[1]*row[1] + coef[2]*row[2] + coef[3]*row[3] + coef[4]*row[4] + coef[5])
    f.write(str(pearsonr(y,pred))) """
    
    
    # Normal training
    
    """ x = sm.add_constant(x)

    model2 = sm.OLS(y,x)
    results = model2.fit()
    # f.write(str(results.summary()))
    # f.write(str(results.aic) + "\t" + str(results.rsquared) + "\n")
    f.write(str(results.rsquared) + "\n") """

train(13,62,81)