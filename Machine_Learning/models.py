import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from scipy.stats import pearsonr
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import pandas as pd
import math


def train(n):

    cwd = os.getcwd()

    # adjust size based on test set that is being used
    x = np.empty((81,n)) 
    y = np.empty(81)
    

    # load data into arrays x and y
    with open(cwd + "\\Machine_Learning\\prodigy_data.txt") as data:
        lines = data.readlines()
        count = 0
        for line in lines:
            line = line.split(' ')
            for i in range(0,n):
                x[count][i] = float(line[i])
            y[count] = float(line[n])
            count += 1

    scaler = MinMaxScaler()
    x = scaler.fit_transform(x)
    # print(x)

    model = RandomForestRegressor()

    '''
    # repeated cross validation
    cv = RepeatedKFold(n_splits=4,n_repeats=10)
    scores = cross_val_score(model,x,y,cv=cv)
    print(scores)
    print("R^2 mean: " + str(np.mean(scores)))
    print("R^2 standard deviation: " + str(scores.std()))
    '''
    
    arr = [[],[],[],[],[],[],[],[],[],[]]

    cwd = os.getcwd()
    f = open(cwd + "\\Machine_Learning\\output.txt", 'a')
    for i in range(0,1000):
        X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.30)
        model.fit(X_train,y_train)
        arr[0].append(model.score(X_test,y_test))
        index = 1
        for val in model.feature_importances_:
            arr[index].append(val)
            index += 1
    for i in range(0,10):
        f.write(str(sum(arr[i])/len(arr[i])) + "\t")    
    f.write("\n")
    f.flush()


    # print("Pearson correlation coefficient (r) and p-value: " + str(pearsonr(predicted, actual)))

    # print()

    

    # model.fit(x,y)
    # print(model.coef_)
    # pred = list(model.predict(x))
    # print("Pearson correlation coefficient (r) and p-value on whole dataset: " + str(pearsonr(pred, list(y))))
    # print("R^2 on entire dataset: " + str(model.score(x,y)))
    # print(model.coef_)

    # f.write(str(model.score(x,y)) + "\n")

    '''
    model.fit(x,y)
    importances = model.feature_importances_
    print(importances)
    print(model.score(x,y))
    '''
    '''
    x = sm.add_constant(x)

    model2 = sm.OLS(y,x)
    results = model2.fit()
    f.write(str(results.aic) + "\t" + str(results.rsquared) + "\t")

    ypred = results.predict(x)
    f.write(str(pearsonr(y,ypred)) + "\t")
    '''


    # print(str(pearsonr(pred, y)))
    # f.write(str(model.coef_[0]) + " " + str(model.coef_[1]) + "\n")
    

    '''
    for i in range(0,81):
        f.write(str(pred[i]) + " " + str(y[i]) + "\n")
    '''

    # print()
    # for i in range(0,81):
    #    print(str(x[i]) + " " + str(y[i]) + " " + str(pred[i]))


    # plt.plot(np.array(pred), np.array(y), 'o')
    # plt.plot(np.array([min(np.min(y),np.min(pred)),max(np.max(y),np.max(pred))]), np.array([min(np.min(y),np.min(pred)),max(np.max(y),np.max(pred))]), color='red')
    # plt.show()

# train(11)

'''
[0.04796388 0.05097372 0.24011638 0.02770147 0.30584535 0.04650581 0.09848382 0.04965073 0.13275883]

ICs_charged-charged	ICs_charged-polar	ICs_charged-apolar	ICs_polar-polar	ICS_polar-apolar	ICs_apolar-apolar	NIS_polar	NIS_apolar	NIS_charged

'''