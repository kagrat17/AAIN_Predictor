import os
import numpy as np
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

    x = np.empty((19,n)) 
    y = np.empty(19)

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

    
    with open(cwd + "\\Machine_Learning\\prodigy_data.txt") as data:
        lines = data.readlines()
        count = 0
        for line in lines:
            line = line.split(' ')
            for i in range(0,n):
                x[count][i] = float(line[i])
            y[count] = float(line[n])
            count += 1
    

    # scaler = MinMaxScaler()
    # x = scaler.fit_transform(x)
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
    
    arr = [[],[],[],[],[],[],[],[],[],[],[],[],[]]

    cwd = os.getcwd()
    f = open(cwd + "\\Machine_Learning\\output.txt", 'a')

    # model.fit(x,y)
    # for val in model.feature_importances_:
    #   f.write(str(val) + "\t")
    
    
    for i in range(0,100):
        X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.30)
        model.fit(X_train,y_train)
        arr[0].append(model.score(X_test,y_test))
        arr[1].append(model.score(X_train,y_train))
        index = 2
        
        for val in model.feature_importances_:
            arr[index].append(val)
            index += 1
        
        predictions = model.predict(X_test)
        arr[index].append(sm.tools.eval_measures.rmse(predictions, y_test, axis=0))

        for j in range(len(arr)):
            f.write(str(arr[j][len(arr[j])-1]) + "\t")
        f.write("\n")    

        if i % 20 == 0:
            print(i)
    
    f.write("\n")
    f.flush()
    

    
    '''
    x1 = sm.add_constant(x)
    model = sm.OLS(y,x1)
    results = model.fit()
    f.write(str(results.summary()) + "\n\n\n")
    '''

    '''
    predictions = results.predict(x1)
    f.write("RMSE: " + str(sm.tools.eval_measures.rmse(predictions, y, axis=0)) + "\n\n\n")

    for pred in predictions:
        f.write(str(pred) + "\n")
    '''
    


    '''
    model.fit(x,y)
    for val in model.feature_importances_:
        f.write(str(val) + "\t")
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

train(10)