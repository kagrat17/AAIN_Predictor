import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, RepeatedKFold, KFold
from sklearn import svm
from scipy.stats import pearsonr
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
import math

def run():

    cwd = os.getcwd()

    # adjust size based on test set that is being used
    x = np.empty((46,1))
    y = np.empty(46)
    

    # load data into arrays x and y
    with open(cwd + "\\Machine_Learning\\test2.txt") as data:
        lines = data.readlines()
        count = 0
        for line in lines:
            line = line.split(' ')
            for i in range(0,1):
                x[count][i] = int(line[i])
                y[count] = float(line[1])
            count += 1

    scaler = MinMaxScaler()
    # scaler = StandardScaler()
    x = scaler.fit_transform(x)
    # print(x)

    model = LinearRegression()

    '''
    # repeated cross validation
    cv = RepeatedKFold(n_splits=4,n_repeats=10)
    scores = cross_val_score(model,x,y,cv=cv)
    print(scores)
    print("R^2 mean: " + str(np.mean(scores)))
    print("R^2 standard deviation: " + str(scores.std()))


    rkf = RepeatedKFold(n_splits=4, n_repeats=10)
    for train_index, test_index in rkf.split(x):
        x_train, x_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]
        model.fit(x_train,y_train)
        actual += list(y_test)
        predicted += list(model.predict(x_test))


    print("Pearson correlation coefficient (r) and p-value: " + str(pearsonr(predicted, actual)))

    print()'''

    model.fit(x,y)
    pred = list(model.predict(x))
    # print("Pearson correlation coefficient (r) and p-value on whole dataset: " + str(pearsonr(pred, list(y))))
    # print("R^2 on entire dataset: " + str(model.score(x,y)))
    # print(model.coef_)

    cwd = os.getcwd()
    f = open(cwd + "\\Machine_Learning\\output.txt", 'a')

    f.write(str(model.score(x,y)) + "\t")
    f.write(str(model.coef_[0]) + "\t")
    # f.write(str(np.mean(x,axis=0)[0]) + "\t")
    f.write("\n")

    # print()
    # for i in range(0,81):
    #    print(str(x[i]) + " " + str(y[i]) + " " + str(pred[i]))


    # plt.plot(np.array(pred), np.array(y), 'o')
    # plt.plot(np.array([min(np.min(y),np.min(pred)),max(np.max(y),np.max(pred))]), np.array([min(np.min(y),np.min(pred)),max(np.max(y),np.max(pred))]), color='red')
    # plt.show()



    '''
    lr = LinearRegression().fit(x,y)
    rfr = RandomForestRegressor().fit(x,y)
    regr = svm.SVR().fit(x,y)

    print(lr.score(x,y))
    print(rfr.score(x,y))
    print(regr.score(x,y))
    '''

run()