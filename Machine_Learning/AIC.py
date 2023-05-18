import os
import numpy as np
import statsmodels.api as sm


def train(n):

    cwd = os.getcwd()

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

    cwd = os.getcwd()
    f = open(cwd + "\\Machine_Learning\\output.txt", 'a')

    x = sm.add_constant(x)
    model = sm.OLS(y,x)
    results = model.fit()
    f.write(str(results.aic) + "\t" + str(results.rsquared) + "\t")