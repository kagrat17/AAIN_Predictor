import os
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr

# train a linear regression model with n features
def train(n):

    cwd = os.getcwd()

    x = np.empty((81,n))  # values of the features
    y = np.empty(81) # experimental deltaG values
    
    # load data into arrays x and y
    with open(cwd + "prodigy_data.txt") as data:
        lines = data.readlines()
        count = 0
        for line in lines:
            line = line.split(' ')
            for i in range(0,n):
                x[count][i] = float(line[i])
            y[count] = float(line[n])
            count += 1

    model = LinearRegression()
    model.fit(x,y)

    # coefficients of each feature
    print(str(model.coef_))
    # constant term in the regression
    print(str(model.intercept_))

    pred = list(model.predict(x))
    print("Pearson correlation coefficient (r) and p-value: " + str(pearsonr(pred, list(y))))

train(6)