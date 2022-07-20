import os
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, RepeatedKFold
from sklearn import svm

cwd = os.getcwd()

# adjust size based on test set that is being used
x = np.empty((81,9))
y = np.empty(81)

# load data into arrays x and y
with open(cwd + "\\Machine_Learning\\prodigy_data_2.txt") as data:
    lines = data.readlines()
    count = 0
    for line in lines:
        line = line.split(' ')
        for i in range(0,9):
            x[count][i] = float(line[i])
            y[count] = float(line[9])
        count += 1


# repeated cross validation
model = RandomForestRegressor()
cv = RepeatedKFold(n_splits=4,n_repeats=10)
scores = cross_val_score(model,x,y,cv=cv)
# print(scores)
print(scores.mean())
print(scores.std())

print()

model.fit(x,y)
print(model.score(x,y))

'''
lr = LinearRegression().fit(x,y)
rfr = RandomForestRegressor().fit(x,y)
regr = svm.SVR().fit(x,y)

print(lr.score(x,y))
print(rfr.score(x,y))
print(regr.score(x,y))
'''