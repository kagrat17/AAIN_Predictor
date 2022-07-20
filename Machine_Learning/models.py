import os
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, RepeatedKFold
from sklearn import svm

cwd = os.getcwd()

# adjust size based on test set that is being used
x = np.empty((81,7))
y = np.empty(81)

# load data into arrays x and y
with open(cwd + "\\Machine_Learning\\prodigy_data.txt") as data:
    lines = data.readlines()
    count = 0
    for line in lines:
        line = line.split(' ')
        x[count] = [float(line[0]), float(line[1]), float(line[2]), float(line[3]), float(line[4]), float(line[5]), float(line[6])]
        y[count] = float(line[7])
        count += 1


# repeated cross validation
model = svm.SVR()
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