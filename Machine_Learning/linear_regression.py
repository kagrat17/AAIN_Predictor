import os
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from sklearn.linear_model import Lasso
from sklearn.linear_model import ElasticNet

cwd = os.getcwd()

x = np.empty((81,7))
y = np.empty(81)

with open(cwd + "\\Machine_Learning\\prodigy_data.txt") as data:
    lines = data.readlines()
    count = 0
    for line in lines:
        line = line.split(' ')
        x[count] = [float(line[0]), float(line[1]), float(line[2]), float(line[3]), float(line[4]), float(line[5]), float(line[6])]
        y[count] = float(line[7])
        count += 1

model = ElasticNet(alpha=0.01).fit(x,y)

r_sq = model.score(x,y)
print(r_sq)
print(model.coef_)