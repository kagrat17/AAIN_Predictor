'''generate histogram graphs from data using matplotlib'''

import matplotlib.pyplot as plt

x = {18: 172, 8: 85, 13: 90, 7: 67, 9: 131, 6: 45, 20: 92, 11: 95, 4: 2, 19: 82, 16: 137, 14: 93, 17: 143, 10: 149, 15: 101, 5: 38, 12: 79}

dict = {}

for i in sorted(x):
    dict[i] = x[i]

keys = list(dict.keys())
keyStrings = [str(key) for key in keys]
vals = list(dict.values())

plt.bar(keyStrings, vals)
plt.title('TRP')
plt.show()