from contactCalc import *;

cutoff = float(input("Enter the cutoff distance to calculate for: "))

files = ["6M0J", "7EKC", "7EKF", "7EKG", "7WBP", "7WBQ"]
for file in files:
    calculate(file, cutoff)