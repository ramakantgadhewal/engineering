import numpy as np

# Sturgis rule for the optimal number of bins in a histogram. Taken from 
# Bechwith, Mechanical Measurements, 5th edition page 65 footnote *
# N = 1 + 3.3 log n, n is the total number of points

n = 31

N = round(1 + 3.3 * np.log10(n)) # 5.92 or 6

# https://www.statology.org/sturges-rule/ equation is log2(n) + 1

N = round(np.log2(n) + 1) # 5.95 or 6