"""
Plot benchmarking results
"""
import numpy as np
import matplotlib.pyplot as plt

# Benchmarking results
t1 = np.array([])
tN = np.array([])
N = np.array([])

# Plot
fig = plt.figure(figsize = (6.5, 4), dpi = 200)
plt.plot(np.linspace(1,8), np.linspace(1,8), 'k--',
        label = "perfect strong scaling")
plt.plot(N, t1/tN, 'bo', label = "benchmarks")
plt.xlabel('Processors')
plt.ylabel('Speedup')
plt.legend(loc = 'upper left')
plt.savefig('../img/shallow_water_benchmarks.png', 
        dpi = fig.dpi, bbox_inches = 'tight')
