#!/usr/bin/env python
import sys
import numpy as np

np.random.seed(200)
spinNum = 500 
spins = np.random.randint(-500, 500, size = (spinNum, 3))

with open(f"bath_{spinNum}.txt", "w") as f:
    line = "{:>12.6f} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(len(spins)+1, 0, 0, 0)
    f.write(line)
    
    for i in range(len(spins)):
        line = "{:>12.6f} {:>12.6f} {:>12.6f} {:>12}\n".format(spins[i][0], spins[i][1], spins[i][2], "13C")
        f.write(line)
