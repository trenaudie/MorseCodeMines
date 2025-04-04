#%%
import numpy as np 
import matplotlib.pyplot as plt
with open("morse_output.txt", "r") as f:
    lines = f.readlines()
lines = [float(line.strip()) for line in lines]
lines = np.array(lines)
plt.plot(lines)
# %%
