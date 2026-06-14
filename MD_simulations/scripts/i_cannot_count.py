import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import sys

file = sys.argv[1]

df = pd.read_csv(file)

# extract base pair (remove trailing _H1 / _H2)
df["pair_base"] = df["pair"].str.rsplit("_", n=1).str[0]

# sum occupancy_percent per base pair
try:
    result = df.groupby("pair_base", as_index=False)["occupancy_percent"].sum()
except:
    result = df.groupby("pair_base", as_index=False)["p"].sum()

print(result)