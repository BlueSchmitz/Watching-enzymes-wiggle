import numpy as np
import matplotlib.pyplot as plt

# Reaction coordinate
x = np.linspace(0, 10, 1000)

# Schematic free-energy landscape
y = (
    -5.5 * np.exp(-(x - 2.0)**2 / 0.4)   # first (deeper) minimum
    -3.0 * np.exp(-(x - 5.0)**2 / 0.4)   # second minimum
    -1.5 * np.exp(-(x - 8.0)**2 / 0.5)   # third shallower minimum
    +1.3 * np.exp(-(x - 3.5)**2 / 0.3)   # barrier between first and second
)

# Shift upward so energies are positive
y = y - y.min() + 0.2

plt.figure(figsize=(6, 4))
plt.plot(x, y, lw=2)

plt.xlabel("Reaction coordinate", fontsize=12)
plt.ylabel("Free energy", fontsize=12)

# Visible axes
plt.xlim(0, 10)
plt.ylim(-3, y.max() + 0.5)

# Remove numerical tick labels but keep axis lines
plt.xticks([])
plt.yticks([])

# Keep left and bottom spines, remove top/right
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig("PMF.svg")
plt.close()