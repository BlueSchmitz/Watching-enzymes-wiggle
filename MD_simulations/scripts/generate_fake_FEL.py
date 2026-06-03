import numpy as np
import matplotlib.pyplot as plt

# Grid
x = np.linspace(-5, 5, 300)
y = np.linspace(-5, 5, 300)
X, Y = np.meshgrid(x, y)

# Gentle funnel
Z = 0.05 * (X**2 + Y**2)

# Deep global minimum
Z += -6.0 * np.exp(
    -((X + 2.0)**2 + (Y + 1.5)**2) / 1.2
)

# Shallow metastable basins
Z += -1.5 * np.exp(
    -((X - 2.0)**2 + (Y - 2.5)**2) / 0.8
)

Z += -1.2 * np.exp(
    -((X - 1.0)**2 + (Y + 2.0)**2) / 0.7
)

Z += -1.0 * np.exp(
    -((X + 3.0)**2 + (Y - 2.0)**2) / 0.6
)

Z += -0.8 * np.exp(
    -((X - 3.0)**2 + (Y + 3.0)**2) / 0.5
)

# Barriers separating basins
Z += 1.5 * np.exp(
    -((X + 0.3)**2 + (Y + 0.3)**2) / 1.5
)

Z += 0.8 * np.exp(
    -((X - 1.5)**2 + (Y - 0.5)**2) / 1.0
)

# Plot
plt.figure(figsize=(6,5))
plt.contourf(X, Y, Z, levels=30, cmap='viridis')
plt.colorbar(label='Free energy')
plt.xlabel('Reaction coordinate 1')
plt.ylabel('Reaction coordinate 2')
plt.tight_layout()
plt.show()