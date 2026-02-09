import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

if len(sys.argv) != 5:
    print("Usage: python vfield_ellipsoid.py <id number> <muscle_length> <Rmax> <Rmin>")
    sys.exit(1)

id = str(sys.argv[1])
L = float(sys.argv[2])
Rmax = float(sys.argv[3])
Rmin = float(sys.argv[4])

zmin, zmax = -L/2, L/2

# Tube radius and its derivative
def R(z):
    return Rmin + (Rmax - Rmin) * (1 + np.cos(np.pi * z / zmax)) / 2.0

def dRdz(z):
    return -(Rmax - Rmin) * (np.pi / 2 / zmax) * np.sin(np.pi * z / zmax)

# Vector field (tangent to tube)
def F(x, y, z):
    Rz = R(z)
    k = dRdz(z) / Rz
    return np.array([k * x, k * y, 1.0])

# ------------------------------
# Structured 3D grid
# ------------------------------
x_vals = np.linspace(-1.1*Rmax, 1.1*Rmax, 50)
y_vals = np.linspace(-1.1*Rmax, 1.1*Rmax, 50)
z_vals = np.linspace(zmin, zmax, 50)
X, Y, Z = np.meshgrid(x_vals, y_vals, z_vals, indexing='ij')

# Compute vector components on grid
U = np.zeros_like(X)
V = np.zeros_like(Y)
W = np.zeros_like(Z)

for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        for k in range(X.shape[2]):
            fx, fy, fz = F(X[i,j,k], Y[i,j,k], Z[i,j,k])
            U[i,j,k], V[i,j,k], W[i,j,k] = fx, fy, fz

# Optionally normalize vectors (for consistent arrow size)
norm = np.sqrt(U**2 + V**2 + W**2)
U /= norm
V /= norm
W /= norm

# ------------------------------
# Save to VTK file
# ------------------------------
import pyvista as pv
# Create a VTK structured grid
grid = pv.StructuredGrid()
grid.points = np.c_[X.ravel(), Y.ravel(), Z.ravel()]
grid.dimensions = X.shape
grid["vectors"] = np.c_[U.ravel(), V.ravel(), W.ravel()]

# Save the structured grid to a VTK file
grid.save("../muscle_meshes/muscle_"+id+"/vector_field_"+id+".vtk")
print("Vector field generated and saved to VTK file vector_field_"+id+".vtk.")

# ------------------------------
# Plot quiver field
# ------------------------------
# fig = plt.figure(figsize=(11, 8))
# ax = fig.add_subplot(111, projection='3d')

# ax.quiver(X, Y, Z, U, V, W, length=0.8, normalize=False, color='crimson', alpha=0.7)

# # Axes and appearance
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# ax.set_xlim([-5, 5])
# ax.set_ylim([-5, 5])
# ax.set_zlim([zmin, zmax])
# ax.set_title('Vector Field F(x, y, z) on Structured Grid')
# ax.view_init(elev=25, azim=45)

# plt.tight_layout()
# plt.show()