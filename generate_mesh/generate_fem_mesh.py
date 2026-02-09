import numpy as np
import sys


def write_structured_vtk(filename, X, Y, Z):
    nx, ny, nz = X.shape
    with open(filename, "w") as f:
        f.write("# vtk DataFile Version 3.0\nStructured hexahedral tube mesh\nASCII\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write(f"DIMENSIONS {nx} {ny} {nz}\n")
        n_points = nx * ny * nz
        f.write(f"POINTS {n_points} float\n")
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    f.write(f"{X[i,j,k]} {Y[i,j,k]} {Z[i,j,k]}\n")


# def read_structured_vtk(filename):
#     with open(filename, "r") as f:
#         lines = f.readlines()

#     # Find DIMENSIONS line
#     dim_line = next(line for line in lines if line.startswith("DIMENSIONS"))
#     _, nx, ny, nz = dim_line.split()
#     nx, ny, nz = int(nx), int(ny), int(nz)

#     # Find POINTS section
#     start_idx = next(i for i, line in enumerate(lines) if line.startswith("POINTS")) + 1

#     # Read all points
#     coords = []
#     for line in lines[start_idx:]:
#         parts = line.strip().split()
#         if len(parts) == 3:
#             coords.append([float(parts[0]), float(parts[1]), float(parts[2])])

#     coords = np.array(coords)
    
#     # Reshape to (nz, ny, nx, 3)
#     coords = coords.reshape((nz, ny, nx, 3))
#     coords = np.transpose(coords, (2, 1, 0, 3))  # -> (nx, ny, nz, 3)

#     # Flip x-axis to go from small → large
#     coords = np.flip(coords, axis=0)

#     # Separate into X, Y, Z
#     X = coords[..., 0]
#     Y = coords[..., 1]
#     Z = coords[..., 2]

#     # Also return flattened vector list
#     vector_list = coords.reshape(-1, 3).tolist()

#     return X, Y, Z, vector_list



if len(sys.argv) != 8:
    print("Usage: python vfield_ellipsoid.py <id number> <muscle_length> <Rmax> <Rmin> <nx> <ny> <nz>")
    sys.exit(1)

id = str(sys.argv[1])
L = float(sys.argv[2])
Rmax = float(sys.argv[3])
Rmin = float(sys.argv[4])
nx = int(sys.argv[5])
ny = int(sys.argv[6])
nz = int(sys.argv[7])


zmin, zmax = -L/2, L/2

# Tube radius function
def R(z):
    """Tube radius (widest at z=0, narrow at z=±10)."""
    return Rmin + (Rmax - Rmin) * (1 + np.cos(np.pi * z / zmax)) / 2.0

# ------------------------------
# Structured grid (cube)
# ------------------------------
x_lin = np.linspace(-1, 1, nx)
y_lin = np.linspace(-1, 1, ny)
z_lin = np.linspace(zmin, zmax, nz)

X, Y, Z = np.meshgrid(x_lin, y_lin, z_lin, indexing="ij")

# ------------------------------
# Map cube → circular tube
# ------------------------------
Rz = R(Z)
r_xy = np.sqrt(X**2 + Y**2)
r_xy = np.maximum(r_xy, 1e-8)  # avoid division by zero

# Normalize (x,y) to unit circle, scale by R(z)
X_mapped = X / np.max(np.abs([X, Y])) * Rz * np.sign(X)
Y_mapped = Y / np.max(np.abs([X, Y])) * Rz * np.sign(Y)

# A smoother, radially symmetric version:
# (instead of using corners, use circular scaling)
Rmax_xy = np.sqrt(1.0)  # cube boundary radius
X_mapped = X / r_xy * Rz * np.minimum(r_xy, Rmax_xy)
Y_mapped = Y / r_xy * Rz * np.minimum(r_xy, Rmax_xy)

# ------------------------------
# Flatten node coordinates
# ------------------------------
points = np.vstack([X_mapped.ravel(), Y_mapped.ravel(), Z.ravel()]).T

# ------------------------------
# Generate hexahedral connectivity
# ------------------------------
hexes = []
for i in range(nx - 1):
    for j in range(ny - 1):
        for k in range(nz - 1):
            n0 = i * ny * nz + j * nz + k
            n1 = (i + 1) * ny * nz + j * nz + k
            n2 = (i + 1) * ny * nz + (j + 1) * nz + k
            n3 = i * ny * nz + (j + 1) * nz + k
            n4 = i * ny * nz + j * nz + (k + 1)
            n5 = (i + 1) * ny * nz + j * nz + (k + 1)
            n6 = (i + 1) * ny * nz + (j + 1) * nz + (k + 1)
            n7 = i * ny * nz + (j + 1) * nz + (k + 1)
            hexes.append([n0, n1, n2, n3, n4, n5, n6, n7])

hexes = np.array(hexes, dtype=int)

write_structured_vtk("../muscle_meshes/muscle_"+id+"/3D_mesh_"+id+".vtk", X_mapped, Y_mapped, Z)
print("3D FEM mesh generated and saved to 3D_mesh_"+id+".vtk")
print(f"Nodes: {len(points)}, Elements: {len(hexes)}")


# FEM mesh generation from .vts file
# -------------------------------------------------------------------

# X,Y,Z, points =read_structured_vtk("structured_muscle"+id+".vtk")

# write_structured_vtk("rewrite_structured_muscle"+id+".vtk", X, Y, Z)
# print(points)