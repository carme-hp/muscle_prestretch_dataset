import numpy as np
import matplotlib.pyplot as plt
import json
import sys


def write_points_to_json(points, filename="coordinates.json"):
    """Writes the list of points to a JSON file with a 'seed_points' root key."""
    data = {"seed_points": points}
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python new_pattern.py <id number> <n_points> <radius> <z_min>")
        sys.exit(1)

    # Parameters
    id = str(sys.argv[1])
    n_points = int(sys.argv[2])
    radius = float(sys.argv[3])
    zmin = float(sys.argv[4])

    # Golden angle in radians
    golden_angle = np.pi * (3 - np.sqrt(5))

    # Generate evenly distributed points (Fermat spiral)
    bias = 0.5  # smaller = more points at center; try values 0.2–0.7
    r = radius * np.linspace(0, 1, n_points) ** bias
    theta = golden_angle * np.arange(n_points)

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    #points = [x, y, zmin * np.ones_like(x)]
    # points = [[float(xi), float(yi), zmin] for xi, yi in zip(x, y)]
    points = [{"x": float(xi), "y": float(yi), "z": float(zmin)} for xi, yi in zip(x, y)]
    write_points_to_json(points, "../muscle_meshes/muscle_"+id+"/fiber_seed_points_"+id + ".json")
    print("Seed points generated and saved to JSON file fiber_seed_points_"+id+".json.")

    # # Plot
    # fig, ax = plt.subplots(figsize=(6, 6))
    # ax.scatter(x, y, s=10)
    # circle = plt.Circle((0, 0), radius, color='r', fill=False, linewidth=2)
    # ax.add_artist(circle)

    # ax.set_aspect('equal')
    # ax.set_xlim(-radius-1, radius+1)
    # ax.set_ylim(-radius-1, radius+1)
    # ax.set_title("Evenly Distributed Points in a Circle (Fermat Spiral)")
    # plt.show()