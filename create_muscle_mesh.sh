#!/bin/bash

# Ask for user input
read -p "Enter dataset id: " id
read -p "Enter number of fibers: " n_points
read -p "Enter number of fem points: " nx ny nz
read -p "Enter length of the muscle: " L
read -p "Enter radius at the center: " Rmax
read -p "Enter radius at the extreme: " Rmin

zmin=$(echo "- $L /2" | bc -l)
zmax=$(echo "$L /2" | bc -l)
seed_radius=$(echo "0.9*$Rmin" | bc -l)

muscle_dir="muscle_meshes/muscle_$id"
mkdir -p "$muscle_dir"

# Save parameters to a text file
cat > "$muscle_dir/configuration.txt" << EOF
id = $id
n_fibers = $n_points
nx = $nx
ny = $ny
nz = $nz
L = $L
Rmax = $Rmax
Rmin = $Rmin
EOF

# Generate FEM mesh
cd generate_mesh
python generate_fem_mesh.py "$id" "$L" "$Rmax" "$Rmin" "$nx" "$ny" "$nz"

# Generate biomesh input
python generate_vector_field.py "$id" "$L" "$Rmax" "$Rmin"
python generate_seed_points.py "$id" "$n_points" "$seed_radius" "$zmin"
cd ..

# Generate fibers with biomesh
cp "$muscle_dir/vector_field_$id.vtk" ~/Software/biomesh/build/examples
cp "$muscle_dir/fiber_seed_points_$id.json" ~/Software/biomesh/build/examples


#cd ~/Software/biomesh/build/examples
cd $BIOMESH_BUILD/examples || exit 1

./ellipsoid fibers_$id vector_field_$id.vtk fiber_seed_points_$id.json
cd -

mv ~/Software/biomesh/build/results/fibers_$id* $muscle_dir

