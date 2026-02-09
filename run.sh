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

# Generate vector field
cd generate_mesh
python generate_vector_field.py "$id" "$L" "$Rmax" "$Rmin"

# Generate seed points
python generate_seed_points.py "$id" "$n_points" "$seed_radius" "$zmin"
# python json_to_vtp.py "seed_points$id.json"

# Generate fibers
echo "Run biomesh"
cd ../results/
cp "vector_fields/vector_field_$id.vtk" ~/Software/biomesh/build/examples
cp "fiber_seed_points/seed_points_$id.json" ~/Software/biomesh/build/examples

cd ~/Software/biomesh/build/examples
./ellipsoid fibers_$id vector_field_$id.vtk seed_points_$id.json
# # Generate fem mesh
# cd ~/Documents/ellipsoid_dataset
# python fem_mesh_generator.py "$id" "$L" "$Rmax" "$Rmin" "$nx" "$ny" "$nz"


# # Clean up and organize files
# echo "Copy files"

# mkdir "muscle$id"

# mv "seed_points$id.json" "muscle$id"
# mv "seed_points$id.vtp" "muscle$id"
# mv "vfield$id.vtk" "muscle$id"
# mv "structured_muscle$id.vtk" "muscle$id"
cd -
cd muscle_meshes
mkdir muscle_$id
cp ~/Software/biomesh/build/results/fibers_$id* muscle_$id
# rm ~/Software/biomesh/build/results/fibersmuscle$id*
