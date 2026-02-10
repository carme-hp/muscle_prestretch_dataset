This repository provides a pipeline from mesh generation to a prestretch + contraction simulation of a muscle within the OpenDiHu framework.  

## Prerequisites

The user must install [OpenDiHu](https://github.com/opendihu/opendihu) (use the branch `develop`) and [BioMesh](https://github.com/opendihu/biomesh) (use branch `filename-as-input`) and define the paths `OPENDIHU_HOME` and `BIOMESH_BUILD`. You can define the paths by running the next lines in every new terminal or adding them to your `.bashrc`.
```
export OPENDIHU_HOME=/path/to/opendihu
export BIOMESH_BUILD=/path/to/biomesh/build
```
Note that your build directory might not be in `biomesh/build` and this needs to be changed accordingly. 

## 1. Generate Meshes

Run `. create_muscle_mesh.sh` to generate all the mesh files required for an OpenDiHu simulation, e.g., `fibers_<muscle_id>.json` and `3D_mesh_<muscle_id>.vtk`. 

The script generates the meshes for a symmetric muscle by 3 geometrical parameters: length `L`, radius at the center `Rmax` and radius at the extremes `Rmin`. In addition, the user must provide an integer `muscle_id` to identify the muscle. All the generated files are saved to `muscle_meshes/muscle_<muscle_id>`.

For example:
```
Enter dataset id: 3
Enter number of fibers: 6
Enter number of fem points: 7 7 21
Enter length of the muscle: 5
Enter radius at the center: 2
Enter radius at the extreme: 1.1

```

Instead of running the script, you can also generate the meshes manually as follows.

### How to generate 3D meshes

- Run `generate_fem_mesh.py`
    It generates a structured 3D mesh and saves it to `3D_mesh_<muscle_id>.vtk`. 

    How-to run: 
    ```
    python generate_fem_mesh.py "$id" "$L" "$Rmax" "$Rmin" "$nx" "$ny" "$nz"
    ```

    **Warning:** Always choose odd numbers for `nx`, `ny` and `nz`.

### How to generate 3D meshes

- Create input files for biomesh
    - Run `generate_vector_field.py`

        It generates a 3D vector field used to define fiber directions in biomesh.

        How-to run:
        ```
        python generate_vector_field.py "$id" "$L" "$Rmax" "$Rmin"
        ```
        
    - Run `generate_seed_points.py`

        It generates a list of points located on a plane normal to the z-axis. The points are used as starting vertices for fiber generation in biomesh.

        How-to run:
        ```
        python generate_vector_field.py "$id" "$L" "$Rmax" "$Rmin"
        ```
    
- Run biomesh

    The biomesh executable is located in `$BIOMESH_BUILD` , inside the folder `examples/`.
    
    How-to run:
    ```
    ./ellipsoid fibers_$id vector_field_$id.vtk fiber_seed_points_$id.json
    ```
    
## 2. Simulate muscle prestretch + contraction

Since the parameters are specified in python scripts, you only need to compile the simulation once. 

### How to compile the simulation

Assuming you have defined `OPENDIHU_HOME`, the easiest way to compile the code is by defining 
```
alias sr='$OPENDIHU_HOME/scripts/shortcuts/sr.sh'
alias mkorn='$OPENDIHU_HOME/scripts/shortcuts/mkorn.sh'
``` 
and then run
```
mkorn && sr
```
This will compile the code and write the executable to the newly generated folder `build_release/`. 

### How to run the simulation

Run `. run_muscle_simulation.sh` to run a prestretch + contraction simulation for a given muscle. All the simulation files are moved to the `results/` folder. For the simulation you require three inputs: `muscle_id`, `prestretch_force` and the activation frequency `specific_states_call_frequency`.

To run it manually:

```
cd simulate_muscle/build_release
./muscle ../settings_contraction_with_prestretch.py muscle.py "$id" "$prestretch_force" "$specific_states_call_frequency"
```

You can also run multiple simulations for the same muscle with different prestretch and activation. Simply run `. run_multiple_simulation.sh` after choosing your desired parameters. For example:

```
# Define parameter sweeps
prestretch_forces=(0.0 1.0 2.0)
frequencies=(0.1 0.01)
```
