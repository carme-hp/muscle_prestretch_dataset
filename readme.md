This repository provides a pipeline from mesh generation to a prestretch + contraction simulation of a muscle within the OpenDiHu framework.  

## Prerequisites

The user must install OpenDiHu and BioMesh and define the paths `OPENDIHU_HOME` and `BIOMESH_HOME`. You can define the paths by running
```
export OPENDIHU_HOME=/path/to/opendihu
export BIOMESH_HOME=/path/to/biomesh
```

### How-to generate meshes

Run `. create_muscle_mesh.sh` to generate all the mesh files required for an OpenDiHu simulation, e.g., `fibers_<muscle_id>.json` and `3D_mesh_<muscle_id>.vtk`. 

The script generates the meshes for a symmetric muscle by 3 geometrical parameters: length `L`, radius at the center `Rmax` and radius at the extremes `Rmin`. In addition, the user must provide an integer `muscle_id` to identify the muscle. All the generated files are saved to `muscle_meshes/muscle_<muscle_id>`.

Instead of running the script, you can also generate the meshes manually as follows.

** How to generate 3D meshes**

- Run `generate_fem_mesh.py`
    It generates a structured 3D mesh and saves it to `3D_mesh_<muscle_id>.vtk`. 

    How-to run: 
    ```
    python generate_fem_mesh.py "$id" "$L" "$Rmax" "$Rmin" "$nx" "$ny" "$nz"
    ```

    **Warning:** Always choose odd numbers for `nx`, `ny` and `nz`.

** How to generate 3D meshes**

- Create input files for biomesh
    - Run `generate_vector_field.py`

        It generates a 3D vector field that is used to define fiber directions in biomesh.

        How-to run:
        ```
        python generate_vector_field.py "$id" "$L" "$Rmax" "$Rmin"
        ```
        
    - Run `generate_seed_points.py`

        It generates a list of points located at a plane normal to the z-axis. The points are used as starting vertices for fiber generation in biomesh.

        How-to run:
        ```
        python generate_vector_field.py "$id" "$L" "$Rmax" "$Rmin"
        ```
    
- Run biomesh

    TODO: version and path 
    
    How-to run:
    ```
    ./ellipsoid fibers_$id vector_field_$id.vtk fiber_seed_points_$id.json
    ```
    
### 2. Simulate muscle prestretch + contraction

**How to compile the simulation**

Assuming you have defined `OPENDIHU_HOME`, the easiest way to compile the code is by defining 
```
alias sr='$OPENDIHU_HOME/scripts/shortcuts/sr.sh'
alias mkorn='$OPENDIHU_HOME/scripts/shortcuts/mkorn.sh'
``` 
and then run
```
mkorn && sr
```
This wil compile the code and write the executable to the newly generated folder `build_release/`. 

**How to run the simulation**

    TODO:


