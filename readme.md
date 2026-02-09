This repository provides a pipeline from mesh generation to a prestretch + contraction simulation of a muscle within the OpenDiHu framework.  

## Prerequisites

The user must install OpenDiHu and BioMesh and define the paths `OPENDIHU_HOME` and `BIOMESH_HOME`. You can define the paths by running
```
export OPENDIHU_HOME=/path/to/opendihu
export BIOMESH_HOME=/path/to/biomesh
```

## Structure

TODO:

### 1. Generate mesh

The `generate_mesh` directory contains scripts to generate the fiber and the mechanics mesh required to run an OpenDiHu simulation. 

TODO:

### 2. Simulate muscle prestretch + contraction

- **Compile the code:**

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

- **Run the code:**

    TODO:


