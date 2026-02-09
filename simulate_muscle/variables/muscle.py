import sys
import json
import read_structured_vtk

muscle_id = int(sys.argv[-5])
prestretch_force = float(sys.argv[-4])
specific_states_call_frequency = float(sys.argv[-3])

rank_no = int(sys.argv[-2])
n_ranks = int(sys.argv[-1])


# Mesh generation inputs
# muscle_id = 5
muscle_dir = "../../muscle_meshes/muscle_" + str(muscle_id)
geometry_name = "muscle_" + str(muscle_id)

# Simulation input
isometric = True
# prestretch_force = 4
specific_states_call_enable_begin = 0.01  # time of first fiber activation
# specific_states_call_frequency = 0.0525    # frequency of fiber activation

scenario_name = geometry_name + "_prestretch_" + str(prestretch_force) + "_frequency_" + str(specific_states_call_frequency)

# -------------------------------------------------------------------
# FEM mesh generation from .vts file
# -------------------------------------------------------------------
vtk_filename = muscle_dir+"/3D_mesh_" + str(muscle_id) + ".vtk"
points, bs_x, bs_y, bs_z = read_structured_vtk.read_structured_vtk(vtk_filename)
el_x, el_y, el_z = int((bs_x-1)/2), int((bs_y-1)/2), int((bs_z-1)/2)

meshes = { # create 3D mechanics mesh
    "mesh3D": {
        "nElements":            [el_x, el_y, el_z],
        "nodePositions":        points,
        "logKey":               "mesh3D",
        "inputMeshIsGlobal":    True,
        "nRanks":               1,
    }
}

# -------------------------------------------------------------------
# fiber mesh generation from .json file
# -------------------------------------------------------------------

fiber_file = muscle_dir+"/fibers_"+str(muscle_id)+".json"
with open(fiber_file,"r") as f:
	fdata = json.load(f)
    
if not isinstance(fdata, dict):
    raise TypeError(f"Expected dict in {fiber_file}, got {type(fdata)}")

fiber_idx = 0
for fiber in fdata:
	fdict = fdata[fiber]
	npos = [[fdict[ii]['x'],fdict[ii]['y'],fdict[ii]['z']] for ii in range(len(fdict)) ]
	meshName = "fiber{}".format(fiber_idx)
	meshes[meshName] = {
			"nElements":		    [len(fdict)-1],
			"nodePositions":	    npos,
			"inputMeshIsGlobal":	True,
			"nRanks":				n_ranks
	}
	fiber_idx += 1
     
n_fibers = fiber_idx

# Boundary conditions: 
prestretch_dirichlet_bc = {}
contraction_dirichlet_bc = {}

# set Dirichlet BC
k = 0
for j in range(bs_y):
  for i in range(bs_x):
    prestretch_dirichlet_bc[k*bs_x*bs_y + j*bs_x + i] = [None,None,0.0,None,None,None]
    contraction_dirichlet_bc[k*bs_x*bs_y + j*bs_x + i] = [None,None,0.0,None,None,None]
for j in range(bs_y):
  prestretch_dirichlet_bc[k*bs_x*bs_y + j*bs_x + 0][0] = 0.0
  
for i in range(bs_x):
  prestretch_dirichlet_bc[k*bs_x*bs_y + 0*bs_x + i][1] = 0.0

if isometric:
    k = bs_z-1
    for j in range(bs_y):
        for i in range(bs_x):
            contraction_dirichlet_bc[k*bs_x*bs_y + j*bs_x + i] = [None,None,0.0,None,None,None]
       
# set Neumann BC
k = el_z-1
prestretch_neumann_bc = [{"element": k*el_x*el_y + j*el_x + i, "constantVector": [0, 0, prestretch_force], "face": "2+"} for j in range(el_y) for i in range(el_x)]
contraction_neumann_bc = [{"element": k*el_x*el_y + j*el_x + i, "constantVector": [0, 0, 0], "face": "2+"} for j in range(el_y) for i in range(el_x)]


