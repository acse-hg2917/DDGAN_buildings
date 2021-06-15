import os,sys
import vtk,vtktools
import u2r
import numpy as np
import matplotlib.pyplot as plt
from nirom_dd_tools import * 


def write_sing_values(s_values):
    f= open('singular_values.dat',"w+")
    f.write('# index, s_values, normalised s_values, cumulative energy \n' )
    for k in range(len(s_values)):
        #f.write('# field: %s\n' % field[k])
        total = 0.0
        s_values = s_values[k]
        for i in range(len(s_values)):
            total = total + s_values[i]*s_values[i]

        running_total = 0.0
        for i in range(len(s_values)):
            running_total = running_total + s_values[i]*s_values[i]
            f.write ('%d %g %g %18.10g \n' % (i, s_values[i], s_values[i]/s_values[0], running_total/total) )
    f.close()
    return


def get_clean_vtk_file(filename):
    "Removes fields and arrays from a vtk file, leaving the coordinates/connectivity information."
    vtu_data = vtktools.vtu(filename)
    clean_vtu = vtktools.vtu()
    clean_vtu.ugrid.DeepCopy(vtu_data.ugrid)
    fieldNames = clean_vtu.GetFieldNames()
# remove all fields and arrays from this vtu
    for field in fieldNames:
        clean_vtu.RemoveField(field)
        fieldNames = clean_vtu.GetFieldNames()
        vtkdata=clean_vtu.ugrid.GetCellData()
        arrayNames = [vtkdata.GetArrayName(i) for i in range(vtkdata.GetNumberOfArrays())]
    for array in arrayNames:
        vtkdata.RemoveArray(array)
    return clean_vtu

#(nNodes, reconstruction_on_mesh[iTime*nScalar:(iTime+1)*nScalar,:], template_vtu, original_data[0][iTime*nDim:(iTime+1)*nDim], iTime)
def create_vtu_file(path, nNodes, value_mesh_twice_interp, filename, orig_vel, iTime):
    velocity_field = np.zeros((nNodes,3))
    velocity_field[:,0:nDim] = np.transpose(value_mesh_twice_interp[0:nDim,:]) # streamwise component only

    difference = np.zeros((nNodes,3))
    difference[:,0:nDim] = np.transpose(value_mesh_twice_interp[0:nDim,:]) - orig_vel # streamwise component only
    difference = difference / np.max(velocity_field)

    clean_vtk = get_clean_vtk_file(filename)
    new_vtu = vtktools.vtu()
    new_vtu.ugrid.DeepCopy(clean_vtk.ugrid)
    new_vtu.filename = path + 'recon_' + str(iTime) + '.vtu'
    new_vtu.AddField('Velocity',velocity_field)
    new_vtu.AddField('Original',orig_vel)
    new_vtu.AddField('Velocity_diff',difference)
    new_vtu.Write()
    return

def create_vtu_file_timelevel(nNodes, value_mesh_twice_interp, template_vtu, iTime):
    velocity_field = np.zeros((nNodes,3))
    velocity_field[:,0:nDim] = np.transpose(value_mesh_twice_interp[0:nDim,:]) # streamwise component only

#    difference = np.zeros((nNodes,3))
#    difference[:,0:nDim] = np.transpose(value_mesh_twice_interp[0:nDim,:]) - orig_vel # streamwise component only

    clean_vtk = get_clean_vtk_file(template_vtu)
    new_vtu = vtktools.vtu()
    new_vtu.ugrid.DeepCopy(clean_vtk.ugrid)
    new_vtu.filename = 'reconstructed_' + str(iTime) + '.vtu'
    new_vtu.AddField('Velocity',velocity_field)
    #new_vtu.AddField('Velocity_diff',difference)
    new_vtu.Write()
    return


def get_grid_end_points(grid_origin,grid_width,iGrid ):
    return np.array(( grid_origin[0]+iGrid*grid_width[0], grid_origin[1] +iGrid*grid_width[1]))



## settings

snapshot_data_location = '/home/cheaney/Results/nirom_test_fpc_nPOD_20_nSnap_100/snapshots_CG/'
snapshot_file_base = 'fpc_2D_Re3900_CG_'

nTime = 200
field_names = ['Velocity']
nFields = len(field_names)

xlength = 2.2  # should work this out from coordinates, below 
ylength = 0.41 # 

# set grid size
nGrids = 1 # 4 or 1
if nGrids==4:
    nx = 55
    ny = 42
    nz = 1 # nz = 1 for 2D problems
elif nGrids==1:
    nx = 221
    ny = 42  
    nz = 1 # nz = 1 for 2D problems
else:
    print('nx, ny, nz not known for ', nGrids, 'grids')

grid_origin = [0.0,0.0] 
grid_width = [xlength/nGrids,0.0] 
#ddx = np.array((0.01,0.01))
ddx = np.array((xlength/(nGrids*(nx-1)),ylength/(ny-1))) 
print ('ddx', ddx)


# get a vtu file (any will do as the mesh is not adapted)
filename = snapshot_data_location + snapshot_file_base + '0.vtu'
representative_vtu = vtktools.vtu(filename)
coordinates = representative_vtu.GetLocations()

nNodes = coordinates.shape[0] # vtu_data.ugrid.GetNumberOfPoints()
nEl = representative_vtu.ugrid.GetNumberOfCells()
print('nEl', nEl, type(nEl), 'nNodes', nNodes) 
#nNodes = 3571
#nEl = 6850
nloc = 3 # number of local nodes, ie three nodes per element (in 2D)
nScalar = 2 # dimension of fields
nDim = 2 # dimension of problem (no need to interpolate in the third dimension)

x_all = np.transpose(coordinates[:,0:nDim]) ### coords n,3  x_all 2,n


# get global node numbers
x_ndgln = np.zeros((nEl*nloc), dtype=int)
for iEl in range(nEl):
    n = representative_vtu.GetCellPoints(iEl) + 1
    x_ndgln[iEl*nloc:(iEl+1)*nloc] = n

# -------------------------------------------------------------------------------------------------
# find node duplications when superposing results
my_field = representative_vtu.GetField(field_names[0])[:,0]
my_field = 1
nScalar_test = 1
# for one timestep
# for one field
value_mesh = np.zeros((nScalar_test,nNodes,1)) # nTime=1
value_mesh[:,:,0] = np.transpose(my_field)
superposed_grids = np.zeros((nNodes))
for iGrid in range(nGrids):
    block_x_start = get_grid_end_points(grid_origin, grid_width, iGrid) 

    zeros_on_mesh = 0
    value_grid = u2r.simple_interpolate_from_mesh_to_grid(value_mesh,x_all,x_ndgln,ddx,block_x_start,nx,ny,nz,zeros_on_mesh, nEl,nloc,nNodes,nScalar_test,nDim,1)

    zeros_on_grid = 1
    value_back_on_mesh = u2r.interpolate_from_grid_to_mesh(value_grid, block_x_start, ddx, x_all, zeros_on_grid, nScalar_test,nx,ny,nz,nNodes,nDim, 1)

    superposed_grids = superposed_grids + np.rint(np.squeeze(value_back_on_mesh))

superposed_grids = np.array(superposed_grids, dtype='int') 
duplicated_nodal_values = []
for iNode in range(nNodes):
    if superposed_grids[iNode] == 0:
        # this is bad news - the node hasn't appeared in any grid
        print ('zero:', iNode) 
    elif superposed_grids[iNode] == 2:
        print ('two:', iNode)
        # the node appears in two grids - deal with this later
        duplicated_nodal_values.append(iNode) 
    elif superposed_grids[iNode] != 1:
        # most of the nodes will appear in one grid
        print ('unknown:', iNode, superposed_grids[iNode]) 


# -------------------------------------------------------------------------------------------------
# build up the snapshots matrix from solutions on each of the grids 
offset = 500 # for time level - at which time level to start taking the snapshots
snapshots_data = []
for iField in range(nFields):
    #nDoF = nNodes # could be different value per field
    snapshots_data.append(np.zeros((nx*ny*nz*nDim, nGrids*nTime)))


#value_mesh = np.zeros((nScalar,nNodes)) # no need to initialise - overwritten 
for iTime in range(nTime):

    #print('')
    #print('time level', iTime)

    filename = snapshot_data_location + snapshot_file_base + str(offset+iTime) + '.vtu'
    vtu_data = vtktools.vtu(filename)

    for iField in range(nFields):

        my_field = vtu_data.GetField(field_names[iField])[:,0:nDim]

        for iGrid in range(nGrids):

            block_x_start = get_grid_end_points(grid_origin, grid_width, iGrid) 
            if iTime==0:
                print('block_x_start', block_x_start)
  

            #value_mesh = np.zeros((nScalar,nNodes,nTime)) # nTime - this must need initialising here
            #value_mesh[:,:,iTime] = np.transpose(my_field)

            value_mesh = np.transpose(my_field) # size nScalar,nNodes

            # interpolate field onto structured mesh
            # feed in one result at t time (no need to store in value_mesh in this case)
            zeros_beyond_mesh = 0 # 0 extrapolate solution (for the cylinder in fpc); 1 gives zeros for nodes outside mesh
            #value_grid = u2r.simple_interpolate_from_mesh_to_grid(value_mesh[:,:,iTime],x_all,x_ndgln,ddx,block_x_start,nx,ny,nz,zeros_beyond_mesh,nEl,nloc,nNodes,nScalar,nDim,1) 
            value_grid = u2r.simple_interpolate_from_mesh_to_grid(value_mesh, x_all, x_ndgln, ddx, block_x_start, nx, ny, nz, zeros_beyond_mesh, nEl, nloc, nNodes, nScalar, nDim,1) 

            snapshots_data[iField][:,iTime*nGrids+iGrid] = value_grid.reshape(-1)

# ---------------------------------------------------------------------------------------
# apply POD to the snapshots
# some POD truncation settings
cumulative_tol = 0.99
nPOD = [nTime] # len(nPOD) = nFields
nPOD = [-2]
nPOD = [10] # 100 50 10

bases = []
singular_values = []

for iField in range(nFields):

    snapshots_matrix = snapshots_data[iField]
    nrows, ncols = snapshots_matrix.shape

    if nrows > ncols:
        SSmatrix = np.dot(snapshots_matrix.T, snapshots_matrix)
    else:
        SSmatrix = np.dot(snapshots_matrix, snapshots_matrix.T)
        print('WARNING - CHECK HOW THE BASIS FUNCTIONS ARE CALCULATED WITH THIS METHOD')

    print('SSmatrix', SSmatrix.shape)
    eigvalues, v = np.linalg.eigh(SSmatrix)
    eigvalues =  eigvalues[::-1]
    # get rid of small negative eigenvalues (there shouldn't be any as the eigenvalues of a real, symmetric
    # matrix are non-negative, but sometimes very small negative values do appear)
    eigvalues[eigvalues<0] = 0
    s_values = np.sqrt(eigvalues)
    #print('s values', s_values[0:20]) 

    singular_values.append(s_values)

    cumulative_info = np.zeros(len(eigvalues))
    for j in range(len(eigvalues)):
        if j==0:
            cumulative_info[j] = eigvalues[j]
        else: 
            cumulative_info[j] = cumulative_info[j-1] + eigvalues[j]

    cumulative_info = cumulative_info / cumulative_info[-1]
    nAll = len(eigvalues)

#if nPOD = -1, use cumulative tolerance
#if nPOD = -2 use all coefficients (or set nPOD = nTime)
#if nPOD > 0 use nPOD coefficients as defined by the user

    if nPOD[iField] == -1:
        # SVD truncation - percentage of information captured or number 
        cumulative_tol = nirom_options.compression.cumulative_tol[iField]
        nPOD_iField = sum(cumulative_info <= cumulative_tol) #tolerance
        nPOD[iField] = nPOD_iField                
    elif nPOD[iField] == -2:
        nPOD_iField = nAll
        nPOD[iField] = nPOD_iField                
    else:
        nPOD_iField = nPOD[iField]

    print "retaining", nPOD_iField, "basis functions of a possible", len(eigvalues) 


    basis_functions = np.zeros((nx*ny*nz*nDim,nPOD_iField)) # nDim should be nScalar?
    for j in reversed(range(nAll-nPOD_iField,nAll)):
        Av = np.dot(snapshots_matrix,v[:,j])
        basis_functions[:,nAll-j-1] = Av/np.linalg.norm(Av)

    bases.append(basis_functions)

write_sing_values(singular_values)

# get reconstructed snapshots
reconstruction_data = []
for iField in range(nFields):

    basis = bases[iField]

    snapshots_matrix = snapshots_data[iField]
    print('snapshots_matrix', snapshots_matrix.shape)

    reconstruction_on_mesh = np.zeros((nScalar*nTime,nNodes))
    #reconstruction_on_mesh_from_one_grid = np.zeros((nScalar,nNodes))

    for iGrid in range(nGrids):

        #:,iTime*nGrids+iGrid
        # want solutions in time for a particular grid
        snapshots_per_grid = np.zeros((nx*ny*nz*nDim,nTime))
        for iTime in range(nTime):
            #print('taking snapshots from', iTime*nGrids+iGrid ) 
            snapshots_per_grid[:,iTime] = snapshots_matrix[:,iTime*nGrids+iGrid]

        reconstruction =  np.dot( basis, np.dot( basis.T, snapshots_per_grid ) )
        #print('reconstruction', reconstruction.shape)         
        #reconstruction_data.append(reconstruction)    

        #print ('recon shape',reconstruction.shape)
        reconstruction_grid = reconstruction.reshape(nScalar,nx,ny,nTime)
        #print ('recon shape just before interpolating back onto mesh',reconstruction.reshape(nScalar,nx,ny,nTime).shape)

        # plot solution on each grid at 4 time steps
        #fig, axs = plt.subplots(2, 2, figsize=(15,15))
        #if iGrid==0:
        #    levels = np.linspace(0, 4, 5)
        #elif iGrid==1:
        #    levels = np.linspace(5, 9, 5)
        #icount = 0
        #for col in range(2):
        #    for row in range(2):
        #        ax = axs[row, col]
        #        ax.set_title('time '+str(icount))
        #        pcm = ax.contourf(reconstruction_grid[0,:,:,icount],levels=levels) 
        #        fig.colorbar(pcm,ax=ax) 
        #        icount += 1 
        #plt.show()


        block_x_start = get_grid_end_points(grid_origin, grid_width, iGrid) 
        if iTime==0:
            print('block_x_start', block_x_start)

        for iTime in range(nTime):

            zeros_beyond_grid = 1 # 0 extrapolate solution; 1 gives zeros for nodes outside grid
            reconstruction_on_mesh_from_one_grid = u2r.interpolate_from_grid_to_mesh(reconstruction_grid[:,:,:,iTime], block_x_start, ddx, x_all, zeros_beyond_grid, nScalar,nx,ny,nz,nNodes,nDim, 1)

            #print('reconstruction_on_mesh_from_one_grid - about to add solutions',reconstruction_on_mesh_from_one_grid.shape)
            reconstruction_on_mesh[nScalar*iTime:nScalar*(iTime+1),:] = reconstruction_on_mesh[nScalar*iTime:nScalar*(iTime+1),:] + np.squeeze(reconstruction_on_mesh_from_one_grid)

    reconstruction_on_mesh[:,duplicated_nodal_values] = 0.5*reconstruction_on_mesh[:,duplicated_nodal_values]        
    reconstruction_data.append(reconstruction_on_mesh)    



original_data = []
#for ifield in range(nFields):
#    nDoF = nNodes # could be different value per field
#    original_data.append(np.zeros((nNodes, nDim*nTime)))
original = np.zeros((nNodes, nDim*nTime))
for iTime in range(nTime):

    #print('')
    #print('time level', iTime)

    filename = snapshot_data_location + snapshot_file_base + str(offset+iTime) + '.vtu'
    vtu_data = vtktools.vtu(filename)

    #original = np.zeros((nNodes, nDim*nTime))

    for iField in range(nFields):

        #vtu_data = vtktools.vtu(filename)
        my_field = vtu_data.GetField(field_names[iField])[:,0:nDim]
        original[:,iTime*nDim:(iTime+1)*nDim] = my_field
        
    #print('original.shape',original.shape)
    #original_data.append(original)

# make diretory for results
path_to_reconstructed_results = 'reconstructed_results/'
if not os.path.isdir(path_to_reconstructed_results):
    os.mkdir(path_to_reconstructed_results)

template_vtu = snapshot_data_location + snapshot_file_base + '0.vtu'
for iTime in range(nTime):
    for iField in range(nFields):
        reconstruction_on_mesh = reconstruction_data[iField]
        #print ('recon shape',reconstruction.shape)

    # for more than one field, will this work? 
    #create_vtu_file_timelevel(nNodes, reconstruction_on_mesh[iTime*nScalar:(iTime+1)*nScalar,:], template_vtu, iTime)
    create_vtu_file(path_to_reconstructed_results, nNodes, reconstruction_on_mesh[iTime*nScalar:(iTime+1)*nScalar,:], template_vtu, original[:,iTime*nDim:(iTime+1)*nDim], iTime)


print('ddx',ddx.shape,ddx)

print('grid lengths',   (nx-1)*ddx[0], (ny-1)*ddx[1] ) #, (nz-1)*ddx[2])

print('Finished.')







