# %%
import os,sys
from numpy.lib.npyio import save
import vtk,vtktools
import u2r
import numpy as np
import matplotlib.pyplot as plt
import random
from nirom_dd_tools import * 
import copy

### MY version of nirom_dd_orig.py for FLOW_PAST_IMPLICIT_BLOCK_84.vtu

# some functions of tools.io copied in
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

#code for full domain case
# def get_grid_end_points(grid_origin,grid_width,iGrid ):
#     return np.array(( grid_origin[0]+iGrid*grid_width[0], grid_origin[1] +iGrid*grid_width[1]))#

def get_grid_end_points(grid_origin,grid_width):
    return np.array((grid_origin[0]+grid_width[0], grid_origin[1] +grid_width[1]))

def plot_grid(grid_origin, grid_width, nx, ny):
    # include plot of entire coordinates with grid
    # plt.figure(figsize=(9,9))
    plt.plot(coordinates[:,0], coordinates[:,1], 'g.', ms = 0.3, label = 'angle = {}˚'.format(random_angle)) # corrdinates

    # code for just the edges
    # plt.plot([grid_origin[0], grid_origin[0]+grid_width[0]], [grid_origin[1], grid_origin[1]], 'ko-') #1
    # plt.plot([grid_origin[0], grid_origin[0]], [grid_origin[1], grid_origin[1]+grid_width[1]], 'ko-') #2
    # plt.plot([grid_origin[0], grid_origin[0]+grid_width[0]], [grid_origin[1]+grid_width[1], grid_origin[1]+grid_width[1]], 'ko-') #3
    # plt.plot([grid_origin[0]+grid_width[0], grid_origin[0]+grid_width[0]], [grid_origin[1], grid_origin[1]+grid_width[1]], 'ko-') #4

    for d in range(ny + 1):
        if d%4 == 0:
            plt.plot([grid_origin[0], grid_origin[0]+grid_width[0]], [grid_origin[1]+d*ddx[1], grid_origin[1]+d*ddx[1]], 'k-', lw = 1.2)   #horizontal
            if ny == nx:
                plt.plot([grid_origin[0]+d*ddx[1], grid_origin[0]+d*ddx[1]], [grid_origin[1], grid_origin[1]+grid_width[1]], 'k-', lw = 1.2) #vertical

    if ny != nx:
        for d in range (nx + 1): #vertical
            if d%4 == 0:
                plt.plot([grid_origin[0]+d*ddx[0], grid_origin[0]+d*ddx[0]], [grid_origin[1], grid_origin[1]+grid_width[1]], 'k-', lw = 1.2) #vertical


    # plt.grid(':')
    # plt.tight_layout()
    # plt.show()

def rotate_mesh(angle):
    theta = np.radians(angle)

    #shift coordinates so that they are centred at (0,0)
    # for i in range(coordinates.shape[0]):
    #     coordinates[i][0] -= 1.5
    #     coordinates[i][1] -= 1.5

    new_mesh = np.zeros(coordinates.shape)

    for i in range(coordinates.shape[0]):
        new_mesh[i][0] = (coordinates[i][0]-1.5)*np.cos(theta) - (coordinates[i][1]-1.5)*np.sin(theta)
        new_mesh[i][1] = (coordinates[i][0]-1.5)*np.sin(theta) + (coordinates[i][1]-1.5)*np.cos(theta)

    #rotate the velocity field as well

    return new_mesh

def rotate_vel(angle):
    theta = np.radians(angle)
    new_mesh = np.zeros(velocity_field.shape)

    for i in range(coordinates.shape[0]):
        new_mesh[i][0] = (velocity_field[i][0])*np.cos(theta) - (velocity_field[i][1])*np.sin(theta)
        new_mesh[i][1] = (velocity_field[i][0])*np.sin(theta) + (velocity_field[i][1])*np.cos(theta)

    return new_mesh

def select_gridpoint():
    min_x = min(coordinates[:,0])
    max_x = max(coordinates[:,0])
    min_y = min(coordinates[:,1])
    max_y = max(coordinates[:,1])

    # plt.plot(min_x+0.5,min_y+0.5, 'ro' )
    # plt.plot(min_x+0.5,max_y-1.0, 'ro' )
    # plt.plot(max_x-1.0,min_y+0.5, 'ro' )
    # plt.plot(max_x-1.0,max_y-1.0, 'ro' )

    grid_origin = [3,3]
    while np.sqrt(grid_origin[0]**2+grid_origin[1]**2) >= 1.3:
        # print("finding point - ", np.sqrt(grid_origin[0]**2+grid_origin[1]**2))
        grid_origin = [random.uniform(min_x+0.5, max_x-1.2), random.uniform(min_y+0.5, max_y-1.2)]

    return grid_origin

def sample_starshape(mesh, grid_origin):
    """
    Returns a snapshot matrix of shape (5,nx*ny) and 
    snapshot_ae of shape (5,nx,ny) with given
    mesh and central grid origin for the starshape grid formation
    """
    grid_point_0 = [grid_origin[0], grid_origin[1]+grid_width[1]]
    grid_point_1 = [grid_origin[0]-grid_width[0], grid_origin[1]]
    grid_point_2 = [grid_origin[0]+grid_width[0], grid_origin[1]]
    grid_point_3 = [grid_origin[0], grid_origin[1]-grid_width[1]]
    #grid_point_4 = grid_origin
    
    grid_list = [grid_point_0,grid_point_1, grid_point_2, grid_point_3, grid_origin]

    s_matrix = np.zeros((nx*ny, 5))
    s_ae = np.zeros((5,nx,ny))

    for iloc in range(5):
        value_grid = u2r.simple_interpolate_from_mesh_to_grid(mesh, x_all, x_ndgln, ddx, grid_list[iloc], nx, ny, nz, zeros_beyond_mesh, nEl, nloc, nNodes, nScalar, nDim,1)
        s_matrix[:,iloc] = value_grid.reshape(-1)
        s_ae[iloc,:,:] = value_grid.reshape((nx,ny))

    return s_matrix, s_ae

def plot_starshape(nSGrids):
    plt.figure(figsize=(8,8))
    plt.subplot(3,3,2)
    plt.imshow(np.rot90(Ssnapshot_ae[5*nSGrids,:,:,2]))

    plt.subplot(3,3,4)
    plt.imshow(np.rot90(Ssnapshot_ae[5*nSGrids+1,:,:,2]))

    plt.subplot(3,3,5)
    plt.imshow(np.rot90(Ssnapshot_ae[5*nSGrids+4,:,:,2]))

    plt.subplot(3,3,6)
    plt.imshow(np.rot90(Ssnapshot_ae[5*nSGrids+2,:,:,2]))

    plt.subplot(3,3,8)
    plt.imshow(np.rot90(Ssnapshot_ae[5*nSGrids+3,:,:,2]))

    plt.tight_layout()
    # plt.show()


## settings

snapshot_data_location = '/Users/gohannaago/Desktop/IC_MSc_ACSE/9_PROJECT/steady_state_data/'
snapshot_file_base = 'Flow_past_buildings_Re1_training_'

nTime = 1
field_names = ['Velocity', 'VelocityAbsorption']
nFields = len(field_names)

xlength = 3.0  
ylength = 3.0

grid_width = [0.48,0.48]
# spacing inside small grid 
nx = int(grid_width[0]*100)
ny = nx
nz = 1
ddx = np.array((0.01,0.01))

# set number of grids - samples/snapshots to take
nGrids = 3000
# Turn on/off snapshots matrix 
save_snapshots = False
save_stargrid = True
# Turn on/off save first 20 images
save_imgs = False

# get a vtu file (any will do as the mesh is not adapted)
filename = snapshot_data_location + snapshot_file_base + '400.vtu'
representative_vtu = vtktools.vtu(filename)
coordinates_org = representative_vtu.GetLocations() #coordinates of the nodes
coordinates = coordinates_org

nNodes = coordinates.shape[0] # vtu_data.ugrid.GetNumberOfPoints()
nEl = representative_vtu.ugrid.GetNumberOfCells()
# print('nEl', nEl, type(nEl), 'nNodes', nNodes) 
#nNodes = 375627
#nEl = 125209
nloc = 3 # number of local nodes, ie three nodes per element (in 2D)
# nScalar = 2 # dimension of fields , 2 = u and v
nScalar = 1 #because I calculate u and v separately
nDim = 2 # dimension of problem (no need to interpolate in the third dimension)

# x_all = np.transpose(coordinates[:,0:nDim]) ### coords n,3  x_all 2,n

# get global node numbers
x_ndgln = np.zeros((nEl*nloc), dtype=int)
for iEl in range(nEl):
    n = representative_vtu.GetCellPoints(iEl) + 1
    x_ndgln[iEl*nloc:(iEl+1)*nloc] = n

# %%
# -------------------------------------------------------------------------------------------------
#SNAPSHOT MATRIX
#snapshot matrix configuration outside loop
snapshots_matrix = np.zeros((nx*ny*3, nGrids))
snapshot_ae = np.zeros((nGrids,nx,ny,3))

#saving coordinates of grid_origins of each grid (for later starshape grid)
origin_save = np.zeros((nGrids,2))
#saving rotation angle
rangle_save = np.zeros(nGrids)



# full run - later to iterate
for iGrid in range (nGrids): #replace with ngrids later
    # plt.subplot(nGrids//4,(nGrids//4+nGrids%4),iGrid+1)
    if iGrid <= 20 and save_imgs:
        plt.figure(figsize=(9,9))

    if iGrid%10 == 0:
        print("Currently Generating Grid ", iGrid+1, " out of ", nGrids, " Grids")
    # rotate mesh and select grid
    coordinates = coordinates_org
    random_angle = random.randint(0,360)
    # print("Rotation angle of mesh = ", random_angle)
    rangle_save[iGrid] = random_angle
    #rotate coordinates
    coordinates = rotate_mesh(random_angle)

    # velocity_field = representative_vtu.GetField(field_names[0])[:,:nDim] #field name 0 is velocity field
    # #rotate velocity
    # # print('velocity_field[0]', velocity_field[0])
    # velocity_field = rotate_vel(random_angle)
    # # print('velocity_field[0]', velocity_field[0])
    # #call velocity absorption field
    # va_field = representative_vtu.GetField(field_names[1])[:,0] #field name 0 is velocity field
    # # print('va_field', va_field[:20])
    # x_all = np.transpose(coordinates[:,0:nDim])
    
    # plot the orientation
    # plt.rcParams.update({'font.size': 18})
    grid_origin = select_gridpoint()
    # print(grid_origin)
    origin_save[iGrid] = grid_origin
    
    # if iGrid <= 20:
    #     plot_grid(grid_origin, grid_width, nx, ny)
    #     plt.tight_layout()
    #     plt.legend(loc = 'best')

    """
    #variables for interpolation
    # block_x_start = get_grid_end_points(grid_origin, grid_width)
    zeros_beyond_mesh = 0

    #interpolate & append stream velocity
    u_mesh = np.zeros((1,nNodes,1))
    u_mesh[:,:,0] = np.transpose(velocity_field[:,0])
    uvalue_grid = u2r.simple_interpolate_from_mesh_to_grid(u_mesh, x_all, x_ndgln, ddx, grid_origin, nx, ny, nz, zeros_beyond_mesh, nEl, nloc, nNodes, nScalar, nDim,1)
    # print('uvalue_grid.shape', uvalue_grid.shape)
    snapshot_ae[iGrid,:,:,0] = uvalue_grid.reshape((nx,ny))
    snapshots_matrix[:nx*ny*1,iGrid] = uvalue_grid.reshape(-1)
    
    #interpolate & append v velocity
    v_mesh = np.zeros((1,nNodes,1))
    v_mesh[:,:,0] = np.transpose(velocity_field[:,1])
    vvalue_grid = u2r.simple_interpolate_from_mesh_to_grid(v_mesh, x_all, x_ndgln, ddx, grid_origin, nx, ny, nz, zeros_beyond_mesh, nEl, nloc, nNodes, nScalar, nDim,1)
    snapshot_ae[iGrid,:,:,1] = vvalue_grid.reshape((nx,ny))
    snapshots_matrix[nx*ny:nx*ny*2,iGrid] = vvalue_grid.reshape(-1)

    #interpolate & append velocity absoprtion
    va_mesh = np.zeros((1,nNodes,1))
    va_mesh[:,:,0] = np.transpose(va_field)
    vavalue_grid = u2r.simple_interpolate_from_mesh_to_grid(va_mesh, x_all, x_ndgln, ddx, grid_origin, nx, ny, nz, zeros_beyond_mesh, nEl, nloc, nNodes, nScalar, nDim,1)
    # print(vavalue_grid.reshape(-1)[:30])
    snapshot_ae[iGrid,:,:,2] = vavalue_grid.reshape((nx,ny))
    snapshots_matrix[nx*ny*2:nx*ny*3,iGrid] = vavalue_grid.reshape(-1)    

    if iGrid <= 20 and save_imgs:
        plt.savefig('plots/grid_selection_{}'.format(iGrid))
        plt.close()

# print('snapshots_matrix.shape: ', snapshots_matrix.shape)
# print('saved grid origins for later - shape:', origin_save.shape)

# Scale snapshot matrix
min_vel = np.amin(snapshots_matrix[:nx*ny*2,:]) #minimum among u and v velocity
max_vel = np.amax(snapshots_matrix[:nx*ny*2,:])
vel_scaling = 1/(max_vel-min_vel)
va_scaling = 1e-5 #maximum of buildings information is 100000

snapshots_matrix[:nx*ny*2,:] = vel_scaling*(snapshots_matrix[:nx*ny*2,:]-min_vel)
snapshots_matrix[nx*ny*2:nx*ny*3,:] = va_scaling*snapshots_matrix[nx*ny*2:nx*ny*3,:]

#scale snapshots for ae
snapshot_ae[:,:,:,:2] = vel_scaling*(snapshot_ae[:,:,:,:2]-min_vel)
snapshot_ae[:,:,:,2] = va_scaling*snapshot_ae[:,:,:,2]

# Export snapshots for training ae
if save_snapshots:
    np.save("ae_data.npy", snapshot_ae)
    np.save("snapshots_pod.npy", snapshots_matrix)

"""
#%%
np.save("grid_origins.npy", origin_save)
np.save("rotation_angles.npy", rangle_save)

# %%
# -------------------------------------------------------------------------------------------------

# origin_save = np.load("grid_origins.npy")
# rangle_save = np.load("rotation_angles.npy")

print("Generating starshape grids for true ss case")
zeros_beyond_mesh = 0

#set number of starshape grids to use
nSGrids = 2000
assert nSGrids <= nGrids, "Cannot make more starshape grids than number of central grids we have"
Ssnapshots_matrix = np.zeros((nx*ny*3, nSGrids*5))
Ssnapshot_ae = np.zeros((nSGrids*5,nx,ny,3))

#saving coordinates of grid_origins of each grid (for later starshape grid)
origin_save = np.zeros((nSGrids,2))
#saving rotation angle
rangle_save = np.zeros(nSGrids)

for iGrid in range(nSGrids):
    if iGrid % 10 == 0:
        print("Sampling starshape grid ", iGrid+1, " out of ", nSGrids)
    
    #Use Generation of new random grid locations
    # rotate mesh and select grid
    coordinates = coordinates_org
    rangle = random.randint(0,360)
    rangle_save[iGrid] = rangle
    #rotate coordinates
    coordinates = rotate_mesh(rangle)
    #selct gridpoint
    grid_origin = select_gridpoint()
    origin_save[iGrid] = grid_origin

    # # call saved values for grid origin and random angle
    # grid_origin = origin_save[iGrid]
    # rangle = rangle_save[iGrid]

    # rotate mesh and velocity
    coordinates = coordinates_org
    coordinates = rotate_mesh(rangle)
    x_all = np.transpose(coordinates[:,0:nDim])
    velocity_field = representative_vtu.GetField(field_names[0])[:,:nDim] #field name 0 is velocity field
    velocity_field = rotate_vel(rangle)
    va_field = representative_vtu.GetField(field_names[1])[:,0]

    #starshape for u_vel
    u_mesh = np.zeros((1,nNodes,1))
    u_mesh[:,:,0] = np.transpose(velocity_field[:,0])
    u_smatrix, u_sae = sample_starshape(u_mesh, grid_origin)
    Ssnapshots_matrix[:nx*ny,5*iGrid:5*iGrid+5] = u_smatrix
    Ssnapshot_ae[5*iGrid:5*iGrid+5,:,:,0] = u_sae

    #starshape for v_vel
    v_mesh = np.zeros((1,nNodes,1))
    v_mesh[:,:,0] = np.transpose(velocity_field[:,1])
    v_smatrix, v_sae = sample_starshape(v_mesh, grid_origin)
    Ssnapshots_matrix[nx*ny:nx*ny*2,5*iGrid:5*iGrid+5] = v_smatrix
    Ssnapshot_ae[5*iGrid:5*iGrid+5,:,:,1] = v_sae

    #starshape for va
    va_mesh = np.zeros((1,nNodes,1))
    va_mesh[:,:,0] = np.transpose(va_field)
    va_smatrix, va_sae = sample_starshape(va_mesh, grid_origin)
    Ssnapshots_matrix[nx*ny*2:nx*ny*3,5*iGrid:5*iGrid+5] = va_smatrix
    Ssnapshot_ae[5*iGrid:5*iGrid+5,:,:,2] = va_sae

#Scale values
min_vel = np.amin(Ssnapshots_matrix[:nx*ny*2,:]) #minimum among u and v velocity
max_vel = np.amax(Ssnapshots_matrix[:nx*ny*2,:])
vel_scaling = 1/(max_vel-min_vel)
min_va = np.amin(Ssnapshots_matrix[nx*ny*2:nx*ny*3,:]) #minimum among u and v velocity
max_va = np.amax(Ssnapshots_matrix[nx*ny*2:nx*ny*3,:])
va_scaling = 1/(max_va-min_va) #maximum of buildings information is 100000

Ssnapshots_matrix[:nx*ny*2,:] = vel_scaling*(Ssnapshots_matrix[:nx*ny*2,:]-min_vel)
Ssnapshots_matrix[nx*ny*2:nx*ny*3,:] = va_scaling*Ssnapshots_matrix[nx*ny*2:nx*ny*3,:]

#scale snapshots for ae
Ssnapshot_ae[:,:,:,:2] = vel_scaling*(Ssnapshot_ae[:,:,:,:2]-min_vel)
Ssnapshot_ae[:,:,:,2] = va_scaling*Ssnapshot_ae[:,:,:,2]


if save_stargrid:
    np.save("starshape_ae_truess817.npy",Ssnapshot_ae)
    np.save("starshape_pod_truess817.npy", Ssnapshots_matrix)

np.save("grid_origins_0817.npy", origin_save)
np.save("rotation_angles_0817.npy", rangle_save)


textfile = open("nirom_dd_truess_817.txt", 'w')
textfile.write("min_vel = "+str(min_vel) +" , max_vel = "+ str(max_vel)+ " , vel_scaling = "+ str(vel_scaling))
textfile.close()

print("min_vel = ", min_vel, " , max_vel = ", max_vel, " , vel_scaling = ", vel_scaling)
#%%

# for i in range(10):
#     plot_starshape(i)
#     plt.savefig('plots/starshape_{}'.format(i))
#     plt.close()
# # -------------------------------------------------------------------------------------------------
#find POD coefficients for starshape grids
#in POD_imports.py




# # -------------------------------------------------------------------------------------------------
# find node duplications when superposing results (don't need this part for my project)
# my_field = representative_vtu.GetField(field_names[0])[:,0] #u(x) velocity only
# my_field = 1
# nScalar_test = 1
# # for one timestep
# # for one field
# value_mesh = np.zeros((nScalar_test,nNodes,1)) # nTime=1
# value_mesh[:,:,0] = np.transpose(my_field)
# superposed_grids = np.zeros((nNodes))
# for iGrid in range(nGrids):
#     block_x_start = get_grid_end_points(grid_origin, grid_width, iGrid) 

#     zeros_on_mesh = 0
#     value_grid = u2r.simple_interpolate_from_mesh_to_grid(value_mesh,x_all,x_ndgln,ddx,block_x_start,nx,ny,nz,zeros_on_mesh, nEl,nloc,nNodes,nScalar_test,nDim,1)

#     zeros_on_grid = 1
#     value_back_on_mesh = u2r.interpolate_from_grid_to_mesh(value_grid, block_x_start, ddx, x_all, zeros_on_grid, nScalar_test,nx,ny,nz,nNodes,nDim, 1)

#     superposed_grids = superposed_grids + np.rint(np.squeeze(value_back_on_mesh))

# superposed_grids = np.array(superposed_grids, dtype='int') 
# duplicated_nodal_values = []
# for iNode in range(nNodes):
#     if superposed_grids[iNode] == 0:
#         # this is bad news - the node hasn't appeared in any grid
#         print ('zero:', iNode) 
#     elif superposed_grids[iNode] == 2:
#         print ('two:', iNode)
#         # the node appears in two grids - deal with this later
#         duplicated_nodal_values.append(iNode) 
#     elif superposed_grids[iNode] != 1:
#         # most of the nodes will appear in one grid
#         print ('unknown:', iNode, superposed_grids[iNode]) 


# -------------------------------------------------------------------------------------------------
# build up the snapshots matrix from solutions on each of the grids 
# offset = 500 # for time level - at which time level to start taking the snapshots
# snapshots_data = []
# for iField in range(nFields):
#     #nDoF = nNodes # could be different value per field
#     snapshots_data.append(np.zeros((nx*ny*nz*nDim, nGrids*nTime))) #create snapshot for each field we have


# #value_mesh = np.zeros((nScalar,nNodes)) # no need to initialise - overwritten 
# for iTime in range(nTime):

#     #print('')
#     #print('time level', iTime)

#     filename = snapshot_data_location + snapshot_file_base + str(offset+iTime) + '.vtu'
#     vtu_data = vtktools.vtu(filename)

#     for iField in range(nFields):

#         my_field = vtu_data.GetField(field_names[iField])[:,0:nDim]

#         for iGrid in range(nGrids):

#             block_x_start = get_grid_end_points(grid_origin, grid_width, iGrid) #randomly locate
#             if iTime==0:
#                 print('block_x_start', block_x_start)
  

#             #value_mesh = np.zeros((nScalar,nNodes,nTime)) # nTime - this must need initialising here
#             #value_mesh[:,:,iTime] = np.transpose(my_field)

#             value_mesh = np.transpose(my_field) # size nScalar,nNodes

#             # interpolate field onto structured mesh
#             # feed in one result at t time (no need to store in value_mesh in this case)
#             zeros_beyond_mesh = 0 # 0 extrapolate solution (for the cylinder in fpc); 1 gives zeros for nodes outside mesh
#             #value_grid = u2r.simple_interpolate_from_mesh_to_grid(value_mesh[:,:,iTime],x_all,x_ndgln,ddx,block_x_start,nx,ny,nz,zeros_beyond_mesh,nEl,nloc,nNodes,nScalar,nDim,1) 
#             value_grid = u2r.simple_interpolate_from_mesh_to_grid(value_mesh, x_all, x_ndgln, ddx, block_x_start, nx, ny, nz, zeros_beyond_mesh, nEl, nloc, nNodes, nScalar, nDim,1) 

#             snapshots_data[iField][:,iTime*nGrids+iGrid] = value_grid.reshape(-1)

# # ---------------------------------------------------------------------------------------
# apply POD to the snapshots
# some POD truncation settings
# cumulative_tol = 0.99
# nPOD = [nTime] # len(nPOD) = nFields
# nPOD = [-2]
# nPOD = [10] # 100 50 10

# bases = []
# singular_values = []

# for iField in range(nFields):

#     snapshots_matrix = snapshots_data[iField]
#     nrows, ncols = snapshots_matrix.shape

#     if nrows > ncols:
#         SSmatrix = np.dot(snapshots_matrix.T, snapshots_matrix)
#     else:
#         SSmatrix = np.dot(snapshots_matrix, snapshots_matrix.T)
#         print('WARNING - CHECK HOW THE BASIS FUNCTIONS ARE CALCULATED WITH THIS METHOD')

#     print('SSmatrix', SSmatrix.shape)
#     eigvalues, v = np.linalg.eigh(SSmatrix)
#     eigvalues =  eigvalues[::-1]
#     # get rid of small negative eigenvalues (there shouldn't be any as the eigenvalues of a real, symmetric
#     # matrix are non-negative, but sometimes very small negative values do appear)
#     eigvalues[eigvalues<0] = 0
#     s_values = np.sqrt(eigvalues)
#     #print('s values', s_values[0:20]) 

#     singular_values.append(s_values)

#     cumulative_info = np.zeros(len(eigvalues))
#     for j in range(len(eigvalues)):
#         if j==0:
#             cumulative_info[j] = eigvalues[j]
#         else: 
#             cumulative_info[j] = cumulative_info[j-1] + eigvalues[j]

#     cumulative_info = cumulative_info / cumulative_info[-1]
#     nAll = len(eigvalues)

# #if nPOD = -1, use cumulative tolerance
# #if nPOD = -2 use all coefficients (or set nPOD = nTime)
# #if nPOD > 0 use nPOD coefficients as defined by the user

#     if nPOD[iField] == -1:
#         # SVD truncation - percentage of information captured or number 
#         cumulative_tol = nirom_options.compression.cumulative_tol[iField]
#         nPOD_iField = sum(cumulative_info <= cumulative_tol) #tolerance
#         nPOD[iField] = nPOD_iField                
#     elif nPOD[iField] == -2:
#         nPOD_iField = nAll
#         nPOD[iField] = nPOD_iField                
#     else:
#         nPOD_iField = nPOD[iField]

#     print("retaining", nPOD_iField, "basis functions of a possible", len(eigvalues))


#     basis_functions = np.zeros((nx*ny*nz*nDim,nPOD_iField)) # nDim should be nScalar?
#     for j in reversed(range(nAll-nPOD_iField,nAll)):
#         Av = np.dot(snapshots_matrix,v[:,j])
#         basis_functions[:,nAll-j-1] = Av/np.linalg.norm(Av)

#     bases.append(basis_functions)

# write_sing_values(singular_values)

# # get reconstructed snapshots
# reconstruction_data = []
# for iField in range(nFields):

#     basis = bases[iField]

#     snapshots_matrix = snapshots_data[iField]
#     print('snapshots_matrix', snapshots_matrix.shape)

#     reconstruction_on_mesh = np.zeros((nScalar*nTime,nNodes))
#     #reconstruction_on_mesh_from_one_grid = np.zeros((nScalar,nNodes))

#     for iGrid in range(nGrids):

#         #:,iTime*nGrids+iGrid
#         # want solutions in time for a particular grid
#         snapshots_per_grid = np.zeros((nx*ny*nz*nDim,nTime))
#         for iTime in range(nTime):
#             #print('taking snapshots from', iTime*nGrids+iGrid ) 
#             snapshots_per_grid[:,iTime] = snapshots_matrix[:,iTime*nGrids+iGrid]

#         reconstruction =  np.dot( basis, np.dot( basis.T, snapshots_per_grid ) )
#         #print('reconstruction', reconstruction.shape)         
#         #reconstruction_data.append(reconstruction)    

#         #print ('recon shape',reconstruction.shape)
#         reconstruction_grid = reconstruction.reshape(nScalar,nx,ny,nTime)
#         #print ('recon shape just before interpolating back onto mesh',reconstruction.reshape(nScalar,nx,ny,nTime).shape)

#         # plot solution on each grid at 4 time steps
#         #fig, axs = plt.subplots(2, 2, figsize=(15,15))
#         #if iGrid==0:
#         #    levels = np.linspace(0, 4, 5)
#         #elif iGrid==1:
#         #    levels = np.linspace(5, 9, 5)
#         #icount = 0
#         #for col in range(2):
#         #    for row in range(2):
#         #        ax = axs[row, col]
#         #        ax.set_title('time '+str(icount))
#         #        pcm = ax.contourf(reconstruction_grid[0,:,:,icount],levels=levels) 
#         #        fig.colorbar(pcm,ax=ax) 
#         #        icount += 1 
#         #plt.show()


#         block_x_start = get_grid_end_points(grid_origin, grid_width, iGrid) 
#         if iTime==0:
#             print('block_x_start', block_x_start)

#         for iTime in range(nTime):

#             zeros_beyond_grid = 1 # 0 extrapolate solution; 1 gives zeros for nodes outside grid
#             reconstruction_on_mesh_from_one_grid = u2r.interpolate_from_grid_to_mesh(reconstruction_grid[:,:,:,iTime], block_x_start, ddx, x_all, zeros_beyond_grid, nScalar,nx,ny,nz,nNodes,nDim, 1)

#             #print('reconstruction_on_mesh_from_one_grid - about to add solutions',reconstruction_on_mesh_from_one_grid.shape)
#             reconstruction_on_mesh[nScalar*iTime:nScalar*(iTime+1),:] = reconstruction_on_mesh[nScalar*iTime:nScalar*(iTime+1),:] + np.squeeze(reconstruction_on_mesh_from_one_grid)

#     reconstruction_on_mesh[:,duplicated_nodal_values] = 0.5*reconstruction_on_mesh[:,duplicated_nodal_values]        
#     reconstruction_data.append(reconstruction_on_mesh)    



# original_data = []
# #for ifield in range(nFields):
# #    nDoF = nNodes # could be different value per field
# #    original_data.append(np.zeros((nNodes, nDim*nTime)))
# original = np.zeros((nNodes, nDim*nTime))
# for iTime in range(nTime):

#     #print('')
#     #print('time level', iTime)

#     filename = snapshot_data_location + snapshot_file_base + str(offset+iTime) + '.vtu'
#     vtu_data = vtktools.vtu(filename)

#     #original = np.zeros((nNodes, nDim*nTime))

#     for iField in range(nFields):

#         #vtu_data = vtktools.vtu(filename)
#         my_field = vtu_data.GetField(field_names[iField])[:,0:nDim]
#         original[:,iTime*nDim:(iTime+1)*nDim] = my_field
        
#     #print('original.shape',original.shape)
#     #original_data.append(original)

# # make diretory for results
# path_to_reconstructed_results = 'reconstructed_results/'
# if not os.path.isdir(path_to_reconstructed_results):
#     os.mkdir(path_to_reconstructed_results)

# template_vtu = snapshot_data_location + snapshot_file_base + '0.vtu'
# for iTime in range(nTime):
#     for iField in range(nFields):
#         reconstruction_on_mesh = reconstruction_data[iField]
#         #print ('recon shape',reconstruction.shape)

#     # for more than one field, will this work? 
#     #create_vtu_file_timelevel(nNodes, reconstruction_on_mesh[iTime*nScalar:(iTime+1)*nScalar,:], template_vtu, iTime)
#     create_vtu_file(path_to_reconstructed_results, nNodes, reconstruction_on_mesh[iTime*nScalar:(iTime+1)*nScalar,:], template_vtu, original[:,iTime*nDim:(iTime+1)*nDim], iTime)


# print('ddx',ddx.shape,ddx)

# print('grid lengths',   (nx-1)*ddx[0], (ny-1)*ddx[1] ) #, (nz-1)*ddx[2])

# print('Finished.')








# %%

# %%
