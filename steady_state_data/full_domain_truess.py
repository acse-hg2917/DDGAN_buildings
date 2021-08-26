#%%
#full domain dataset generation notebook for steady state case
import vtk,vtktools
import u2r
import numpy as np
import matplotlib.pyplot as plt

nGrids = 36 # full domain is represented by 6*6 square subdomains
grid_origins = np.zeros((nGrids,2)) # index of subdomain and x, y coordinates of lower-left corner

for i in range(6): #row
    for j in range(6):#col
        grid_origins[(6*j)+i] = [(0.06)+(0.48*(i)), (3-0.06)-(0.48*(j+1))]
#%%
def plot_square(grid_index):
    x = grid_origins[grid_index,0]
    y = grid_origins[grid_index,1]
    plt.plot([x,x], [y,y+0.48],'g-')
    plt.plot([x,x+0.48], [y,y],'g-')
    plt.plot([x+0.48,x+0.48], [y,y+0.48],'g-')
    plt.plot([x,x+0.48], [y+0.48,y+0.48],'g-')

plt.figure(figsize = (6,6)) #make a square plot
for i in range(len(grid_origins)):
    plot_square(i)
plt.plot(grid_origins[:,0], grid_origins[:,1], 'o')
plt.plot([0,0], [0,3],'k-')
plt.plot([0,3], [0,0],'k-')
plt.plot([3,3], [0,3],'k-')
plt.plot([0,3], [3,3],'k-')
plt.show()
#%%
# import vtu file and make interpolations
snapshot_data_location = '/Users/gohannaago/Desktop/IC_MSc_ACSE/9_PROJECT/steady_state_data/'
snapshot_file_base = 'Flow_past_buildings_Re1_training_'
nTime = 1
nDim = 2
field_names = ['Velocity', 'VelocityAbsorption']
nFields = len(field_names)
grid_width = [0.48,0.48]
# spacing inside small grid 
nx = int(grid_width[0]*100)
ny = nx
nz = 1
ddx = np.array((0.01,0.01))


filename = snapshot_data_location + snapshot_file_base + '400.vtu'
representative_vtu = vtktools.vtu(filename)
coordinates_org = representative_vtu.GetLocations() #coordinates of the nodes
coordinates = coordinates_org
nNodes = coordinates.shape[0] # vtu_data.ugrid.GetNumberOfPoints()
nEl = representative_vtu.ugrid.GetNumberOfCells()
nloc = 3
nScalar = 1

#snapshot matrix in shape (nGrid, 48, 48, 3)
coordinates = coordinates_org
velocity_field = representative_vtu.GetField(field_names[0])[:,:nDim] #field name 0 is velocity field
#call velocity absorption field
va_field = representative_vtu.GetField(field_names[1])[:,0] #field name 0 is velocity field
# print('va_field', va_field[:20])
x_all = np.transpose(coordinates[:,0:nDim])

# get global node numbers
x_ndgln = np.zeros((nEl*nloc), dtype=int)
for iEl in range(nEl):
    n = representative_vtu.GetCellPoints(iEl) + 1
    x_ndgln[iEl*nloc:(iEl+1)*nloc] = n

#variables for interpolation
zeros_beyond_mesh = 0
#%%
snapshot_ae = np.zeros((nGrids,48,48,3))

for iGrid in range(nGrids):
    print("Interpolating Grid ", iGrid, "among ", nGrids, 'Grids')
    #interpolate & append stream velocity
    u_mesh = np.zeros((1,nNodes,1))
    u_mesh[:,:,0] = np.transpose(velocity_field[:,0])
    uvalue_grid = u2r.simple_interpolate_from_mesh_to_grid(u_mesh, x_all, x_ndgln, ddx, grid_origins[iGrid], nx, ny, nz, zeros_beyond_mesh, nEl, nloc, nNodes, nScalar, nDim,1)
    # print('uvalue_grid.shape', uvalue_grid.shape)
    snapshot_ae[iGrid,:,:,0] = np.rot90(uvalue_grid.reshape((nx,ny)))

    #interpolate & append v velocity
    v_mesh = np.zeros((1,nNodes,1))
    v_mesh[:,:,0] = np.transpose(velocity_field[:,1])
    vvalue_grid = u2r.simple_interpolate_from_mesh_to_grid(v_mesh, x_all, x_ndgln, ddx, grid_origins[iGrid], nx, ny, nz, zeros_beyond_mesh, nEl, nloc, nNodes, nScalar, nDim,1)
    snapshot_ae[iGrid,:,:,1] = np.rot90(vvalue_grid.reshape((nx,ny)))
    # snapshots_matrix[nx*ny:nx*ny*2,iGrid] = vvalue_grid.reshape(-1)

    #interpolate & append velocity absoprtion
    va_mesh = np.zeros((1,nNodes,1))
    va_mesh[:,:,0] = np.transpose(va_field)
    vavalue_grid = u2r.simple_interpolate_from_mesh_to_grid(va_mesh, x_all, x_ndgln, ddx, grid_origins[iGrid], nx, ny, nz, zeros_beyond_mesh, nEl, nloc, nNodes, nScalar, nDim,1)
    # print(vavalue_grid.reshape(-1)[:30])
    snapshot_ae[iGrid,:,:,2] = np.rot90(vavalue_grid.reshape((nx,ny)))
    # snapshots_matrix[nx*ny*2:nx*ny*3,iGrid] = vavalue_grid.reshape(-1)    

# Scale snapshot matrix
min_vel = np.amin(snapshot_ae[:,:,:,:2]) #minimum among u and v velocity
max_vel = np.amax(snapshot_ae[:,:,:,:2])
min_va = np.amin(snapshot_ae[:,:,:,2])
max_va = np.amax(snapshot_ae[:,:,:,2])
# vel_scaling = 1/(max_vel-min_vel)
va_scaling = 1/(max_va - min_va) #maximum of buildings information is 100000

# #maybe should try using that of 08153000 rather than averaging
# min_vel =  -6.329753160476685
# max_vel = 6.316451787948608
vel_scaling = 1/(max_vel-min_vel)

print("min_vel = ", min_vel, " , max_vel = ", max_vel, " , vel_scaling = ", vel_scaling)

#scale snapshots for ae
snapshot_ae[:,:,:,:2] = vel_scaling*(snapshot_ae[:,:,:,:2]-min_vel)
snapshot_ae[:,:,:,2] = va_scaling*snapshot_ae[:,:,:,2]

print("Scaled_min = ", np.amin(snapshot_ae[:,:,:,:2]), " scaled max = ", np.amax(snapshot_ae[:,:,:,:2]))

#save grid origins for the TD case dataset
np.save("reg_snapshots_truess_0824.npy", snapshot_ae)
# np.save("reg_grid_origins.npy", grid_origins)

#%%
textfile = open("full_domain_truess_0824.txt", 'w')
textfile.write("min_vel = "+str(min_vel) +" , max_vel = "+ str(max_vel)+ " , vel_scaling = "+ str(vel_scaling))
textfile.close()

print("min_vel = ", min_vel, " , max_vel = ", max_vel, " , vel_scaling = ", vel_scaling)
#%%