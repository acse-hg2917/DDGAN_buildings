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
snapshot_data_location = '/Users/gohannaago/Desktop/IC_MSc_ACSE/9_PROJECT/transient_data/transient_data/'
snapshot_file_base = 'Flow_past_buildings_'
nTime = 480
offset = 10
t_step = 3
field_names = ['Velocity', 'VelocityAbsorption']
nFields = len(field_names)
grid_width = [0.48,0.48]
# spacing inside small grid 
nx = int(grid_width[0]*100)
ny = nx
nz = 1
ddx = np.array((0.01,0.01))
nDim = 2


filename = snapshot_data_location + snapshot_file_base + '84.vtu'
representative_vtu = vtktools.vtu(filename)
coordinates_org = representative_vtu.GetLocations() #coordinates of the nodes
# nNodes = coordinates.shape[0] # vtu_data.ugrid.GetNumberOfPoints()
nEl = representative_vtu.ugrid.GetNumberOfCells()
nloc = 3
nScalar = 1

assert nTime <= (500 - offset), "Time range nTime is larger than what comprises full dataset"

nSnapshots = nTime // t_step +1
print("Number of time level Snapshots = ", nSnapshots)
xlength = 3.0  
ylength = 3.0

# get global node numbers
# x_ndgln = np.zeros((nEl*nloc), dtype=int)
# for iEl in range(nEl):
#     n = representative_vtu.GetCellPoints(iEl) + 1
#     x_ndgln[iEl*nloc:(iEl+1)*nloc] = n

#variables for interpolation
zeros_beyond_mesh = 0
#%%
snapshot_ae = np.zeros((nSnapshots,nGrids,nx,ny,3))

for iTime in range(0,nTime,t_step):
    print("Timestep ", offset + iTime, "out of ", 500)
    filename = snapshot_data_location + snapshot_file_base + str(offset+iTime) + '.vtu'
    vtu_data = vtktools.vtu(filename)
    iSnapshot = iTime // t_step
    print("At snapshot ", iSnapshot+1, " out of ", nSnapshots, " snapshots")

    coordinates_org = vtu_data.GetLocations() #coordinates of the nodes
    coordinates = coordinates_org
    nNodes = coordinates.shape[0]# vtu_data.ugrid.GetNumberOfPoints()
    nEl = vtu_data.ugrid.GetNumberOfCells()
    x_all = np.transpose(coordinates[:,0:nDim])
    x_ndgln = np.zeros((nEl*nloc), dtype=int)
    for iEl in range(nEl):
        n = vtu_data.GetCellPoints(iEl) + 1
        x_ndgln[iEl*nloc:(iEl+1)*nloc] = n
    print("nEl for time level ", iTime, " = ", nEl)

    velocity_field = vtu_data.GetField(field_names[0])[:,:nDim] #field name 0 is velocity field
    va_field = vtu_data.GetField(field_names[1])[:,0]
    for iGrid in range(nGrids):
        print("Interpolating Grid ", iGrid, "among ", nGrids, 'Grids')
        #interpolate & append stream velocity
        u_mesh = np.zeros((1,nNodes,1))
        u_mesh[:,:,0] = np.transpose(velocity_field[:,0])
        uvalue_grid = u2r.simple_interpolate_from_mesh_to_grid(u_mesh, x_all, x_ndgln, ddx, grid_origins[iGrid], nx, ny, nz, zeros_beyond_mesh, nEl, nloc, nNodes, nScalar, nDim,1)
        # print('uvalue_grid.shape', uvalue_grid.shape)
        snapshot_ae[iSnapshot,iGrid,:,:,0] = np.rot90(uvalue_grid.reshape((nx,ny)))

        #interpolate & append v velocity
        v_mesh = np.zeros((1,nNodes,1))
        v_mesh[:,:,0] = np.transpose(velocity_field[:,1])
        vvalue_grid = u2r.simple_interpolate_from_mesh_to_grid(v_mesh, x_all, x_ndgln, ddx, grid_origins[iGrid], nx, ny, nz, zeros_beyond_mesh, nEl, nloc, nNodes, nScalar, nDim,1)
        snapshot_ae[iSnapshot,iGrid,:,:,1] = np.rot90(vvalue_grid.reshape((nx,ny)))
        # snapshots_matrix[nx*ny:nx*ny*2,iGrid] = vvalue_grid.reshape(-1)

        #interpolate & append velocity absoprtion
        va_mesh = np.zeros((1,nNodes,1))
        va_mesh[:,:,0] = np.transpose(va_field)
        vavalue_grid = u2r.simple_interpolate_from_mesh_to_grid(va_mesh, x_all, x_ndgln, ddx, grid_origins[iGrid], nx, ny, nz, zeros_beyond_mesh, nEl, nloc, nNodes, nScalar, nDim,1)
        # print(vavalue_grid.reshape(-1)[:30])
        snapshot_ae[iSnapshot,iGrid,:,:,2] = np.rot90(vavalue_grid.reshape((nx,ny)))
        # snapshots_matrix[nx*ny*2:nx*ny*3,iGrid] = vavalue_grid.reshape(-1)    

# Scale snapshot matrix
# min_vel = np.amin(snapshot_ae[:,:,:,:2]) #minimum among u and v velocity
# max_vel = np.amax(snapshot_ae[:,:,:,:2])
min_va = np.amin(snapshot_ae[:,:,:,:,2])
max_va = np.amax(snapshot_ae[:,:,:,:,2])
# vel_scaling = 1/(max_vel-min_vel)
va_scaling = 1/(max_va - min_va) #maximum of buildings information is 100000

#hardcode values to stay consistent with training data
min_vel =  -6.545574188232422
max_vel =  6.623080253601074
vel_scaling =  0.0759379027232465

# print("min_vel = ", min_vel, " , max_vel = ", max_vel, " , vel_scaling = ", vel_scaling)

#scale snapshots for ae
snapshot_ae[:,:,:,:,:2] = vel_scaling*(snapshot_ae[:,:,:,:,:2]-min_vel)
snapshot_ae[:,:,:,:,2] = va_scaling*snapshot_ae[:,:,:,:,2]

print("Scaled minimum = ", np.amin(snapshot_ae[:,:,:,:,2]))
print("Scaled maximum = ", np.amax(snapshot_ae[:,:,:,:,2]))

#save grid origins for the TD case dataset
np.save("reg_snapshots_td.npy", snapshot_ae)
# np.save("reg_grid_origins.npy", grid_origins)

#%%
textfile = open("full_domain_td.txt", 'w')
textfile.write("min_vel = "+str(min_vel) +" , max_vel = "+ str(max_vel)+ " , vel_scaling = "+ str(vel_scaling))
textfile.close()

print("min_vel = ", min_vel, " , max_vel = ", max_vel, " , vel_scaling = ", vel_scaling)
#%%