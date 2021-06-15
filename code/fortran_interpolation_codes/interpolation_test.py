import u2r
import numpy as np
import sys, os
sys.path.append('/usr/lib/python2.7/dist-packages/')
import vtk, vtktools

#import matplotlib.pyplot as plt
#import matplotlib.ticker as mticker
#from matplotlib.patches import Circle, PathPatch
#import matplotlib.cm as cm

from collections import OrderedDict

# f2py -c unstruc_mesh_2_regular_grid_new.f90 -m u2r
# help(u2r)
print(u2r.simple_interpolate_from_mesh_to_grid.__doc__)

# test for u2r.simple_interpolate_from_mesh_to_grid

# info from vtu file - has DG velocities
i=100
filename = 'LSBU_' +str(i)+ '.vtu' #'/home/caq13/Inhale/data/LSBU_100.vtu'
vtu_data =  vtktools.vtu(filename)
tracer = vtu_data.GetField('Tracer')
coordinates = vtu_data.GetLocations()

print('tracer.shape', tracer.shape)

nNodes = coordinates.shape[0] # vtu_data.ugrid.GetNumberOfPoints()
print('nNodes', nNodes) 
nEl = vtu_data.ugrid.GetNumberOfCells()
print('nEl', nEl, type(nEl)) # 6850

x_all = np.transpose(coordinates[:,0:3])

# hardwire for lazyness
# only checking one time level and x component of velocity here 
nscalar = 1
ndim = 3 # 2D problem (nothing to do with the dimension of the field!)
nTime = 1
nloc = 4 #  

value_mesh = np.zeros((nscalar,nNodes,nTime)) #value_mesh(nscalar,nonods,ntime)
value_mesh[0,:,0] = tracer[:,0] # streamwise velocity 
# interpolate coordinate
#value_mesh[0,:,0] = x_all[2,:] # streamwise velocity 

###################### other tests ####################################
# overwrite velocity field
#value_mesh[0,:,0] = 2
#value_mesh[0,:,0] = coordinates[:,1] #+ coordinates[:,0]
#######################################################################

# rectangular domain so
x0 = min(coordinates[:,0])
y0 = min(coordinates[:,1])
z0 = min(coordinates[:,2])
xN = max(coordinates[:,0])
yN = max(coordinates[:,1])
zN = max(coordinates[:,2])
print('(x0,y0,z0)',x0, y0, z0)
print('(xN,yN,zN)',xN, yN, zN)
block_x_start = np.array(( x0, y0, z0 )) 
#block_x_start = np.array(( 0, 0, 10 )) 

# domain for fpc 2.2 by  0.41 #########################################
#2201  #221  #23 for 0.1 - length 2.2
#411   #42   #5  for 0.1 - length 0.4

# get global node numbers
x_ndgln = np.zeros((nEl*nloc), dtype=int)
for iEl in range(nEl):
    n = vtu_data.GetCellPoints(iEl) + 1
    x_ndgln[iEl*nloc:(iEl+1)*nloc] = n


# set grid size
# change BOTH nx and ddx[0] to change resolution in x, etc (leave z fixed as 2D results)
ddx_size = 3
nx = 128#I chose these values because it's easier for the CNN later.
ny = 128
nz = 32
#nx = 98#I chose these values because it's easier for the CNN later.
#ny = 98
#nz = 16
#nx = 98#I chose these values because it's easier for the CNN later.
#ny = 98
#nz = 32
#nx = 328#I chose these values because it's easier for the CNN later.
#ny = 328
#nz = 102
#nx = 528#I chose these values because it's easier for the CNN later.
#ny = 528
#nz = 202
#nx = 1028#I chose these values because it's easier for the CNN later.
#ny = 1028
#nz = 402
#nx = 328#I chose these values because it's easier for the CNN later.
#ny = 328
#nz = 602
#ddx = np.array((ddx_size, ddx_size, ddx_size)) 
#ddx = np.array((5.625,5.28125,7.81875)) # nx=128=ny, nz=32 domain sizes should be almost equal
ddx = np.zeros((3))
print('ddx.shape',ddx.shape)
ddx[0] = (xN-x0)/(nx-1) 
ddx[1] = (yN-y0)/(ny-1) 
ddx[2] = (zN-z0)/(nz-1) 
#ddx = np.array((2,2,7)) 

print('ddx',ddx.shape,ddx)

print('domain lengths', xN-x0, yN-y0, zN-z0)
print('grid lengths',   (nx-1)*ddx[0], (ny-1)*ddx[1], (nz-1)*ddx[2])
#######################################################################
print('np.max(value_mesh)', np.max(value_mesh))
print('np.min(value_mesh)', np.min(value_mesh))

# interpolate from (unstructured) mesh to (structured) grid
zeros_on_mesh = 0
value_grid = u2r.simple_interpolate_from_mesh_to_grid(value_mesh,x_all,x_ndgln,ddx,block_x_start,nx,ny,nz,zeros_on_mesh,nEl,nloc,nNodes,nscalar,ndim,nTime) 

print('np.max(value_grid)', np.max(value_grid))
print('np.min(value_grid)', np.min(value_grid))
zeros_on_grid = 1

value_remesh = u2r.interpolate_from_grid_to_mesh(value_grid, block_x_start, ddx, x_all, zeros_on_grid, nscalar,nx,ny,nz,nNodes,ndim,nTime)

print('np.max(value_remesh)', np.max(value_remesh))
print('np.min(value_remesh)', np.min(value_remesh))


#directory_Model = '/home/caq13/Inhale/data/'
filename = 'LSBU_' + str(i) + '.vtu'
ug=vtktools.vtu(filename)
ug.AddScalarField('Claire', np.squeeze(value_remesh))
print('np.squeeze(value_remesh).shape', np.squeeze(value_remesh).shape)
ug.Write()


