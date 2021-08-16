import os,sys
import vtk,vtktools
import numpy as np
import u2r


class Mesh_Information():
    def __init__(self):
        self.nNodes = 0
        self.nEl = 0
        self.nloc = 0
        self.field_names = []
        self.subtract_mean = False 
        self.nDim = 0

class Grid_Information():
    def __init__(self):
        self.nx = 0
        self.ny = 0
        self.nz = 1
        self.nGrids = 0
        self.ddx = []
        self.grid_origin = [] 
        self.grid_width = [] 

def set_mesh_info(nNodes, nEl, nloc, nDim, nFields, field_names):
    mesh_info = Mesh_Information()
    mesh_info.nNodes = nNodes
    mesh_info.nEl = nEl
    mesh_info.nloc = nloc
    mesh_info.nDim = nDim
    mesh_info.nFields = nFields
    mesh_info.field_names = field_names

    return mesh_info

def set_grid_info(nx, ny, nz, nGrids, ddx, grid_origin, grid_width):
    grid_info = Grid_Information()
    grid_info.nx = nx
    grid_info.ny = ny
    grid_info.nz = nz
    grid_info.nGrids = nGrids
    grid_info.ddx = ddx
    grid_info.grid_origin = grid_origin
    grid_info.grid_width = grid_width

    return grid_info

def get_mesh_info(mesh_info):
    return mesh_info.nNodes, mesh_info.nEl, mesh_info.nloc, mesh_info.nDim, mesh_info.nFields, mesh_info.field_names

def get_grid_info(grid_info):
    return grid_info.nx, grid_info.ny, grid_info.nz, grid_info.nGrids, grid_info.ddx, grid_info.grid_origin, grid_info.grid_width


def get_global_node_numbers(nEl, nloc, represnetative_vtu):
    x_ndgln = np.zeros((nEl*nloc), dtype=int)
    for iEl in range(nEl):
        n = represnetative_vtu.GetCellPoints(iEl) + 1
        x_ndgln[iEl*nloc:(iEl+1)*nloc] = n
    return x_ndgln

def get_block_origin(grid_origin, grid_width, iGrid): 
    return np.array(( grid_origin[0]+iGrid*grid_width[0], grid_origin[1] +iGrid*grid_width[1]))

def find_node_duplications_from_overlapping_grids(representative_vtu, mesh_info, grid_info, x_all, x_ndgln):

    nNodes, nEl, nloc, nDim, nFields, field_names = get_mesh_info(mesh_info)
    nx, ny, nz, nGrids, ddx, grid_origin, grid_width  = get_grid_info(grid_info)
  
    # obtain a field of ones the size of the first field
    # assumes all fields are defined over the same number of nodes
    my_field = representative_vtu.GetField(field_names[0])[:,0]
    my_field = 1
    nScalar_test = 1

    # for one timestep and for one field, map the solution (field of value one) from the mesh to the grid and back to the mesh again
    nTime = 1
    value_mesh = np.zeros((nScalar_test,nNodes,nTime))
    value_mesh[:,:,0] = np.transpose(my_field)
    superposed_grids = np.zeros((nNodes))
    for iGrid in range(nGrids):
        block_x_start = get_block_origin(grid_origin, grid_width, iGrid) 

        zeros_beyond_mesh = 0
        value_grid = u2r.simple_interpolate_from_mesh_to_grid(value_mesh,x_all,x_ndgln,ddx,block_x_start,nx,ny,nz,zeros_beyond_mesh, nEl,nloc,nNodes,nScalar_test,nDim,1)

        zeros_beyond_grid = 1
        value_back_on_mesh = u2r.interpolate_from_grid_to_mesh(value_grid, block_x_start, ddx, x_all, zeros_beyond_grid, nScalar_test,nx,ny,nz,nNodes,nDim, 1)

        superposed_grids = superposed_grids + np.rint(np.squeeze(value_back_on_mesh))

    # superpose the solutions on the mesh and detect how many 2s are present, indicating a node at which the solution is duplicated
    superposed_grids = np.array(superposed_grids, dtype='int') 
    duplicated_nodal_values = []
    for iNode in range(nNodes):
        if superposed_grids[iNode] == 0:
            print('zero:', iNode)
        elif superposed_grids[iNode] == 2:
            print('two:', iNode)
            duplicated_nodal_values.append(iNode)
        elif superposed_grids[iNode] != 1:
            print('unknown:', iNode, superposed_grids[iNode])
            print('error - currently can handle a node being on only two grids, not more') 
            sys.exit()

    return duplicated_nodal_values

def read_in_snapshots_interpolate_to_grids(snapshot_data_location, snapshot_file_base, mesh_info, grid_info, nTime, offset, nScalar, x_all, x_ndgln):

    nNodes, nEl, nloc, nDim, nFields, field_names = get_mesh_info(mesh_info)
    nx, ny, nz, nGrids, ddx, grid_origin, grid_width  = get_grid_info(grid_info)

    snapshots_data = []
    for iField in range(nFields):
        nDoF = nNodes # could be different value per field
        snapshots_data.append(np.zeros((nx*ny*nz*nDim, nGrids*nTime)))


    for iTime in range(nTime):

        #print('')
        #print('time level', iTime)

        filename = snapshot_data_location + snapshot_file_base + str(offset+iTime) + '.vtu'
        vtu_data = vtktools.vtu(filename)

        for iField in range(nFields):
    
            my_field = vtu_data.GetField(field_names[iField])[:,0:nDim]

            for iGrid in range(nGrids):

                block_x_start = get_block_origin(grid_origin, grid_width, iGrid) 
                if iTime==0:
                    print('block_x_start', block_x_start)
  
                # shuld this be up with other grid settings
                #ddx = np.array((0.01,0.01))  #np.array((0.002,0.002)) 

                value_mesh = np.zeros((nScalar,nNodes,nTime)) # nTime
                value_mesh[:,:,iTime] = np.transpose(my_field)

                # interpolate field onto structured mesh
                # feed in one result at t time (no need to store in value_mesh in this case)
                zeros_beyond_mesh = 0 # 0 extrapolate solution (for the cylinder in fpc); 1 gives zeros for nodes outside mesh
                value_grid = u2r.simple_interpolate_from_mesh_to_grid(value_mesh[:,:,iTime],x_all,x_ndgln,ddx,block_x_start,nx,ny,nz,zeros_beyond_mesh,nEl,nloc,nNodes,nScalar,nDim,1) 

                snapshots_data[iField][:,iTime*nGrids+iGrid] = value_grid.reshape(-1)

    return snapshots_data

def write_sing_values(singular_values, field_names):
    
    for iField in range(len(field_names)):
        f= open('singular_values_' + field_names[iField] + '.dat',"w+")
        f.write('# index, s_values, normalised s_values, cumulative energy \n' )
        s_values = singular_values[iField] 
        total = 0.0
        for k in range(len(s_values)):
            total = total + s_values[k]*s_values[k]

        running_total = 0.0
        for i in range(len(s_values)):
            running_total = running_total + s_values[i]*s_values[i]
            f.write ('%d %g %g %18.10g \n' % (i, s_values[i], s_values[i]/s_values[0], running_total/total) )
        f.close()
    return

def get_POD_bases(mesh_info, grid_info, snapshots_data, nPOD):
    
    #nNodes, nEl, nloc, nDim, nFields, field_names = get_mesh_info(mesh_info)
    nFields = mesh_info.nFields
    nDim = mesh_info.nDim
    field_names = mesh_info.field_names

    nx = grid_info.nx; ny = grid_info.ny; nz = grid_info.nz


    bases = []
    singular_values = []

    for iField in range(len(field_names)):

        snapshots_matrix = snapshots_data[iField]
        nrows, ncols = snapshots_matrix.shape

        if nrows > ncols:
            SSmatrix = np.dot(snapshots_matrix.T, snapshots_matrix)
        else:
            SSmatrix = np.dot(snapshots_matrix, snapshots_matrix.T)
            print('WARNING - CHECK HOW THE BASIS FUNCTIONS ARE CALCULATED WITH THIS METHOD')

        print('SSmatrix', SSmatrix.shape)
        #print('SSmatrix', SSmatrix)
        eigvalues, v = np.linalg.eigh(SSmatrix)
        eigvalues =  eigvalues[::-1]
        # get rid of small negative eigenvalues
        eigvalues[eigvalues<0] = 0
        s_values = np.sqrt(eigvalues)

        singular_values.append(s_values)

        cumulative_info = np.zeros(len(eigvalues))
        for j in range(len(eigvalues)):
            if j==0:
                cumulative_info[j] = eigvalues[j]
            else: 
                cumulative_info[j] = cumulative_info[j-1] + eigvalues[j]

        cumulative_info = cumulative_info / cumulative_info[-1]
        nAll = len(eigvalues)

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

        print("retaining", nPOD_iField, "basis functions of a possible", len(eigvalues))


        basis_functions = np.zeros((nx*ny*nz*nDim,nPOD_iField))
        for j in reversed(range(nAll-nPOD_iField,nAll)):
            Av = np.dot(snapshots_matrix,v[:,j])
            basis_functions[:,nAll-j-1] = Av/np.linalg.norm(Av)

        bases.append(basis_functions)

    write_sing_values(singular_values, field_names)

    return bases


def reconstruct_data_on_mesh(snapshots_data, mesh_info, grid_info, bases, nScalar, nTime, x_all, duplicated_nodal_values):

    nNodes, nEl, nloc, nDim, nFields, field_names = get_mesh_info(mesh_info)
    nx, ny, nz, nGrids, ddx, grid_origin, grid_width  = get_grid_info(grid_info)

    reconstructed_data = []
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


            block_x_start = get_block_origin(grid_origin, grid_width, iGrid) 
            if iTime==0:
                print('block_x_start', block_x_start)
    
            for iTime in range(nTime):

                zeros_beyond_grid = 1 # 0 extrapolate solution; 1 gives zeros for nodes outside grid
                reconstruction_on_mesh_from_one_grid = u2r.interpolate_from_grid_to_mesh(reconstruction_grid[:,:,:,iTime], block_x_start, ddx, x_all, zeros_beyond_grid, nScalar,nx,ny,nz,nNodes,nDim, 1)

                #print('reconstruction_on_mesh_from_one_grid - about to add solutions',reconstruction_on_mesh_from_one_grid.shape)
                reconstruction_on_mesh[nScalar*iTime:nScalar*(iTime+1),:] = reconstruction_on_mesh[nScalar*iTime:nScalar*(iTime+1),:] + np.squeeze(reconstruction_on_mesh_from_one_grid)

        reconstruction_on_mesh[:,duplicated_nodal_values] = 0.5*reconstruction_on_mesh[:,duplicated_nodal_values]        
        reconstructed_data.append(reconstruction_on_mesh)    
    return  reconstructed_data


def get_original_data_from_vtu_files(snapshot_data_location, snapshot_file_base, offset, mesh_info, nTime):

    nNodes, nEl, nloc, nDim, nFields, field_names = get_mesh_info(mesh_info)

    original_data = []
    #    nDoF = nNodes # could be different value per field
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
        
        original_data.append(original)

    return original_data

