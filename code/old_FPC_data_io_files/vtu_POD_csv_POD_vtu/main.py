import numpy as np
import os, sys
import time

# PART 1 - READING IN SNAPSHOTS AND WRITING POD COEFFICIENTS ----------------------------------
def read_in_snapshots_and_write_out_POD_coeffs(snapshot_data_location, snapshot_file_base, nTime, nDim, field_name, G, cumulative_tol):

    # read in snapshots from vtu files ------------------------------------------------------------
    print('reading in snapshots from vtu files')

    nNodes = get_nNodes_from_vtu(snapshot_data_location, snapshot_file_base )
    snapshots_matrix = np.zeros((nDim*nNodes, nTime))
    velocity = np.zeros((nNodes, nDim))

    for iTime in range(nTime):
        # iTime+1 excludes the initial condition (sometimes zero)
        filename = snapshot_data_location + snapshot_file_base + str(iTime+1) + '.vtu'
        vtu_data = vtktools.vtu(filename)
        velocity = vtu_data.GetField(field_name)[:,0:nDim] # as 2D data appears as 3D data in vtu file
        #snapshots_matrix[:nNodes,iTime] = velocity[:,0]
        #snapshots_matrix[nNodes:2*nNodes,iTime] = velocity[:,1]
        #if nDim==3:
        #    snapshots_matrix[2*nNodes:,iTime] = velocity[:,2]
        snapshots_matrix[:,iTime] = velocity.reshape((nDim*nNodes),order='F')

    # take SVD and output POD coeffs -----------------------------------------------------------
    print('finding the POD coefficients of the snapshots')

    # get basis functions and singular values (sing values probably not needed here)
    s_values, basis_functions = get_POD_functions(snapshots_matrix, nPOD, cumulative_tol, nNodes)

    # calculate POD coefficients
    print('shape of basis_functions and snapshots', basis_functions.shape, snapshots_matrix.shape)
    POD_coeffs = np.dot(np.transpose(basis_functions), snapshots_matrix)
    print( 'shape of POD coeffs', POD_coeffs.shape)

    # output POD coefficients to file
    np.savetxt('POD_coeffs.csv', POD_coeffs , delimiter=',')
    np.savetxt('basis_functions.csv', basis_functions , delimiter=',')

    return



# PART 2 - READING IN PREDICTIONS OF POD COEFFICIENTS AND WRITING RESULTS --------------------
def read_in_ML_predictions_and_write_results(snapshot_data_location, snapshot_file_base):

    # read in, apply inverse SVD and print out results -------------------------------------------
    POD_coeffs_prediction = np.loadtxt('POD_coeffs.csv', delimiter=',') # _prediction
    basis_functions = np.loadtxt('basis_functions.csv', delimiter=',')

    prediction = np.dot(basis_functions, POD_coeffs_prediction)

    # get clean vtu file
    filename = snapshot_data_location + snapshot_file_base + '0.vtu'
    clean_vtu = get_clean_vtu_file(filename) #delete previous timesteps in file

    # write results to vtu
    nNodes = get_nNodes_from_vtu(snapshot_data_location, snapshot_file_base )
    velocity = np.zeros((nNodes,3))

    cwd = os.getcwd()
    if not os.path.isdir('vtu_results'):
        os.mkdir('vtu_results')  
    os.chdir('vtu_results') # will overwrite files in results

    for i in range(prediction.shape[1]):

        new_vtu = vtktools.vtu()
        new_vtu.ugrid.DeepCopy(clean_vtu.ugrid) #make deepcopy of mesh
        new_vtu.filename = 'prediction_' + str(i) + '.vtu'

        #velocity[:,0] = prediction[:nNodes,i] 
        #velocity[:,1] = prediction[nNodes:2*nNodes,i]
        #if nDim==3:
        #    velocity[:,2] = prediction[2*nNodes:,i]
        velocity[:,0:nDim] = prediction[:,i].reshape((nNodes,nDim),order='F') 

        new_vtu.AddField('Velocity_CAE',velocity) #Adds velocity to the vtu file
        new_vtu.Write()

    return
    
    os.chdir(cwd)

#----------------------------------------------------------------

snapshot_data_location = '../data/FPC_Re3900_2D_CG_old/' #might need to change this name
snapshot_file_base = 'fpc_2D_Re3900_CG_'
#snapshot_data_location = '/home/cheaney/Results/nirom_test_fpc_nPOD_20_nSnap_100/snapshots_CG/'
#snapshot_file_base = 'fpc_2D_Re3900_CG_'
#path = '/home/cheaney/Results/nirom_test_fpc_nPOD_20_nSnap_100/snapshots/'
#filename = path + 'Flowpast_2d_Re3900_0.vtu'
nTime = 204 # number of snapshots to read in -- the max value will be 999
nDim = 2 # physical dimension of the field -- can be 2 or 3
field_name = 'Velocity'


#if nPOD = -1, use cumulative tolerance
#if nPOD = -2 use all coefficients (or set nPOD = nTime)
#if nPOD > 0 use nPOD coefficients as defined by the user
nPOD = 17 #number of POD coefficients
cumulative_tol = 0.9999 # ignored unless nPOD = -1 #percentage of snapshot to be captured

apply_POD_to_vtu_data = True
POD_coeffs_to_vtu = False


if apply_POD_to_vtu_data:
    # reads in vtu files, applies POD, writes POD basis functions, POD coeffs and singular values to csv file
    import vtk, vtktools
    from tools_io import get_nNodes_from_vtu, get_POD_functions

    t0 = time.time()
    read_in_snapshots_and_write_out_POD_coeffs(snapshot_data_location, snapshot_file_base, nTime, nDim, field_name, nPOD, cumulative_tol)
    t_get_POD_coeffs = time.time() - t0

    print('Time taken to read in snapshots and write POD coeffs to file', t_get_POD_coeffs)



if POD_coeffs_to_vtu:
    # reads in prediction csv and writes the results to vtu files
    import vtk, vtktools
    from tools_io import get_nNodes_from_vtu, get_clean_vtu_file
    t0 = time.time()
    read_in_ML_predictions_and_write_results(snapshot_data_location, snapshot_file_base)
    t_generate_vtu = time.time() - t0

    print('Time taken to read in POD coeff predictions and write to vtu', t_generate_vtu)

print('Finished')



