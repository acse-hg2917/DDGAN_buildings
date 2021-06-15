import numpy as np
import vtk, vtktools


# inputs --------------------------------------
# filenames constructed from this
filebase = 'Flowpast_2d_Re3900_'
# how many files do we want to read in?
nFiles = 3 #we have nFiles+1 time level, with 0 as initial state
# how many physical dimensiona does our problem have?
nDim = 2
# which fields do you want to read in? (these strings MUST match the name in the vtu file. Open the file in paraview to see which fields there are. You could also check this out from this python script - see vtktools.py)
field_name = 'Velocity' # you can read in more than one field

#-------------------------------------------------------------
for iFile in range(nFiles):

    filename = filebase + str(iFile) + '.vtu'

    vtu_data = vtktools.vtu(filename)

    if iFile==0:
        nNodes = vtu_data.ugrid.GetNumberOfPoints() 
        field_data = np.zeros((nNodes*nDim, nFiles))

    # you could use reshape command for this but I wrote it out so you can get a feel for how the data is structured
    # this makes sure you only get 2D data if you have a 2D problem, as, in the vtk data, 3D arrays are used for 2D with zeros in the third dimension - we don't need the zeros

    for iDim in range(nDim):
        ia = iDim*nNodes
        ib = (iDim+1)*nNodes 
        field_data[ia:ib,iFile] = vtu_data.GetField(field_name)[:,iDim]

norm_of_all_data_known = 8.403021824816781 # pre-calculated
norm_of_all_data = np.linalg.norm(field_data)
assert abs(norm_of_all_data - norm_of_all_data_known)< 1e-10, "There has been an error in reading in these files."

# max value of y-component of velocity for time level 2
# if you open 'Flowpast_2d_Re3900_2.vtu' (i.e. the results at time level 2) in Paraview, does the maximum value that you see there agree with the value output to screen below?
print('The maximum value of the y-component of velocity for time level 2 is', np.max(field_data[nNodes:2*nNodes,2]))
print('nNodes = ', nNodes)
# make sure in paraview you are looking at the y-component and not the "magnitude" (speed)

print('Finished.')




