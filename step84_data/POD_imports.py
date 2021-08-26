#%%
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
# from nirom_dd_orig import plot_starshape

save_bf = False

snapshots_matrix_1 = np.load('snapshots/snapshots_pod_0709.npy')
snapshots_matrix_2 = np.load('snapshots/snapshots_pod.npy')
snapshots_matrix = np.concatenate((snapshots_matrix_1, snapshots_matrix_2), axis = 1)

print("Number of snapshots = ",snapshots_matrix.shape[1])

xlength = 3.0  
ylength = 3.0
grid_width = [0.48,0.48]
# spacing inside small grid 
# %%
snapshots_matrix_vel = snapshots_matrix[:nx*ny*2,:]
snapshots_matrix_b = snapshots_matrix[nx*ny*2:nx*ny*3,:]
# %%
nx = int(grid_width[0]*100)
ny = nx
nz = 1
ddx = np.array((0.01,0.01))

nloc = 3
nScalar = 1 #because I calculate u and v separately
nDim = 2 # dimension of problem (no need to interpolate in the third dimension)

#%%
# cumulative_tol = 0.95
# nPOD = [nTime] # len(nPOD) = nFields
# nPOD = [-2]
nPOD = [-1]
nPOD = [100] # 100 50 10
bases = []
singular_values = []

#snapshots matrix done!
nrows, ncols = snapshots_matrix.shape
if nrows > ncols:
    SSmatrix = np.dot(snapshots_matrix.T, snapshots_matrix)
else:
    SSmatrix = np.dot(snapshots_matrix, snapshots_matrix.T)
    print('WARNING - CHECK HOW THE BASIS FUNCTIONS ARE CALCULATED WITH THIS METHOD')

# print('SSmatrix', SSmatrix.shape)
eigvalues, v = np.linalg.eigh(SSmatrix)
# print('eigenvalues', eigvalues)
# print('normalized eigenvector ([:j]) v = ', v)
eigvalues =  eigvalues[::-1]
# get rid of small negative eigenvalues (there shouldn't be any as the eigenvalues of a real, symmetric
# matrix are non-negative, but sometimes very small negative values do appear)
eigvalues[eigvalues<0] = 0
s_values = np.sqrt(eigvalues)
# print('s values', s_values[0:20]) 

singular_values.append(s_values)

cumulative_info = np.zeros(len(eigvalues))
for j in range(len(eigvalues)):
    if j==0:
        cumulative_info[j] = eigvalues[j]
    else: 
        cumulative_info[j] = cumulative_info[j-1] + eigvalues[j]

cumulative_info = cumulative_info / cumulative_info[-1]
nAll = len(eigvalues)

# Apply POD
#if nPOD = -1, use cumulative tolerance
#if nPOD = -2 use all coefficients (or set nPOD = nTime)
#if nPOD > 0 use nPOD coefficients as defined by the user

if nPOD[0] == -1:  
    # SVD truncation - percentage of information captured or number 
    # cumulative_tol = nirom_options.compression.cumulative_tol[iField]
    nPOD_iField = sum(cumulative_info <= cumulative_tol) #tolerance
    nPOD[0] = nPOD_iField
elif nPOD[0] == -2:
    nPOD_iField = nAll
    nPOD[0] = nPOD_iField
else:
    nPOD_iField = nPOD[0]

print("retaining", nPOD_iField, "basis functions of a possible", len(eigvalues))


basis_functions = np.zeros((nx*ny*nz*(nDim+1),nPOD_iField)) # nDim should be nScalar?
for j in reversed(range(nAll-nPOD_iField,nAll)):
    Av = np.dot(snapshots_matrix,v[:,j])
    basis_functions[:,nAll-j-1] = Av/np.linalg.norm(Av)
# %%
plt.plot(range(len(singular_values[0])), np.log(singular_values[0]), 'k.-')
plt.xlabel("Log of singular values")
plt.ylabel("Index of basis function")
plt.grid(':')
# plt.savefig("singular_plot.png")
plt.show()
# %%
plt.plot(range(len(cumulative_info)), cumulative_info, 'g.-')
plt.xlabel("Log of singular values")
plt.ylabel("Index of basis function")
plt.grid(':')
# plt.savefig("cumulative_info.png")
plt.show()
#%%
"""
Compress and reconstruct the snapshots and graph the MSE reconstruction error
R- basis_functions
alpha - pod coeff
phi - reconstruction
"""
alphas = np.zeros((snapshots_matrix.shape[1], nPOD_iField)) #array that stores compressions
for i in range(snapshots_matrix.shape[1]):
    alphas[i] = basis_functions.T @ snapshots_matrix[:,i]

#reconstruct
reconstructions = np.zeros(snapshots_matrix.shape)
mse_errors = []
for i in range(snapshots_matrix.shape[1]):
    reconstructions[:,i] = basis_functions @ alphas[i]
    mse_errors.append(mean_squared_error(snapshots_matrix[:,i], reconstructions[:,i]))
#%%
# with 278 basis functions, np.mean(mse_errors) = 0.0024700645895766376

#save basis functions
if save_bf:
    np.save("basis_functions.npy", basis_functions)
    
# %%
"""
Import the starshape grid and compress that information
"""
starshape = np.load("starshape_pod.npy")

alphas_s = np.zeros((starshape.shape[1], nPOD_iField)) #array that stores compressions
for i in range(starshape.shape[1]):
    alphas_s[i] = basis_functions.T @ starshape[:,i]

# reconstruct
reconstructions_s = np.zeros(starshape.shape)
mse_errors_s = []
for i in range(starshape.shape[1]):
    reconstructions_s[:,i] = basis_functions @ alphas_s[i]
    mse_errors_s.append(mean_squared_error(starshape[:,i], reconstructions_s[:,i]))

# Export the compressions
np.save("ss_pod_coeffs.npy",alphas_s)

# %%
def create_basis_func(snapshot_matrix, nPOD):
    nrows, ncols = snapshot_matrix.shape
    # if nrows > ncols:
    #     SSmatrix = np.dot(snapshot_matrix.T, snapshot_matrix)
    # else:
    #     SSmatrix = np.dot(snapshot_matrix, snapshot_matrix.T)
    #     print('WARNING - CHECK HOW THE BASIS FUNCTIONS ARE CALCULATED WITH THIS METHOD')

    SSmatrix = np.dot(snapshot_matrix.T, snapshot_matrix)

    # print('SSmatrix', SSmatrix.shape)
    eigvalues, v = np.linalg.eigh(SSmatrix)
    # print('eigenvalues', eigvalues)
    # print('normalized eigenvector ([:j]) v = ', v)
    eigvalues =  eigvalues[::-1]
    # get rid of small negative eigenvalues (there shouldn't be any as the eigenvalues of a real, symmetric
    # matrix are non-negative, but sometimes very small negative values do appear)
    eigvalues[eigvalues<0] = 0
    # s_values = np.sqrt(eigvalues)
    # singular_values.append(s_values)

    cumulative_info = np.zeros(len(eigvalues))
    for j in range(len(eigvalues)):
        if j==0:
            cumulative_info[j] = eigvalues[j]
        else: 
            cumulative_info[j] = cumulative_info[j-1] + eigvalues[j]

    cumulative_info = cumulative_info / cumulative_info[-1]
    nAll = len(eigvalues)

    # Apply POD
    #if nPOD = -1, use cumulative tolerance
    #if nPOD = -2 use all coefficients (or set nPOD = nTime)
    #if nPOD > 0 use nPOD coefficients as defined by the user

    # if nPOD[0] == -1:  
    #     # SVD truncation - percentage of information captured or number 
    #     # cumulative_tol = nirom_options.compression.cumulative_tol[iField]
    #     nPOD_iField = sum(cumulative_info <= cumulative_tol) #tolerance
    #     nPOD[0] = nPOD_iField
    # elif nPOD[0] == -2:
    #     nPOD_iField = nAll
    #     nPOD[0] = nPOD_iField
    # else:
    #     nPOD_iField = nPOD[0]

    print("retaining", nPOD, "basis functions of a possible", len(eigvalues))

    basis = np.zeros((snapshot_matrix.shape[0],nPOD)) # nDim should be nScalar?
    for j in reversed(range(nAll-nPOD,nAll)):
        Av = np.dot(snapshot_matrix,v[:,j])
        basis[:,nAll-j-1] = Av/np.linalg.norm(Av)

    return basis
#%%    
def extract_POD(basis_func, starshape_sm):
    nPOD = basis_func.shape[1]
    alpha_s = np.zeros((starshape_sm.shape[1], nPOD)) #array that stores compressions
    for i in range(starshape.shape[1]):
        alpha_s[i] = basis_func.T @ starshape_sm[:,i]

    return alpha_s
# %%
basis_vel = create_basis_func(snapshots_matrix_vel, 20)
basis_buildings = create_basis_func(snapshots_matrix_b, 10)
starshape_vel = starshape[:nx*ny*2,:]
starshape_b = starshape[nx*ny*2:nx*ny*3,:]
# %%
ss_pod_vel = extract_POD(basis_vel, starshape_vel)
ss_pod_buildings = extract_POD(basis_buildings, starshape_b)

np.save("ss_pod_vel_3000.npy", ss_pod_vel)
np.save("ss_pod_b_3000.npy", ss_pod_buildings)

np.save("basis_b_3000.npy", basis_buildings)
np.save("basis_vel_3000.npy", basis_vel)
