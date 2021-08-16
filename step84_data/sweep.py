import numpy as np
import matplotlib.pyplot as plt

def patch_together_snapshots(reg_snapshot_ae):
    rows = []
    for i in range(6):
        row = np.concatenate(reg_snapshot_ae[i*6:(i+1)*6,:,:,:2], axis = 1)
        rows.append(row)
    full_recon = np.concatenate(rows, axis = 0)

    return full_recon

def plot_pred(grid_v, decoder_v_40, field_index=0, title=None):
    assert field_index == 0 or field_index == 1, "Specify field_index as 0 (u) or 1 (v)"
    reg_snapshot = decoder_v_40(grid_v)
    full_recon_array = patch_together_snapshots(reg_snapshot)
    
    plt.figure(figsize=(6,6))
    plt.imshow(full_recon_array[:,:,field_index])
    if title is not None:
        plt.title(title)
    plt.show()
    plt.close()

def pred_central_grid_ss(grid_index, grid_v, grid_b, predictor):
    inputs = np.concatenate((grid_v[grid_index-6], grid_v[grid_index-1], grid_v[grid_index+1], grid_v[grid_index+6], np.zeros(40),grid_b[grid_index])).reshape(1,220)
    prediction = predictor.predict(inputs).numpy()

    return prediction[0] #outputs are shape of (nSamples, 40)

def single_sweep(grid,grid_b,predictor):
    for i in range(1,5): #rows
        for j in range(1,5):
            I = (i*6)+j
            grid[I] = pred_central_grid_ss(I, grid, grid_b, predictor)
    
    return grid

def sweep_full_domain(predictor, grid_v, grid_b, decoder_v_40, iterations = 3):
    grid_prev = grid_v.numpy() #shape (36,40)
    for row in range(1,5):
        grid_prev[row*6+1:row*6+5] = 0 #leave the outer domain, zeros for prediction area

    it = 0
    while it < iterations:
        print("Predicting Full Domain - sweep number: ", it+1, " out of ", iterations)
        grid_next = single_sweep(grid_prev, grid_b, predictor)
        plot_pred(grid_next, decoder_v_40, title = 'Prediction plot of sweep {}'.format(it))
        grid_prev = grid_next
        it += 1

    return grid_next