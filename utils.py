import numpy as np
import pandas as pd


'''
Loads data from the diol-solvent dataset and 
create input and output variables.
'''

def make_data(path, task, train=True, test_temps=[343, 393, 443, 493]):
    if train:
        data = pd.read_csv(path + task + '-all.csv', index_col=0)
    else:
        data = pd.read_csv(path + task + '-stat.csv', index_col=0)

    X = data.values[:, :4]
    Y = data.values[:, 4:]

    # Split train/test
    mask_test = np.any(X[:, 0:1] == test_temps, axis=1)

    # Variable transformation
    X[:, 0] = 1000 / X[:, 0]
    X[:, 1] = -np.log10(X[:, 1])
    if train:
        Y /= X[:, 2:4]
    else:
        Y[:, :2] /= X[:, 2:3]
        Y[:, 2:] /= X[:, 3:4]

    X_train = X[np.logical_not(mask_test), :]
    X_test = X[mask_test, :]
    Y_train = Y[np.logical_not(mask_test), :]
    Y_test = Y[mask_test, :]
    return X_train, Y_train, X_test, Y_test

def convert_pressure(data):
    X, Y = data
    if Y.shape[1] == 2:
        nvap = np.sum(X[:, 2:] * (1 - Y), axis=1)
    else:
        nvap = np.sum(X[:, 2:] * (1 - Y[:, ::2]), axis=1)
    X_npt = X[nvap > 0.1, :]
    Y_npt = Y[nvap > 0.1, :]
    V_npt = X[nvap > 0.1, 1:2]
   
    nvap = nvap[nvap > 0.1]
    logp = np.log10(nvap) + np.log10(13806.5) + np.log10(1000 / X_npt[:, 0]) + X_npt[:, 1]
    X_npt = np.vstack([X_npt[:, 0], logp, X_npt[:, 2], X_npt[:, 3]]).T
    return X_npt, Y_npt, V_npt

'''
Stores the Keras history losses as separate numpy array columns.
'''
def list_history(history):
    loss_list = [s for s in history.history.keys() if 'loss' in s and 'val' not in s]
    val_loss_list = [s for s in history.history.keys() if 'loss' in s and 'val' in s]
    
    if len(loss_list) == 0:
        print('Loss is missing in history')
        return 
    
    ## As loss always exists
    epochs = range(1,len(history.history[loss_list[0]]) + 1)
    # Check whether ther is validation data
    if val_loss_list: 
        return np.array([history.history[loss_list[0]]]).T.ravel(), np.array([history.history[val_loss_list[0]]]).T.ravel()
    else:
        return np.array([history.history[loss_list[0]]]).T.ravel()