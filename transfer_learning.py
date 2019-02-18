import keras
import numpy as np

def copy_layers(source, layerid, trainable=True):
    model = keras.models.clone_model(source)
    for i in layerid:
        model.layers[i].set_weights(source.layers[i].get_weights())
        model.layers[i].trainable = trainable
    return model

def make_trainable(source, layerid=None):
    model = keras.models.clone_model(source)
    model.set_weights(source.get_weights())
    for l in model.layers:
        l.trainable = False
    if not layerid:
        layerid = range(len(model.layers))
    for i in layerid:
        model.layers[i].trainable = True
    return model

def swap_sorbates(data):
    X_train, Y_train, X_test, Y_test = data
    for X in [X_train, X_test]:
        temp_col = X[:, 3].copy()
        X[:, 3] = X[:, 2].copy()
        X[:, 2] = temp_col.copy()
    for Y in [Y_train, Y_test]:
        if Y.shape[1] == 2:
            temp_col = Y[:, 1].copy()
            Y[:, 1] = Y[:, 0].copy()
            Y[:, 0] = temp_col.copy()
        else:
            temp_col = Y[:, 2:].copy()
            Y[:, 2:] = Y[:, :2].copy()
            Y[:, :2] = temp_col.copy()
            
def evaluate_by_temp(model, data, temps=np.linspace(343, 493, 16)):
    X, Y = np.vstack([data[0], data[2]]), np.vstack([data[1], data[3]])
    mse_list = []
    for T in temps:
        mask = X[:, 0] == 1000 / T
        Y_pred = model.predict(X[mask, :])
        mse_y1 = np.mean((Y[mask, 0] - Y_pred[:, 0]) ** 2)
        mse_y2 = np.mean((Y[mask, 0] - Y_pred[:, 0]) ** 2)
        mse_list.append((mse_y1 + mse_y2) / 2)
    return mse_list
    