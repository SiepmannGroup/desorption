from keras import optimizers
import numpy as np

def calc_mse(data):
    mse_y1 = np.mean((data[:, 0] - data[:, 1]) ** 2)
    mse_y2 = np.mean((data[:, 3] - data[:, 4]) ** 2)
    return (mse_y1 + mse_y2) / 2

'''
Train a model using the given training set data.
Returns the loss history.
'''
def train_model(model, data_train, data_test=None, loss='binary_crossentropy', epochs=500, opt=optimizers.adam, lr=0.002, decay=0.002):
    X_train, Y_train = data_train
    opt_obj = opt(lr=lr, decay=decay)
    model.compile(opt_obj, loss=loss)
    hist = model.fit(X_train, Y_train, epochs=epochs, 
          batch_size=1024, verbose=0, validation_data=data_test, shuffle=True)
    return hist

'''
Evaluate a model using the given training and test data.
Note that evaluation is done on the average of 32 simulation samples on each input
so Y_true will have error bars.
Return array is organized by column as Y1_true, Y1_pred, Y1_err, Y2_true, Y2_pred, Y2_err
'''
def eval_model(model, data_eval):
    X_eval, Y_eval, X_test, Y_test = data_eval
    Y_pred = model.predict(X_eval)
    Y_pred_test = model.predict(X_test)
    train_res = np.vstack([Y_eval[:, 0], Y_pred[:, 0], Y_eval[:, 1],
                             Y_eval[:, 2], Y_pred[:, 1], Y_eval[:, 3]
                            ]).T
    test_res = np.vstack([Y_test[:, 0], Y_pred_test[:, 0], Y_test[:, 1],
                             Y_test[:, 2], Y_pred_test[:, 1], Y_test[:, 3]
                            ]).T
    return train_res, test_res

'''
Interpolate adsorption isotherms using the neural network.
TEMP specifies the temperatures used.
INIT specifies the initial loading value for Y1 (diol) used.
EXT specifies the range of extrapolation.
'''
def interpolate(model, data, temps, init, npoints=200):
    X_all, Y_all = data
    data_true = np.zeros((16, 20))
    data_pred = np.zeros((200, 12))
    ext = 0.5
    for i, T in enumerate(temps):
        mask = np.logical_and(X_all[:, 2] == init, X_all[:, 0] == 1000 / T)
        X_plot = X_all[mask, :]
        Y_plot = Y_all[mask, :]
        Y_plot[:, :2] *= X_plot[:, 2:3]
        Y_plot[:, 2:] *= X_plot[:, 3:4]
        npoints = 200
        x_cur = np.zeros((npoints, 4))
        x_cur[:, 0] += 1000 / T
        x_cur[:, 1] = np.linspace(np.min(X_plot[:, 1]) - ext, np.max(X_plot[:, 1]) + ext, npoints)
        x_cur[:, 2] += X_plot[0, 2]
        x_cur[:, 3] += X_plot[0, 3]
        y_cur = model.predict(x_cur) * X_plot[0, 2:]
        data_true[:, 5*i : 5*(i+1)] = np.hstack([X_plot[:, 1:2], Y_plot])
        data_pred[:, 3*i : 3*(i+1)] = np.hstack([x_cur[:, 1:2], y_cur])
    return data_true, data_pred
    