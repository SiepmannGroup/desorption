import numpy as np

def evaluate_pvnet(model, data_eval, data_test):
    X_eval, V_eval = data_eval
    X_test, V_test = data_eval
    V_pred = model.predict(X_eval)
    V_pred_test = model.predict(X_test)
    train_res = np.hstack([V_eval, V_pred])
    test_res = np.hstack([V_test, V_pred_test])
    return train_res, test_res

# Returns a function which predicts adsorption loading from pressure
# using p-v network.
def get_npt_func(loading_model, pv_model):
    def npt_func(X):
        V = pv_model.predict(X)
        X[:, 1:2] = V
        Y = loading_model.predict(X)
        return Y
    return npt_func

def isobaric_curve(npt_func, lg_pressures, temp_range, inits, npoints=200):
    t_min, t_max = temp_range
    data_pred = np.zeros((200, 12))
    for i, p_input in enumerate(lg_pressures):
        x_cur = np.zeros((npoints, 4))
        x_cur[:, 0] = np.linspace(1000 / t_max, 1000 / t_min, npoints)
        x_cur[:, 1] = p_input                
        x_cur[:, 2] += inits[0]
        x_cur[:, 3] += inits[1]
        y_cur = npt_func(x_cur) * inits
        data_pred[:, i * 3] = 1000 / x_cur[:, 0]
        data_pred[:, i * 3 + 1: i * 3 + 3] = y_cur
    return data_pred

def loading_image(npt_func, logp_range, temp_range, inits, npoints=200):
    t_min, t_max = temp_range
    p_min, p_max = logp_range
    ts = np.linspace(t_min, t_max, npoints)
    ps = np.linspace(p_min, p_max, npoints)
    ts, ps = np.meshgrid(ts, ps)
    x_cur = np.zeros((ts.shape[0] * ts.shape[1], 4))
    x_cur[:, 0] = 1000 / ts.ravel()
    x_cur[:, 1] = ps.ravel()
    x_cur[:, 2] += inits[0]
    x_cur[:, 3] += inits[1]
    y_cur = npt_func(x_cur) * inits
    q1 = y_cur[:, 0].reshape((npoints, npoints))
    q2 = y_cur[:, 1].reshape((npoints, npoints))
    return q1, q2

# Returns the temperature and adsorbed molar ratio for optimal separation. Uses log pressure.
#    Also applies the constraint of product yield no less than MIN_YIELD.
def find_optimal_t(npt_func, t_min, t_max, p_tot, inits, min_yield=0.99, tspace=0.2):
    ts = np.arange(t_min, t_max + tspace, tspace)
    x_cur = np.zeros((ts.shape[0], 4))
    x_cur[:, 0] = 1000 / ts
    x_cur[:, 1] = p_tot
    x_cur[:, 2] += inits[0]
    x_cur[:, 3] += inits[1]
    y_cur = npt_func(x_cur)
    y_ratio = y_cur[:, 0] / y_cur[:, 1] * inits[0] / inits[1]
    y_ratio = y_ratio[y_cur[:, 0] > min_yield]
    id_max = np.argmax(y_ratio)
    return ts[id_max], y_ratio[id_max]

def optimal_t_curve(npt_func, logp_range, inits_lst, t_range=(293, 543), npoints=100):
    t_min, t_max = t_range
    p_min, p_max = logp_range
    nptot = 100
    data_pred = np.zeros((nptot, 3 * len(inits_lst)))
    for j, inits in enumerate(inits_lst):
        ps = np.linspace(p_min, p_max, nptot)
        t_ops = np.zeros((nptot))
        s_ops = np.zeros((nptot))
        for i in range(nptot):
            t_ops[i], s_ops[i] = find_optimal_t(npt_func, t_min, t_max, ps[i], inits)
        data_pred[:, 3*j] = ps
        data_pred[:, 3*j + 1] = t_ops
        data_pred[:, 3*j + 2] = s_ops
    return data_pred

