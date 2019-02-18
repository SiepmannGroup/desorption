import matplotlib.pyplot as plt
import matplotlib
import numpy as np

# Plot style
plt.style.use('default')
font = {'family' : 'arial',
        'size'   : 8}

matplotlib.rc('font', **font)
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['axes.titlesize'] = 8
matplotlib.rcParams['axes.labelsize'] = 8


#colors
BG = np.array([245, 245, 245]) / 255
B1 = np.array([27, 117, 187]) / 255
B2 = np.array([23, 154, 228]) / 255
G1 = np.array([8, 147, 71]) / 255
G2 = np.array([71, 209, 93]) / 255
P = np.array([144, 39, 142]) / 255

T = [np.array([245, 41, 39]) / 255,
     np.array([255, 208, 38]) / 255,
     np.array([45, 195, 14]) / 255,
     np.array([29, 179, 204]) / 255][::-1]

B3 = np.array([0, 176, 240]) / 255
Y = np.array([255, 192, 0]) / 255


'''
Set up the figure template
'''
def setup_figure(fig, ax):
    fig.set_dpi(600)
    if isinstance(ax, list) or isinstance(ax, np.ndarray):
        for x in ax:
            x.set_facecolor(BG)
            x.tick_params(labelsize=8)
            x.grid(color="0.9", linewidth=1, zorder=-10)
    else:
        ax.set_facecolor(BG)
        ax.tick_params(labelsize=8)
        ax.grid(color="0.9", linewidth=1, zorder=-10)


# Parity plot of training and test errors; Figure 2b
def parity_plot(train_res, test_res):
    fig = plt.figure(figsize=(1.8, 3.6))
    plt.subplots_adjust(hspace=0)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    setup_figure(fig, [ax1, ax2])
    l1 = ax1.scatter(train_res[:, 0], train_res[:, 1], s=22.5, marker='o', edgecolors="1", linewidths=0.75, c=[B3], zorder=20)
    ax1.errorbar(train_res[:, 0], train_res[:, 1], xerr=train_res[:, 2], ls='', c=[B3], zorder=10)
    l2 = ax1.scatter(test_res[:, 0], test_res[:, 1], s=22.5, marker='D', edgecolors="1", linewidths=0.75, c=[Y], zorder=20)
    ax1.errorbar(test_res[:, 0], test_res[:, 1], xerr=test_res[:, 2], ls='', c=[Y], zorder=10)
    ax1.plot([-1, 2], [-1, 2], ls='--', c="0.5", zorder=20)
    ax1.legend([l1,l2], ['C5 Training','C5 Test'], fontsize=6)
    ax1.set_xlim([-0.05, 1.05])
    ax1.set_ylim([-0.05, 1.05])
    ax1.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax1.set_yticks([0, 0.25, 0.5, 0.75, 1])
    l1 = ax2.scatter(train_res[:, 3], train_res[:, 4], s=22.5, marker='o', edgecolors="1", linewidths=0.75, c=[B3], zorder=20)
    ax2.errorbar(train_res[:, 3], train_res[:, 4], xerr=train_res[:, 5], ls='', c=[B3], zorder=10)
    l2 = ax2.scatter(test_res[:, 3], test_res[:, 4], s=22.5, marker='D', edgecolors="1", linewidths=0.75, c=[Y], zorder=20)
    ax2.errorbar(test_res[:, 3], test_res[:, 4], xerr=test_res[:, 5], ls='', c=[Y], zorder=10)
    ax2.plot([-1, 2], [-1, 2], ls='--', c="0.5", zorder=20)
    ax2.legend([l1,l2], ['W Training','W Test'], fontsize=6)
    ax2.set_xlabel("$y_i$ Simulation")
    ax2.set_ylabel("$\hat{y_i}$ Neural Network Prediction")
    ax2.set_xlim([-0.05, 1.05])
    ax2.set_ylim([-0.05, 1.05])
    ax2.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax2.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax1.set_xticklabels([])
    ax2.yaxis.set_label_coords(-0.2, 1)
    plt.show()

# Plot NN-interpolated adsorption isotherms; Figure 2d
def isotherm_plot(data_true, data_pred):
    fig, axes = plt.subplots(2, 1)
    fig.set_size_inches(1.8, 3.6)
    setup_figure(fig, [axes[0], axes[1]])
    temps = ['343', '393', '443', '493']
    ylabels = ['C5 loading [molec/uc]', 'W loading [molec/uc]']
    markers = ['D', '^', 'v', 'o']
    sizes = np.array([25, 40, 40, 30])*0.75
    dots = [[],[]]
    xlabel = '$-\log_{10}{[v/(\mathrm{\AA^{3}uc^{-1}})]}$'
    plt.subplots_adjust(hspace=0)
    for i in range(4):
        for axid, ax in enumerate(axes):
            ax.plot(data_pred[:, 0], data_pred[:, 3 * i + 1 +  axid], c=T[i])
            ax.errorbar(data_true[:, 5 * i], data_true[:, 5 * i + 1 + 2 * axid], 
                        yerr=data_true[:, 5 * i + 2 + 2 * axid], ls='', c=T[i], zorder=10)
            dots[axid].append(axes[axid].scatter(data_true[:, 5 * i], data_true[:, 5 * i + 1 + 2 * axid], 
                                    s=sizes[i], marker=markers[i], edgecolors="1", linewidths=0.75, c=T[i], zorder=20))
            ax.set_xlim([-16.5, -2.5])
            ax.set_ylim([-0.2, 3])
            ax.set_ylabel(ylabels[axid])
            ax.set_xlabel(xlabel)
            ax.set_xticks([-16, -12, -8,-4])
    axes[1].legend(dots[1], ['343 K','393 K', '443 K', '493 K'], fontsize=6)
    plt.show()
    
# Parity plot of Figure 3c inset    
def parity_plot_pvnet(train_res, test_res):
    data_all = np.vstack([train_res, test_res])
    fig, ax = plt.subplots()
    fig.set_size_inches(1.5, 1.25)
    fig.set_dpi(600)
    plt.setp(ax.spines.values(), linewidth=0.6)
    ax.set_facecolor(BG)
    ax.tick_params(labelsize=6, width=0.6)
    ax.grid(color="0.9", linewidth=0.6, zorder=-10)
    ax.scatter(data_all[:, 0], data_all[:, 1], s=10, marker='o', edgecolors="1", linewidths=0.75, c=P, zorder=20)
    plt.plot([-18,0],[-18,0], c="0.5", zorder=10, ls='--', linewidth=0.75)
    ax.set_xlim([-17.5, -2.5])
    ax.set_ylim([-17.5, -2.5])
    ax.set_xticks([-16, -10, -4])
    ax.set_yticks([-16, -10, -4])
    ax.set_xlabel('Simulation log($v$)')
    ax.set_ylabel('NN Prediction log($v$)')
    plt.show()

# Molar ratio surface in Figure 3a
# Returns the colorbar limit to align with Figure 3b
def plot_molar_ratio(q1, q2, pressure_line):
    T_light = np.array(T) / 2 + 0.5
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(3.5, 3)
    fig.set_dpi(600)
    ax.tick_params(labelsize=8)
    for i, p in enumerate(pressure_line):
        ax.plot([343, 493], [p, p], c=T_light[i], ls='--')
    z = q1 / q2
    cs = ax.imshow(z, extent=[343, 493, 2, -4], aspect=24, norm=matplotlib.colors.LogNorm())
    cbar = fig.colorbar(cs, pad=0.05, shrink=0.91, aspect=15)
    cbar.set_label('$q_\mathrm{C5}(p,T)\ /\ q_\mathrm{W}(p,T)$')
    cbar.ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator())
    cbar.set_ticks(cbar.ax.yaxis.get_major_locator().tick_values(z.min(), z.max()))
    ax.set_xlabel('$T$ [K]')
    ax.set_xticks([343, 393, 443, 493])
    ax.set_ylabel('$\log_{10}({p/\mathrm{kPa}})$')
    plt.show()
    return cbar.get_clim()

# Optimal temperature curve in Figure 3d
def plot_temp_curve(data):
    fig, ax = plt.subplots()
    fig.set_size_inches(3, 2.5)
    fig.set_dpi(600)
    setup_figure(fig, ax)
    ls_list = ['-', '--', '--']
    ps = ['7.92:0.92', '2.67:1.92']
    ylabel = '$T_o$ [K]'
    xlabel = '$\log_{10}{(p/\mathrm{kPa})}$'
    plt.subplots_adjust(wspace=0.25)
    for i in [0, 1]:
        ax.plot(data[:, 3 * i], data[:, 3 * i + 1], c=T[2 * i + 1], ls=ls_list[i])
    ax.set_xlim([-3.5, 0.5])
    ax.set_ylim([320, 440])
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.legend(ps, fontsize=8)
    plt.show()

# Molar ratio isobaric curves in Figure 3b
# Also plots the constraint boundary
def plot_isobaric_ops(data_isobaric, optimal_t_curve, clim):
    fig, ax = plt.subplots()
    fig.set_size_inches(3, 2.5)
    fig.set_dpi(600)
    setup_figure(fig, ax)
    ps = ['$10^{-3}$ kPa', '$10^{-2}$ kPa', '$10^{-1}$ kPa', '$1$ kPa']
    ylabel = '$q_\mathrm{C5}(p,T)\ /\ q_\mathrm{W}(p,T)$'
    markers = ['D', '^', 'v', 'o']
    sizes = [25, 40, 40, 30]
    xlabel = '$T$ [K]'
    plt.subplots_adjust(wspace=0.25)
    for i in range(4):
        ax.plot(data_isobaric[:, 3 * i], data_isobaric[:, 3 * i + 1] / data_isobaric[:, 3 * i + 2], c=T[i], zorder=20)
    ax.set_xlim([343, 493])
    ax.set_ylim(clim)
    ax.set_yscale('log')
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_xticks([343, 393, 443, 493])
    ax.legend(ps, fontsize=8)
    ax.plot(optimal_t_curve[:, 1], optimal_t_curve[:, 2], ls='--', c='0.5')
    plt.show()
    
# MSE vs temperature plot in Figure 4
def plot_by_temperature(data):
    temps = np.linspace(343, 493, 16)
    zorders = [30, 10, 20]
    markers = ['o', 'D', '<']
    names = list(data.keys())
    sizes = [4, 4, 4]
    colors = [T[1], T[2], T[3]]
    lines = []
    fig, axes = plt.subplots(4, 1)
    axes = axes.ravel()
    fig.set_size_inches(2.5, 7)
    fig.set_dpi(600)
    fig.subplots_adjust(wspace=0, hspace=0)
    setup_figure(fig, axes)
    for i, key in enumerate(names):
        if names[i] == 'MFI-C5-W':
            names[i] = 'MFI-W-C5'
        for j in range(data[key].shape[1]):
            lines.append(axes[i].errorbar(temps, data[key][:, j], zorder=zorders[j], marker=markers[j],
                             ms=sizes[j], ls='-', c=colors[j]))
            axes[i].plot([423, 423], [1e-6, 1], ls='--', c="0.4", linewidth=1)
            axes[i].set_yscale('log')
            axes[i].set_ylim([5e-5, 0.2])
            axes[i].set_xticks([343, 383, 423, 453, 493])
            axes[i].set_xlabel('Temperature [K]', fontsize=8)
            axes[i].text(495, 8e-5, names[i],
                         color="0.96",
                         ha='right',
                         va='baseline',
                         weight='bold',
                         #backgroundcolor="0.5"
                         bbox=dict(boxstyle="square", fc=(0.5, 0.5, 0.5), ec=(0.5, 0.5, 0.5, 0))
                        )
    axes[-1].legend(lines[-3:], ['SorbNet', 'dense network', 'new training'], fontsize=8) 
    axes[-1].set_ylabel('Mean Square Error', fontsize=8)
    axes[-1].yaxis.set_label_coords(-0.2, 2)
    plt.show()    
    
