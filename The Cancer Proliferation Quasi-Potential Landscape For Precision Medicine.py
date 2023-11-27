# -*- coding = utf-8 -*-
# @Time : 2023/9/27 17:30
# @Author : Yourui Han
# @File : The Cancer Proliferation Quasi-Potential Landscape For Precision Medicine.py
# @Software : PyCharm


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict, Counter
import matplotlib.gridspec as gridspec
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy import optimize
import xlwt
from scipy.interpolate import griddata


# 利用高斯距离法计算临近点的权重
# X,Y 模板大小，c 中心点的位置， sigma 影响半径
def gaussion_neighborhood(X, Y, c, sigma):
    xx, yy = np.meshgrid(np.arange(X), np.arange(Y))
    d = 2 * sigma * sigma
    ax = np.exp(-np.power(xx - xx.T[c], 2) / d)
    ay = np.exp(-np.power(yy - yy.T[c], 2) / d)
    return (ax * ay).T


# 利用bubble距离法计算临近点的权重
# X,Y 模板大小，c 中心点的位置， sigma 影响半径
def bubble_neighborhood(X, Y, c, sigma):
    neigx = np.arange(X)
    neigy = np.arange(Y)

    ax = np.logical_and(neigx > c[0] - sigma,
                        neigx < c[0] + sigma)
    ay = np.logical_and(neigy > c[1] - sigma,
                        neigy < c[1] + sigma)
    return np.outer(ax, ay) * 1.


# 计算学习率
def get_learning_rate(lr, t, max_steps):
    return lr / (1 + t / (max_steps / 2))


# 计算欧式距离
def euclidean_distance(x, w):
    dis = np.expand_dims(x, axis=(0, 1)) - w
    return np.linalg.norm(dis, axis=-1)


# 特征标准化 (x-mu)/std
def feature_normalization(data):
    mu = np.mean(data, axis=0, keepdims=True)
    sigma = np.std(data, axis=0, keepdims=True)
    return (data - mu) / sigma


# 获取激活节点的位置
def get_winner_index(x, w, dis_fun=euclidean_distance):
    # 计算输入样本和各个节点的距离
    dis = dis_fun(x, w)

    # 找到距离最小的位置
    index = np.where(dis == np.min(dis))
    return (index[0][0], index[1][0])


def weights_PCA(X, Y, data):
    N, D = np.shape(data)
    weights = np.zeros([X, Y, D])

    pc_length, pc = np.linalg.eig(np.cov(np.transpose(data)))
    pc_order = np.argsort(-pc_length)
    for i, c1 in enumerate(np.linspace(-1, 1, X)):
        for j, c2 in enumerate(np.linspace(-1, 1, Y)):
            weights[i, j] = c1 * pc[pc_order[0]] + c2 * pc[pc_order[1]]
    return weights


# 计算量化误差
def get_quantization_error(datas, weights):
    w_x, w_y = zip(*[get_winner_index(d, weights) for d in datas])
    error = datas - weights[w_x, w_y]
    error = np.linalg.norm(error, axis=-1)
    return np.mean(error)


def train_SOM(X,
              Y,
              N_epoch,
              datas,
              init_lr=0.5,
              sigma=0.5,
              dis_fun=euclidean_distance,
              neighborhood_fun=gaussion_neighborhood,
              init_weight_fun=None,
              seed=20):
    # 获取输入特征的维度
    N, D = np.shape(datas)

    # 训练的步数
    N_steps = N_epoch * N

    # 对权重进行初始化
    rng = np.random.RandomState(seed)
    if init_weight_fun is None:
        weights = rng.rand(X, Y, D) * 2 - 1
        weights /= np.linalg.norm(weights, axis=-1, keepdims=True)
    else:
        weights = init_weight_fun(X, Y, datas)

    for n_epoch in range(N_epoch):
        print("Epoch %d" % (n_epoch + 1))
        # 打乱次序
        index = rng.permutation(np.arange(N))
        for n_step, _id in enumerate(index):
            # 取一个样本
            x = datas[_id]

            # 计算learning rate(eta)
            t = N * n_epoch + n_step
            eta = get_learning_rate(init_lr, t, N_steps)

            # 计算样本距离每个顶点的距离,并获得激活点的位置
            winner = get_winner_index(x, weights, dis_fun)

            # 根据激活点的位置计算临近点的权重
            new_sigma = get_learning_rate(sigma, t, N_steps)
            g = neighborhood_fun(X, Y, winner, new_sigma)
            g = g * eta

            # 进行权重的更新
            weights = weights + np.expand_dims(g, -1) * (x - weights)

        # 打印量化误差
        print("quantization_error= %.4f" % (get_quantization_error(datas, weights)))

    return weights


def get_U_Matrix(weights):
    X, Y, D = np.shape(weights)
    um = np.nan * np.zeros((X, Y, 8))  # 8邻域

    ii = [0, -1, -1, -1, 0, 1, 1, 1]
    jj = [-1, -1, 0, 1, 1, 1, 0, -1]

    for x in range(X):
        for y in range(Y):
            w_2 = weights[x, y]

            for k, (i, j) in enumerate(zip(ii, jj)):
                if (x + i >= 0 and x + i < X and y + j >= 0 and y + j < Y):
                    w_1 = weights[x + i, y + j]
                    um[x, y, k] = np.linalg.norm(w_1 - w_2)

    um = np.nansum(um, axis=2)
    return um / um.max()


def stage_normalization_weight_calculate(datas, labels):
    stage_weights = np.zeros(5)
    counts = np.zeros(5)
    for i in range(5):
        for j in range(len(labels)):
            if j == i:
                stage_weights[i] = stage_weights[i] + np.sum(datas[j])
                counts[i] = counts[i] + 1
    for i in range(5):
        stage_weights[i] = stage_weights[i] / counts[i]

    return stage_weights


def stage_normalization_entropy_calculate(datas, stage_weights, labels):
    # 计算每个病人的熵
    entropy_list = []
    for i in range(len(labels)):
        entropy = 0
        sum = np.sum(datas[i])
        for xi in datas[i]:
            if xi == 0:
                entropy = entropy + 0
            else:
                entropy = entropy + (sum / stage_weights[labels[i]]) * xi * np.log2(
                    (sum / stage_weights[labels[i]]) * xi)
        entropy_list.append(900 - entropy)
    return entropy_list


def SurfacePlot(func, data, fittedParameters):
    graphWidth = 800  # units are pixels
    graphHeight = 600  # units are pixels

    f = plt.figure(figsize=(graphWidth/100.0, graphHeight/100.0), dpi=100)

    plt.grid(True)
    axes = Axes3D(f)

    x_data = data[0]
    y_data = data[1]
    z_data = data[2]

    xModel = np.linspace(min(x_data), max(x_data), 20)
    yModel = np.linspace(min(y_data), max(y_data), 20)
    X, Y = np.meshgrid(xModel, yModel)

    Z = func(np.array([X, Y]), *fittedParameters)

    axes.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=1, antialiased=True)

    axes.scatter(x_data, y_data, z_data) # show data along with plotted surface

    axes.set_title('Surface Plot (click-drag with mouse)') # add a title for surface plot
    axes.set_xlabel('X Data') # X axis data label
    axes.set_ylabel('Y Data') # Y axis data label
    axes.set_zlabel('Z Data') # Z axis data label

    plt.show()
    plt.close('all') # clean up after using pyplot or else thaere can be memory and process problems


if __name__ == '__main__':
    data = pd.read_excel('TCGA_LUAD_all_deg_sur_expr_sub_smoke.xlsx')
    labs = data['stage'].values
    # print(labs)
    label_names = {0: 'Normal', 1: 'Stage I', 2: 'Stage II', 3: 'Stage III', 4: 'Stage IV', 5: 'interpolation'}
    datas = data[data.columns[1:-7]].values
    N, D = np.shape(datas)
    subtype = data['subtype'].values
    years_smoked = data['years_smoked'].values
    pack_years_smoked = data['pack_years_smoked'].values

    # 计算每个病人的熵
    stage_weights = stage_normalization_weight_calculate(datas, labs)
    entropy_list = stage_normalization_entropy_calculate(datas, stage_weights, labs)
    entropy_list = np.array(entropy_list)

    # SOM的训练 %LUAD数据集在30*20的grid下满足3-限制激活，故直接使用原始som算法
    X = 30
    Y = 20
    weights = train_SOM(X=X, Y=Y, N_epoch=500, datas=datas, sigma=1.5, init_weight_fun=weights_PCA)

    # 获取UMAP
    UM = get_U_Matrix(weights)

    plt.figure(figsize=(30, 30))
    plt.pcolor(UM.T, cmap='bone_r')  # plotting the distance map as background
    plt.colorbar()

    markers = ['o', 's', 'D', 'x', '+']
    colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']

    for i in range(N):
        x = datas[i]
        w = get_winner_index(x, weights)
        i_lab = labs[i] - 1

        plt.plot(w[0] + .5, w[1] + .5, markers[i_lab], markerfacecolor='None',
                 markeredgecolor=colors[i_lab], markersize=10, markeredgewidth=2)

    plt.show()

    UM = get_U_Matrix(weights)

    # print(UM)
    # '''画散点图'''

    # 显示UMAP
    plt.figure(1, figsize=(30, 30))
    plt.pcolor(UM.T, cmap='bone_r')  # plotting the distance map as background
    plt.colorbar()

    # 计算每个样本点投射后的坐标
    w_x, w_y = zip(*[get_winner_index(d, weights) for d in datas])
    w_x = np.array(w_x)
    w_y = np.array(w_y)

    # 分别把每一类的散点在响应的方格内进行打印（+随机位置偏移）
    for c in np.unique(labs):
        idx_target = (labs == c)
        plt.scatter(w_x[idx_target] + .5 + (np.random.rand(np.sum(idx_target)) - .5) * .8,
                    w_y[idx_target] + .5 + (np.random.rand(np.sum(idx_target)) - .5) * .8,
                    s=50, c=colors[c], label=label_names[c])
    plt.legend(loc='upper right')
    plt.grid()
    # plt.show()

    fig = plt.figure(figsize=(30, 30))
    ax = plt.axes(projection="3d")
    for c in np.unique(labs):
        idx_target = (labs == c)
        ax.scatter3D(w_x[idx_target] + .5 + (np.random.rand(np.sum(idx_target)) - .5) * .8,
                     w_y[idx_target] + .5 + (np.random.rand(np.sum(idx_target)) - .5) * .8,
                     entropy_list[idx_target], s=75, c=colors[c], label=label_names[c])
    plt.legend(loc='upper right')
    plt.show()


#  ''' 画饼图'''
    # 计算输出层的每个节点上映射了哪些数据
    win_map = defaultdict(list)
    for x, lab in zip(datas, labs):
        win_map[get_winner_index(x, weights)].append(lab)

    # 统计每个输出节点上，映射了各类数据、各多少个
    for pos in win_map:
        win_map[pos] = Counter(win_map[pos])

    fig = plt.figure(2, figsize=(30, 30))
    # 按照 X,Y对画面进行分格
    the_grid = gridspec.GridSpec(Y, X, fig)

    # 在每个格子里面画饼图
    for pos in win_map.keys():
        label_fracs = [win_map[pos][l] for l in label_names.keys()]

        plt.subplot(the_grid[Y - 1 - pos[1],
                             pos[0]], aspect=1)
        patches, texts = plt.pie(label_fracs)

    plt.legend(labels=label_names.values(), loc='best')
    # plt.savefig('resulting_images/som_seed_pies.png')
    plt.show()

    w_x_new = []
    w_y_new = []
    new_labs = []
    new_subtypes = []
    new_years_smoked = []
    new_pack_years_smoked = []
    node_list = []
    new_node_list = []
    new_entropy_list = []
    for i in range(len(w_x)):
        node_list.append([w_x[i], w_y[i]])

    for node in node_list:
        idx_target = []
        for i in range(len(node_list)):
            if node_list[i] == node:
                idx_target.append(i)
        if len(idx_target) == 1:
            new_node_list.append(node)
            w_x_new.append(node[0])
            w_y_new.append(node[1])
            new_labs.append(labs[idx_target[0]])
            new_subtypes.append(subtype[idx_target[0]])
            new_years_smoked.append(years_smoked[idx_target[0]])
            new_pack_years_smoked.append(pack_years_smoked[idx_target[0]])
            new_entropy_list.append(entropy_list[idx_target[0]])
        if len(idx_target) >= 2 and node not in new_node_list:
            new_node_list.append(node)
            w_x_new.append(node[0])
            w_y_new.append(node[1])
            min_lab = labs[idx_target[0]]
            min_idx = idx_target[0]
            for idx in idx_target:
                if labs[idx] > min_lab:
                    min_lab = labs[idx]
                    min_idx = idx
            new_labs.append(min_lab)
            new_subtypes.append(subtype[idx])
            new_years_smoked.append(years_smoked[idx])
            new_pack_years_smoked.append(pack_years_smoked[idx])
            new_entropy_list.append(entropy_list[idx])

    bottom = []
    width = []
    depth = []
    for i in range(len(w_x_new)):
        bottom.append(0)
        width.append(0.5)
        depth.append(0.5)

    # k递减最近邻插值
    X_som = np.arange(0, X, 1)
    Y_som = np.arange(0, Y, 1)
    X_som, Y_som = np.meshgrid(X_som, Y_som)
    # print(np.shape(X_som), np.shape(Y_som))
    Z_entropy = np.zeros([Y, X], dtype=float)
    # print(np.shape(Z_entropy))
    x_index = 0
    for x in w_x_new:
        Z_entropy[w_y_new[x_index], x] = new_entropy_list[x_index]
        x_index = x_index + 1
    # print(Z_entropy)
    cnt_array = np.where(Z_entropy, 0, 1)
    flag = np.sum(cnt_array)
    count = 0
    while flag >= 1:
        for i in range(Y):
            for j in range(X):
                if Z_entropy[i, j] == 0:
                    neighbor = []
                    if i == 0 and j == 0:
                        for m in range(i, i + 2):
                            for n in range(j, j + 2):
                                if Z_entropy[m, n] != 0:
                                    neighbor.append(Z_entropy[m, n])
                        if len(neighbor) >= 2:
                            Z_entropy[i, j] = sum(neighbor)/len(neighbor)
                    elif i == 0 and j == X-1:
                        for m in range(i, i + 2):
                            for n in range(j - 1, j + 1):
                                if Z_entropy[m, n] != 0:
                                    neighbor.append(Z_entropy[m, n])
                        if len(neighbor) >= 2:
                            Z_entropy[i, j] = sum(neighbor)/len(neighbor)
                    elif i == Y-1 and j == 0:
                        for m in range(i-1, i + 1):
                            for n in range(j, j + 2):
                                if Z_entropy[m, n] != 0:
                                    neighbor.append(Z_entropy[m, n])
                        if len(neighbor) >= 2:
                            Z_entropy[i, j] = sum(neighbor)/len(neighbor)
                    elif i == Y-1 and j == X-1:
                        for m in range(i-1, i + 1):
                            for n in range(j-1, j + 1):
                                if Z_entropy[m, n] != 0:
                                    neighbor.append(Z_entropy[m, n])
                        if len(neighbor) >= 2:
                            Z_entropy[i, j] = sum(neighbor)/len(neighbor)
                    elif i == 0 and j != 0 and j != X-1:
                        for m in range(i, i + 2):
                            for n in range(j-1, j + 2):
                                if Z_entropy[m, n] != 0:
                                    neighbor.append(Z_entropy[m, n])
                        if len(neighbor) >= 3:
                            Z_entropy[i, j] = sum(neighbor)/len(neighbor)
                    elif i == Y-1 and j != 0 and j != X-1:
                        for m in range(i - 1, i + 1):
                            for n in range(j-1, j + 2):
                                if Z_entropy[m, n] != 0:
                                    neighbor.append(Z_entropy[m, n])
                        if len(neighbor) >= 3:
                            Z_entropy[i, j] = sum(neighbor)/len(neighbor)
                    elif j == 0 and i != Y-1 and i != 0:
                        for m in range(i - 1, i + 2):
                            for n in range(j, j + 2):
                                if Z_entropy[m, n] != 0:
                                    neighbor.append(Z_entropy[m, n])
                        if len(neighbor) >= 3:
                            Z_entropy[i, j] = sum(neighbor) / len(neighbor)
                    elif j == X - 1 and i != Y-1 and i != 0:
                        for m in range(i - 1, i + 2):
                            for n in range(j - 1, j + 1):
                                if Z_entropy[m, n] != 0:
                                    neighbor.append(Z_entropy[m, n])
                        if len(neighbor) >= 3:
                            Z_entropy[i, j] = sum(neighbor)/len(neighbor)
                    else:
                        for m in range(i-1, i+2):
                            for n in range(j-1, j+2):
                                if Z_entropy[m, n] != 0:
                                    neighbor.append(Z_entropy[m, n])
                        if len(neighbor) >= 6:
                            Z_entropy[i, j] = sum(neighbor)/len(neighbor)
        cnt_array = np.where(Z_entropy, 0, 1)
        flag = np.sum(cnt_array)
        count = count + 1
        if count >= 100:
            flag = 0

    Z_labs = np.zeros([Y, X], dtype=int)
    for i in range(Y):
        for j in range(X):
            Z_labs[i, j] = 5
    lab_index = 0
    for x in w_x_new:
        Z_labs[w_y_new[lab_index], x] = new_labs[lab_index]
        lab_index = lab_index + 1

    Z_subtypes = np.zeros([Y, X], dtype=int)
    for i in range(Y):
        for j in range(X):
            Z_subtypes[i, j] = 0
    lab_index = 0
    for x in w_x_new:
        Z_subtypes[w_y_new[lab_index], x] = new_subtypes[lab_index]
        lab_index = lab_index + 1

    Z_years_smoked = np.zeros([Y, X], dtype=int)
    for i in range(Y):
        for j in range(X):
            Z_years_smoked[i, j] = 5
    lab_index = 0
    for x in w_x_new:
        Z_years_smoked[w_y_new[lab_index], x] = new_years_smoked[lab_index]
        lab_index = lab_index + 1

    Z_pack_years_smoked = np.zeros([Y, X], dtype=int)
    for i in range(Y):
        for j in range(X):
            Z_pack_years_smoked[i, j] = 5
    lab_index = 0
    for x in w_x_new:
        Z_pack_years_smoked[w_y_new[lab_index], x] = new_pack_years_smoked[lab_index]
        lab_index = lab_index + 1

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(X_som, Y_som, Z_entropy)
    plt.show()

    angles = [30, 45, 60]
    # view_angles 是 angles 的笛卡尔积
    view_angles = [(elev, azim) for elev in angles for azim in angles]

    # 创建一个 3x3 的子图
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(10, 10), subplot_kw={'projection': '3d'})

    # 在每个子图中绘制一个不同视角的 3D 散点图
    for ax, view_angle in zip(axes.flat, view_angles):
        ax.plot_surface(X_som, Y_som, Z_entropy)
        for i in range(Y):
            for j in range(X):
                ax.scatter3D(X_som[i, j], Y_som[i, j], Z_entropy[i, j] + 10, s=5, c=colors[Z_labs[i, j]])
        ax.view_init(elev=view_angle[0], azim=view_angle[1])
        # ax.set_title(f"View angle: {view_angle}")

    plt.tight_layout()
    plt.show()

    # workbook = xlwt.Workbook(encoding="utf-8")  # 创建workbook对象
    # worksheet0 = workbook.add_sheet("X")  # 创建工作表
    # worksheet1 = workbook.add_sheet("Y")
    # worksheet2 = workbook.add_sheet("Z")
    # worksheet3 = workbook.add_sheet("lab")
    # for i in range(Y):
    #     for j in range(X):
    #         worksheet0.write(i, j, str(X_som[i, j]))
    #         worksheet1.write(i, j, str(Y_som[i, j]))
    #         worksheet2.write(i, j, str(Z_entropy[i, j]))
    #         worksheet3.write(i, j, str(Z_labs[i, j]))
    # workbook.save("LUAD_deg_sur_som.xls")
    #
    # workbook = xlwt.Workbook(encoding="utf-8")  # 创建workbook对象
    # worksheet0 = workbook.add_sheet("X")  # 创建工作表
    # worksheet1 = workbook.add_sheet("Y")
    # worksheet2 = workbook.add_sheet("Z")
    # worksheet3 = workbook.add_sheet("lab")
    # worksheet4 = workbook.add_sheet("subtype")
    # worksheet5 = workbook.add_sheet("years_smoked")
    # worksheet6 = workbook.add_sheet("pack_years_smoked")
    # X_som_flatten = X_som.flatten()
    # Y_som_flatten = Y_som.flatten()
    # Z_entropy_flatten = Z_entropy.flatten()
    # Z_labs_flatten = Z_labs.flatten()
    # Z_subtypes_flatten = Z_subtypes.flatten()
    # Z_years_smoked_flatten = Z_years_smoked.flatten()
    # Z_pack_years_smoked_flatten = Z_pack_years_smoked.flatten()
    # for i in range(len(X_som_flatten)):
    #         worksheet0.write(i, 0, str(X_som_flatten[i]))
    #         worksheet1.write(i, 0, str(Y_som_flatten[i]))
    #         worksheet2.write(i, 0, str(Z_entropy_flatten[i]))
    #         worksheet3.write(i, 0, str(Z_labs_flatten[i]))
    #         worksheet4.write(i, 0, str(Z_subtypes_flatten[i]))
    #         worksheet5.write(i, 0, str(Z_years_smoked_flatten[i]))
    #         worksheet6.write(i, 0, str(Z_pack_years_smoked_flatten[i]))
    # workbook.save("LUAD_deg_sur_som_flatten_subtype_smoked.xls")

    # 双调和样条插值曲面
    xi = np.arange(0, X - 0.5, 0.5)
    yi = np.arange(0, Y - 0.5, 0.5)
    xi, yi = np.meshgrid(xi, yi)
    zi = griddata(np.hstack((X_som.flatten()[:, None], Y_som.flatten()[:, None])), Z_entropy.flatten(), (xi, yi),
                  method='cubic')
    Zi_entropy = np.nan_to_num(zi)
    print(Zi_entropy, np.shape(Zi_entropy))
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(xi, yi, Zi_entropy)
    for i in range(Y):
        for j in range(X):
            ax.bar3d(X_som[i, j], Y_som[i, j], 0, 0.5, 0.5, Z_entropy[i, j] - 10,
                     color=colors[Z_labs[i, j]], label=label_names[Z_labs[i, j]])
    plt.show()

    angles = [30, 45, 60]
    # view_angles 是 angles 的笛卡尔积
    view_angles = [(elev, azim) for elev in angles for azim in angles]

    # 创建一个 3x3 的子图
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(10, 10), subplot_kw={'projection': '3d'})

    # 在每个子图中绘制一个不同视角的 3D 散点图
    for ax, view_angle in zip(axes.flat, view_angles):
        ax.plot_surface(xi, yi, Zi_entropy)
        for i in range(Y):
            for j in range(X):
                ax.bar3d(X_som[i, j], Y_som[i, j], 0, 0.5, 0.5, Z_entropy[i, j] - 10,
                         color=colors[Z_labs[i, j]], label=label_names[Z_labs[i, j]])
        ax.view_init(elev=view_angle[0], azim=view_angle[1])
        # ax.set_title(f"View angle: {view_angle}")
    plt.tight_layout()
    plt.show()

    # 3D旋转图像及信息染色图形见matlab程序
