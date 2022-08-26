#encoding=UTF-8
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from sklearn.linear_model import Lasso
from sklearn.linear_model import SGDRegressor

def xxtofxx(x,N):
    list=[]
    for i in range(N):
        list.append(x**i)
    return list
def xxtofxxsz(X1,N):
    XX=[]
    for i in X1:
        XX.append(xxtofxx(i,N))
    return np.array(XX)
def onetotwo(y,n):
    return np.array(y).reshape(( n,1))
def qvxiannihe(X,Y,x,N,tu_name='曲线拟合图',X_name='X轴',Y_name='Y轴'):
    '''通过给出的X一维列表和Y一维列表组成的点进行曲线拟合，曲线拟合幂数为N，通过横坐标x返回该曲线的纵坐标'''
    #通过给定点进行多项式回归
    X1=xxtofxxsz(X,N)
    Y_length=len(Y)
    Y1=onetotwo( Y,Y_length)
    lin_reg = LinearRegression()
    lin_reg.fit(X1, Y1)

    #在图上绘制给定点
    plt.plot(X, Y1, 'ro')

    #求出需要点的纵坐标并绘制在图上
    X_new = np.array(xxtofxxsz([x],N))
    x_list=[]
    x_list.append(x)
    x1=onetotwo(x,1)
    Y_new = lin_reg.predict(X_new)
    # Y_new = sgd_reg.predict(X_new)
    Y_new1=Y_new[0][0]
    aa='({0},{1})'.format(x,Y_new1)
    plt.text(x, Y_new1, aa, ha='left', va='center', fontsize=10.5)
    plt.plot([x, x], [min(Y)-0.5, Y_new1], c='b', linestyle='--')
    plt.plot([min(X)-10, x], [Y_new1, Y_new1], c='b', linestyle='--')
    plt.plot(x1, Y_new, 'bp')



    #绘制拟合曲线
    x_max=1000*int(max(X))
    x_min=1000*int(min(X))
    #
    X2 = [0.001 * i for i in range(x_min,x_max)]
    X22 = xxtofxxsz(X2,N)
    Y2 = lin_reg.predict(X22)
    plt.plot(onetotwo(X2, x_max-x_min), Y2, 'k-', linewidth=1)

    # 设置中文乱码问题
    plt.rcParams['font.sans-serif'] = ['SimHei']
    # 设置图标标题，并在坐标轴上添加标签
    plt.title(tu_name, fontsize=24)
    plt.xlabel(X_name, fontsize=14)
    plt.ylabel(Y_name, fontsize=14)
    plt.axis([min(X)-10,max(X)+10,  min(Y)-0.1, max(Y)+0.2])

    # 保存图
    # 调用savefig将拟合曲线保存为表名
    plt.savefig(tu_name)
    plt.show()
    # print(Y_new1)
    return Y_new1


# qvxiannihe([400,700,800,950,1200,1450,1600,2000,2400,2800],[0.67,1.07,1.19,1.37,1.66,1.92,2.07,2.44,2.74,2.98],1440,6,tu_name='A-125基本额定功率Po图',X_name='小带轮转速(r/min)',Y_name='基本额定功率Po(KW)')
# qvxiannihe([400,700,800,950,1200,1450,1600,2000,2400,2800],[0.05,0.09,0.10,0.11,0.15,0.17,0.19,0.24,0.29,0.34],1440,5,tu_name='A-2额定功率增量图',X_name='小带轮转速(r/min)',Y_name='额定功率增量(KW)')

# qvxiannihe([180,175,170,165,160,155,150,145,140,135,130,125,120],[1.00,0.99,0.98,0.96,0.95,0.93,0.92,0.91,0.89,0.88,0.86,0.84,0.82],166.27,5,tu_name='包角修正系数图',X_name='小带轮包角',Y_name='包角修正系数')
# qvxiannihe([17,18,19,20,21,22,23,24,25,26,27,28,29,30,35,40,45,50,60,70,80,90,100,125,150,175,200],
#                            [2.97,2.91,2.85,2.76,2.72,2.69,2.65,2.62,2.60,2.57,2.55,2.53,2.52,2.45,2.40,2.35,2.32,2.28,2.24,2.22,2.20,2.20,2.18,2.16,2.14,2.13,2.12], 20, 7 )
# qvxiannihe([17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 150, 200],
#            [1.52, 1.53, 1.54, 1.55, 1.56, 1.57, 1.575, 1.58, 1.59, 1.595, 1.60, 1.61, 1.62, 1.625, 1.65, 1.67, 1.68,
#             1.70, 1.73, 1.75, 1.77, 1.78, 1.79, 1.83, 1.865], 20, 6)
# qvxiannihe([40,80,120,160,200],
#                    [1.201,1.210,1.219,1.229,1.238],10, 5 )
# kvlist_x=[0,5,10,15,20,25,30,40]
# kvlist_y=[1,1.12,1.12,1.22,1.26,1.28,1.31,1.34]
# qvxiannihe([0,5,10,15,20,25,30,40],
#                         [1,1.08,1.12,1.22,1.26,1.28,1.31,1.34], 10, 3,tu_name='动载系数Kv',X_name='v/(m/s)',Y_name='Kv')