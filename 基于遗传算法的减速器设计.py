#encoding=UTF-8
from math import *
import numpy as np
# import 曲线拟合
from 一对齿轮传动设计计算 import *
import random
import matplotlib.pyplot as plt


DNA_SIZE = 11#碱基对数
RST_SIZE = 250#染色体数
CROSSOVER_RATE = 0.8#交叉率
MUTATION_RATE = 0.005#突变率
N_GENERATIONS = 5000#迭代数目
isum=17.905#减速器总传动比
n1=720#输入轴转速
P1=7.128#输入轴传递效率
xiaolv=0.9603#输入轴到中间轴的传递效率
i1_BOUND = [0.001, 17.905]#i1限定范围
Pc1=0.7#交叉概率常数
Pc2=0.7#交叉概率常数
Pm1=0.08#突变概率常数
Pm2=0.08#突变概率常数
fm=0.5
fmax=fzj=1
NN=1
zbb1=[]
zbb2=[]
zbb3=[]
zbb4=[]
zbb5=[]
zbb6=[]
zbb7=[]
zbbn=[]
fig, ax = plt.subplots()
# _BOUND = [, ]#限定范围
# _BOUND = [, ]#限定范围
# _BOUND = [, ]#限定范围
# _BOUND = [, ]#限定范围
# _BOUND = [, ]#限定范围
# _BOUND = [, ]#限定范围
# _BOUND = [, ]#限定范围
# _BOUND = [, ]#限定范围
# _BOUND = [, ]#限定范围
# _BOUND = [, ]#限定范围
#11维向量
# BOUND=[[0.001,17.904],[20,40],[20,40],[8,20],[8,20],[0.8,1.0],[0.8,1.0],[1,1.2],[1,1.2],[16,18],[16,18]]
#对应于第一级传动比，第一级主动齿轮齿数，第二级主动齿轮齿数，第一级螺旋角，
# 第二级螺旋角，第一级齿宽系数，第二级齿宽系数，第一级齿顶高系数，第二级齿顶高系数，第一级压力角，第二级压力角
#5维向量
BOUND=[[4.5,5.5],[20,40],[20,40],[10,15],[10,15],[0.8,1.0],[0.8,1.0],[1,1],[1,1],[20,20],[20,20]]
#自适应遗传算法

#生成初始种群
POP=[]
for j in range(RST_SIZE):
    RST=[]
    i=0
    for list in BOUND:
        # print(i)
        if i==1 or i==2:
            RST.append(random.randint(list[0], list[1]))
        elif i==5 or i==6:
            RST.append(float(np.random.choice(list, size=1, replace=True)))
        else:
            RST.append(random.uniform(list[0],list[1]))
        i=i+1
    POP.append(RST)

#目标函数
def F(i1,z11,z21,bt1,bt2,sdd1,sdd2,ha1,ha2,aef1,aef2):
    f = open(r"记录last1.txt", "a", encoding="utf-8")
    print('-------个体--------')
    f.write('-------个体--------\n')
    print('<初选参数>')
    f.write('<初选参数>\n')
    print('高速级传动比：{0}，低速级传动比：{1}，高速级小齿轮齿数：{2},低速级小齿轮齿数：{3}，高速级齿宽系数：{4}，低速级齿宽系数：{5}'
          '高速级齿顶高系数：{6}，低速级齿顶高系数：{7}，'
          '高速级压力角：{8},低速级压力角：{9},'
          '高速级螺旋角：{10}，低速级螺旋角：{11}'.format(i1,isum/i1,z11,z21,sdd1,sdd2,ha1,ha2,aef1,aef2,bt1,bt2))
    f.write('高速级传动比：{0}，低速级传动比：{1}，高速级小齿轮齿数：{2},低速级小齿轮齿数：{3}，高速级齿宽系数：{4}，低速级齿宽系数：{5}'
          '高速级齿顶高系数：{6}，低速级齿顶高系数：{7}，'
          '高速级压力角：{8},低速级压力角：{9},'
          '高速级螺旋角：{10}，低速级螺旋角：{11}\n'.format(i1,isum/i1,z11,z21,sdd1,sdd2,ha1,ha2,aef1,aef2,bt1,bt2))
    P2=xiaolv*P1
    n2=n1/i1
    i2=isum/i1
    sz1=sheji(i1,z11,bt1,sdd1,ha1,aef1,n1,0.5*P1)
    z11, z12, mn1, aef1, bt1, a1, b11, b12=sz1
    sz2=sheji(i2, z21, bt2, sdd2, ha2, aef2, n2, P2)
    z21, z22, mn2, aef2, bt2, a2, b21, b22=sz2
    # print('<测试>:',sz2)
    #设计目标
    #指标一：总中心距  越小越好
    zb1=a1+a2
    aeft1 = degrees(atan(tan(hudu(aef1)) / cos(hudu(bt1))))
    aefat1 = degrees(acos(z11 * cos(hudu(aeft1)) / (z11 + 2 * ha1 * cos(hudu(bt1)))))
    aefat2 = degrees(acos(z12 * cos(hudu(aeft1)) / (z12 + 2 * ha1 * cos(hudu(bt1)))))
    yfxnaef = (z11 * (tan(hudu(aefat1)) - tan(hudu(aeft1))) + z12 * (tan(hudu(aefat2)) - tan(hudu(aeft1)))) / (2 * pi)
    yfxnbt = sdd1 * z11 * tan(hudu(bt1)) / pi
    aeft2 = degrees(atan(tan(hudu(aef2)) / cos(hudu(bt2))))
    aefat1 = degrees(acos(z21 * cos(hudu(aeft2)) / (z21 + 2 * ha2 * cos(hudu(bt2)))))
    aefat2 = degrees(acos(z22 * cos(hudu(aeft2)) / (z22 + 2 * ha2 * cos(hudu(bt2)))))
    yfxnaef = (z21 * (tan(hudu(aefat1)) - tan(hudu(aeft2))) + z22 * (tan(hudu(aefat2)) - tan(hudu(aeft2)))) / (2 * pi)
    yfxnbt = sdd2 * z21 * tan(hudu(bt2)) / pi
    #指标二：重合度的倒数 越小越好
    zb2=1/(yfxnaef+yfxnbt)
    #指标三：两大齿轮浸油深度相近程度度 越大越好（指标越小越好）
    zb3=abs(mn1*z12/cos(pi*bt1/180)-mn2*z22/cos(pi*bt2/180))
    #指标四：  保证高速级齿轮大齿轮与低速级轴相近（采用梯度法）
    if mn1*z12/cos(bt1)/2>=a2:
        zb4=100
    elif mn1*z12/cos(bt1)/2+20>=a2:
        zb4=50
    elif mn1*z12/cos(bt1)/2+30>=a2:
        zb4=25
    elif mn1*z12/cos(bt1)/2+40>=a2:
        zb4=12
    elif mn1*z12/cos(bt1)/2+50>=a2:
        zb4=5
    elif mn1*z12/cos(bt1)/2+52>=a2:
        zb4=4
    elif mn1*z12/cos(bt1)/2+54>=a2:
        zb4=3
    elif mn1*z12/cos(bt1)/2+56>=a2:
        zb4=2
    elif mn1*z12/cos(bt1)/2+58>=a2:
        zb4=1
    elif mn1*z12/cos(bt1)/2+60>=a2:
        zb4=0.5
    elif mn1*z12/cos(bt1)/2+62>=a2:
        zb4=0.25
    elif mn1*z12/cos(bt1)/2+64>=a2:
        zb4=0.125
    elif mn1*z12/cos(bt1)/2+65>=a2:
        zb4=0.01
    else:
        zb4=0
    #7900*pi/4(b11*(mn1*z11/cos(bt1))**2+b21*(mn2*z21/cos(bt2))**2)*10**(-9)+7850*pi/4(b12*(mn1*z12/cos(bt1))**2+b22*(mn2*z22/cos(bt2))**2)*10**(-9)
    #指标五：4个齿轮的总质量 越小越好
    zb5=7900*pi/4*(b11*(mn1*z11/cos(bt1))**2+b21*(mn2*z21/cos(bt2))**2)*10**(-9)+7850*pi/4*(b12*(mn1*z12/cos(bt1))**2+b22*(mn2*z22/cos(bt2))**2)*10**(-9)
    #指标六 模数 越小越好
    # zb6=mn1+mn2
    #指标七 是否满足接触强度和弯曲强度
    #jiaohe(z1,z2,mn,aef,bt,a,b1,b2,n1,P1):
    zb7=jiaohe(z11, z12, mn1, aef1, bt1, a1, b11, b12,n1, 0.5*P1)+jiaohe(z21, z22, mn2, aef2, bt2, a2, b21, b22,n2, P2)

    #权重
    w1=1/392.15*3
    w2=1/0.288877
    w3=1/8
    w4=1/50
    w5=1/5002.366
    # w6=1
    zbb1.append(zb1 * w1)
    zbb2.append(zb2)
    zbb3.append(zb3)
    zbb4.append(zb4)
    zbb5.append(zb5)
    # zbb6.append(zb6)
    zbb7.append(zb7)
    # NN = len(zbb7)
    # NNN = [i for i in range(NN)]
    #总指标(适应度），权重包含了对量纲的考量
    zbsum=w1*zb1+w2*zb2+w3*zb3+w4*zb4+w5*zb5+zb7
    print('<结果>')
    f.write('<结果>\n')
    print('高速级：齿数为z11={0},z12={1},模数mn1为{2}mm，压力角aef1={3}度,螺旋角bt1={4}度，中心距a1={5}mm,齿宽b11={6}mm,b2={7}mm.'
          '低速级：齿数为z21={8},z22={9},模数mn2为{10}mm，压力角aef2={11}度,螺旋角bt2={12}度，中心距a2={13}mm,齿宽b21={14}mm,'
          'b2={15}mm'.format(z11,z12,mn1, aef1, bt1, a1, b11, b12,z21, z22, mn2,aef2, bt2,  a2, b21, b22))
    f.write('高速级：齿数为z11={0},z12={1},模数mn1为{2}mm，压力角aef1={3}度,螺旋角bt1={4}度，中心距a1={5}mm,齿宽b11={6}mm,b2={7}mm.'
          '低速级：齿数为z21={8},z22={9},模数mn2为{10}mm，压力角aef2={11}度,螺旋角bt2={12}度，中心距a2={13}mm,齿宽b21={14}mm,'
          'b2={15}mm\n'.format(z11,z12,mn1, aef1, bt1, a1, b11, b12,z21, z22, mn2, aef2, bt2, a2, b21, b22))

    # print(sz1,sz2)
    print('指标1：{0}，指标2：{1}，指标3：{2}，指标4：{3}，指标5：{4},指标6：{5}'.format(w1*zb1,w2*zb2,w3*zb3,w4*zb4,w5*zb5,zb7))
    f.write('指标1：{0}，指标2：{1}，指标3：{2}，指标4：{3},，指标5：{4},指标6：{5}\n'.format(w1*zb1,w2*zb2,w3*zb3,w4*zb4,w5*zb5,zb7))
    print('不适应度:', zbsum)
    f.write('不适应度:'+str(zbsum)+'\n')
    # if 'sca' in locals():
    #     sca.remove()
    # sca = plt.scatter([1,2,3,4,5,6,7],[w1*zb1,w2*zb2,w3*zb3,w4*zb4,w5*zb5,zb7,zbsum] , s=200, lw=0, c='red', alpha=0.5)
    # plt.pause(0.0000000001)
    # sca = plt.scatter(NNN, zbb7, s=200,
    #                   lw=0, c='red', alpha=0.5)
    # plt.pause(0.0000000001)
    # ax.cla()  # clear plot
    # ax.plot(NNN, zbb1, 'r', lw=1)  # draw line chart
    # # ax.bar(zbb1, height=zbb1, width=0.3) # draw bar chart
    # # do_something()
    # plt.pause(0.0000000000000000000000000000000000000001)
    f.close()
    return zbsum
    # print(y)
# F(4.92,970,11,0.99,3.2)







#获取适应度
def get_fitness(POP):
    PRED=[]
    for list in POP:
        #F(i1,z11,z21,bt1,bt2,sdd1,sdd2,ha1,ha2,aef1,aef2)
        i1, z11, z21, bt1, bt2, sdd1, sdd2, ha1, ha2, aef1, aef2 = list
        pred = 1/F(i1,z11,z21,bt1,bt2,sdd1,sdd2,ha1,ha2,aef1,aef2)
        PRED.append(pred)
    return np.array(PRED)
        # shiyingdu=(pred - np.min(pred)) + 1e-3  # 减去最小的适应度是为了防止适应度出现负数，通过这一步fitness的范围为[0, np.max(pred)-np.min(pred)],最后在加上一个很小的数防止出现为0的适应度



#交叉与突变运算
def crossover_and_mutation(POP, CROSSOVER_RATE=0.8):
    new_POP = []

    for father in POP:  # 遍历种群中的每一个个体，将该个体作为父亲
        i1, z11, z21, bt1, bt2, sdd1, sdd2, ha1, ha2, aef1, aef2 = father
        fzj = F(i1, z11, z21, bt1, bt2, sdd1, sdd2, ha1, ha2, aef1, aef2)
        #自适应遗传算法
        if fzj>fm:
            CROSSOVER_RATE=Pc1*(Pc1*Pc2)*(fzj-fm)/(fmax-fm)
        else:
            CROSSOVER_RATE =Pc1
        child = father  # 孩子先得到父亲的全部基因（这里我把一串二进制串的那些0，1称为基因）
        if np.random.rand() < CROSSOVER_RATE:  # 产生子代时不是必然发生交叉，而是以一定的概率发生交叉
            mother = POP[np.random.randint(RST_SIZE)]  # 再种群中选择另一个个体，并将该个体作为母亲
            cross_points = np.random.randint(low=0, high=DNA_SIZE )  # 随机产生交叉的点
            child[cross_points:] = mother[cross_points:]  # 孩子得到位于交叉点后的母亲的基因
        mutation(child)  # 每个后代有一定的机率发生变异
        new_POP.append(child)
    #
    return new_POP


#突变函数
def tubian(i):
    if i == 1 or i == 2:
        return random.randint(BOUND[i][0], BOUND[i][1])
    elif i == 5 or i == 6:
        return float(np.random.choice(BOUND[i], size=1, replace=True))
    else:
        return random.uniform(BOUND[i][0], BOUND[i][1])

#子类突变
def mutation(child, MUTATION_RATE=0.003):
    if fzj> fm:
        MUTATION_RATE = Pm1 * (Pm1 * Pm2)*(fzj- fm) / (fmax - fm)
    else:
        MUTATION_RATE = Pc1
    if np.random.rand() < MUTATION_RATE:  # 以MUTATION_RATE的概率进行变异
        mutate_point = np.random.randint(0, DNA_SIZE)  # 随机产生一个实数，代表要变异基因的位置
        child[mutate_point] =tubian( mutate_point)# 将变异点的二进制为反转

#自然选择 有放缩效果
def select(POP, fitness):  # nature selection wrt POP's fitness
    #replace等于True表示放回，从[0,POP_SIZE)中对应的概率p抽取POP_SIZE个数组成一维序列
    idx = np.random.choice(np.arange(RST_SIZE), size=RST_SIZE, replace=True,
                           p=(fitness) / (fitness.sum()))
    zbb7.append(fitness.sum()/len(fitness))
    return POP[idx]

if __name__ == "__main__":


    for i in range(N_GENERATIONS):  # 迭代N代

        print('----------种群{}开始进化----------'.format(i))
        f = open(r"记录last1.txt", "a", encoding="utf-8")
        f.write('----------种群{}开始进化----------\n'.format(i))
        f.close()
        POP = np.array(crossover_and_mutation(POP, CROSSOVER_RATE))
        # F_values = F(translateDNA(POP)[0], translateDNA(POP)[1])#x, y --> Z matrix
        fitness = get_fitness(POP)
        fmax=fitness.max()+0.01
        fm=fitness.mean()

        # CROSSOVER_RATE = 0.8  # 交叉率
        # MUTATION_RATE = 0.005  # 突变率
        print('----------种群{}进化完成----------'.format(i))
        f = open(r"记录last1.txt", "a", encoding="utf-8")
        f.write('----------种群{}进化完成----------\n'.format(i))
        f.close()
        # ax.cla()  # clear plot
        NN = len(zbb7)
        NNN = [i for i in range(NN)]
        print(fitness.sum())
        # ax.plot(NNN, zbb7, 'r', lw=1)  # draw line chart
        # ax.bar(zbb1, height=zbb1, width=0.3) # draw bar chart
        # do_something()
        # plt.pause(0.0000000000000000000000000000000000000001)
        plt.plot(NNN, zbb7, 'bp')
        plt.show()


        POP = select(POP, fitness)  # 选择生成新的种群
    # ax.cla()  # clear plot
    NN = len(zbb7)
    NNN = [i for i in range(NN)]
    print(fitness.sum())
    # ax.plot(NNN, zbb7, 'r', lw=1)  # draw line chart
    # ax.bar(zbb1, height=zbb1, width=0.3) # draw bar chart
    # do_something()
    # plt.pause(0.0000000000000000000000000000000000000001)
    plt.plot(NNN, zbb7, 'bp')
