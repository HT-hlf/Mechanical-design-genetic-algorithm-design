#encoding=UTF-8
from math import *
import numpy as np
import 曲线拟合
import random
import matplotlib.pyplot as plt

DNA_SIZE = 11#碱基对数
POP_SIZE = 200#染色体数
CROSSOVER_RATE = 0.8#交叉率
MUTATION_RATE = 0.005#突变率
N_GENERATIONS = 500#迭代数目
i1_BOUND = [0.001, 17.905]#i1限定范围
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
BOUND=[[0.001,17.904],[20,40],[20,40],[8,20],[8,20],[0.8,1.0],[0.8,1.0],[1,1.2],[1,1.2],[16,18],[16,18]]
POP=[]
i=0
# for list in BOUND:
#     i = i + 1
#     # print(i)
#     if i==2 or i==3:
#         POP.append(random.randint(list[0], list[1]))
#     elif i==6 or i==7:
#         POP.append(np.random.choice(list, size=1, replace=True))
#     else:
#         POP.append(random.uniform(list[0],list[1]))
# print(POP)
# print(random.uniform(1, 10)) #1到10之间取一个随机小数
# print(random.randint(1, 10)) #1到10之间取一个随机整数
# print(random.random()) #0.0到1.0之间取一个随机小数
def hudu(jiaodu):
    return jiaodu/180*pi
def cKhbt(sdd,b):
    if sdd==0.8:
        Khbt=曲线拟合.qvxiannihe([40,80,120,160,200],
                   [1.201,1.210,1.219,1.229,1.238],b , 5 )
    else:
        Khbt=曲线拟合.qvxiannihe([40, 80, 120, 160, 200],
                        [1.315,1.342,1.334,1.343,1.352], b, 5,)
    return Khbt

def cKhbt(sdd,b):
    if sdd==0.8:
        Khbt=曲线拟合.qvxiannihe([40,80,120,160,200],
                   [1.201,1.210,1.219,1.229,1.238],b , 5 )
    else:
        Khbt=曲线拟合.qvxiannihe([40, 80, 120, 160, 200],
                        [1.315,1.342,1.334,1.343,1.352], b, 5,)
    return Khbt

def cYfa(zv):
    Yfa = 曲线拟合.qvxiannihe([17,18,19,20,21,22,23,24,25,26,27,28,29,30,35,40,45,50,60,70,80,90,100,125,150,175,200],
                           [2.97,2.91,2.85,2.76,2.72,2.69,2.65,2.62,2.60,2.57,2.55,2.53,2.52,2.45,2.40,2.35,2.32,2.28,2.24,2.22,2.20,2.20,2.18,2.16,2.14,2.13,2.12], zv, 7 )
    return Yfa
def cYsa(zv):
    Ysa=曲线拟合.qvxiannihe([17,18,19,20,21,22,23,24,25,26,27,28,29,30,35,40,45,50,60,70,80,90,100,150,200],
                        [1.52,1.53,1.54,1.55,1.56,1.57,1.575,1.58,1.59,1.595,1.60,1.61,1.62,1.625,1.65,1.67,1.68,1.70,1.73,1.75,1.77,1.78,1.79,1.83,1.865], zv, 6)
    return Ysa
def shortlength_m(mnf):
    standard_m = [1,1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 16, 20, 25, 32, 40, 50]
    for i in range(18):
        if mnf<=standard_m[i] and i==0:
            return 1
        elif mnf<=standard_m[i]:
            if abs(standard_m[i-1]-mnf)>abs(standard_m[i]-mnf):
                return standard_m[i]
            else:
                return standard_m[i-1]
        else:
            pass




#目标函数
def F(i1,z1,bt,sdd,ha,aef,n1=970,P1=11):
    z2=round(i1*z1)
    print('z2:',z2)
    u1=(z2/z1)

    #接触疲劳强度设计
    Kht=1.3
    aeft=degrees(atan(tan(hudu(aef))/cos(hudu(bt))))
    btb=degrees(atan(tan(hudu(bt))*cos(hudu(aeft))))
    Zh=(2*cos(hudu(btb))/(cos(hudu(aeft))*sin(hudu(aeft))))**0.5
    print('Zh:', Zh)
    aefat1=degrees(acos(z1*cos(hudu(aeft))/(z1+2*ha*cos(hudu(bt)))))
    aefat2=degrees(acos(z2*cos(hudu(aeft))/(z2+2*ha*cos(hudu(bt)))))
    yfxnaef=(z1*(tan(hudu(aefat1))-tan(hudu(aeft)))+z2*(tan(hudu(aefat2))-tan(hudu(aeft))))/(2*pi)
    yfxnbt=sdd*z1*tan(hudu(bt))/pi
    Zyfxn=((4-yfxnaef)*(1-yfxnbt)/3+yfxnbt/yfxnaef)**0.5
    print('Zyfxn:',Zyfxn )
    Zbt=(cos(hudu(bt)))**0.5
    print('Zbt:', Zbt)
    #分流式齿轮
    # T1=0.5*9.55*10**6*P1/n1
    # 展开式齿轮
    T1 = 9.55 * 10 ** 6 * P1 / n1
    Ze=189.8
    sdhlim1=600
    sdhlim2=550
    N1=60*n1*1*16*350*10
    print('N1:', N1)
    N2=N1/u1
    print('N2:', N2)
    S=1
    Khn1=0.88
    Khn2=0.91
    xysdh1=Khn1*sdhlim1/S
    xysdh2 = Khn2 * sdhlim2 / S
    xysdh=min(xysdh1,xysdh2)
    print('xysdh1:', xysdh1)
    print('xysdh2:',xysdh2 )
    print('xysdh:',xysdh )
    d1t=(2*Kht*T1/sdd*(u1+1)/u1*(Zh*Ze*Zyfxn*Zbt/xysdh)**2)**(1/3)
    print('d1t:', d1t)

    # print('n1:', n1)
    v=pi*d1t*n1/60000
    print('v:', v)
    b=sdd*d1t
    print('b:', b)
    KA=1.00
    Kvh=1.10
    Ft1=2*T1/d1t
    pd1=KA*Ft1/b
    if pd1>=100:
        Khaef=1.2
    else:
        Khaef=1.4
    Khbt=cKhbt(sdd,b)
    print('Khbt:',Khbt )
    # print(Khbt)
    Kh=KA*Kvh*Khaef*Khbt
    print('Kh:', Kh)
    d1h=d1t*(Kh/Kht)**(1/3)
    print('d1h:', d1h)
    mnh=d1h*cos(hudu(bt))/z1
    print('mnh:', mnh)

    Kft=1.3
    yfxnaefv=yfxnaef/cos(hudu(btb))**2
    Yyfxn=0.25+0.75/yfxnaefv
    print('Yyfxn',Yyfxn )
    Ybt=1-yfxnbt*bt/120
    print('Ybt', Ybt)
    zv1=z1/cos(hudu(bt))**3
    zv2 = z2 / cos(hudu(bt)) ** 3
    Yfa1=cYfa(zv1)
    Yfa2=cYfa(zv2)
    Ysa1=cYsa(zv1)
    Ysa2=cYsa(zv2)
    print('Yfa1', Yfa1)
    print('Yfa2:', Yfa2)
    print('Ysa1:', Ysa1)
    print('Ysa2:', Ysa2)

    Kfn1=0.85
    Kfn2=0.88
    sdflim1=500
    sdflim2=320
    S=1.4
    xyf1=Kfn1*sdflim1/S
    xyf2 = Kfn2 * sdflim2 / S
    YYsd1=Yfa1*Ysa1/xyf1
    YYsd2 = Yfa2 * Ysa2 / xyf2
    YYsdqvda=max(YYsd1,YYsd2)
    print('YYsd1', YYsd1)
    print('YYsd2', YYsd2)
    print('YYsdqvda', YYsdqvda)
    mnt=(2*Kft*T1*Yyfxn*Ybt*cos(hudu(bt))**2/sdd/z1**2*YYsdqvda)**(1/3)
    print('mnt:', mnt)

    #调整模数
    d1=mnt*z1/cos(hudu(bt))
    v=pi*d1*n1/60
    b=sdd*d1
    cn=0.25
    h=(2*ha+cn)*mnt
    bh=b/h
    print('bh',bh )
    Kvf=1.08
    Ft1=2*T1/d1
    pd2=KA*Ft1/b
    if pd2 >= 100:
        Kfaef = 1.2
    else:
        Kfaef = 1.4
    Kfbt =1.28
    Kf=KA*Kvf*Kfaef*Kfbt
    print('Kfaef', Kfaef)
    print('Kf', Kf)

    mnf=mnt*(Kf/Kft)**(1/3)
    print('mnf', mnf)
    d1f=mnf*z1/cos(hudu(bt))
    print('d1f', d1f)
    # standard_m=[1,1.25,1.5,2,2.5,3,4,5,6,8,10,12,16,20,25,32,40,50]
    mn=shortlength_m(mnf)
    print('mn', mn)
    z1=round(d1h*cos(hudu(bt))/mn)
    print('z1', z1)
    z2=round(u1*z1)
    print('z2',z2)


    a=(z1+z2)*mn/2/cos(hudu(bt))
    print('a', a)
    bt=degrees(acos((z1+z2)*mn/2/a))
    print('bt', bt)

    d1=z1*mn/cos(hudu(bt))
    d2=z2*mn/cos(hudu(bt))
    print('d1', d1)
    print('d2',d2)

    b=sdd*d1
    b1=b+5
    b2=b
    print('b1', b1)
    print('b2',b2 )
    #(z1,z2,mn,aef,bt,a,b1,b2)
    return (z1,z2,mn,aef,bt,a,b1,b2)
# F(i1,z1,bt,sdd,ha,aef,n1=720,P1=7.128):
x=F(3.2,24,14,1,1,20,n1=970,P1=11)
print(x)
