#encoding=UTF-8
from math import *
import numpy as np
import 曲线拟合
#动载系数
# kvlist_x=[0,5,10,15,20,25,30,40]
# kvlist_y=[1,1.08,1.12,1.22,1.26,1.28,1.31,1.34]

def cKv(v):
    C = 7
    B = 0.25 * (C - 5.0) ** 0.667
    A = 50 + 56 * (1.0 - B)
    return (A / (A + (200 * v) ** 0.5)) ** (-B)
def cKfn(N):
    x=log(N, 10)
    if 2<=x<3:
        return 2.5
    elif 3<=x<6.4:
        return (1-2.5)/(6.4-3)*(x-3)+2.5
    elif 6.4<=x<10:
        return (0.8-1)/(10-6.4)*(x-10)+0.8
    else:
        print('Kfn出错')
        return 1

def cKhn(N):
    x=log(N, 10)
    if 4<=x<5:
        return 1.6
    elif 5<=x<7.4:
        return (1-1.6)/(7.67-5)*(x-5)+1.6
    elif 7.4<=x<10:
        return (0.86-1)/(10-7.67)*(x-10)+0.86
    else:
        print('Khn出错')
        return 1


def hudu(jiaodu):
    return jiaodu/180*pi
def cKhbt(sdd,b):
    if sdd==0.8:
        Khbt=曲线拟合.qvxiannihe([40,80,120,160,200],
                   [1.201,1.210,1.219,1.229,1.238],b , 5 ,tu_name='Khbt1')
    else:
        Khbt=曲线拟合.qvxiannihe([40, 80, 120, 160, 200],
                        [1.315,1.342,1.334,1.343,1.352], b, 5,tu_name='Khbt1')
    return Khbt
def cKhbt(sdd,b):
    if sdd==0.8:
        Khbt=曲线拟合.qvxiannihe([40,80,120,160,200],
                   [1.201,1.210,1.219,1.229,1.238],b , 5 ,tu_name='Khbt2')
    else:
        Khbt=曲线拟合.qvxiannihe([40, 80, 120, 160, 200],
                        [1.315,1.342,1.334,1.343,1.352], b, 5,tu_name='Khbt2')
    return Khbt



def cYfa1(zv):
    Yfa = 曲线拟合.qvxiannihe([17,18,19,20,21,22,23,24,25,26,27,28,29,30,35,40,45,50,60,70,80,90,100,125,150,175,200],
                           [2.97,2.91,2.85,2.76,2.72,2.69,2.65,2.62,2.60,2.57,2.55,2.53,2.52,2.45,2.40,2.35,2.32,2.28,2.24,2.22,2.20,2.20,2.18,2.16,2.14,2.13,2.12], zv, 7 ,tu_name='Yfa1')
    return Yfa
def cYfa2(zv):
    Yfa = 曲线拟合.qvxiannihe([17,18,19,20,21,22,23,24,25,26,27,28,29,30,35,40,45,50,60,70,80,90,100,125,150,175,200],
                           [2.97,2.91,2.85,2.76,2.72,2.69,2.65,2.62,2.60,2.57,2.55,2.53,2.52,2.45,2.40,2.35,2.32,2.28,2.24,2.22,2.20,2.20,2.18,2.16,2.14,2.13,2.12], zv, 7 ,tu_name='Yfa2')
    return Yfa
def cYsa1(zv):
    Ysa=曲线拟合.qvxiannihe([17,18,19,20,21,22,23,24,25,26,27,28,29,30,35,40,45,50,60,70,80,90,100,150,200],
                        [1.52,1.53,1.54,1.55,1.56,1.57,1.575,1.58,1.59,1.595,1.60,1.61,1.62,1.625,1.65,1.67,1.68,1.70,1.73,1.75,1.77,1.78,1.79,1.83,1.865], zv, 6,tu_name='Ysa1')
    return Ysa
def cYsa2(zv):
    Ysa=曲线拟合.qvxiannihe([17,18,19,20,21,22,23,24,25,26,27,28,29,30,35,40,45,50,60,70,80,90,100,150,200],
                        [1.52,1.53,1.54,1.55,1.56,1.57,1.575,1.58,1.59,1.595,1.60,1.61,1.62,1.625,1.65,1.67,1.68,1.70,1.73,1.75,1.77,1.78,1.79,1.83,1.865], zv, 6,tu_name='Ysa2')
    return Ysa
def shortlength_m(mnf):
    standard_m = [1,1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 16, 20, 25, 32, 40, 50]
    for i in range(18):
        if mnf<=standard_m[i] and i==0:
            return 1
        elif mnf<=standard_m[i]:
            return standard_m[i]
        else:
            pass
    #坏数据
    return 1
#一对齿轮强度校核函数
#z1,z2,mn,aef,bt,a,b1,b2
#i1,z1,bt,sdd,ha,aef,n1=970,P1=11
def jiaohe(z1,z2,mn,aef,bt,a,b1,b2,n1,P1):
    print('--------强度校核---------')
    print('--------齿面接触疲劳强度校核---------')
    # 齿面接触疲劳强度设计
    d1t = mn * z1 / cos(hudu(bt))
    sdd=b2/d1t
    print('d1t:',d1t)
    print('sdd:',sdd )
    u=z2/z1
    print('u:', u)
    # d1t=mn*z1/cos(hudu(bt))
    v = pi * d1t * n1 / 60000
    print('v:', v)
    KA = 1.00
    print('KA', KA)
    # Kvh=1.10
    Kvh = cKv(v)
    print('Kvh:', Kvh)
    T1 = 9.55 * 10 ** 6 * P1 / n1
    print('T1',T1)
    Ft1 = 2 * T1 / d1t
    print('Ft1', Ft1)
    pd1 = KA * Ft1 / b2
    print('pd1', pd1)
    if pd1 >= 100:
        Khaef = 1.2
    else:
        Khaef = 1.4
    print('Khaef', Khaef)
    Khbt = cKhbt(sdd, b2)
    # print('Khbt', Khbt)
    print('Khbt:', Khbt)
    # print(Khbt)
    Kh = KA * Kvh * Khaef * Khbt
    print('Kh:', Kh)
    aeft = degrees(atan(tan(hudu(aef)) / cos(hudu(bt))))
    btb = degrees(atan(tan(hudu(bt)) * cos(hudu(aeft))))
    Zh = (2 * cos(hudu(btb)) / (cos(hudu(aeft)) * sin(hudu(aeft)))) ** 0.5
    # print('Zh:', Zh)
    ha=1
    aefat1 = degrees(acos(z1 * cos(hudu(aeft)) / (z1 + 2 * ha * cos(hudu(bt)))))
    aefat2 = degrees(acos(z2 * cos(hudu(aeft)) / (z2 + 2 * ha * cos(hudu(bt)))))
    yfxnaef = (z1 * (tan(hudu(aefat1)) - tan(hudu(aeft))) + z2 * (tan(hudu(aefat2)) - tan(hudu(aeft)))) / (2 * pi)
    yfxnbt = sdd * z1 * tan(hudu(bt)) / pi
    Zyfxn = ((4 - yfxnaef) * (1 - yfxnbt) / 3 + yfxnbt / yfxnaef) ** 0.5
    # print('Zyfxn:',Zyfxn )
    Zbt = (cos(hudu(bt))) ** 0.5
    # print('Zbt:', Zbt)
    print('aeft', aeft)
    print('btb', btb)
    print('Zh', Zh)
    print('aefat1', aefat1)
    print('aefat2', aefat2)
    print('yfxnaef', yfxnaef)
    print('yfxnbt', yfxnbt)
    print('Zyfxn', Zyfxn)
    print('Zbt', Zbt)
    # 分流式齿轮
    # T1=0.5*9.55*10**6*P1/n1
    # 展开式齿轮
    T1 = 9.55 * 10 ** 6* P1 / n1
    print('T1', T1)
    Ze = 189.8
    print('Ze', Ze)
    sdhlim1 = 600
    print('sdhlim1', sdhlim1)
    sdhlim2 = 550
    print('sdhlim2', sdhlim2)
    print('n1:',n1)
    N1 = 60 * n1 * 1 * 16 * 350 * 10
    print('N1:', N1)
    N2 = N1 / u
    print('N2:', N2)
    S = 1
    print('S', S)
    Khn1 = cKhn(N1)
    print('Khn1', Khn1)
    Khn2 = cKhn(N2)
    print('Khn2', Khn2)
    xysdh1 = Khn1 * sdhlim1 / S
    xysdh2 = Khn2 * sdhlim2 / S
    xysdh = min(xysdh1, xysdh2)
    print('xysdh1:', xysdh1)
    print('xysdh2:', xysdh2)
    print('xysdh:', xysdh)
    sdhsj=(2*Kh*T1*(u+1)/(sdd*d1t**3*u))**0.5*Zh*Ze*Zyfxn*Zbt
    print('sdhsj:',sdhsj)
    if(sdhsj<=xysdh):
        print('满足齿面接触疲劳强度极限')
    else:
        print('不满足齿面接触疲劳强度极限')


    print('--------齿根弯曲疲劳强度校核---------')
    # 齿根弯曲疲劳强度校核
    cn = 0.25
    print('cn', cn)
    h = (2 * ha + cn) * mn
    print('h', h)
    bh = b2/ h
    print('bh', bh)
    # Kvf=1.08
    Kvf = cKv(v)
    print('Kvf:', Kvf)
    Ft1 = 2 * T1 / d1t
    print('Ft1', Ft1)
    pd2 = KA * Ft1 / b2
    print('pd2', pd2)
    if pd2 >= 100:
        Kfaef = 1.2
    else:
        Kfaef = 1.4
    print('Kfaef', Kfaef)
    Kfbt = 1.28
    print('Kfbt', Kfbt)
    Kf = KA * Kvf * Kfaef * Kfbt
    # print('Kfaef', Kfaef)
    print('Kf', Kf)
    yfxnaefv = yfxnaef / cos(hudu(btb)) ** 2
    print('yfxnaefv', yfxnaefv)
    Yyfxn = 0.25 + 0.75 / yfxnaefv
    # print('Yyfxn', Yyfxn)
    print('Yyfxn', Yyfxn)
    Ybt = 1 - yfxnbt * bt / 120
    print('Ybt', Ybt)
    zv1 = z1 / cos(hudu(bt)) ** 3
    zv2 = z2 / cos(hudu(bt)) ** 3
    print('zv1', zv1)
    print('zv2', zv2)
    Yfa1 = cYfa1(zv1)
    Yfa2 = cYfa2(zv2)
    Ysa1 = cYsa1(zv1)
    Ysa2 = cYsa2(zv2)
    print('Yfa1', Yfa1)
    print('Yfa2:', Yfa2)
    print('Ysa1:', Ysa1)
    print('Ysa2:', Ysa2)

    Kfn1 = cKfn(N1)
    Kfn2 = cKfn(N2)
    sdflim1 = 500
    sdflim2 = 320
    S = 1.4
    xyf1 = Kfn1 * sdflim1 / S
    xyf2 = Kfn2 * sdflim2 / S
    print('Kfn1', Kfn1)
    print('Kfn2', Kfn2)
    print('sdflim1', sdflim1)
    print('sdflim2', sdflim2)
    print('S', S)
    print('xyf1', xyf1)
    print('xyf2', xyf2)

    YYsd1 = Yfa1 * Ysa1 / xyf1
    YYsd2 = Yfa2 * Ysa2 / xyf2
    YYsdqvda = max(YYsd1, YYsd2)
    print('YYsd1', YYsd1)
    print('YYsd2', YYsd2)
    print('YYsdqvda', YYsdqvda)
    sdfsj1=2*Kf*T1*Yfa1*Ysa1*Yyfxn*Ybt*cos(hudu(bt))**2/(sdd*mn**3*z1**2)
    sdfsj2 = 2 * Kf * T1* Yfa2 * Ysa2 * Yyfxn * Ybt * cos(hudu(bt)) ** 2 / (sdd * mn ** 3 * z1 ** 2)
    sdfsjmax=max(sdfsj1,sdfsj2)
    print('sdfsj1:',sdfsj1)
    print('sdfsj2:', sdfsj2)
    print('sdfsjmax:', sdfsjmax)
    if (sdfsjmax <= xysdh):
        print('满足齿根弯曲强度极限')
    else:
        print('不满足齿根弯曲强度极限')





#一对齿轮设计函数
def sheji(i1,z1,bt,sdd,ha,aef,n1=970,P1=11):
    z2=round(i1*z1)
    print('z2:',z2)
    u1=(z2/z1)

    #接触疲劳强度设计
    Kht=1.3
    print('Kht', Kht)
    aeft=degrees(atan(tan(hudu(aef))/cos(hudu(bt))))
    btb=degrees(atan(tan(hudu(bt))*cos(hudu(aeft))))
    Zh=(2*cos(hudu(btb))/(cos(hudu(aeft))*sin(hudu(aeft))))**0.5
    # print('Zh:', Zh)
    aefat1=degrees(acos(z1*cos(hudu(aeft))/(z1+2*ha*cos(hudu(bt)))))
    aefat2=degrees(acos(z2*cos(hudu(aeft))/(z2+2*ha*cos(hudu(bt)))))
    yfxnaef=(z1*(tan(hudu(aefat1))-tan(hudu(aeft)))+z2*(tan(hudu(aefat2))-tan(hudu(aeft))))/(2*pi)
    yfxnbt=sdd*z1*tan(hudu(bt))/pi
    Zyfxn=((4-yfxnaef)*(1-yfxnbt)/3+yfxnbt/yfxnaef)**0.5
    # print('Zyfxn:',Zyfxn )
    Zbt=(cos(hudu(bt)))**0.5
    # print('Zbt:', Zbt)
    print('aeft', aeft)
    print('btb', btb)
    print('Zh', Zh)
    print('aefat1', aefat1)
    print('aefat2', aefat2)
    print('yfxnaef', yfxnaef)
    print('yfxnbt', yfxnbt)
    print('Zyfxn', Zyfxn)
    print('Zbt', Zbt)
    #分流式齿轮
    # T1=0.5*9.55*10**6*P1/n1
    # 展开式齿轮
    T1 = 9.55 * 10 ** 6 * P1 / n1
    print('T1',T1 )
    Ze=189.8
    print('Ze',Ze )
    sdhlim1=600
    print('sdhlim1',sdhlim1 )
    sdhlim2=550
    print('sdhlim2',sdhlim2 )
    print('n1:',n1)
    N1=60*n1*1*16*350*10
    print('N1:', N1)
    N2=N1/u1
    print('N2:', N2)
    S=1
    print('S', S)
    Khn1=cKhn(N1)
    print('Khn1',Khn1 )
    Khn2=cKhn(N2)
    print('Khn2',Khn2 )
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
    print('KA', KA)
    # Kvh=1.10
    Kvh=cKv(v)
    print('Kvh:', Kvh)
    Ft1=2*T1/d1t
    print('Ft1', Ft1)
    pd1=KA*Ft1/b
    print('pd1', pd1)
    if pd1>=100:
        Khaef=1.2
    else:
        Khaef=1.4
    print('Khaef',Khaef )
    Khbt=cKhbt(sdd,b)
    # print('Khbt', Khbt)
    print('Khbt:',Khbt )
    # print(Khbt)
    Kh=KA*Kvh*Khaef*Khbt
    print('Kh:', Kh)
    d1h=d1t*(Kh/Kht)**(1/3)
    print('d1h:', d1h)
    mnh=d1h*cos(hudu(bt))/z1
    print('mnh:', mnh)

    Kft=1.3
    print('Kft', Kft)
    print('btb',btb)
    yfxnaefv=yfxnaef/cos(hudu(btb))**2
    print('yfxnaefv',yfxnaefv )
    Yyfxn=0.25+0.75/yfxnaefv
    # print('Yyfxn', Yyfxn)
    print('Yyfxn',Yyfxn )
    Ybt=1-yfxnbt*bt/120
    print('Ybt', Ybt)
    zv1=z1/cos(hudu(bt))**3
    zv2 = z2 / cos(hudu(bt)) ** 3
    print('zv1', zv1)
    print('zv2', zv2)
    Yfa1=cYfa1(zv1)
    Yfa2=cYfa2(zv2)
    Ysa1=cYsa1(zv1)
    Ysa2=cYsa2(zv2)
    print('Yfa1', Yfa1)
    print('Yfa2:', Yfa2)
    print('Ysa1:', Ysa1)
    print('Ysa2:', Ysa2)

    Kfn1=cKfn(N1)
    Kfn2=cKfn(N2)
    sdflim1=500
    sdflim2=320
    S=1.4
    xyf1=Kfn1*sdflim1/S
    xyf2 = Kfn2 * sdflim2 / S
    print('Kfn1',Kfn1 )
    print('Kfn2', Kfn2)
    print('sdflim1', sdflim1)
    print('sdflim2', sdflim2)
    print('S', S)
    print('xyf1',xyf1 )
    print('xyf2', xyf2)

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
    print('d1', d1)
    v=pi*d1*n1/60000
    print('v', v)
    b=sdd*d1
    print('b', b)
    cn=0.25
    print('cn',cn )
    h=(2*ha+cn)*mnt
    print('h', h)
    bh=b/h
    print('bh',bh )
    # Kvf=1.08
    Kvf=cKv(v)
    print('Kvf:',Kvf)
    Ft1=2*T1/d1
    print('Ft1', Ft1)
    pd2=KA*Ft1/b
    print('pd2', pd2)
    if pd2 >= 100:
        Kfaef = 1.2
    else:
        Kfaef = 1.4
    print('Khaef:',Khaef)
    print('Kfaef',Kfaef )
    Kfbt =1.28
    print('Kfbt', Kfbt)
    Kf=KA*Kvf*Kfaef*Kfbt
    # print('Kfaef', Kfaef)
    print('Kf', Kf)

    mnf=mnt*(Kf/Kft)**(1/3)
    print('mnf', mnf)
    d1f=mnf*z1/cos(hudu(bt))
    print('d1f', d1f)
    # standard_m=[1,1.25,1.5,2,2.5,3,4,5,6,8,10,12,16,20,25,32,40,50]
    mn=shortlength_m(mnf)
    print('mn', mn)
    z1y=d1h*cos(hudu(bt))/mn
    print('z1y:',z1y)
    z1=round(z1y)
    print('z1', z1)
    z2y=u1*z1
    print('z2y:', z2y)
    z2 = round(z2y)
    print('z2',z2)


    ayv=(z1+z2)*mn/2/cos(hudu(bt))
    print('ayv:', ayv)
    a=round(ayv)
    print('a', a)
    bt=degrees(acos((z1+z2)*mn/2/a))
    print('bt', bt)

    d1=z1*mn/cos(hudu(bt))
    d2=z2*mn/cos(hudu(bt))
    print('d1', d1)
    print('d2',d2)

    byv=sdd*d1
    print('byv:', byv)
    b=round(byv)
    print('b', b)
    b1=b+5
    b2=b
    print('b1', b1)
    print('b2',b2 )
    xfdyzj=mn*z1/cos(pi*bt/180)
    dfdyzj=mn*z2/cos(pi*bt/180)
    hf=mn*(ha+cn)
    haa=mn*ha
    xcdyzj = xfdyzj + 2 * haa
    dcdyzj = dfdyzj + 2 * haa
    xcgyzj=xfdyzj-2*hf
    dcgyzj = dfdyzj - 2 * hf
    #(z1,z2,mn,aef,bt,a,b1,b2)
    print('模数：{0}，中心距：{1}，小齿数：{2}，大齿数：{3}，小齿宽：{4}，大齿宽：{5}\n'
          '小分度圆直径：{6}，大分度圆直径：{7}，小齿顶圆直径：{8}，大齿顶圆直径：{9}，\n小齿根圆直径：{10}，'
          '大齿根圆直径：{11}，螺旋角：{12}，\n齿顶高系数：{13}，顶隙系数：{14}，'
          '齿顶高：{15}，齿根高：{16}，全齿高：{17}'.format(mn,
                                        a,z1,z2,b1,b2,xfdyzj,dfdyzj,
                                        xcdyzj,dcdyzj,xcgyzj,dcgyzj,bt,ha,cn,haa,hf,haa+hf
                                        ))





    return (z1,z2,mn,aef,bt,a,b1,b2)
# x=sheji(3.2,24,14,1,1,20,n1=970,P1=11)
# z1,z2,mn,aef,bt,a,b1,b2=x
# print(x)
# jiaohe(z1,z2,mn,aef,bt,a,b1,b2,n1=970,P1=11)

#def sheji(i1,z1,bt,sdd,ha,aef,n1=970,P1=11)
#sheji(i1,z11,bt1,sdd1,ha1,aef1,n1,0.5*P1)
'''-------个体--------
<初选参数>
高速级传动比：5.49861076894544，低速级传动比：3.256277040215731，
高速级小齿轮齿数：20.0,低速级小齿轮齿数：30.0，
高速级齿宽系数：1.0，低速级齿宽系数：1.0
高速级齿顶高系数：1.0，低速级齿顶高系数：1.0，
高速级压力角：20.0,低速级压力角：20.0,
高速级螺旋角：10.302336090323761，低速级螺旋角：14.7461742953384
<结果>
高速级：齿数为z11=23.0,z12=126.0,模数mn1为2mm，压力角aef1=20.0度,螺旋角bt1=9.335652473769246度，中心距a1=151.0mm,齿宽b11=52.0mm,b2=47.0mm.
低速级：齿数为z21=31.0,z22=101.0,模数mn2为2.5mm，压力角aef2=20.0度,螺旋角bt2=15.222757034189314度，中心距a2=171.0mm,齿宽b21=85.0mm,b2=80.0mm
指标1：2.4633431085043993，指标2：0.7954405814620085，指标3：0.7874084807809645，指标4：0.0,，指标5：0.01268866325682789,指标6：0
不适应度:4.0588808340042'''
# x=sheji(4.93,39,13.21,1,1,20,n1=720,P1=7.128*0.5)
# z1,z2,mn,aef,bt,a,b1,b2=x
# print(x)
#最优的
x=sheji(5.4986,20,10.3023,1,1,20,n1=720,P1=7.128*0.5)
z1,z2,mn,aef,bt,a,b1,b2=x
print(x)
jiaohe(z1,z2,mn,aef,bt,a,b1,b2,n1=720,P1=7.128*0.5)
#
# x=sheji(17.905/5.4986,30,14.7462,1,1,20,n1=720/5.4986,P1=6.845)
# z1,z2,mn,aef,bt,a,b1,b2=x
# print(x)
# jiaohe(z1,z2,mn,aef,bt,a,b1,b2,n1=720/5.4986,P1=6.845)
#推荐的
# x=sheji(4.92,24,14,1,1,20,n1=720,P1=7.128*0.5)

# z1,z2,mn,aef,bt,a,b1,b2=x
# print(x)
# jiaohe(z1,z2,mn,aef,bt,a,b1,b2,n1=720,P1=7.128*0.5)
#
# x=sheji(17.905/4.92,24,14,1,1,20,n1=720/4.92,P1=6.845)
# z1,z2,mn,aef,bt,a,b1,b2=x
# print(x)
# jiaohe(z1,z2,mn,aef,bt,a,b1,b2,n1=720/5.4986,P1=6.845)