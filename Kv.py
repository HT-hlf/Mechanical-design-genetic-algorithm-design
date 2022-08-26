#encoding=UTF-8
import matplotlib.pyplot as plt
# Cv1=0.32
# def Cv2_value(yfxny=2):
#     if 1<yfxny<=2:
#         return 0.34
#     else:
#         return 0.57/(yfxny-1.56)
# def Cv3_value(yfxny=2):
#     if 1<yfxny<=2:
#         return 0.23
#     else:
#         return 0.096/(yfxny-1.56)
#第一种计算公式
# def Kvv(n1,C=7,Cy=23):
#     Bn = (n1 / 1200) ** 2
#     Bc = ((C - 5) / 6) ** 2
#     Cv1 = 0.32
#     Cv2 = Cv2_value()
#     Cv3 = Cv3_value()
#     C1 = Cv1 + Cv2
#     C2 = Cv3
#     By=Cy/24
#     Kv = 1 + Bn*(C1 * Bc + C2 * By)
#     return Kv

C=7
B=0.25*(C-5.0)**0.667
A=50+56*(1.0-B)
Vmax=(A+(14-C))**2/200
X=[0.001*i for i in range(1000*round(Vmax))]
#第二种计算公式
# Y=(A/(A+(200*v)**0.5))
# Bn=(n1/1200)**2
# Bc=((C-5)/6)**2
# Cv1=0.32
# Cv2=Cv2_value()
# Cv3=Cv3_value()
# C1=Cv1+Cv2
# C2=Cv3
#Kv=1+Bn(C1*Bc+C2*By)
Y=[(A/(A+(200*v)**0.5))**(-B) for v in X]
plt.plot(X,Y, 'k-', linewidth=1)
# v=25
# Kv=(A/(A+(200*v)**0.5))**(-B)
# v1=30
# print((A/(A+(200*v)**0.5))**(-B))
# # 设置中文乱码问题
plt.rcParams['font.sans-serif'] = ['SimHei']
# plt.text(v, Kv, 'C=7', ha='left', va='center', fontsize=10.5)
# 设置图标标题，并在坐标轴上添加标签
plt.title('精度等级为7的动载系数图', fontsize=24)
plt.xlabel('v/(m/s)', fontsize=14)
plt.ylabel('Kv', fontsize=14)
plt.axis([0,50,1.0,1.8])

# 保存图
# 调用savefig将拟合曲线保存为表名
plt.savefig('动载系数')

plt.show()