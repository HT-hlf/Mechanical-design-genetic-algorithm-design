#encoding=UTF-8

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
reg = LinearRegression()
poly_features = PolynomialFeatures(degree=3, include_bias=False)
X_poly = poly_features.fit_transform([17,18,19,20,21,22,23,24,25,26,27,28,29,30,35,40,45,50,60,70,80,90,100,150,200])
reg.fit(X_poly, [2.97,2.91,2.85,2.76,2.72,2.69,2.65,2.62,2.60,2.57,2.55,2.53,2.52,2.45,2.40,2.35,2.32,2.28,2.24,2.22,2.20,2.20,2.18,2.14,2.12])
print(reg.coef_)
print(reg.intercept_)
y_predict = reg.predict(X_poly)