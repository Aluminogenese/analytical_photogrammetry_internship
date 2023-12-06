# 解析摄影测量实习
## 空间后方交会
利用影响覆盖范围内一定数量的控制点的空间坐标与影响坐标，基于共线方程，反求该影像的的六个外方位元素： $X_S$, $Y_S$, $Z_S$, $\varphi$, $\omega$, $\kappa$ 
### 共线条件方程
$$
\begin{equation}
\begin{split}
x-x_0=-f\dfrac{a_1(X-X_S)+b_1(Y-Y_S)+c_1(Z-Z_S)}{a_3(X-X_S)+b_3(Y-Y_S)+c_3(Z-Z_S)} \\
y-y_0=-f\dfrac{a_2(X-X_S)+b_2(Y-Y_S)+c_2(Z-Z_S)}{a_3(X-X_S)+b_3(Y-Y_S)+c_3(Z-Z_S)}
\end{split}
\end{equation}
$$
### 求解过程
1. 获取已知数据。摄影比例尺m，摄影航高H，内方位元素$x_0$, $y_0$, f,控制点空间坐标$X_t$, $Y_t$, $Z_t$
2. 量测控制点对应像点坐标$x$, $y$
3. 确定未知数初始值$X_S^0=\dfrac{1}{n}\sum_{i=1}^{n}X_{ti}$, $Y_S^0=\dfrac{1}{n}\sum_{i=1}^{n}Y_{ti}$, $\varphi^0=\omega^0=\kappa^0=0$
4. 计算旋转矩阵（以Y轴为主轴的$\varphi$-$\omega$-$\kappa$转角系统
$$
R=\begin{bmatrix}
 cos\varphi cos\kappa -sin\varphi sin\omega sin\kappa & -cos\varphi sin\kappa -sin\varphi sin\omega cos\kappa & -sin\varphi sin\omega \\
 cos\omega sin\kappa  & cos\omega cos\kappa & -sin\omega \\
 sin\varphi cos\kappa +cos\varphi sin\omega sin\kappa & -sin\varphi sin\kappa +cos\varphi sin\omega cos\kappa & cos\varphi cos\omega 
\end{bmatrix}
$$
5. 利用未知数近似值按共线条件方程$(1)$计算像点坐标近似值$(x), (y)$
6. 逐点计算误差方程的系数项和常数项，组成误差方程
$$
A=\begin{bmatrix}
 a_{11} & a_{12} & a_{13} & a_{14} & a_{15} & a_{16} \\
a_{21} & a_{22} & a_{23} & a_{24} & a_{25} & a_{26} 
\end{bmatrix}
$$
$$
\begin{equation}
\begin{split}
&a_{11}=\frac{\partial x}{\partial X_S} =\frac{1}{\bar{Z} }[a_1f+a_3(x-x_0)] \\
&a_{12}=\frac{\partial x}{\partial Y_S} =\frac{1}{\bar{Z} }[b_1f+b_3(x-x_0)] \\
&a_{13}=\frac{\partial x}{\partial Z_S} =\frac{1}{\bar{Z} }[c_1f+c_3(x-x_0)] \\
&a_{14}=\frac{\partial x}{\partial \varphi} =(y-y_0)sin\omega -{\frac{x-x_0}{f}[(x-x_0)cos\kappa -(y-y_0)sin\kappa]+fcos\kappa}cos\omega \\
&a_{15}=\frac{\partial x}{\partial \omega }=-fsin\kappa -\frac{x-x_0}{f}[(x-x_0)sin\kappa +(y-y_0)cos\kappa ] \\
&a_{16}=\frac{\partial x}{\partial \kappa }=(y-y_0) \\
&a_{21}=\frac{\partial y}{\partial X_S} =\frac{1}{\bar{Z} }[a_2f+a_3(y-y_0)] \\
&a_{22}=\frac{\partial y}{\partial Y_S} =\frac{1}{\bar{Z} }[b_2f+b_3(y-y_0)] \\
&a_{23}=\frac{\partial y}{\partial Z_S} =\frac{1}{\bar{Z} }[c_2f+c_3(y-y_0)] \\
&a_{24}=\frac{\partial y}{\partial \varphi} =-(x-x_0)sin\omega -{\frac{y-y_0}{f}[(x-x_0)cos\kappa -(y-y_0)sin\kappa]-fsin\kappa}cos\omega \\
&a_{25}=\frac{\partial y}{\partial \omega }=-fcos\kappa -\frac{y-y_0}{f}[(x-x_0)sin\kappa +(y-y_0)cos\kappa ] \\
&a_{26}=\frac{\partial y}{\partial \kappa }=-(x-x_0)
\end{split}
\end{equation}
$$
## 空间前方交会
利用点投影系数法根据两张像片的外方位元素计算模型点坐标。
## 解析内定向
扫描坐标转为相片坐标。
## 相对定向
恢复两张像片间的相互关系，构建立体模型。
## 绝对定向
将相对定向的模型用控制点纳入大地测量坐标系中，恢复模型的绝对位置。
