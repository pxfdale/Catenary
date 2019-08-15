# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import math

""" 参数 """
""" 
悬链线左端点：(0, 0)
悬链线右端点：(1, h)
悬链线长度: L
等分段数: N
""" 
h = 1
L = 5
N = 50

# 其他参数
a = 1 # 悬链线右端点横坐标
Delta = L / N # 每一段长度

def f(alpha, beta):
    s = 0
    for k in range(1, N+1):
        s = s + 1/np.sqrt(1 + ((0.5 + N - k)*alpha + beta)**2)
    s = s - a/Delta
    return s

delta = 1.e-6

def fx(x, y):
    return (f(x + delta, y) - f(x, y))/delta

def fy(x, y):
    return (f(x, y + delta) - f(x, y))/delta

def g(alpha, beta):
    s = 0
    for k in range(1, N+1):
        s = s + ((0.5 + N - k)*alpha + beta)/np.sqrt(1 + ((0.5 + N - k)*alpha + beta)**2)
    s = s + (h)/Delta
    return s

def gx(x, y):
    return (g(x + delta, y) - g(x, y))/delta

def gy(x, y):
    return (g(x, y+delta) - g(x, y))/delta

def newton(x0, y0):
    count = 0
    iterationMax = 20
    eps = 1.0e-16
    vector = np.zeros(2)
    J = np.zeros((2,2))
    vector[0] = x0
    vector[1] = y0
    fg = np.zeros(2)
    while (count < iterationMax):
        count = count + 1
        old = vector
        x = vector[0]
        y = vector[1]
        J[0,0] = fx(x, y)
        J[0,1] = fy(x, y)
        J[1,0] = gx(x, y)
        J[1,1] = gy(x, y)
        fg[0] = f(x, y)
        fg[1] = g(x, y)
        det = np.linalg.det(J)
        if (abs(det) < 1.e-16):
            break
        vector = vector - np.linalg.inv(J).dot(fg)
        diff = vector - old
        error = diff.dot(diff)
        #print "count = ", count, ", error = ", error, ", det = ", det
        if (error < eps):
            break
    return vector

def createCurve(v):
    x = []
    y = []
    theta = []
    for k in range(1, N+1):
        theta.append(np.arctan((0.5 + N - k)*v[0] + v[1]))
    sx = 0
    sy = 0
    x.append(sx)
    y.append(sy)
    for i in range(len(theta)):
        sx = sx + Delta*np.cos(theta[i])
        sy = sy + Delta*np.sin(theta[i])
        x.append(sx)
        y.append(-(sy))
    return x, y

if __name__=='__main__':
    v = newton(0.6, -5)
    x, y = createCurve(v)
    plt.plot(x, y, "b-o")
    # y_star = x
    # N = len(x)
    # for i in range(N):
    # y_star[i] = math.cosh(x[i] - 1/2) - math.cosh(1/2)
    # plt.plot(x, y_star, "o")
    plt.grid()
    plt.xlim(min(x) - 0.05, max(x) + 0.05)
    plt.ylim(min(y) - 0.05, max(y) + 0.05)
    plt.show()
