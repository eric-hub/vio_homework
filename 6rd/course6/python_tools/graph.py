#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import matplotlib.pyplot as plt
import numpy as np


# y_ = [232534,58093.2,25801.3,14503.3,9275.77,6437.13,4726.13,3616.01,2855.19,2311.17]
# y = np.array(y_)
# x = range(1,11)

y_ = [7812.66,8114.69,9087.38,11899.2,11319.1,14503.3,21219,29305.9,26469.7]
y = np.array(y_)
x = range(2,11)
# 设置画布
plt.figure()


# 标题
plt.title("proportional change")
#数据
plt.plot(x,y, '-p', color='grey',
        marker = 'o',
        markersize=8, linewidth=2,
        markerfacecolor='red',
        markeredgecolor='grey',
        markeredgewidth=2)
# 坐标描述
plt.xlabel('number')
plt.ylabel('sigema3/sigma4')

# 设置数字标签
for a,b in zip(x,y):
    plt.text(a,b,'%.2f'%b,ha='center', va='bottom', fontsize=10)


plt.legend()
plt.show()
