#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import matplotlib.pyplot as plt
import numpy as np


#获取文本内容
f = open("../log/lambda.txt",encoding = "utf-8")
lambdas_txt = f.read()
f.close()
#分割字符串，删除''
lambdas = lambdas_txt.split(',')
while '' in lambdas:
    lambdas.remove('')
#字符串转float
lambdas = [float(x) for x in lambdas]
print(lambdas)

x = range(len(lambdas))

# 设置画布
plt.figure()


# 标题
plt.title("lambdas")
#数据
plt.plot(x,lambdas, '-p', color='grey',
        marker = 'o',
        markersize=8, linewidth=2,
        markerfacecolor='red',
        markeredgecolor='grey',
        markeredgewidth=2,label='lambdas')
# 坐标描述
plt.xlabel('index')
plt.ylabel('lambda')

# 设置数字标签
for a,b in zip(x,lambdas):
    plt.text(a,b,'%.2f'%b,ha='center', va='bottom', fontsize=10)


plt.legend()
plt.show()
