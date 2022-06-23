# Загрузка необходимых модулей:
from sympy import *
import time
import matplotlib.pyplot as plt
import numpy as np

start = time.time()
# Объявляем символьные переменные:
var('s Kp Ti Kd Td')

# Накладываем ограничение на символьную переменную времени:
var('t', positive=True)
Kp = 2
Ti = 2
Kd = 4
Td = Rational(1, 2)

# Передаточная функция ПИД-регулятора с оператором s:
fp = (1 + (Kd * Td * s) / (1 + Td*s)) * Kp * (1 + 1/(Ti*s)) * (1/s)

# Переходная характеристика ПИД-регулятора,
# получаемая методом обратного преобразования Лапласа:
ht = inverse_laplace_transform(fp, s, t)
Kd = 0

# Передаточная функция ПИ-регулятора (Kd = 0) с оператором s:
fpp = (1 + (Kd * Td * s) / (1 + Td*s)) * Kp * (1 + 1/(Ti*s)) * (1/s)

# Переходная характеристика ПИ-регулятора,
# получаемая методом обратного преобразования Лапласа:
htt = inverse_laplace_transform(fpp, s, t)
stop = time.time()
print ('Время на обратное визуальное преобразование Лапласа: %s s' % N((stop-start), 3))

# Переходим из символьной области в численную:
f = lambdify(t, ht, 'numpy')
F = lambdify(t, htt, 'numpy')
tt = np.arange(0.01, 20, 0.05)

# Построение графика:
plt.title('Переходные характеристики регуляторов \n с передаточными функциями: \n ПИД - W(s)=%s \n ПИ - W(s)=%s' % (fp, fpp))
plt.plot(tt, f(tt), color='r', linewidth=2, label='ПИД-регулятор: h(t)=%s' % ht)
plt.plot(tt, F(tt), color='b', linewidth=2, label='ПИ-регулятор: h(t)=%s' % htt)
plt.grid(True)
plt.legend(loc='best')
plt.show()
