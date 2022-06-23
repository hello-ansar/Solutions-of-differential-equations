from sympy import *
import time
import numpy as np
import matplotlib.pyplot as plt
start = time.time()
var('s')
var('t', positive=True)
var('X', cls=Function)
# Начальные условия:
x0 = 0
x01 = 0
x02 = 0
x03 = 0
# Запись левой части дифференциального уравнения:
g = 4*t*exp(t)
# Прямое преобразование Лапласа:
Lg = laplace_transform(g, t, s, noconds=True)
d4 = s**4*X(s) - s**3*x0 - s**2*x01 - s*x02 - x03
d3 = s**3*X(s) - s**2*x0 - s*x01 - x02
d2 = s**2*X(s) - s*x0 - x01
d1 = s*X(s) - x0
d0 = X(s)
# Запись правой части дифференциального уравнения:
d = factor(d4 - 3*d3 + 3*d2 - d1)
de = Eq(d, Lg)
# Решение алгебраического уравнения:
rez = solve(de, X(s))[0]
# Обратное преобразование Лапласа:
soln = collect(inverse_laplace_transform(rez, s, t), t)
f = lambdify(t, soln, 'numpy')
x = np.linspace(0, 6*np.pi, 100)
stop = time.time()
print ('Время решения уравнения с использованием преобразования Лапласа: %s s' % round((stop-start), 3))
plt.title('Решение с использованием преобразования Лапласа:\n х (t)=%s\n' % soln, fontsize=11)
plt.grid(True)
plt.xlabel('t', fontsize=12)
plt.ylabel('x(t)', fontsize=12)
plt.plot(x, f(x), 'g', linewidth=2)
plt.show()