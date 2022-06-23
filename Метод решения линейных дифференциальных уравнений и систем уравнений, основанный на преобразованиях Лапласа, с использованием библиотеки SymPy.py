from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import sympy
var('s')
var('t', positive=True)
var('X', cls=Function)
print ("Версия библиотеки sympy – %s" % (sympy.__version__))
# Приведенное начальное положение материальной точки заданной массы:
x0 = Rational(6, 5)
# Приведенная начальная скорость материальной точки заданной массы:
x01 = Rational(1, 1)
g = sin(3*t)
# Прямое преобразование Лапласа:
Lg = laplace_transform(g, t, s, noconds=True)
d2 = s**2*X(s) - s*x0 - x01
d0 = X(s)
d = d2 + 4*d0
de = Eq(d, Lg)
# Решение алгебраического уравнения:
rez = solve(de, X(s))[0]
# Обратное преобразование Лапласа:
soln = inverse_laplace_transform(rez, s, t)
f = lambdify(t, soln, "numpy")
x = np.linspace(0, 6*np.pi, 100)
plt.title('Функция, дающая положение материальной точки \n заданной массы:\n х (t)=%s' % soln)
plt.grid(True)
plt.xlabel('t', fontsize=12)
plt.ylabel('x(t)', fontsize=12)
plt.plot(x, f(x), 'g', linewidth=2)
plt.show()