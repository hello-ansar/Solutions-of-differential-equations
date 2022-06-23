from sympy import *
import time
import numpy as np
import matplotlib.pyplot as plt
start = time.time()
var('t C1 C2 C3 C4')
u = Function("u")(t)
# Запись дифференциального уравнения:
de = Eq(u.diff(t, t, t, t) - 3*u.diff(t, t, t) + 3*u.diff(t, t) - u.diff(t), 4*t*exp(t))
# Решение дифференциального уравнения:
des = dsolve(de, u)
# Начальные условия:
eq1 = des.rhs.subs(t, 0)
eq2 = des.rhs.diff(t).subs(t, 0)
eq3 = des.rhs.diff(t, t).subs(t, 0)
eq4 = des.rhs.diff(t, t, t).subs(t, 0)
# Решение системы алгебраических уравнений для начальных условий:
seq = solve([eq1, eq2-1, eq3-2, eq4-3], C1, C2, C3, C4)
rez = des.rhs.subs([(C1, seq[C1]), (C2, seq[C2]), (C3, seq[C3]), (C4, seq[C4])])


def F(t): return rez
f = lambdify(t, rez, 'numpy')
x = np.linspace(0, 6*np.pi, 100)
stop = time.time()
print ('Время решения уравнения с использованием функции dsolve(): %s s' % round((stop-start), 3))
plt.title('Решение с использованием функции dsolve():\n х (t)=%s\n' % rez, fontsize=11)
plt.grid(True)
plt.xlabel('Time t seconds', fontsize=12)
plt.ylabel('f(t)', fontsize=16)
plt.plot(x, f(x), color='#008000', linewidth=3)
plt.show()