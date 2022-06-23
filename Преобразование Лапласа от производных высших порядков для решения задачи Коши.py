from sympy import*
var('s a b')
var('t', positive=True)
a = 5
f = sin(3*t)
b = 7
q = cos(7*t)
# Прямое преобразование a*L{f(t)}:
L1 = a * laplace_transform(f, t, s, noconds=True)
# Прямое преобразование b*L{q(t)}:
L2 = b*laplace_transform(q, t, s, noconds=True)
# Сумма a*L{f(t)}+b*L{q(t)}:
L = factor(L1 + L2)
print (L)
# Прямое преобразование L{a*f(t)+b*q(t)}:
LS = factor(laplace_transform(a*f + b*q, t, s, noconds=True))
print (LS)
print (LS == L)
# Обратное преобразование a* L^-1{f(t)}:
L_1 = a * inverse_laplace_transform(L1/a, s, t)
# Обратное преобразование b* L^-1{q(t)}
L_2 = b * inverse_laplace_transform(L2/b, s, t)
# a* L^-1{f(t)}+b* L^-1{q(t)}:
L_S = L_1 + L_2
print (L_S)
# Обратное преобразование L^-1{a*f(t)+b*q(t)}:
L_1_2 = inverse_laplace_transform(L1 + L2, s, t)
print (L_1_2)
print (L_1_2 == L_S)