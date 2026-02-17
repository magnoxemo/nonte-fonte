import matplotlib.pyplot as plt
import numpy as np

co_effs = [23.5603, 8.43377, -1.73564, -0.249877, 0.0317756, 0.00620425, 0.00122843, 0.00107179, 0.000440832, -0.000711827]


fig, ax1 = plt.subplots()

x = np.linspace(-1, 1, 100)
order = len(co_effs)
n = np.arange(order)
a_n = (2 * n + 1) / 2 * co_effs
phi = np.polynomial.Legendre(a_n / (order - 1), domain=(-1, 1))

ax1.plot(x, phi(x), label = 'FET')
ax1.plot(x, x*x-3*x, label='actual')
ax1.set_xlabel('x position [cm]')
ax1.set_ylabel('Y')
ax1.legend()

plt.show()

