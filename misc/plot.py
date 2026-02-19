import matplotlib.pyplot as plt
import numpy as np

co_effs = [[6.0004, -4.00024, 1.99825, 0.0014865, -0.00060914, 0.00148858, 0.000890668, 0.00132551, -0.0031603]]
fig, ax1 = plt.subplots()

x = np.linspace(-1, 1, 1000)


for co_eff in co_effs:
    order = len(co_eff) - 1  # max polynomial degree
    n = np.arange(order + 1)
    a_n = np.array(co_eff)
    # np.polynomial.Legendre(c) evaluates sum_n c[n] * P_n(x)
    phi = np.polynomial.Legendre(a_n, domain=[-1, 1])
    ax1.plot(x, phi(x), label=f"FET order = {order}", linestyle="--")

ax1.plot(x, (3*x*x-4*x+5) , label="actual", color="black", linewidth=1, alpha=0.8)
ax1.set_xlabel("x position [cm]")
ax1.set_ylabel("Y")
ax1.legend()
plt.tight_layout()
plt.show()

plt.show()
