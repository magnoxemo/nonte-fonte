import matplotlib.pyplot as plt
import numpy as np

co_effs = [
    [
        -0.000216194,
        0.34569,
        4.29324e-05,
        -0.152019,
        7.97447e-05,
        0.0163828,
        -7.11651e-05,
        -0.000822597,
        2.65694e-05,
        -1.95318e-05,
        -1.40579e-06,
        3.54206e-05,
        1.95068e-05,
    ]
]
fig, ax1 = plt.subplots()

x = np.linspace(-1, 1, 1000)


for co_eff in co_effs:
    order = len(co_eff) - 1  # max polynomial degree
    n = np.arange(order + 1)

    # Each coefficient scaled by (2n+1)/2 gives the Legendre series weights
    a_n = (2 * n + 1) * np.array(co_eff)

    # np.polynomial.Legendre(c) evaluates sum_n c[n] * P_n(x)
    phi = np.polynomial.Legendre(a_n, domain=[-1, 1])
    ax1.plot(x, phi(x), label=f"FET order = {order}", linestyle="--")

ax1.plot(x, np.sin(3 * x), label="actual", color="black", linewidth=1, alpha=0.2)
ax1.set_xlabel("x position [cm]")
ax1.set_ylabel("Y")
ax1.legend()
plt.tight_layout()
plt.show()

plt.show()
