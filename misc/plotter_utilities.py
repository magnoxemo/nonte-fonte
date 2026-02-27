import numpy as np
from scipy.special import eval_legendre
import matplotlib.pyplot as plt
from matplotlib import cm


def evaluate_fet_1D(co_efficients,x_position):
    soln = 0
    for i in range(len(co_efficients)):
        soln+=co_efficients[i]*eval_legendre(i, x_position)

    return soln


def evaluate_fet_2D(co_efficients,x_position, y_position):

    val = 0
    for i in range(len(co_efficients[0])):
        for j in range(len(co_efficients[1])):
            val += (
                co_efficients[i][j]
                * eval_legendre(i, x_position)
                * eval_legendre(j, y_position)
            )

    return val


def plot_pdf( domain, co_efficients, number_of_sampled_points=1000):

    if len(domain) == 1:
        x = np.linspace(domain[0][0], domain[0][1], number_of_sampled_points)
        solution = np.zeros_like(x)

        for index, pos in enumerate(x):
            solution[index] = evaluate_fet_1D(co_efficients, pos)

        plt.plot(x, solution, label="FET 1D")
        plt.xlabel("x axis")
        plt.grid(True)
        plt.legend()

    if len(domain) == 2:

        x = np.linspace(domain[0][0], domain[0][1], number_of_sampled_points)
        y = np.linspace(domain[1][0], domain[1][1], number_of_sampled_points)

        X, Y = np.meshgrid(x, y)
        Z = np.zeros_like(X)

        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                Z[i, j] = evaluate_fet_2D(co_efficients, x[i],y[j])

        plt.contourf(X, Y, Z, levels=50)
        plt.xlabel("x axis")
        plt.ylabel("y axis")
        plt.colorbar(label="PDF value")
        plt.title("FET 2D PDF")


if __name__ =="__main__":
    pass