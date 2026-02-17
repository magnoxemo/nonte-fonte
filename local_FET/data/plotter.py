import openmc
import numpy as np
import matplotlib.pyplot as plt

sp_file = "/home/ebny-walid-ahammed/Documents/github/nonte-fonte/local_FET/statepoint.200.h5"
xmin, xmax = -1, 1


def plot_flux(data_frames, xmin, xmax):
    fig, ax1 = plt.subplots()

    for index, data_frame in enumerate(data_frames):
        df = data_frame.get_pandas_dataframe()

        x = np.linspace(xmin[index], xmax[index], 100)
        order = len(df['mean'])
        n = np.arange(order)
        a_n = (2 * n + 1) / 2 * df['mean']
        phi = np.polynomial.Legendre(a_n / (order - 1), domain=(xmin[index], xmax[index]))

        ax1.plot(x, phi(x))
        ax1.set_xlabel('Z position [cm]')
        ax1.set_ylabel('Flux [n/src]')

    plt.show()


def main():
    with openmc.StatePoint(sp_file) as sp:
        df = sp.tallies
    plot_flux(data_frames=[df[1], df[2]], xmin=[0, 30], xmax=[50, 100])


if __name__ == '__main__':
    main()

