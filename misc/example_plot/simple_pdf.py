import sys
sys.path.append("../")
from misc.plotter_utilities import *



co_efficients = [6.00026, -4.00084, 1.99313, 0.0140851, -0.012394, 0.00601481, 0.0138056]

x = np.linspace(-1, 1, 100)
plt.plot(x,3*x*x-4*x+5,marker='o')

plot_pdf(domain=[[-1,1]],co_efficients=co_efficients)
plt.show()
