from scipy import signal
from pylab import *
x1 = linspace(1,1000, num=1000000)

x2 = signal.resample(x1,100)