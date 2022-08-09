# Chin Ho Kua, NYIT 1267789, Euler Method Implementation
# -------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

# requires installation of numpy and matplotlib (which depends on numpy, so install numpy first) using 'pip install numpy' and 
# 'pip install matplotlib' respectively.

def odeEuler(f,y0,t):
    '''Approximate the solution of y'=f(y,t) by Euler's method.

    Parameters
    ----------
    f : function
        Right-hand side of the differential equation y'=f(t,y), y(t_0)=y_0
    y0 : number
        Initial value y(t0)=y0 wher t0 is the entry at index 0 in the array t
    t : array
        1D NumPy array of t values where we approximate y values. Time step
        at each iteration is given by t[n+1] - t[n].

    Returns
    -------
    y : 1D NumPy array
        Approximation y[n] of the solution y(t_n) computed by Euler's method.
    '''
    y = np.zeros(len(t))
    y[0] = y0
    for n in range(0,len(t)-1):
        y[n+1] = y[n] + f(y[n],t[n])*(t[n+1] - t[n])
    return y


# function arguments
t = np.linspace(0,2,21)
y0 = 0
f = lambda y,t: np.sin(t)


# plotting shenanigans
y = odeEuler(f,y0,t)
y_true = np.sin(t)
plt.plot(t,y,'b.-',t,y_true,'r-')
plt.legend(['Euler','True'])
plt.axis([0,2,0,9])
plt.grid(True)
plt.title("Solution of $y'=y , y(0)=1$")
plt.show()