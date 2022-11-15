# Chin Ho Kua, NYIT 1267789, Runge-Kutta method (4th order) implementation
# -------------------------------------------------------------------
# Python program to implement Runge Kutta method for coupled equations
import numpy as np
import matplotlib.pyplot as plt
import math
import tkinter
from matplotlib.widgets import Slider, Button

# INITIAL VALUES
x0 = 0     # initial value of primary independent variable
y0 = 0     # where y(x0) = y0
z0 = 0     # where z(x0) = z0
dx = 0.1   # step size
x_end = 30 # final accepted value for x before calculation terminates
factor = 1 # THIS IS A/B!!!!!

# Initial variables, and building the array to place calculated values of every iteration of RK4 of both equations
x = x0
y = y0 
z = z0
x_values = [x0]
y_values = [y0]
z_values = [x0] 

# derivative functions definition
def dydx(x, y, z): # dy/dx, where y and z are functions of x
    return x + math.pow(0.01, y)*y

def dzdx(x, y, z): # dz/dx, where y and z are functions of x
    return 47 - x + (math.pow(0.02, z)*z)

# main RK4 functiom
def RungeKuttaCoupled(x, y, z, dx, dydx, dzdx):
    
    # calculating RK4 variables as before, but this time alternating between the two equations
    # k values are for RK4 of dydx, and m values are for RK4 of dzdx
    k1 = dx * dydx(x, y, z) 
    m1 = dx * dzdx(x, y, z)
    k2 = dx * dydx(x + dx/2., y + k1/2., z + m1/2.)
    m2 = dx * dzdx(x + dx/2., y + k1/2., z + m1/2.)
    k3 = dx * dydx(x + dx/2., y + k2/2., z + m2/2.)
    m3 = dx * dzdx(x + dx/2., y + k2/2., z + m2/2.)
    k4 = dx * dydx(x + dx, y + k3, z + m3)
    m4 = dx * dzdx(x + dx, y + k3, z + m3)

    # summative calculation
    y = y + 1./6.*(k1+2*k2+2*k3+k4)
    z = z + 1./6.*(m1+2*m2+2*m3+m4)
    x = x + dx
    return x, y, z # outputted so they can be appended as values to a table (as well as update the x value)

while x <= x_end: # while we haven't reached intended final calculation area (x)
    
    # perform one iteration of the RK4
    x, y, z = RungeKuttaCoupled(x, y, z, dx, dydx, dzdx)
    
    #append outputted values into the value lists so they can be plotted by matplotlib later
    x_values.append(x)
    y_values.append(y)
    z_values.append(z)

# Plotting the graphs 
plt.autoscale(enable=True, axis='y', tight=None) # forces a y-axis autofit 
plt.plot(x_values, y_values, 'bo') 
plt.plot(x_values, z_values,'ro')
plt.legend(['RK4 dy/dx','dz/dx'])
plt.grid(True)
plt.title("Solution")
plt.show()

# NOT FUNCTIONAL YET
# little fun section for rendering latex in python using matplotlib
# -----------------------------------------------------------------------------
# import matplotlib.pyplot as plt
# a = '\\frac{a}{b}'  #notice escaped slash
# plt.plot()
# plt.text(0.5, 0.5,'$%s$'%a)
# plt.show()