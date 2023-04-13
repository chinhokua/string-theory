# Chin Ho Kua, NYIT 1267789, Runge-Kutta method (4th order) implementation
# -------------------------------------------------------------------
# Python program to implement Runge Kutta method for coupled equations
import numpy as np
import matplotlib.pyplot as plt
import math
import tkinter
from matplotlib.widgets import Slider, Button

# additional metric variables
mu = 1 
p = 4
alpha = 2 / ( 5 - p )
z_0 = 1       # different from function z0 
b1 = 1 / 3   

# INITIAL VALUES
x0 = 0.0001     # initial value of primary independent variable
y0 =  x0 * x0 * b1 / 2   # where y(x0) = y0
z0 = x0 * b1      # where z(x0) = z0 [1/3 for A=B=1, might change otherwise]
dx = x0   # step size
x_end = 100  # final accepted value for x before calculation terminates
# a = 1      # the value of a 
# b2 = 1      # the value of b 


# Initial variables, and building the array to place calculated values of every iteration of RK4 of both equations
x = x0
y = y0 
z = z0
x_values = [x0]
y_values = [y0]
z_values = [x0] 

#a wit no z
# a = np.sqrt( mu * np.power(x * x + 1, (p-3)/ (10 - 2*p)) * (x * x + 1 / (alpha * alpha)  ) /(x*x + 1 ) ) 

#a wit z
a = np.sqrt( mu * np.power(z_0,  (p-3) / (5-p))  * np.power(x * x + 1, (p-3)/ (10 - 2*p)) * (x * x + 1 / (alpha * alpha)  ) /(x*x + 1 ) )

#b wit no z
# b2 = np.power( (x*x + 1) , (p-3) / (10-2*p) ) * 1/(alpha * alpha)

#b wit z
b2 = np.power(z_0, (p-3) / (5 - p)) * np.power( (x*x + 1) , (p-3) / (10-2*p) ) * 1/(alpha * alpha)

# derivative functions definition
def dydx(x, y, z): # dy/dx, where y and z are functions of x
    return   a * z * (1-y)
    # return x * z

def dzdx(x, y, z): # dz/dx, where y and z are functions of x
    return  ((a * y * (2-y))/(b2 * x * x))
    # return x / y


def pvalueplot(pstart): #plots p to final values, incrementing by 1.
    step = 1

def recalculate_initial(): # this recalculates the initial vairables when b1 is varied
    global x0, y0, z0, dx, x_end, a, b2
    x0 = 0.0001     # initial value of primary independent variable
    y0 =  x0 * x0 * b1 / 2   # where y(x0) = y0
    z0 = x0 * b1      # where z(x0) = z0 [1/3 for A=B=1, might change otherwise]
    dx = x0   # step size
    x_end = 100  # final accepted value for x before calculation terminates
    # a = 1      # the value of a 
    # b2 = 1      # the value of b

    #a wit z
    a = np.sqrt( mu * np.power(z_0,  (p-3) / (5-p))  * np.power(x * x + 1, (p-3)/ (10 - 2*p)) * (x * x + 1 / (alpha * alpha)  ) /(x*x + 1 ) )

    #b wit no z
    # b2 = np.power( (x*x + 1) , (p-3) / (10-2*p) ) * 1/(alpha * alpha)

    #b wit z
    b2 = np.power(z_0, (p-3) / (5 - p)) * np.power( (x*x + 1) , (p-3) / (10-2*p) ) * 1/(alpha * alpha)  
   
# main RK4 function
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
    # print(x, y, z, k1, m1, k2, m2, k3, m3, k4, m4)
    return x, y, z # outputted so they can be appended as values to a table (as well as update the x value)

def main1():

    while x <= x_end: # while we haven't reached intended final calculation area (x)
        
        # perform one iteration of the RK4
        x, y, z = RungeKuttaCoupled(x, y, z, dx, dydx, dzdx)
        #append outputted values into the value lists so they can be plotted by matplotlib later
    
    # Plotting the graphs 
    print(x_values[-1], y_values[-1], z_values[-1])
    plt.autoscale(enable=True, axis='y', tight=None) # forces a y-axis autofit 
    plt.plot(x_values, y_values, 'bo') 
    plt.plot(x_values, z_values,'ro')
    plt.legend(['F','H'])
    plt.grid(True)
    plt.title("Solution")
    plt.show()
    return x,y,z # returns the final value of the calculation

def main2():
    global b1
    b1step = 0.05
    x = x0
    while x <= x_end:
        x,y,z = RungeKuttaCoupled(x, y, z, dx, dydx, dzdx)

    
# Plotting the graphs 
print(x_values[-1], y_values[-1], z_values[-1])
plt.autoscale(enable=True, axis='y', tight=None) # forces a y-axis autofit 
plt.plot(x_values, y_values, 'bo') 
plt.plot(x_values, z_values,'ro')
plt.legend(['F','H'])
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