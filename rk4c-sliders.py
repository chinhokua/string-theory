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
a, b = 1 # THIS IS A/B!!!!!

# Initial variables, and building the array to place calculated values of every iteration of RK4 of both equations
x = x0
y = y0 
z = z0
x_values = [x0]
y_values = [y0]
z_values = [x0] 

# derivative functions definition
def dydx(x, y, z, a, b): # dy/dx, where y and z are functions of x
    return 

def dzdx(x, y, z, a, b): # dz/dx, where y and z are functions of x
    return 47 - x + (math.pow(0.02, z)*dydx(x, y, z, a, b))

# main RK4 functiom
def RungeKuttaCoupled(x, y, z, dx, dydx, dzdx, a, b):
    
    # calculating RK4 variables as before, but this time alternating between the two equations
    # k values are for RK4 of dydx, and m values are for RK4 of dzdx
    k1 = dx * dydx(x, y, z, a, b) 
    m1 = dx * dzdx(x, y, z, a, b)
    k2 = dx * dydx(x + dx/2., y + k1/2., z + m1/2., a, b)
    m2 = dx * dzdx(x + dx/2., y + k1/2., z + m1/2., a, b)
    k3 = dx * dydx(x + dx/2., y + k2/2., z + m2/2., a, b)
    m3 = dx * dzdx(x + dx/2., y + k2/2., z + m2/2., a, b)
    k4 = dx * dydx(x + dx, y + k3, z + m3, a, b)
    m4 = dx * dzdx(x + dx, y + k3, z + m3, a, b)

    # summative calculation
    y = y + 1./6.*(k1+2*k2+2*k3+k4)
    z = z + 1./6.*(m1+2*m2+2*m3+m4)
    x = x + dx
    return x, y, z # outputted so they can be appended as values to a table (as well as update the x value)

while x <= x_end: # while we haven't reached intended final calculation area (x)
    
    # perform one iteration of the RK4
    x, y, z = RungeKuttaCoupled(x, y, z, dx, dydx, dzdx, a, b)
    
    #append outputted values into the value lists so they can be plotted by matplotlib later
    x_values.append(x)
    y_values.append(y)
    z_values.append(z)

fig, ax = plt.subplots()
line, = ax.plot(x_values, y_values, 'b')
line2, = ax.plot(x_values, z_values, 'g')
plt.title("RK4 of dy/dx and dz/dx")
ax.legend(["dy/dx", "dz/dx"])

# adjust the main plot to make room for the sliders
fig.subplots_adjust(left=0.25)

# display settings 
plt.autoscale(enable=True, axis='y', tight=None) # forces a y-axis autofit 
plt.grid(True)

# Make a vertically oriented slider to control the amplitude
axamp = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
ab_slider = Slider(
    ax=axamp,
    label="A/B",
    valmin=0,
    valmax=30,
    valinit=1,
    orientation="vertical"
)

def update(val):
    while x <= x_end: # while we haven't reached intended final calculation area (x)
        # perform one iteration of the RK4
        x, y, z = RungeKuttaCoupled(x, y, z, dx, dydx, dzdx, a, b)
        
        #append outputted values into the value lists so they can be plotted by matplotlib later
        x_values.append(x)
        y_values.append(y)
        z_values.append(z)
        fig.canvas.draw_idle()

# making a reset button to reset the sliders
resetax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')
def reset(event):
    ab_slider.reset()
button.on_clicked(reset)
plt.show()

# NOT FUNCTIONAL YET
# little fun section for rendering latex in python using matplotlib
# -----------------------------------------------------------------------------
# import matplotlib.pyplot as plt
# a = '\\frac{a}{b}'  #notice escaped slash
# plt.plot()
# plt.text(0.5, 0.5,'$%s$'%a)
# plt.show()