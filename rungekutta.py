# Chin Ho Kua, NYIT 1267789, Runge-Kutta method (4th order) implementation
# -------------------------------------------------------------------
# Python program to implement Runge Kutta method


# A sample differential equation "dy / dx = (x - y)/2"
# place differential equation after return on line 9
 
# Finds value of y for a given x using step size h
# and initial value y0 at x0.
def rungeKutta(f, x0, y0, x, h):
    # Count number of iterations using step size or
    # step height h
    n = (int)((x - x0)/h)
    # Iterate for number of iterations
    y = y0
    for i in range(1, n + 1):
        # four factors involved in Runge Kutta 4th Order
        k1 = h * f(x0, y)
        k2 = h * f(x0 + 0.5 * h, y + 0.5 * k1)
        k3 = h * f(x0 + 0.5 * h, y + 0.5 * k2)
        k4 = h * f(x0 + h, y + k3)
 
        # Update next value of y
        y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
 
        # Update next value of x
        x0 = x0 + h
    return y
 
# Driver method
f = lambda x, y: (x-y)/2 # insert DERIVATIVE function here (this equals f = dy/dx = (x-y/2) at present)
x0 = 0 # boundary x-value
y = 1  # boundary y-value
x = 2  # instruct program to find y at this x-value using Runge Kutta 
h = 0.2 # step size - NOTE: If [(x-x0) mod h != 0] then the terms summed will be mildly skewed 
print ('The value of y at x is:', rungeKutta(f, x0, y, x, h))