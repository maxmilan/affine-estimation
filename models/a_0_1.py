from math import *
import numpy as np
from models.univariate_model import *

# Parameters binding
# theta[0] - b_11
# theta[1] - delra_1
# theta[2] - lambda_1

class A01(UnivariateModel):
  def __init__(self, without_error_maturity=1.0/12, with_error_maturity=2.0):
    UnivariateModel.__init__(self, 0, without_error_maturity, with_error_maturity)

  def g_0(self, tau, theta):
    return (theta[1] * ((3 - 4 * exp(tau * theta[0]) + exp(2 * tau * theta[0]) + 2 * tau * theta[0]) * theta[1] - 4 * theta[0] * (1 - exp(tau * theta[0]) + tau * theta[0]) * theta[2])) / (4 * tau * theta[0]**3)

  def g(self, tau, theta):
    return ((exp(tau * theta[0]) - 1) * theta[1]) / (theta[0] * tau)

  def l_x_true(self, delta, x, x_0, theta):
    res = sqrt(- theta[0] / (pi * (1 - exp(2 * theta[0] * delta)))) * exp((theta[0] * (x - theta[2] / theta[0] - (x_0 - theta[2] / theta[0]) * exp(theta[0] * delta))**2)/(1 - exp(2 * theta[0] * delta)))
    if res > 1e-20:
      return log(res)
    else:
      return 0

  def l_x_euler(self, delta, x, x_0, theta):
    return log(sqrt(1 / (2 * pi * delta)) * exp(-(x - x_0 - (-theta[2] + theta[0] * x_0) * delta)**2 / (2 * delta)))

  def p_x_0(self, delta, x, x_0, theta):
    return 1.0 / sqrt(2 * pi * delta) * exp(- (x - x_0)**2 / (2 * delta) + theta[0] * x**2 / 2 - theta[0] * x_0**2 / 2 - theta[2] * x + theta[2] * x_0)

  def c1(self, x, x_0, theta):
    return theta[0] / 6 * (-3 * (theta[2])**2 / theta[0] + 3 * theta[2] * (x + x_0) - 3 - theta[0] * x**2 - theta[0] * x * x_0 - theta[0] * x_0**2)

  def c2(self, x, x_0, theta):
    return theta[0]**2 / 36 * ((9 * theta[2]**4) / theta[0]**2 - 18 * theta[2]**3 / theta[0] * x + 3 * theta[2]**2 / theta[0] * (6 + 5 * theta[0] * x**2) - 6 * theta[2] * x * (3 + theta[0] * x**2) + 3 + 6 * theta[0] * x**2 + theta[0]**2 * x**4 + 2 * theta[0] * x_0 * (x - 3 * theta[2] / theta[0]) * (3 * theta[2]**2 / theta[0] - 3 * theta[2] * x + 3 + theta[0] * x**2) + 3 * theta[0] * x_0**2 * (5 * theta[2]**2 / theta[0] - 4 * theta[2] * x + 2 + theta[0] * x**2) + 2 * theta[0]**2 * x_0**3 * (x - 3 * theta[2] / theta[0]) + theta[0]**2 * x_0**4)

  def ex(self, delta, x, x_0, theta):
    return x_0 * exp(theta[0] * delta) - theta[2] / theta[0] * (exp(theta[0] * delta) - 1)

  def vx(self, delta, x, x_0, theta):
    return 1 / (2 * theta[0]) * (exp(2 * theta[0] * delta) - 1)

  def is_valid(self, theta):
    return theta[0] < 0

  def print_params(self, theta):
    print("| b_11 = " + str(theta[0]) + ", delta_1 = " + str(theta[1]) + ", lambda_1 = " + str(theta[2]) + " |")
