from math import *
import numpy as np
from lib.net_function import NetFunction
from lib.enumerable import equal_enumerables
from models.univariate_model import *

# Parameters binding
# theta[0] - b_11
# theta[1] - delta_1
# theta[2] - lambda_1
# theta[3] - a_1

MAX_MATURITY = 30
MATURITY_STEP = 1 / (12 * 32)

class A11(UnivariateModel):
  def __init__(self, without_error_maturity=1.0/12, with_error_maturity=2.0):
    UnivariateModel.__init__(self, 1, without_error_maturity, with_error_maturity)
    self.theta0 = []
    self.theta = []
    self._gam0 = None
    self._gam = None
    self._nodes = None

  @property
  def nodes(self):
    if self._nodes is None:
      values = []
      for i in range(int(MAX_MATURITY / MATURITY_STEP) + 1):
        values.append(i * MATURITY_STEP)
      self._nodes = values

    return self._nodes

  def gamma0(self, theta):
    if not(equal_enumerables(self.theta0, theta)):
      values = [0]
      for i in range(1, len(self.nodes)):
        values.append(values[i - 1] - MATURITY_STEP * theta[3] * self.gamma(theta).y[i - 1])

      self._gam0 = NetFunction(self.nodes, values)
      self.theta0 = theta

    return self._gam0

  def gamma(self, theta):
    if not(equal_enumerables(self.theta, theta)):
      values = [0]
      for i in range(1, len(self.nodes)):
        values.append(values[i - 1] + MATURITY_STEP * ((theta[0] - theta[2]) * values[i - 1] - 0.5 * values[i - 1]**2 + theta[1]))

      self._gam = NetFunction(self.nodes, values)
      self.theta = theta

    return self._gam

  def g_0(self, tau, theta):
    return self.gamma0(theta).call(tau) / tau

  def g(self, tau, theta):
    return self.gamma(theta).call(tau) / tau

  def l_x_euler(self, delta, x, x_0, theta):
    return log(sqrt(1 / (2 * pi * delta * abs(x_0))) * exp(-(x - x_0 - (theta[3] + (theta[0] - theta[2]) * x_0) * delta)**2 / (2 * delta * abs(x_0))))

  def p_x_0(self, delta, x, x_0, theta):
    return x**(2 * theta[3] - 0.5) * x_0**(0.5 - 2 * theta[3]) / sqrt(2 * pi * delta) * exp(- (x - x_0)**2 / (2 * delta) - (theta[2] - theta[0]) * x**2 / 4 + (theta[2] - theta[0]) * x_0**2 / 4)

  def c1(self, x, x_0, theta):
    return - 1 / (24 * x * x_0) * (48 * theta[3]**2 - 48 * theta[3] + 9 + x * x_0 * (theta[2] - theta[0])**2 * (x**2 - 24 * theta[3] / (theta[2] - theta[0])) + x**2 * x_0**2 * (theta[2] - theta[0])**2 + x * x_0**3 * (theta[2] - theta[0])**2)

  def c2(self, x, x_0, theta):
    return 1 / (576 * x**2 * x_0**2) * (9 * (256 * theta[3]**4 - 512 * theta[3]**3 + 224 * theta[3]**2 + 32 * theta[3] - 15) + 6 * x * x_0 * (theta[0] - theta[2])**2 * (x**2 + (24 * theta[3]) / (theta[0] - theta[2])) * (16 * theta[3]**2 - 16 * theta[3] + 3) + x**2 * x_0**2 * (theta[0] - theta[2])**2 * (672 * theta[3]**2 - 48 * theta[3] * (2 - (theta[0] - theta[2]) * x**2) - 6 + (theta[0] - theta[2])**2 * x**4) + 2 * x * x_0**3 * (theta[0] - theta[2])**2 * (48 * theta[3]**2 - 24 * theta[3] * (2 - (theta[0] - theta[2]) * x**2) + 9 + (theta[0] - theta[2])**2 * x**4) + 3 * x**2 * x_0**4 * (theta[0] - theta[2])**4 * (x**2 + (16 * theta[3]) / (theta[0] - theta[2])) + 2 * x**3 * x_0**5 * (theta[0] - theta[2])**4 + x**2 * x_0**6 * (theta[0] - theta[2])**4)

  def ex(self, delta, x, x_0, theta):
    return exp((theta[0] - theta[2]) * delta) * (x_0 + theta[3] / (theta[0] - theta[2])) - theta[3] / (theta[0] - theta[2])

  def vx(self, delta, x, x_0, theta):
    return (exp((theta[0] - theta[2]) * delta) - 1) * (-theta[3] + exp((theta[0] - theta[2]) * delta) * (theta[3] + 2 * x_0 * (theta[0] - theta[2]))) / (2 * (theta[0] - theta[2])**2)

  def is_valid(self, theta):
    return (theta[0] < 0) & (theta[1] >= 0) & (theta[3] >= 0.5)

  def print_params(self, theta):
    print("| b_11 = " + str(theta[0]) + ", delta_1 = " + str(theta[1]) + ", lambda_1 = " + str(theta[2]) + ", a_1 = " + str(theta[3]) + " |")
