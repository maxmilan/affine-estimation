from math import *

class UnivariateModel:
  def __init__(self, m, without_error_maturity, with_error_maturity):
    self._n = 1
    self.without_error_maturity = without_error_maturity
    self.with_error_maturity = with_error_maturity

  @property
  def m(self):
    return self._m

  @property
  def n(self):
    return self._n

  def G_0(self, theta):
    return self.g_0(self.without_error_maturity, theta)

  def G_0_ext(self, theta):
    return self.g_0(self.with_error_maturity, theta)

  def G(self, theta):
    return self.g(self.without_error_maturity, theta)

  def G_ext(self, theta):
    return self.g(self.with_error_maturity, theta)

  def g(self, g_0, tau_0, theta, tau):
    x = (g_0 + self.g_0(tau_0, theta)) / self.g(tau_0, theta)

    return -self.g_0(tau, theta) + self.g(tau, theta) * x

  def l_g(self, k, delta, g, g_0, theta):
    x = (g + self.G_0(theta)) / self.G(theta)
    x_0 = (g_0 + self.G_0(theta)) / self.G(theta)

    if isinstance(k, int):
      return log(abs(1 / self.G(theta))) + self.l_x(k, delta, x, x_0, theta)
    elif k == "true":
      return log(abs(1 / self.G(theta))) + self.l_x_true(delta, x, x_0, theta)
    elif k == "euler":
      return log(abs(1 / self.G(theta))) + self.l_x_euler(delta, x, x_0, theta)
    elif k == "qml":
      return log(abs(1 / self.G(theta))) + self.l_x_qml(delta, x, x_0, theta)
    else:
      raise BaseException("Wrong argument k = " + str(k) + " in l_g!")

  def c(self, k, x, x_0, theta):
    def coefficient_not_found(x, x_0, theta):
      raise BaseException("No coefficient " + str(k) + " found!")

    function_name = 'c' + str(k)
    function = getattr(self, function_name, coefficient_not_found)
    return function(x, x_0, theta)

  def l_x(self, k, delta, x, x_0, theta):
    sum = 1
    for i in range(1, k + 1):
      sum = sum + self.c(i, x, x_0, theta) * delta**i / factorial(i)
    return log(self.p_x_0(delta, x, x_0, theta) * sum)

  def l_x_qml(self, delta, x, x_0, theta):
    value = sqrt(1 / (2 * pi * self.vx(delta, x, x_0, theta))) * exp(- (x - self.ex(delta, x, x_0, theta))**2 / (2 * self.vx(delta, x, x_0, theta)))
    return log(value)
