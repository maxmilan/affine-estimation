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

  def Γ_0(self, θ):
    return self.γ_0(self.without_error_maturity, θ)

  def Γ_0_ext(self, θ):
    return self.γ_0(self.with_error_maturity, θ)

  def Γ(self, θ):
    return self.γ(self.without_error_maturity, θ)

  def Γ_ext(self, θ):
    return self.γ(self.with_error_maturity, θ)

  def l_g(self, k, Δ, g, g_0, θ):
    x = (g + self.Γ_0(θ)) / self.Γ(θ)
    x_0 = (g_0 + self.Γ_0(θ)) / self.Γ(θ)

    if isinstance(k, int):
      return log(abs(1 / self.Γ(θ))) + self.l_x(k, Δ, x, x_0, θ)
    elif k == "true":
      return log(abs(1 / self.Γ(θ))) + self.l_x_true(Δ, x, x_0, θ)
    elif k == "euler":
      return log(abs(1 / self.Γ(θ))) + self.l_x_euler(Δ, x, x_0, θ)
    elif k == "qml":
      return log(abs(1 / self.Γ(θ))) + self.l_x_qml(Δ, x, x_0, θ)
    else:
      raise BaseException("Wrong argument k = " + str(k) + " in l_g!")

  def c(self, k, x, x_0, θ):
    def coefficient_not_found(x, x_0, θ):
      raise BaseException("No coefficient " + str(k) + " found!")

    function_name = 'c' + str(k)
    function = getattr(self, function_name, coefficient_not_found)

    return function(x, x_0, θ)

  def l_x(self, k, Δ, x, x_0, θ):
    sum = 1

    for i in range(1, k + 1):
      sum = sum + self.c(i, x, x_0, θ) * Δ**i / factorial(i)

    return log(self.p_x_0(Δ, x, x_0, θ) * sum)

  def l_x_qml(self, Δ, x, x_0, θ):
    value = sqrt(1 / (2 * pi * self.vx(Δ, x, x_0, θ))) * exp(- (x - self.ex(Δ, x, x_0, θ))**2 / (2 * self.vx(Δ, x, x_0, θ)))

    return log(value)
