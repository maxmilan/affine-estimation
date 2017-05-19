from math import *
import numpy as np
from numpy.linalg import det
from numpy.linalg import inv
from lib.matrix import *
import code

class MultivariateModel:
  def __init__(self, m, n, without_error_maturities, with_error_maturities):
    self._n = n
    self._m = m
    self.without_error_maturities = without_error_maturities
    self.with_error_maturities = with_error_maturities

  @property
  def m(self):
    return self._m

  @property
  def n(self):
    return self._n

  def Γ_0(self, θ):
    return np.array(list(map(lambda x: self.γ_0(x, θ), self.without_error_maturities)))

  def Γ_0_ext(self, θ):
    return np.array(list(map(lambda x: self.γ_0(x, θ), self.with_error_maturities)))

  def Γ(self, θ):
    return np.array(list(map(lambda x: self.γ(x, θ), self.without_error_maturities))).transpose()

  def Γ_ext(self, θ):
    return np.array(list(map(lambda x: self.γ(x, θ), self.with_error_maturities))).transpose()

  def l_g(self, k, Δ, g, g_0, θ):
    x = inv(self.Γ(θ).transpose()).dot(sum_vectors(g, self.Γ_0(θ)))
    x_0 = inv(self.Γ(θ).transpose()).dot(sum_vectors(g_0, self.Γ_0(θ)))

    if isinstance(k, int):
      return log(abs(det(inv(self.Γ(θ).transpose())))) + self.l_x(k, Δ, x, x_0, θ)
    elif k == "true":
      return log(abs(det(inv(self.Γ(θ).transpose())))) + self.l_x_true(Δ, x, x_0, θ)
    elif k == "euler":
      return log(abs(det(inv(self.Γ(θ).transpose())))) + self.l_x_euler(Δ, x, x_0, θ)
    elif k == "qml":
      return log(abs(det(inv(self.Γ(θ).transpose())))) + self.l_x_qml(Δ, x, x_0, θ)
    else:
      raise BaseException("Wrong argument k = " + str(k) + " in l_g!")

  def c(self, k, x, x_0, θ):
    def coefficient_not_found(x, x_0, θ):
      raise BaseException("No coefficient " + str(k) + " found!")

    function_name = 'c' + str(k)
    function = getattr(self, function_name, coefficient_not_found)

    return function(x, x_0, θ)

  def l_x(self, k, Δ, x, x_0, θ):
    sum = -(self.n / float(2)) * log(2 * pi * Δ) + self.c("inv", x, x_0, θ) / Δ

    for i in range(0, k + 1):
      sum = sum + self.c(i, x, x_0, θ) * Δ**i / factorial(i)

    return sum

  def l_x_qml(self, Δ, x, x_0, θ):
    value = sqrt(1 / (2 * pi * self.vx(Δ, x, x_0, θ))) * exp(- (x - self.ex(Δ, x, x_0, θ))**2 / (2 * self.vx(Δ, x, x_0, θ)))

    return log(value)
