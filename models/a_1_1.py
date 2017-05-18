from math import *
import numpy as np
from lib.net_function import NetFunction
from lib.enumerable import equal_enumerables
from models.univariate_model import *

# Parameters binding
# θ[0] - b_11
# θ[1] - δ_1
# θ[2] - λ_1
# θ[3] - a_1

MAX_MATURITY = 30
MATURITY_STEP = 1 / (12 * 32)

class A11(UnivariateModel):
  def __init__(self, without_error_maturity=1.0/12, with_error_maturity=2.0):
    UnivariateModel.__init__(self, 1, without_error_maturity, with_error_maturity)
    self.theta0 = []
    self.theta = []
    self._g0 = None
    self._g = None
    self._nodes = None

  @property
  def nodes(self):
    if self._nodes is None:
      values = []
      for i in range(int(MAX_MATURITY / MATURITY_STEP) + 1):
        values.append(i * MATURITY_STEP)
      self._nodes = values

    return self._nodes

  def gamma0(self, θ):
    if not(equal_enumerables(self.theta0, θ)):
      values = [0]
      for i in range(1, len(self.nodes)):
        values.append(values[i - 1] - MATURITY_STEP * θ[3] * self.gamma(θ).y[i - 1])

      self._g0 = NetFunction(self.nodes, values)
      self.theta0 = θ

    return self._g0

  def gamma(self, θ):
    if not(equal_enumerables(self.theta, θ)):
      values = [0]
      for i in range(1, len(self.nodes)):
        values.append(values[i - 1] + MATURITY_STEP * ((θ[0] - θ[2]) * values[i - 1] - 0.5 * values[i - 1]**2 + θ[1]))

      self._g = NetFunction(self.nodes, values)
      self.theta = θ

    return self._g

  def γ_0(self, τ, θ):
    return self.gamma0(θ).call(τ) / τ

  def γ(self, τ, θ):
    return self.gamma(θ).call(τ) / τ

  def l_x_euler(self, Δ, x, x_0, θ):
    return log(sqrt(1 / (2 * pi * Δ * abs(x_0))) * exp(-(x - x_0 - (θ[3] + (θ[0] - θ[2]) * x_0) * Δ)**2 / (2 * Δ * abs(x_0))))

  def p_x_0(self, Δ, x, x_0, θ):
    return x**(2 * θ[3] - 0.5) * x_0**(0.5 - 2 * θ[3]) / sqrt(2 * pi * Δ) * exp(- (x - x_0)**2 / (2 * Δ) - (θ[2] - θ[0]) * x**2 / 4 + (θ[2] - θ[0]) * x_0**2 / 4)

  def c1(self, x, x_0, θ):
    return - 1 / (24 * x * x_0) * (48 * θ[3]**2 - 48 * θ[3] + 9 + x * x_0 * (θ[2] - θ[0])**2 * (x**2 - 24 * θ[3] / (θ[2] - θ[0])) + x**2 * x_0**2 * (θ[2] - θ[0])**2 + x * x_0**3 * (θ[2] - θ[0])**2)

  # def c2(self, x, x_0, θ):
  #   return θ[0]**2 / 36 * ((9 * θ[2]**4) / θ[0]**2 - 18 * θ[2]**3 / θ[0] * x + 3 * θ[2]**2 / θ[0] * (6 + 5 * θ[0] * x**2) - 6 * θ[2] * x * (3 + θ[0] * x**2) + 3 + 6 * θ[0] * x**2 + θ[0]**2 * x**4 + 2 * θ[0] * x_0 * (x - 3 * θ[2] / θ[0]) * (3 * θ[2]**2 / θ[0] - 3 * θ[2] * x + 3 + θ[0] * x**2) + 3 * θ[0] * x_0**2 * (5 * θ[2]**2 / θ[0] - 4 * θ[2] * x + 2 + θ[0] * x**2) + 2 * θ[0]**2 * x_0**3 * (x - 3 * θ[2] / θ[0]) + θ[0]**2 * x_0**4)

  def ex(self, Δ, x, x_0, θ):
    return exp((θ[0] - θ[2]) * Δ) * (x_0 + θ[3] / (θ[0] - θ[2])) - θ[3] / (θ[0] - θ[2])

  def vx(self, Δ, x, x_0, θ):
    return (exp((θ[0] - θ[2]) * Δ) - 1) * (-θ[3] + exp((θ[0] - θ[2]) * Δ) * (θ[3] + 2 * x_0 * (θ[0] - θ[2]))) / (2 * (θ[0] - θ[2])**2)

  def is_valid(self):
    return (θ[0] < 0) & (θ[1] >= 0) & (θ[3] >= 0.5)
