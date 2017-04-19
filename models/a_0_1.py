from sympy import *
from math import *

class A01:
  def __init__(self, without_error_maturity=1.0/12, with_error_maturity=2.0):
    self._m = 0
    self._n = 1
    self.b_11, self.λ_1, self.δ_1, self.τ = var("b_11 λ_1 δ_1 τ")
    self._A = None
    self._B = None
    self.without_error_maturity = without_error_maturity
    self.with_error_maturity = with_error_maturity

  @property
  def m(self):
    return self._m

  @property
  def n(self):
    return self._n

  @property
  def A(self):
    if self._A is None:
      expression = (self.δ_1 * ((3 - 4 * exp(self.τ * self.b_11) + exp(2 * self.τ * self.b_11) + 2 * self.τ * self.b_11) * self.δ_1 - 4 * self.b_11 * (1 - exp(self.τ * self.b_11) + self.τ * self.b_11) * self.λ_1)) / (4 * self.b_11**3)
      self._A = expression

    return self._A

  @property
  def B(self):
    if self._B is None:
      expression = ((exp(self.τ * self.b_11)) * self.δ_1) / self.b_11
      self._B = expression

    return self._B

  def A(self, τ, θ):
    return (θ[1] * ((3 - 4 * exp(τ * θ[0]) + exp(2 * τ * θ[0]) + 2 * τ * θ[0]) * θ[1] - 4 * θ[0] * (1 - exp(τ * θ[0]) + τ * θ[0]) * θ[2])) / (4 * θ[0]**3)

  def B(self, τ, θ):
    return ((exp(τ * θ[0])) * θ[1]) / θ[0]

  def Γ0(self, θ):
    # f = lambdify(self.τ, -self.A.subs([(self.b_11, θ[0]), (self.λ_1, θ[2]), (self.δ_1, θ[1])]))
    # return f(self.without_error_maturity)
    return -self.A(self.without_error_maturity, θ)

  def Γ0_ext(self, θ):
    # f = lambdify(self.τ, -self.A.subs([(self.b_11, θ[0]), (self.λ_1, θ[2]), (self.δ_1, θ[1])]))
    # return f(self.with_error_maturity)
    return -self.A(self.with_error_maturity, θ)

  def Γ(self, θ):
    return self.B(self.without_error_maturity, θ)

  def Γ_ext(self, θ):
    return self.B(self.with_error_maturity, θ)

  def l_x(self, k, Δ, x, x_0, θ):
    sum = 1

    for i in range(1, k + 1):
      sum = sum + self.c(i, x, x_0, θ) * Δ**i / factorial(i)

    return self.l_x_0(Δ, x, x_0, θ) + log(sum)

  def l_x_0(self, Δ, x, x_0, θ):
    return -log(sqrt(2 * pi * Δ)) - (x - x_0)**2 / (2 * Δ) + θ[0] * x**2 / 2 - θ[0] * x_0**2 / 2

  def c(self, k, x, x_0, θ):
    def coefficient_not_found(x, x_0, θ):
      raise BaseException("No coefficient " + str(k) + " found!")

    function_name = 'c' + str(k)
    function = getattr(self, function_name, coefficient_not_found)

    return function(x, x_0, θ)

  def c1(self, x, x_0, θ):
    return θ[0] / 6 * (3 + θ[0] * x**2 + θ[0] * x * x_0 + θ[0] * x_0**2)

  def c2(self, x, x_0, θ):
    return θ[0]**2 / 36 * (3 + 6 * θ[0] * x**2 + θ[0]**2 * x**4 + 2 * θ[0] * x * x_0 * (3 + θ[0] * x**2) + 3 * θ[0] * x_0**2 * (2 + θ[0] * x**2) - 2 * θ[0]**2 * x * x_0**3 + θ[0]**2 * x_0**4)

  def l_g(self, k, Δ, g, g_0, θ):
    x = (g - self.Γ0(θ)) / self.Γ(θ)
    x_0 = (g_0 - self.Γ0(θ)) / self.Γ(θ)
    return log(abs(1 / self.Γ(θ))) + self.l_x(k, Δ, x, x_0, θ)
