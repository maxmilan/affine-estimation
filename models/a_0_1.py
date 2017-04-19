from sympy import *
from math import *

class A01:
  def __init__(self, without_error_maturity=1.0/12, with_error_maturity=2.0):
    self._m = 0
    self._n = 1
    self.b11, self.λ1, self.δ0, self.τ = var("b11 λ1 δ0 τ")
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
      expression = (self.δ0 * ((3 - 4 * exp(self.τ * self.b11) + exp(2 * self.τ * self.b11) + 2 * self.τ * self.b11) * self.δ0 - 4 * self.b11 * (1 - exp(self.τ * self.b11) + self.τ * self.b11) * self.λ1)) / (4 * self.b11**3)
      self._A = expression
    
    return self._A

  @property
  def B(self):
    if self._B is None:
      expression = ((exp(self.τ * self.b11)) * self.δ0) / self.b11
      self._B = expression

    return self._B  

  def A(self, τ, θ):
    return (θ[1] * ((3 - 4 * exp(τ * θ[0]) + exp(2 * τ * θ[0]) + 2 * τ * θ[0]) * θ[1] - 4 * θ[0] * (1 - exp(τ * θ[0]) + τ * θ[0]) * θ[2])) / (4. * θ[0]**3)

  def B(self, τ, θ):
    return ((exp(τ * θ[0])) * θ[1]) / θ[0]

  def Γ0(self, θ):
    # f = lambdify(self.τ, -self.A.subs([(self.b11, θ[0]), (self.λ1, θ[2]), (self.δ0, θ[1])]))
    # return f(self.without_error_maturity)
    return -self.A(self.without_error_maturity, θ)

  def Γ0_ext(self, θ):
    # f = lambdify(self.τ, -self.A.subs([(self.b11, θ[0]), (self.λ1, θ[2]), (self.δ0, θ[1])]))
    # return f(self.with_error_maturity)
    return -self.A(self.with_error_maturity, θ)

  def Γ(self, θ):
    return self.B(self.without_error_maturity, θ)

  def Γ_ext(self, θ):
    return self.B(self.with_error_maturity, θ)


