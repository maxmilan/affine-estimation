from math import *

# Parameters binding
# θ[0] - b_11
# θ[1] - δ_1
# θ[2] - λ_1

class A02():
  def __init__(self, without_error_maturity=1.0/12, with_error_maturity=2.0):
    UnivariateModel.__init__(self, 0, without_error_maturity, with_error_maturity)

    def Γ(θ):
  τ = [1 / 12, 2]
  return np.array(list(map(lambda x: [- ((math.exp(θ[0]) - math.exp(θ[2])) * x * θ[1]**2) / (θ[0] * θ[2] * (θ[0] - θ[2])), - (1 - x * math.exp(θ[2])) / θ[2]], τ)))

def Γ_ext(θ):
  τ = [4, 6]
  return np.array(list(map(lambda x: [- ((math.exp(θ[0]) - math.exp(θ[2])) * x * θ[1]**2) / (θ[0] * θ[2] * (θ[0] - θ[2])), - (1 - x * math.exp(θ[2])) / θ[2]], τ)))

def Γ0(θ):
  τ = [1 / 12, 2]
  return np.array(list(map(lambda x: - x * (x**2 * (float(1) / 6 * math.exp(2 * θ[0]) + float(1) / 6 * math.exp(2 * θ[2]) - float(1) / 3 * math.exp(θ[0] + θ[2])) * θ[1]**4 + x * (0.5 * math.exp(θ[0]) - 0.5 * math.exp(θ[2])) * θ[0] * θ[1]**2 * θ[2]**2 * θ[3] + θ[0]**3 * θ[2] * (-1 + x * math.exp(θ[2]) - float(1) / 3 * x**2 * math.exp(2 * θ[2]) + (2 - x * math.exp(θ[2])) * θ[2] * θ[4]) + θ[0]**4 * (0.5 - 0.5 * x * math.exp(θ[2]) + float(1) / 6 * x**2 * math.exp(2 * θ[2]) + (-1 + 0.5 * x * math.exp(θ[2])) * θ[2] * θ[4]) + θ[0]**2 * θ[2] * ((0.5 - 0.5 * x * math.exp(θ[2]) + float(1) / 6 * x**2 * math.exp(2 * θ[2])) * θ[2] + x * (-0.5 * math.exp(θ[0]) + 0.5 * math.exp(θ[2])) * θ[1]**2 * θ[3] + (-1 + 0.5 * x * math.exp(θ[2])) * θ[2]**2 * θ[4])) / (θ[0]**2 * θ[2]**2 * (θ[0]**2 - 2 * θ[0] * θ[2] + θ[2]**2)), τ)))

def Γ0_ext(θ):
  τ = [4, 6]
  return np.array(list(map(lambda x: - x * (x**2 * (float(1) / 6 * math.exp(2 * θ[0]) + float(1) / 6 * math.exp(2 * θ[2]) - float(1) / 3 * math.exp(θ[0] + θ[2])) * θ[1]**4 + x * (0.5 * math.exp(θ[0]) - 0.5 * math.exp(θ[2])) * θ[0] * θ[1]**2 * θ[2]**2 * θ[3] + θ[0]**3 * θ[2] * (-1 + x * math.exp(θ[2]) - float(1) / 3 * x**2 * math.exp(2 * θ[2]) + (2 - x * math.exp(θ[2])) * θ[2] * θ[4]) + θ[0]**4 * (0.5 - 0.5 * x * math.exp(θ[2]) + float(1) / 6 * x**2 * math.exp(2 * θ[2]) + (-1 + 0.5 * x * math.exp(θ[2])) * θ[2] * θ[4]) + θ[0]**2 * θ[2] * ((0.5 - 0.5 * x * math.exp(θ[2]) + float(1) / 6 * x**2 * math.exp(2 * θ[2])) * θ[2] + x * (-0.5 * math.exp(θ[0]) + 0.5 * math.exp(θ[2])) * θ[1]**2 * θ[3] + (-1 + 0.5 * x * math.exp(θ[2])) * θ[2]**2 * θ[4])) / (θ[0]**2 * θ[2]**2 * (θ[0]**2 - 2 * θ[0] * θ[2] + θ[2]**2)), τ)))

  def γ_0(self, τ, θ):
    return (θ[1] * ((3 - 4 * exp(τ * θ[0]) + exp(2 * τ * θ[0]) + 2 * τ * θ[0]) * θ[1] - 4 * θ[0] * (1 - exp(τ * θ[0]) + τ * θ[0]) * θ[2])) / (4 * τ * θ[0]**3)

  def γ(self, τ, θ):
    return ((exp(τ * θ[0]) - 1) * θ[1]) / (θ[0] * τ)

  def l_x_true(self, Δ, x, x_0, θ):
    res = sqrt(- θ[0] / (pi * (1 - exp(2 * θ[0] * Δ)))) * exp((θ[0] * (x - θ[2] / θ[0] - (x_0 - θ[2] / θ[0]) * exp(θ[0] * Δ))**2)/(1 - exp(2 * θ[0] * Δ)))
    if res > 1e-20:
      return log(res)
    else:
      return 0

  def l_x_euler(self, Δ, x, x_0, θ):
    return log(sqrt(1 / (2 * pi * Δ)) * exp(-(x - x_0 - (-θ[2] + θ[0] * x_0) * Δ)**2 / (2 * Δ)))

  def p_x_0(self, Δ, x, x_0, θ):
    return 1.0 / sqrt(2 * pi * Δ) * exp(- (x - x_0)**2 / (2 * Δ) + θ[0] * x**2 / 2 - θ[0] * x_0**2 / 2 - θ[2] * x + θ[2] * x_0)

  def c1(self, x, x_0, θ):
    return θ[0] / 6 * (-3 * (θ[2])**2 / θ[0] + 3 * θ[2] * (x + x_0) - 3 - θ[0] * x**2 - θ[0] * x * x_0 - θ[0] * x_0**2)

  def c2(self, x, x_0, θ):
    return θ[0]**2 / 36 * ((9 * θ[2]**4) / θ[0]**2 - 18 * θ[2]**3 / θ[0] * x + 3 * θ[2]**2 / θ[0] * (6 + 5 * θ[0] * x**2) - 6 * θ[2] * x * (3 + θ[0] * x**2) + 3 + 6 * θ[0] * x**2 + θ[0]**2 * x**4 + 2 * θ[0] * x_0 * (x - 3 * θ[2] / θ[0]) * (3 * θ[2]**2 / θ[0] - 3 * θ[2] * x + 3 + θ[0] * x**2) + 3 * θ[0] * x_0**2 * (5 * θ[2]**2 / θ[0] - 4 * θ[2] * x + 2 + θ[0] * x**2) + 2 * θ[0]**2 * x_0**3 * (x - 3 * θ[2] / θ[0]) + θ[0]**2 * x_0**4)

  def ex(self, Δ, x, x_0, θ):
    return x_0 * exp(θ[0] * Δ) - θ[2] / θ[0] * (exp(θ[0] * Δ) - 1)

  def vx(self, Δ, x, x_0, θ):
    return 1 / (2 * θ[0]) * (exp(2 * θ[0] * Δ) - 1)

  def print_params(self, θ):
    print("| b_11 = " + str(θ[0]) + ", δ_1 = " + str(θ[1]) + ", λ_1 = " + str(θ[2]) + " |")

  def is_valid(self, θ):
    return θ[0] < 0
