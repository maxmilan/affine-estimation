from math import *
from models.multivariate_model import *
import code

# Parameters binding
# theta[0] - b_11
# theta[1] - delta_1
# theta[2] - lambda_1
# theta[3] - b_21
# theta[4] - b_22
# theta[5] - delta_2
# theta[6] - lambda_2

class A02(MultivariateModel):
  def __init__(self, without_error_maturities=[1.0/12, 2.0], with_error_maturities=[4.0, 6.0]):
    MultivariateModel.__init__(self, 0, 2, without_error_maturities, with_error_maturities)

  def g_0(self, tau, theta):
    return (theta[4]**5 * ((0.75 - exp(tau * theta[0]) + 0.25 * exp(2 * tau* theta[0])) * theta[4]**4 * theta[1]**2 + (0.5 - exp(tau * theta[0]) + 0.5 * exp(2 * tau * theta[0])) * theta[3]**2 * theta[4]**2 * theta[1] * theta[5] +(-0.25 + 0.25 * exp(2 * tau * theta[0])) * theta[3]**4 * theta[5]**2) + theta[0]**2 * theta[4]**3 * ((-0.75 + exp(tau * theta[0]) - 0.25 * exp(2 * tau * theta[0])) * theta[4]**4 * theta[1]**2 + (-2 + 2 * exp(tau * theta[0]) - exp(2 * tau * theta[0]) - exp(tau * theta[0]) + 2 * exp(tau * (theta[0] + theta[4]))) * theta[3]**2 * theta[4]**2 * theta[1] * theta[5] + (-0.5 - 0.5 * exp(2 * tau * theta[0]) + exp(tau * (theta[0] + theta[4]))) * theta[3]**4 * theta[5]**2 - tau * theta[4]**6 * theta[1] * theta[2] + (2 - exp(tau * theta[0]) - exp(tau * theta[4])) * theta[3]**2 * theta[4]**3 * theta[5] * theta[2] + theta[4]**5 * theta[1]* (-tau * theta[1] + (2 - 2 * exp(tau * theta[0])) * theta[2])) + theta[0] * theta[4]**4 * ((-1.5 + 2 * exp(tau * theta[0]) - 0.5 * exp(2 * tau * theta[0])) * theta[4]**4 * theta[1]**2 + (-0.5 + exp(tau * theta[0])- 0.5 * exp(2 * tau * theta[0]) + exp(tau * theta[4]) - exp(tau * (theta[0] + theta[4]))) * theta[3]**2 * theta[4]**2 * theta[1] * theta[5] + (0.75 + 0.25 * exp(2 * tau * theta[4]) - exp(tau * (theta[0] + theta[4]))) * theta[3]**4 * theta[5]**2 + (exp(tau * theta[0]) - 1) * theta[3]**2 * theta[4]**3 * theta[5] * theta[2] + theta[4]**5 * theta[1] * (0.5 * tau * theta[1] + (exp(tau * theta[0]) - 1) * theta[2])) + theta[0]**9 * theta[5] * ((0.75 - exp(tau * theta[4]) + 0.25 * exp(2 * tau * theta[4]) + 0.5 * tau * theta[4]) * theta[5] + theta[4] * (-1 + exp(tau * theta[4]) - tau * theta[4]) * theta[6]) + theta[0]**8 * theta[4] * ((-1.5 + 2 * exp(tau * theta[4]) - 0.5 * exp(2 * tau * theta[4]) - tau * theta[4]) * theta[5]**2 -tau * theta[4]**2 * theta[1] * theta[2] + theta[4] * (2 - 2 * exp(tau * theta[4]) + 2 * tau * theta[4]) * theta[5] * theta[6]) + theta[0]**4 * theta[4] * ((0.5 - exp(tau * theta[0]) + 0.5 * exp(2 * tau * theta[0]) + 2 * exp(tau * theta[4]) - 2 * exp(tau * (theta[0] + theta[4]))) * theta[3]**2 * theta[4]**2 * theta[1] * theta[5] + (0.75 + 0.25 * exp(2 * tau * theta[0]) - exp(tau * (theta[0] + theta[4]))) * theta[3]**4 * theta[5]**2 + theta[4]**4 * ((-0.75 + exp(tau * theta[0]) - 0.25 * exp(2 * tau * theta[0])) * theta[1]**2 + (-1.5 + 2 * exp(tau * theta[4]) - 0.5 * exp(2 * tau * theta[4])) * theta[5]**2) + (-4 + 2 * exp(tau * theta[0]) + 2 * exp(tau * theta[4])) * theta[3]**2 * theta[4]**3 * theta[5] * theta[2] + theta[4]**6 * (tau * theta[1] * theta[2] + 2 * tau * theta[5] * theta[6]) + theta[4]**5 * (2 * tau * theta[1]**2 + (-4 + 4 * exp(tau * theta[0])) * theta[1] * theta[2] + theta[5] * (-tau * theta[5] + (2 - 2 * exp(tau * theta[4])) * theta[6]))) + theta[0]**5 * ((-2.5 + exp(tau * theta[0]) - 0.5 * exp(2 * tau * theta[0]) + exp(tau * theta[4]) + exp(tau * (theta[0] + theta[4]))) * theta[3]**2 * theta[4]**2 * theta[1] * theta[5] + (-0.25 + 0.25 * exp(2 * tau * theta[4])) * theta[3]**4 * theta[5]**2 + theta[4]**4 * ((-1.5 + 2 * exp(tau * theta[0]) - 0.5 * exp(2 * tau * theta[0])) * theta[1]**2 + (-0.75 + exp(tau * theta[4]) - 0.25 * exp(2 * tau * theta[4])) * theta[5]**2) + (1 + exp(tau * theta[0]) - 2 * exp(tau * theta[4])) * theta[3]**2 * theta[4]**3 * theta[5] * theta[2] + theta[4]**6 * (-4 * tau * theta[1] * theta[2] + tau * theta[5] * theta[6]) + theta[4]**5 * (-0.5 * tau * theta[1]**2 + (exp(tau * theta[0]) - 1) * theta[1] * theta[2] + theta[5] * (-0.5 * tau * theta[5] + (exp(tau * theta[4]) - 1) * theta[6]))) + theta[0]**7 * theta[4] * ((-0.75 + exp(tau * theta[4]) - 0.25 * exp(2 * tau * theta[4])) * theta[4] * theta[5]**2 + (exp(tau * theta[4]) - 1) * theta[3]**2 * theta[5] * theta[2] + theta[4]**3 * (2 * tau * theta[1] * theta[2] + tau * theta[5] * theta[6]) + theta[4]**2 * (0.5 * tau * theta[1]**2 + (exp(tau * theta[0]) - 1) * theta[1] * theta[2] + theta[5] * (-0.5 * tau * theta[5] + (exp(tau * theta[4]) - 1) * theta[6]))) + theta[0]**3 * theta[4]**2 * ((3 - 3 * exp(tau * theta[0]) + exp(2 * tau * theta[0]) - 2 * exp(tau * theta[4])) * theta[3]**2 * theta[4]**2 * theta[1] * theta[5] + (-0.5 - 0.5 * exp(2 * tau * theta[4]) + exp(tau * (theta[0] + theta[4]))) * theta[3]**4 * theta[5]**2 + theta[4]**4 * ((3 - 4 * exp(tau * theta[0]) + exp(2 * tau * theta[0])) * theta[1]**2 + (0.75 - exp(tau * theta[4]) + 0.25 * exp(2 * tau * theta[4])) * theta[5]**2) + (1 - 2 * exp(tau * theta[0]) + exp(tau * theta[4])) * theta[3]**2 * theta[4]**3 * theta[5] * theta[2] + theta[4]**6 * (2 * tau * theta[1] * theta[2] - tau * theta[5] * theta[6]) + theta[4]**5 * (-0.5 * tau * theta[1]**2 + (exp(tau * theta[0]) - 1) * theta[1] * theta[2] + theta[5] * (0.5 * tau * theta[5] + (exp(tau * theta[4]) - 1) * theta[6]))) + theta[0]**6 * theta[4] * ((exp(tau * theta[4]) - 1) * theta[3]**2 * theta[1] * theta[5] + theta[4]**2 * ((0.75 - exp(tau * theta[0]) + 0.25 * exp(2 * tau * theta[0])) * theta[1]**2 + (3 - 4 * exp(tau * theta[4]) + exp(2 * tau * theta[4])) * theta[5]**2) + (2 - exp(tau * theta[0]) - exp(tau * theta[4]) * theta[3]**2 * theta[4] * theta[5] * theta[2] + theta[4]**4 * (tau * theta[1] * theta[2] - 4 * tau * theta[5] * theta[6]) + theta[4]**3 * (-tau * theta[1]**2 + (2 - 2 * exp(tau * theta[0])) * theta[1] * theta[2] + theta[5] * (2 * tau * theta[5] + (-4 + 4 * exp(tau * theta[4])) * theta[6]))))) / (tau * theta[0]**3 * theta[4]**3 * (theta[0] + theta[4]) * (theta[0]**2 - 2 * theta[0] * theta[4] + theta[4]**2) * (theta[0]**3 - theta[0]**2 * theta[4] - theta[0] * theta[4]**2 + theta[4]**3))

  def g(self, tau, theta):
    # code.interact(local=dict(globals(), **locals()))
    return [(theta[1] * (exp(tau * theta[0]) - 1)) / (theta[0] * tau) + ((exp(tau * theta[0]) - exp(tau * theta[4])) * theta[3]**2 * theta[5]) / (tau * theta[0] * theta[4] * (theta[4] - theta[0])), (theta[5] * (exp(tau * theta[4]) - 1)) / (theta[4] * tau)]

  def cinv(self, x, x_0, theta):
    return -0.5 * (x[0] - x_0[0])**2 - 0.5 * (x[1] - x_0[1])**2

  def c0(self, x, x_0, theta):
    return 0.5 * theta[0] * (x[0] - x_0[0]) * (x[0] + x_0[0]) + 0.5 * (x[1] - x_0[1]) * (theta[3] * (x[0] + x_0[0]) + theta[4] * (x[1]  +  x_0[1]))

  def c1(self, x, x_0, theta):
    return -0.5 * (theta[0] + theta[4]) - theta[0]**2 / 6 * (x[0]**2 + x[0] * x_0[0] + x_0[0]**2) - theta[4]**2 / 6 * (x[1]**2 + x[1] * x_0[1] + x_0[1]**2) - theta[3]**2 / 8 * (x[0] + x_0[0])**2 + theta[3]**2 / 24 * (x[1] - x_0[1])**2 - theta[0] * theta[4] / 6 * (2 * x[0] * x[1] + x[0] * x_0[1] + x_0[0] * x[1] + 2 * x_0[0] * x_0[1])

  def print_params(self, theta):
    print("| b_11 = " + str(theta[0]) + ", delta_1 = " + str(theta[1]) + ", lambda_1 = " + str(theta[2]) + ", lambda_1 = " + str(theta[3]) + ", b_21 = " + str(theta[4]) + ", b_22 = " + str(theta[5]) + ", delta_2 = " + str(theta[6]) + " |")

  def is_valid(self, theta):
    return (theta[0] < 0) & (theta[4] < 0)
