import numpy as np
from numpy.linalg import inv
from scipy.optimize import minimize
from scipy.stats import norm
from lib.yield_series import *
from lib.matrix import *
from lib.statistics import *
from lib.enumerable import *
from models.a_0_1 import A01
from models.a_1_1 import A11
from models.a_0_2 import A02
from datetime import datetime
import time

delta = 1/12

def univariate_likelihood(θ, model, observations, likelihood_method):
  yield_errors = []
  n = observations.length()

  if not(model.is_valid(θ)):
    return 0

  # for i in range(0, n):
  #   gt = observations[i, 1, model.n]
  #   gt_ext = observations[i, model.n + 1, 2 * model.n]
  #   xt = (gt + model.Γ_0(θ)) / model.Γ(θ)
  #   gt_calculated = -model.Γ_0_ext(θ) + model.Γ_ext(θ) * xt
  #   yield_errors.append(gt_calculated - gt_ext)

  # σ = np.std(yield_errors)
  # error_distribution = norm(0, σ)

  # joint_errors_likelihood = inject(lambda memo, x: memo + log(x), 0.0, list(map(lambda x: error_distribution.pdf(x), yield_errors)))

  observations_likelihood = 0.0
  for i in range(1, n):
    g = observations[i, 1, model.n]
    g_0 = observations[i - 1, 1, model.n]
    observations_likelihood += model.l_g(likelihood_method, delta, g, g_0, θ)

  return -observations_likelihood / n


def multivariate_likelihood(θ, model, observations, likelihood_method):
  yield_errors = []
  n = observations.length()

  if not(model.is_valid(θ)):
    return 0

  # for i in range(0, y_s.length()):
  #   gt = y_s[i, 1, 2]
  #   gt_ext = y_s[i, 3, 4]
  #   xt = inv(Γ(θ).transpose()).dot(subtract_vectors(gt, Γ0(θ)))
  #   gt_calculated = Γ0_ext(θ) + Γ_ext(θ).dot(xt)
  #   yield_errors.append(subtract(gt_calculated, gt_ext))

  # σ = multi_std(yield_errors)
  # error_distribution = multivariate_normal(mean = np.zeros(2), cov = np.diag(σ))

  # joint_errors_likelihood = inject(lambda memo, x: memo * x, 1.0, list(map(lambda x: error_distribution.pdf(x), yield_errors)))

  observations_likelihood = 0.0
  for i in range(1, n):
    g = observations[i, 1, model.n]
    g_0 = observations[i - 1, 1, model.n]
    observations_likelihood += model.l_g(likelihood_method, delta, g, g_0, θ)

  return -observations_likelihood / n

def run_estimation(model, initial_value, data, estimation_method):
  start_time = time.time()
  optimization_result = minimize(univariate_likelihood, initial_value, args=(model, data, estimation_method,), method='nelder-mead', options= { 'disp': False, 'maxiter': 1000 })
  end_time = time.time()
  model.print_params(optimization_result.x)
  print("Duration: " + "{:6.3f}".format(end_time - start_time) + "s, iterations: " + str(optimization_result.nit))

data_storage = DataStorage("data/yields.xls", "ZEROYLD", datetime(1965, 1, 1))

y_s_1_factor = import_data(data_storage, 1)
print("Number of observations: " + str(y_s_1_factor.length()) + "\n")

print("-------------------------- A01 Model --------------------------------")
theta0 = [-0.03, 0.05, 1]
model01 = A01()
print("__________________________ TRUE likelihood __________________________")
run_estimation(model01, theta0, y_s_1_factor, "true")
print("___________________________ Euler method ____________________________")
run_estimation(model01, theta0, y_s_1_factor, "euler")
print("____________________________ QML method _____________________________")
run_estimation(model01, theta0, y_s_1_factor, "qml")
print("______________________ Approximations (k = 1) _______________________")
run_estimation(model01, theta0, y_s_1_factor, 1)
print("______________________ Approximations (k = 2) _______________________")
run_estimation(model01, theta0, y_s_1_factor, 2)


print("-------------------------- A11 Model --------------------------------")
theta0 = [-0.03, 0.05, 1, 0.5]
model11 = A11()
print("___________________________ Euler method ____________________________")
run_estimation(model11, theta0, y_s_1_factor, "euler")
print("____________________________ QML method _____________________________")
run_estimation(model11, theta0, y_s_1_factor, "qml")
print("______________________ Approximations (k = 1) _______________________")
run_estimation(model11, theta0, y_s_1_factor, 1)
print("______________________ Approximations (k = 2) _______________________")
run_estimation(model11, theta0, y_s_1_factor, 2)

y_s_2_factors = import_data(data_storage, 2)
