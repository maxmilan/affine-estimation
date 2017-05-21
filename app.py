import numpy as np
from numpy.linalg import inv
from scipy.optimize import minimize
from scipy.stats import norm
from lib.data_storage import *
from lib.matrix import *
from lib.statistics import *
from lib.enumerable import *
from models.a_0_1 import A01
from models.a_1_1 import A11
from models.a_0_2 import A02
from datetime import datetime
from math import *
import time

delta = 1/12

def univariate_likelihood(theta, model, observations, likelihood_method):
  yield_errors = []
  n = observations.length()

  if not(model.is_valid(theta)):
    return 0

  for i in range(0, n):
    gt = observations[i, 1, model.n]
    gt_ext = observations[i, model.n + 1, 2 * model.n]
    xt = (gt + model.G_0(theta)) / model.G(theta)
    gt_calculated = -model.G_0_ext(theta) + model.G_ext(theta) * xt
    yield_errors.append(gt_calculated - gt_ext)

  sigma = np.std(yield_errors)
  error_distribution = norm(0, sigma)

  joint_errors_likelihood = inject(lambda memo, x: memo + log(x), 0.0, list(map(lambda x: error_distribution.pdf(x), yield_errors)))

  observations_likelihood = 0.0
  for i in range(1, n):
    g = observations[i, 1, model.n]
    g_0 = observations[i - 1, 1, model.n]
    observations_likelihood += model.l_g(likelihood_method, delta, g, g_0, theta)

  return -observations_likelihood / n


def multivariate_likelihood(theta, model, observations, likelihood_method):
  yield_errors = []
  n = observations.length()

  if not(model.is_valid(theta)):
    return 0

  for i in range(0, observations.length()):
    gt = observations[i, 1, 2]
    gt_ext = observations[i, 3, 4]
    xt = inv(G(theta).transpose()).dot(subtract_vectors(gt, G0(theta)))
    gt_calculated = G0_ext(theta) + G_ext(theta).dot(xt)
    yield_errors.append(subtract(gt_calculated, gt_ext))

  sigma = multi_std(yield_errors)
  error_distribution = multivariate_normal(mean = np.zeros(2), cov = np.diag(sigma))

  joint_errors_likelihood = inject(lambda memo, x: memo * x, 1.0, list(map(lambda x: error_distribution.pdf(x), yield_errors)))

  observations_likelihood = 0.0
  for i in range(1, n):
    g = observations[i, 1, model.n]
    g_0 = observations[i - 1, 1, model.n]
    observations_likelihood += model.l_g(likelihood_method, delta, g, g_0, theta)

  return -observations_likelihood / n

def run_estimation(model, initial_value, data, estimation_method, true_estimates = None):
  start_time = time.time()
  optimization_result = minimize(univariate_likelihood, initial_value, args=(model, data, estimation_method,), method='nelder-mead', options= { 'disp': False, 'maxiter': 1000 })
  end_time = time.time()
  model.print_params(optimization_result.x)

  if not(true_estimates is None):
    deviations = []
    for i in range(len(true_estimates)):
      deviations.append(abs((optimization_result.x[i] - true_estimates[i])/(true_estimates[i])))
    print("Deviation:")
    print(deviations)

  print("Duration: " + "{:6.3f}".format(end_time - start_time) + "s, iterations: " + str(optimization_result.nit) + ", iteration duration: " + str() + "{:6.3f}".format((end_time - start_time) / optimization_result.nit))

  return optimization_result.x

data_storage = DataStorage("data/yields.xls", "ZEROYLD", datetime(1965, 1, 1))

y_s_1_factor = import_data(data_storage, 1)
print("Number of observations: " + str(y_s_1_factor.length()) + "\n")

print("-------------------------- A01 Model --------------------------------")
theta0 = [-0.03, 0.05, 1]
model01 = A01()
print("__________________________ TRUE likelihood __________________________")
true_estimates = run_estimation(model01, theta0, y_s_1_factor, "true")

print("___________________________ Euler method ____________________________")
run_estimation(model01, theta0, y_s_1_factor, "euler", true_estimates)
print("____________________________ QML method _____________________________")
run_estimation(model01, theta0, y_s_1_factor, "qml", true_estimates)
print("______________________ Approximations (k = 1) _______________________")
run_estimation(model01, theta0, y_s_1_factor, 1, true_estimates)
print("______________________ Approximations (k = 2) _______________________")
run_estimation(model01, theta0, y_s_1_factor, 2, true_estimates)


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
