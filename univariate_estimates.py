import xlrd
import numpy as np
import math
from numpy.linalg import det
from numpy.linalg import inv
from scipy.optimize import minimize
from scipy.stats import norm
from sympy import *
from datetime import datetime
from lib.yield_series import YieldSeries
from lib.matrix import *
from lib.statistics import *
from lib.enumerable import *
from models.a_0_1 import A01
from models.a_1_1 import A11
import code
import time

FILENAME = "data/yields.xls"

book = xlrd.open_workbook(FILENAME)
sheet_index = book.sheet_names().index('ZEROYLD')
sheet = book.sheet_by_index(sheet_index)
delta = 1/12

def prepare_data():
  first_row = sheet.row_values(0)
  columns_indices = range(1, len(first_row))
  data = []

  for row_index in range(sheet.nrows):
    row_values = sheet.row_values(row_index)
    if row_values[0] and datetime.strptime(row_values[0], "%m/%Y") >= datetime(1965, 1, 1):
      item = {}

      for i in range(len(columns_indices)):
        if not row_values[columns_indices[i]]:
          value = None
        else:
          value = float(row_values[columns_indices[i]]) / 100
        item[first_row[columns_indices[i]]] = value

      data.append(item)

  return data

def likelihood(θ, model, observations, likelihood_method):
  yield_errors = []
  n = observations.length()
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

  # return -(joint_errors_likelihood + observations_likelihood) / n
  return -observations_likelihood / n

# theta0 = [-0.185793266805, 0.0264870240577, -0.866720304618]
# theta1 = [-0.0999041150336, -0.106014406172, 0.242671357612]
# -model01.γ_0(2, theta0) + model01.γ(2, theta0) * xt0
# yields = y_s.table[0]
# gt = yields[1]
# xt0 = (gt + model01.γ_0(1/12, theta0)) / model01.γ(1/12, theta0)
# xt1 = (gt + model01.γ_0(1/12, theta1)) / model01.γ(1/12, theta1)

# code.interact(local=dict(globals(), **locals()))

y_s = YieldSeries(table = prepare_data(), nfactors = 1)
print("Number of observations: " + str(y_s.length()) + "\n")

print("-------------------------- A01 Model --------------------------------")
theta0 = [-0.03, 0.05, 1]
model01 = A01()

print("__________________________ TRUE likelihood __________________________")
start_time = time.time()
a01_true = minimize(likelihood, theta0, args=(model01, y_s, "true",), method='nelder-mead', options= { 'xtol': 1e-6, 'disp': False, 'maxiter': 1000 })
end_time = time.time()
model01.print_params(a01_true.x)
print("Duration: " + "{:6.3f}".format(end_time - start_time) + "s, iterations: " + str(a01_true.nit))
print("___________________________ Euler method ____________________________")
start_time = time.time()
a01_euler = minimize(likelihood, theta0, args=(model01, y_s, "euler",), method='nelder-mead', options= { 'xtol': 1e-6, 'disp': False, 'maxiter': 1000 })
end_time = time.time()
model01.print_params(a01_euler.x)
print("Duration: " + "{:6.3f}".format(end_time - start_time) + "s, iterations: " + str(a01_euler.nit))
print("____________________________ QML method _____________________________")
start_time = time.time()
a01_qml = minimize(likelihood, theta0, args=(model01, y_s, "qml",), method='nelder-mead', options= { 'xtol': 1e-6, 'disp': False, 'maxiter': 1000 })
end_time = time.time()
model01.print_params(a01_qml.x)
print("Duration: " + "{:6.3f}".format(end_time - start_time) + "s, iterations: " + str(a01_qml.nit))
print("______________________ Approximations (k = 1) _______________________")
start_time = time.time()
a01_approx_1 = minimize(likelihood, theta0, args=(model01, y_s, 1,), method='nelder-mead', options= { 'xtol': 1e-6, 'disp': False, 'maxiter': 1000 })
end_time = time.time()
model01.print_params(a01_approx_1.x)
print("Duration: " + "{:6.3f}".format(end_time - start_time) + "s, iterations: " + str(a01_approx_1.nit))
print("______________________ Approximations (k = 2) _______________________")
start_time = time.time()
a01_approx_2 = minimize(likelihood, theta0, args=(model01, y_s, 2,), method='nelder-mead', options= { 'xtol': 1e-6, 'disp': False, 'maxiter': 1000 })
end_time = time.time()
model01.print_params(a01_approx_2.x)
print("Duration: " + "{:6.3f}".format(end_time - start_time) + "s, iterations: " + str(a01_approx_2.nit))

# theta0 = [-0.03, 0.05, 1, 0.5]
# model11 = A11()
# a01_approx_1 = minimize(likelihood, theta0, args=(model11, y_s, 1,), method='nelder-mead', options= { 'xtol': 1e-6, 'disp': True, 'maxiter': 1000 })
# print(a01_approx_1.x)
