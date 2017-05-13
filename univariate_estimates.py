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
    if row_values[0] and datetime.strptime(row_values[0], "%m/%Y") >= datetime(1972, 1, 1):
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

  # joint_errors_likelihood = inject(lambda memo, x: memo + math.log(x), 0.0, list(map(lambda x: error_distribution.pdf(x), yield_errors)))
  observations_likelihood = 0.0
  for i in range(1, n):
    g = observations[i, 1, model.n]
    g_0 = observations[i - 1, 1, model.n]
    observations_likelihood += model.l_g(likelihood_method, delta, g, g_0, θ)

  # return -(joint_errors_likelihood + observations_likelihood) / n
  # return -joint_errors_likelihood / n
  return -observations_likelihood / n

theta0 = [-0.03, 0.05, 1]
model = A01()
y_s = YieldSeries(table = prepare_data(), nfactors = 1)
# OK!!!
a01_true = minimize(likelihood, theta0, args=(model, y_s, "true",), method='nelder-mead', options= { 'xtol': 1e-6, 'disp': True, 'maxiter': 1000 })
print(a01_true.x)
a01_euler = minimize(likelihood, theta0, args=(model, y_s, "euler",), method='nelder-mead', options= { 'xtol': 1e-6, 'disp': True, 'maxiter': 1000 })
print(a01_euler.x)
a01_approx_1 = minimize(likelihood, theta0, args=(model, y_s, 1,), method='nelder-mead', options= { 'xtol': 1e-6, 'disp': True, 'maxiter': 1000 })
print(a01_approx_1.x)
a01_approx_2 = minimize(likelihood, theta0, args=(model, y_s, 2,), method='nelder-mead', options= { 'xtol': 1e-6, 'disp': True, 'maxiter': 1000 })
print(a01_approx_2.x)
