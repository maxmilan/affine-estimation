import xlrd
import numpy as np
import math
from numpy.linalg import det
from numpy.linalg import inv
from scipy.optimize import minimize
from scipy.stats import multivariate_normal
from datetime import datetime
from lib.yield_series import YieldSeries
from lib.matrix import *

FILENAME = "data/yields.xls"

book = xlrd.open_workbook(FILENAME)
sheet_index = book.sheet_names().index('ZEROYLD')
sheet = book.sheet_by_index(sheet_index)

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

def Γ(θ):
  τ = [1 / 12, 2]
  return np.array(list(map(lambda x: [- ((math.exp(θ[0]) - math.exp(θ[2])) * x * θ[1]**2) / (θ[0] * θ[2] * (θ[0] - θ[2])), - (1 - x * math.exp(θ[2])) / θ[2]], τ)))

def Γ0(θ):
  τ = [1 / 12, 2]
  return np.array(list(map(lambda x: x * (x**2 * (float(1) / 6 * math.exp(2 * θ[0]) + float(1) / 6 * math.exp(2 * θ[2]) - float(1) / 3 * math.exp(θ[0] + θ[2])) * θ[1]**4 + x * (0.5 * math.exp(θ[0]) - 0.5 * math.exp(θ[2])) * θ[0] * θ[1]**2 * θ[2]**2 * θ[3] + θ[0]**3 * θ[2] * (-1 + x * math.exp(θ[2]) - float(1) / 3 * x**2 * math.exp(2 * θ[2]) + (2 - x * math.exp(θ[2])) * θ[2] * θ[4]) + θ[0]**4 * (0.5 - 0.5 * x * math.exp(θ[2]) + float(1) / 6 * x**2 * math.exp(2 * θ[2]) + (-1 + 0.5 * x * math.exp(θ[2])) * θ[2] * θ[4]) + θ[0]**2 * θ[2] * ((0.5 - 0.5 * x * math.exp(θ[2]) + float(1) / 6 * x**2 * math.exp(2 * θ[2])) * θ[2] + x * (-0.5 * math.exp(θ[0]) + 0.5 * math.exp(θ[2])) * θ[1]**2 * θ[3] + (-1 + 0.5 * x * math.exp(θ[2])) * θ[2]**2 * θ[4])) / (θ[0]**2 * θ[2]**2 * (θ[0]**2 - 2 * θ[0] * θ[2] + θ[2]**2)), τ)))

y_s = YieldSeries(table = prepare_data(), nfactors = 2)

# error_distribution = multivariate_normal(mean = np.zeros(2), cov = np.diag([1, 1]))
# print(error_distribution.pdf([1,0]))
 
θ = [0.1, 0.2, 0.3, 0.03, 0.04]
gt = y_s[0, 1, 2]
xt = inv(Γ(θ).transpose()).dot(subtract(gt, Γ0(θ)))
print(xt)
