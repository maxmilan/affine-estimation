import numpy as np
from numpy.linalg import det
from numpy.linalg import inv
from scipy.optimize import minimize
import math

def likelihood(theta):
  tau = [1 / 12, 2]
  g = np.array(list(map(lambda x: [- ((math.exp(theta[0]) - math.exp(theta[2])) * x * theta[1]**2) / (theta[0] * theta[2] * (theta[0] - theta[2])), - (1 - x * math.exp(theta[2])) / theta[2]], tau)))
  result = det(inv(g.transpose()))
  return result

theta0 = np.array([-0.6, -1.2, 0.135])

estimate = minimize(likelihood, theta0, method='nelder-mead', options={'xtol': 1e-8, 'disp': True})  
print(estimate.x)
