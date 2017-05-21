# theta0 = [-0.185793266805, 0.0264870240577, -0.866720304618]
# theta1 = [-0.0999041150336, -0.106014406172, 0.242671357612]
# -model01.γ_0(2, theta0) + model01.γ(2, theta0) * xt0
# yields = y_s.table[0]
# gt = yields[1]
# xt0 = (gt + model01.γ_0(1/12, theta0)) / model01.γ(1/12, theta0)
# xt1 = (gt + model01.γ_0(1/12, theta1)) / model01.γ(1/12, theta1)

# code.interact(local=dict(globals(), **locals()))


# print("-------------------------- A02 Model --------------------------------")
# theta0 = [-0.1, 0.8, 0.14, 0.15, -0.01, 0.16, 0.17]
# model02 = A02()
# # a11_qml = minimize(univariate_likelihood, theta0, args=(model11, y_s, "qml",), bounds = model11.bounds(), method='L-BFGS-B', options= { 'disp': True, 'maxiter': 1000 })
# print(a11_qml.x)
# a11_euler = minimize(univariate_likelihood, theta0, args=(model11, y_s, "euler",), bounds = model11.bounds(), method='L-BFGS-B', options= { 'disp': False, 'maxiter': 1000 })
# print(a11_euler.x)
# a11_approx_1 = minimize(multivariate_likelihood, theta0, args=(model02, y_s, 2,), method='nelder-mead', options= { 'xtol': 1e-6, 'disp': True, 'maxiter': 1000 })
# print(a11_approx_1.x)
# a11_approx_2 = minimize(univariate_likelihood, theta0, args=(model11, y_s, 2,), method='nelder-mead', options= { 'xtol': 1e-6, 'disp': True, 'maxiter': 1000 })
# print(a11_approx_2.x)


# Number of observations: 314

# -------------------------- A01 Model --------------------------------
# __________________________ TRUE likelihood __________________________
# | b_11 = -0.541283023755, δ_1 = 0.0270055331701, λ_1 = -1.39693665894 |
# Duration:  5.767s, iterations: 130
# ___________________________ Euler method ____________________________
# | b_11 = -0.529246987694, δ_1 = 0.0263945240161, λ_1 = -1.39752847953 |
# Duration:  5.135s, iterations: 119
# ____________________________ QML method _____________________________
# | b_11 = -0.541283023755, δ_1 = 0.0270055331701, λ_1 = -1.39693665894 |
# Duration:  5.934s, iterations: 130
# ______________________ Approximations (k = 1) _______________________
# | b_11 = -0.529966393226, δ_1 = 0.0269878628418, λ_1 = -1.37579241587 |
# Duration:  9.516s, iterations: 129
# ______________________ Approximations (k = 2) _______________________
# | b_11 = -0.541704337033, δ_1 = 0.0270060932315, λ_1 = -1.39779907979 |
# Duration: 18.224s, iterations: 137
# -------------------------- A11 Model --------------------------------
# ___________________________ Euler method ____________________________
# | b_11 = -3.61754739434e-08, δ_1 = 0.00738654079711, λ_1 = 0.310469387671, a_1 = 2.99428906425 |
# Duration: 31.301s, iterations: 191
# ____________________________ QML method _____________________________
# | b_11 = -1.00100032412e-07, δ_1 = 0.00752280440205, λ_1 = 0.236037909638, a_1 = 2.24951340998 |
# Duration: 34.456s, iterations: 193
# ______________________ Approximations (k = 1) _______________________
# | b_11 = -1.09569001229e-08, δ_1 = 0.0264426145877, λ_1 = 0.599075403724, a_1 = 1.17410989886 |
# Duration: 27.664s, iterations: 147
# ______________________ Approximations (k = 2) _______________________
# | b_11 = -2.04720362824e-08, δ_1 = 0.0264955200624, λ_1 = 0.58712463531, a_1 = 1.14777104854 |
# Duration: 46.645s, iterations: 190


yield_errors = []

for i in range(y_s_1_factor.length()):
  gt = y_s_1_factor[i, 1, model01.n]
  gt_ext = y_s_1_factor[i, model01.n + 1, 2 * model01.n]
  xt = (gt + model01.Γ_0(true_estimates)) / model01.Γ(true_estimates)
  gt_calculated = -model01.Γ_0_ext(true_estimates) + model01.Γ_ext(true_estimates) * xt
  yield_errors.append([abs(gt_calculated - gt_ext), i])

yield_errors.sort()
print(yield_errors)
