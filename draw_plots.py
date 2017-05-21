import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from lib.data_storage import *
from models.a_0_1 import A01
from datetime import *

PATH = "/home/max/Projects/diploma/pics/"

data_storage = DataStorage("data/yields.xls", "ZEROYLD", datetime(1965, 1, 1))

plot_list = np.array(data_storage.yields_by_maturity(0)).transpose()

# plt.plot_date(x=plot_list[0], y=plot_list[1], fmt="b-")
# plt.style.use('ggplot')
# plt.grid(True)
# # plt.show()
# plt.savefig(PATH + 'ir_series.png', transparent=False, bbox_inches='tight', pad_inches=0.1)

theta01_opt = [-0.541283023755, 0.0270055331701, -1.39693665894]
plot_list = np.array(data_storage.yields_by_date(datetime(1966, 10, 1))).transpose()
model01 = A01()

maturities = list(map(lambda x: x / 12, plot_list[0]))[1:30]
yields = plot_list[1][1:30]
calculated_yields = list(map(lambda x: model01.g(yields[0], 1 / 12, theta01_opt, x), maturities))

plt.plot(maturities, yields, 'b-', maturities, calculated_yields, 'r--')
plt.style.use('ggplot')
plt.grid(True)
plt.savefig(PATH + '1966_10_yields.png', transparent=False, bbox_inches='tight', pad_inches=0.1)
