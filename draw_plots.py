import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from lib.yield_series import *
from datetime import *

PATH = "/home/max/Projects/diploma/pics/"

data_storage = DataStorage("data/yields.xls", "ZEROYLD", datetime(1965, 1, 1))

plot_list = np.array(data_storage.yields_by_maturity(0)).transpose()

plt.plot_date(x=plot_list[0], y=plot_list[1], fmt="b-")
plt.style.use('ggplot')
plt.grid(True)
# plt.show()
plt.savefig(PATH + 'ir_series.png', transparent=False, bbox_inches='tight', pad_inches=0.1)
