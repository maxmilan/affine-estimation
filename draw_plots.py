import xlrd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import *

FILENAME = "data/yields.xls"

book = xlrd.open_workbook(FILENAME)
sheet_index = book.sheet_names().index('ZEROYLD')
sheet = book.sheet_by_index(sheet_index)
delta = 1/12

first_row = sheet.row_values(0)
columns_indices = range(1, len(first_row))

data = []
dates = []
values = []

for row_index in range(sheet.nrows):
  row_values = sheet.row_values(row_index)
  if row_values[0] and datetime.strptime(row_values[0], "%m/%Y") >= datetime(1946, 12, 1):
    item = {}

    for i in range(len(columns_indices)):
      if not row_values[columns_indices[i]]:
        value = None
      else:
        value = float(row_values[columns_indices[i]]) / 100
      item[first_row[columns_indices[i]]] = value

    data.append(item)
    dates.append(datetime.strptime(row_values[0], "%m/%Y").date())
    values.append(item[0])

plt.plot_date(x=dates, y=values, fmt="b-")
plt.style.use('ggplot')
plt.grid(True)
# plt.show()
plt.savefig('ir_series.png', transparent=False, bbox_inches='tight', pad_inches=0.1)
