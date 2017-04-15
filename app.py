import xlrd
import code
from scipy.optimize import minimize
from datetime import datetime

FILENAME = "data/yields.xls"

book = xlrd.open_workbook(FILENAME)
sheet_index = book.sheet_names().index('ZEROYLD')
sheet = book.sheet_by_index(sheet_index)

def prepare_data(nfactors = 1):
  maturities_in_months = list(map(lambda x: x * 24, range(0, 2 * nfactors)))
  first_row = sheet.row_values(0)
  columns_indices = list(map(lambda x: first_row.index(x), maturities_in_months))
  data = []

  for row_index in range(sheet.nrows):
    row_values = sheet.row_values(row_index)
    if row_values[0] and datetime.strptime(row_values[0], "%m/%Y") >= datetime(1972, 1, 1):
      item = {}
   
      for i in range(len(columns_indices)):
        item[maturities_in_months[i]] = row_values[columns_indices[i]]
      data.append(item)

  return data

result = prepare_data(nfactors = 2)

# code.interact(local=dict(globals(), **locals()))
