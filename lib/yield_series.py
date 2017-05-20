import xlrd
from datetime import *
from lib.hash import *
from lib.enumerable import *

def import_data(data_storage, nfactors = 1):
  return YieldSeries(table = list(map(lambda x: x[1], data_storage.table)), nfactors = nfactors)


class DataStorage():
  def __init__(self, filename = "data/yields.xls", sheet_name = "ZEROYLD", from_date = datetime(1965, 1, 1)):
    self.book = xlrd.open_workbook(filename)
    self.sheet = self.book.sheet_by_index(self.book.sheet_names().index(sheet_name))
    self.from_date = from_date
    self._table = None

  @property
  def table(self):
    if self._table is None:
      first_row = self.sheet.row_values(0)
      columns_indices = range(1, len(first_row))
      data = []

      for row_index in range(self.sheet.nrows):
        row_values = self.sheet.row_values(row_index)
        if row_values[0] and datetime.strptime(row_values[0], "%m/%Y") >= self.from_date:
          item = {}

          for i in range(len(columns_indices)):
            if not row_values[columns_indices[i]]:
              value = None
            else:
              value = float(row_values[columns_indices[i]]) / 100
            item[first_row[columns_indices[i]]] = value

          data.append([row_values[0], item])

      self._table = data


    return self._table

  def yields_by_date(self, date):
    index = find_first(lambda x: x[0] == date.strftime("%m/%Y"), self.table)
    plot_list = hash_to_a(self.table[index][1])
    plot_list.sort()

    return plot_list

  def yields_by_maturity(self, maturity):
    data = {}
    for x in self.table:
      data[datetime.strptime(x[0], "%m/%Y")] = x[1][maturity]

    plot_list = hash_to_a(data)
    plot_list.sort()

    return plot_list

class YieldSeries():
  def __init__(self, table=[], nfactors=1):
    self.table = table
    self.nfactors = nfactors
    self._series = None
    self._maturities = None

  @property
  def series(self):
    if self._series is None:
      self._series = list(map(lambda x: { k:v for (k, v) in x.items() if k in self.maturities }, self.table))

    return self._series

  @property
  def maturities(self):
    if self._maturities is None:
      memo = list(map(lambda x: x * 24, range(1, 2 * self.nfactors)))
      memo.insert(0, 1.0)
      self._maturities = memo

    return self._maturities

  def length(self):
    return len(self.series)

  def __getitem__(self, keys):
    if isinstance(keys, int):
      time = keys
      maturity = None
    else:
      time = keys[0]
      if len(keys) < 3:
        last_item = None
      else:
        last_item = keys[2]

      maturity = slice(keys[1] - 1, last_item)

    section = self.series[time]

    if maturity is None:
      return sorted_dict_values(section)
    else:
      return sorted_dict_values(section)[maturity]
