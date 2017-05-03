from lib.hash import *

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
