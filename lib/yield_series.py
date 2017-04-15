import time

class YieldSeries():
  def __init__(self, table=[]):
    self.table = table
    self._bar = None

  @property
  def bar(self):
    if self._bar is None:
      print("starting long calculation")
      time.sleep(1)
      self._bar = self.table
      print("finished long caclulation")
    
    return self._bar
