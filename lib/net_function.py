from lib.enumerable import find_first

class NetFunction():
  def __init__(self, x, y):
    self.x = x
    self.y = y

  def call(self, value):
    if (value < self.x[0]) | (value > self.x[-1]):
      raise BaseException("Net function domain error: " + str(value))

    prev_index = find_first(lambda item: item <= value , self.x)

    return self.y[prev_index] + (value - self.x[prev_index]) * (self.y[prev_index + 1] - self.y[prev_index]) / (self.x[prev_index + 1] - self.x[prev_index])
