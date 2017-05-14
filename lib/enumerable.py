def inject(f, memo, object):
  for x in object:
    memo = f(memo, x)
  return memo

def find_first(f, object):
  for i in range(len(object)):
    if f(object[i]):
      return i

  return -1
