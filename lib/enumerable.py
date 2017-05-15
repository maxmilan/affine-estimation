def inject(f, memo, object):
  for x in object:
    memo = f(memo, x)
  return memo

def find_first(f, object):
  for i in range(len(object)):
    if f(object[i]):
      return i

  return -1

def equal_enumerables(a, b):
  if (a is None) or (b is None) or (len(a) != len(b)):
    return False

  for i in range(len(a)):
    if a[i] != b[i]:
      return False

  return True
