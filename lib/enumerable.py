def inject(f, memo, object):
  for x in object:
    memo = f(memo, x)
  return memo
