import numpy

def multi_std(vectors):
  vectors = numpy.array(vectors)
  return list(map(lambda x: numpy.std(x), vectors.transpose()))
