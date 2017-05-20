def sorted_dict_values(dict):
  sorted_keys = sorted(dict.keys())
  return list(map(lambda x: dict[x], sorted_keys))

def hash_to_a(dict):
  dictlist = []
  for key, value in dict.items():
    dictlist.append([key,value])

  return dictlist
