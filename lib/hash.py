def sorted_dict_values(dict):
  sorted_keys = sorted(dict.keys())
  return list(map(lambda x: dict[x], sorted_keys))
