class StatEnum:
  stat_name = ["average", "min", "max", "first", "maxratio"]
  AVG = stat_name.index("average")
  AVERAGE = stat_name.index("average")
  MIN = stat_name.index("min")
  MAX = stat_name.index("max")
  FIRST = stat_name.index("first")
  RATIO = stat_name.index("maxratio")
  MAXRATIO = stat_name.index("maxratio")

  @classmethod
  def get_name_list(cls, stat_index):
    if (type(stat_index) is int):
      stat_index = [stat_index]

    curr_stat_name = []
    for i in range(len(stat_index)):
      curr_stat_name.append(cls.stat_name[stat_index[i]])

    return curr_stat_name

  @classmethod
  def get_name(cls, stat_index):
    if (type(stat_index) is list):
      return cls.get_name_list(stat_index)
    else:
      return cls.stat_name[stat_index]

