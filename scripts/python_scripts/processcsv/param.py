import copy
from utility import index_of

class Param:
  """
  Class Param:
    Information about parameters used
  """
  def __init__(self):
    """ Constructor """
    ## Is it better to just read this from the first num_params cols of the header?
    self.param_names = [
      "INSTANCE",
      "SICS",
      "PHA",
      "HPLANE_SCORING_FN",
      "NUM_ALG2_ROUNDS",
      "NUM_RAYS_CUT",
      "MAX_FRAC_VAR",
      "CUT_LIMIT",
      "USE_SPLIT_SHARE",
      "NUM_CUTS_ITER_BILINEAR",
      "USE_UNIT_VECTORS_HEUR",
      "USE_CUT_VERT_HEUR",
      "USE_TIGHT_POINTS_HEUR",
      "CUT_PRESOLVE",
      "ROUNDS",
      "MIN_ORTHOGONALITY",
      "EPS",
      "RAYEPS",
      "TIMELIMIT"
    ]
    self.num_params = len(self.param_names) # Counting EPS as a parameter
    self.type_param = [int] * (self.num_params - 4)
    self.type_param.extend([float] * 4)
    self.type_param[0] = str
    self.param_default = [[] for i in range(self.num_params)]
    self._param_val = copy.deepcopy(self.param_default)

  def index(self, param_name):
    return self.param_names.index(param_name)

  def get_param(self, which_param):
    """ Get param value """
    if (type(which_param) is int):
      return copy.deepcopy(self._param_val[which_param])
    elif (type(which_param) is str):
      return copy.deepcopy(self._param_val[self.index(which_param)])
    else:
      raise TypeError(
        "Trying to get param, but is not a string or integer. Type is %s." 
        % type(which_param)
      )

  def set_param(self, param_name, param, set_val = None, type_of_param = None):
    """ Set param to either be itself, a list of the value given, or the current value """
    if (set_val is None):
      set_val = True
    param_index = self.index(param_name);
    
    if (type_of_param is None):
      type_of_param = self.type_param[param_index]
    if (param is None):
      param = copy.deepcopy(self._param_val[param_index])
    elif (type(param) is type_of_param):
      param = [param]
    elif (type(param) is not list):
      error_str = (
        """ 
        Trying to set parameter %s to non-integer value (or string if instance name, float if eps+rayeps value), and not a list.
        Type of param: %s
        """ % (param_name, str(type(param)))
      )
      raise TypeError(error_str)

    if (set_val):
      self._param_val[param_index] = copy.deepcopy(param)
    return param

  def set_default_param(self, param_name, param, set_val = None, type_of_param = None):
    """ Set default for param to either be itself or a list of the value given """
    if (set_val is None):
      set_val = False
    param_index = self.index(param_name);

    if (param is None):
      return self.param_default[param_index]

    if (type_of_param is None):
      type_of_param = self.type_param[param_index]
    if (type(param) is type_of_param):
      param = [param]
    elif (type(param) is not list):
      error_str = (
        """
        Trying to set parameter to non-integer value (or string if instance name, float if eps+rayeps value), and not a list.
        Type of param: %s
        """ % str(type(param))
      )
      raise TypeError(error_str)
    self.param_default[param_index] = copy.deepcopy(param)
    
    if (set_val):
      self._param_val[param_index] = copy.deepcopy(param)
    return param

  def get_inst(self, inst):
    """ Return name or index of inst from self._param_val[0] """
    if (type(inst) is str):
      return index_of(inst, self._param_val[0])
    elif (type(inst) is int):
      return self._param_val[0][inst]
    else:
      raise TypeError("Trying to get_inst for inst of type %s (not int or str)." % type(inst))
