### Shell type ###
UNAME := $(shell uname)

ifeq ($(UNAME),Linux)
  CC = g++
endif
ifeq ($(UNAME),Darwin)
  CC = clang++
endif

### Build type ###
# Choose 'debug' or 'release'
# Can also be chosen through make "BUILD_CONFIG=XX" from command line 
# Or one can call make debug or make release directly
#BUILD_CONFIG = release
BUILD_CONFIG = debug

### Variables user should set ###
EXECUTABLE_STUB = PHA
SRC_DIR = src
USE_GUROBI=0

SOURCES = \
		$(EXECUTABLE_STUB).cpp \
		classes/AdvCut.cpp \
		classes/Hplane.cpp \
		classes/Point.cpp \
		classes/Ray.cpp \
		classes/SolutionInfo.cpp \
		classes/Vertex.cpp \
		cuts/CutHelper.cpp \
		cuts/CglPHA/CglPHA.cpp \
		cuts/CglPHA/hplaneActivation.cpp \
		cuts/CglSIC/CglSIC.cpp \
		utility/CglGICParam.cpp \
		utility/GlobalConstants.cpp \
		utility/optionhelper.cpp \
		utility/Output.cpp \
		utility/Utility.cpp
ifeq ($(BUILD_CONFIG),debug)
  SOURCES += utility/Debug.cpp
endif
ifeq ($(USE_GUROBI),1)
  SOURCES += cuts/BBHelper.cpp
endif

DIR_LIST = \
		$(SRC_DIR) \
		$(SRC_DIR)/classes \
		$(SRC_DIR)/cuts \
		$(SRC_DIR)/cuts/CglPHA \
		$(SRC_DIR)/cuts/CglSIC \
		$(SRC_DIR)/utility

ifeq ($(USE_GUROBI),1)
  GUROBI_DEFS = -DSHOULD_USE_GUROBI
endif

### Set build values based on user variables ###
ifeq ($(BUILD_CONFIG),debug)
  # "Debug" build - no optimization, include debugging symbols, and keep inline functions
  OUT_DIR = Debug
  DEBUG_FLAG = -g3
  OPT_FLAG = -O0
  DEFS = -DTRACE -DSHOULD_CALC_ACTIVE_FINAL_POINTS $(GUROBI_DEFS)
  EXTRA_FLAGS = -fmessage-length=0
  ifeq ($(CC),g++)
    EXTRA_FLAGS += -fkeep-inline-functions
  endif
endif

ifeq ($(BUILD_CONFIG),release)
  # "Release" build - maximum optimization, no debug symbols
  OUT_DIR = Release
  DEBUG_FLAG = 
  OPT_FLAG = -O0
  DEFS = $(GUROBI_DEFS)
	DEFS += -DSHOULD_CALC_ACTIVE_FINAL_POINTS
  EXTRA_FLAGS = -fmessage-length=0 -ffast-math
  ifeq ($(USE_GUROBI),1)
    OUT_DIR = ReleaseGurobi
    #DEFS += -DSHOULD_CALC_ACTIVE_FINAL_POINTS
  endif
endif

EXECUTABLE = $(OUT_DIR)/$(EXECUTABLE_STUB)

# It is important that the only thing that changes about these directories
# is OUT_DIR depending on whether it is the debug or release build
# This is because later (building dependencies, archive file, cleaning, etc.)
# depends on this fact (when doing *_debug targets)
OBJ_DIR = $(OUT_DIR)/$(SRC_DIR)

OUT_DIR_LIST = $(addprefix $(OUT_DIR)/,$(DIR_LIST))

OBJECTS = $(SOURCES:.cpp=.o)

OUT_OBJECTS = $(addprefix $(OBJ_DIR)/,$(OBJECTS))

# Make sure that environment variables are properly defined
# If not defined for the environment, define COIN-OR variables here
COIN_OR_DIR = ${PHA_DIR}/lib
COIN_VERSION = Cgl-0.59
COIN_DEBUG = $(COIN_OR_DIR)/$(COIN_VERSION)/buildg
COIN_RELEASE = $(COIN_OR_DIR)/$(COIN_VERSION)/build
ifeq ($(BUILD_CONFIG),debug)
  COIN = ${COIN_DEBUG}
endif
ifeq ($(BUILD_CONFIG),release)
  COIN = ${COIN_RELEASE}
endif
COINlib = ${COIN}/lib
COINinc = ${COIN}/include/coin

INCL_SRC_DIRS = $(addprefix -I,$(DIR_LIST))

APPLINCLS = $(INCL_SRC_DIRS)
APPLINCLS += -isystem $(COINinc)
APPLLIB = -L$(COINlib)

ifeq ($(USE_GUROBI),1)
  APPLINCLS += -I${GUROBI_INC}
endif

APPLLIB += \
		-lOsiClp \
    -lCgl \
		-lOsi \
		-lClp \
		-lCoinUtils

ifneq (${ENV_LAPACK_LIB},)
  APPLLIB += -L${ENV_LAPACK_LIB} -llapack
endif

ifneq (${ENV_BLAS_LIB},)
  APPLLIB += -L${ENV_BLAS_LIB} -lblas
endif

APPLLIB += -lm -lz -lbz2 -lreadline

ifeq ($(USE_GUROBI),1)
  APPLLIB += -L${GUROBI_LIB} -lgurobi_c++ -l${GUROBI_LINK} -lm
endif

CXXLINKFLAGS = -Wl,-rpath $(COINlib)

CXXFLAGS = -Wall -MMD -MP 
#CXXFLAGS = -Wall -Wextra -Wpedantic -MMD -MP 
ifeq ($(CC),clang++)
  CXXFLAGS += -Wno-gnu-zero-variadic-macro-arguments
endif
CXXFLAGS += -std=gnu++11
ifeq ($(CC),clang++)
  CXXFLAGS += -stdlib=libc++
  CXXLINKFLAGS += -stdlib=libc++
  ifeq ($(BUILD_CONFIG),debug)
    CXXLINKFLAGS += -fsanitize=address
  endif
  #CXXFLAGS += -stdlib=libstdc++
  #CXXLINKFLAGS += -stdlib=libstdc++
  APPLLIB += -framework Accelerate
endif
CXXFLAGS += -m64 $(DEBUG_FLAG) $(OPT_FLAG) $(EXTRA_FLAGS)

### Targets ###
all: | directories	$(EXECUTABLE)
debug: FORCE
	@$(MAKE) "BUILD_CONFIG=debug"
release: FORCE
	@$(MAKE) "BUILD_CONFIG=release"

$(EXECUTABLE): $(OUT_OBJECTS)
		@echo ' '
		@echo 'Building target: $@'
		@echo 'Invoking' $(CC) 'linker'
		$(CC) $(DEFS) $(CXXLINKFLAGS) $(APPLINCLS) -o $@ $^ $(APPLLIB)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
		@echo ' '
		@echo 'Building target: $@'
		@echo 'Invoking' $(CC) 'compiler'
		$(CC) $(CXXFLAGS) $(DEFS) $(APPLINCLS) -c $< -o $@ 
		@echo 'Finished building: $@'

### Dependencies ###
# Dependencies (the -include says to ignore errors)
DEPENDENCIES = $(OUT_OBJECTS:.o=.d)
-include $(DEPENDENCIES)

### Phony ###
#.PHONY = \
#		all \
#		clean clean_debug distclean_debug clean_release distclean_release \
#		directories dir_debug dir_release \
#		print print_dep

.PHONY = all clean directories print

### Cleaning ###
RM = rm -f
clean_%: FORCE
	@$(MAKE) clean "BUILD_CONFIG=$*"
distclean_%: FORCE
	@$(MAKE) distclean "BUILD_CONFIG=$*"

clean: FORCE
	@$(RM) $(OUT_OBJECTS) $(EXECUTABLE)

distclean: FORCE
	@$(RM) $(OUT_OBJECTS) $(EXECUTABLE) $(DEPENDENCIES)

### Making directories that you need ###
MKDIR_P = mkdir -p

dir_%: FORCE
	@$(MAKE) directories "BUILD_CONFIG=$*"
directories: $(OUT_DIR_LIST)
$(OUT_DIR_LIST):
	$(MKDIR_P) $(OUT_DIR_LIST)

print: FORCE
	@echo 'EXECUTABLE: $(EXECUTABLE)'
	@echo 'OUT_DIR: $(OUT_DIR)'
	@echo 'DEBUG_FLAG: $(DEBUG_FLAG)'
	@echo 'OPT_FLAG: $(OPT_FLAG)'
	@echo 'DEFS: $(DEFS)'
	@echo 'EXTRA_FLAGS: $(EXTRA_FLAGS)'
	@echo 'SOURCES: $(SOURCES)'
	@echo 'OUT_OBJECTS: $(OUT_OBJECTS)'

print_dep: FORCE
	@echo 'DEPENDENCIES: $(DEPENDENCIES)'

FORCE: 
