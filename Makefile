BUILD_MODE = RELEASE
OPENMP = NO
JOBS=1
CXXFLAGS = -Wall -static -fmessage-length=0 -std=gnu++0x
LNKFLAGS = -Wl,--wrap,memcpy

# Change CXX to MPIC++
ifneq ($(OS),Windows_NT)
	CXX:= mpic++
endif

# Needs gperftools to be installed and configured before use
PROFILE = NO

PROJECT_ROOT = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SRC = $(PROJECT_ROOT)/src
MSTOOLKIT = $(PROJECT_ROOT)/mstoolkit

DEPS = src/dslim_query.cpp src/dslim.cpp src/lbe.cpp src/lbe_internal.cpp src/utils.cpp src/mods.cpp src/msquery.cpp include/common.h include/config.h include/keyval.h include/utils.h include/mods.h include/msquery.h include/slm_dsts.h include/dslim.h include/slmerr.h include/lbe.h

LIBS = -ldslim -lmstoolkitlite -lstdc++ -lc
LIBPATHS = -L$(SRC) -L$(MSTOOLKIT)
EXECUTEABLE = LoadBalancer.exe

ifeq ($(PROFILE),YES)
	LIBS += -lprofiler
	CXXFLAGS += -D_PROFILE -g
endif

ifeq ($(OPENMP),YES)
	CXXFLAGS += -fopenmp -D_GLIBCXX_PARALLEL
	LNKFLAGS += -fopenmp
endif

# Setup build mode
ifeq ($(BUILD_MODE),DEBUG)
	CXXFLAGS += -g3 -O0 
else ifeq ($(BUILD_MODE),RELEASE)
	CXXFLAGS += -O3
else
	$(error Build mode $(BUILD_MODE) not supported by this Makefile)
endif

# Setup parallel build jobs
ifeq ($(OS),Windows_NT)
	JOBS := $(NUMBER_OF_PROCESSORS)
else
	JOBS := $(shell grep -c ^processor /proc/cpuinfo)
endif

# Export Build Variables
export BUILD_MODE
export OPENMP
export JOBS
export CXXFLAGS

all: mstoolkit lbe
	$(CXX) $(LNKFLAGS) $(LIBPATHS) -Wl,--start-group $(LIBS) -Wl,--end-group -o $(EXECUTEABLE) -Wl,-Map=$(EXECUTEABLE).map

mstoolkit:
	$(MAKE) -j$(JOBS) -C $(MSTOOLKIT)

lbe: $(DEPS)
	$(MAKE) -j$(JOBS) -C $(SRC)

clean:
	$(MAKE) -C $(SRC) clean
	rm -rf $(EXECUTEABLE) $(EXECUTEABLE).map

allclean:
	$(MAKE) -C $(SRC) clean
	rm -rf $(EXECUTEABLE) $(EXECUTEABLE).map
	$(MAKE) -C $(MSTOOLKIT) clean

.PHONY: all mstoolkit lbe clean allclean
