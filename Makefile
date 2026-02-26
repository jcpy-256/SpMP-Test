LOADIMBA=0
LOADIMBA=${loadimba}

LIBRARY		= libspmp.a
CUSTOM_CXX = g++
# 用户自定义编译器
ifdef CUSTOM_CXX
  CXX   = $(CUSTOM_CXX)
  # Assume gcc
  CXXFLAGS		= -fopenmp -std=c++11 -Wno-deprecated -Wall
else
  CXX		= icpc
  CXXFLAGS		= -qopenmp -std=c++11 -Wno-deprecated -Wall
endif

# 调试模式
ifeq (yes, $(DBG))
  CXXFLAGS += -O0 -g
else
# 发布模式
  CXXFLAGS += -O3 -DNDEBUG
endif

ifeq (yes, $(XEON_PHI))
  CXXFLAGS	+= -mmic
else
  ifeq (yes, $(KNL))
    CXXFLAGS += -xMIC-AVX512
  else
    ifndef CUSTOM_CXX
      CXXFLAGS += -xHost
    endif
  endif
endif

ifeq (yes, $(MKL))
  ifndef CUSTOM_CXX
    CXXFLAGS += -mkl
  endif
  CXXFLAGS += -DMKL
endif

ifeq (${LOADIMBA}, 1)
  yesnolist += LOADIMBA
endif

INCLUDE = -I./sparseMtxOp/
# LDFLAGS = -L$(MKLROOT)/lib -lmkl_rt -lpthread -lm -ldl

DEFS += $(strip $(foreach var, $(yesnolist), $(if $(filter 1, $($(var))), -D$(var))))
DEFS += $(strip $(foreach var, $(deflist), $(if $($(var)), -D$(var)=$($(var)))))

CXXFLAGS 	+= ${DEFS}
CXXFLAGS  += ${INCLUDE}

CXXFLAGS += -g -w
SRCS = $(wildcard *.cpp) $(wildcard reordering/*.cpp)
SYNK_SRCS = $(wildcard synk/*.cpp)
OBJS = $(SRCS:.cpp=.o) $(addsuffix .o, $(addprefix synk/, $(basename $(notdir $(SYNK_SRCS)))))

# out = test/sptrsv_test
out = test/gs_test test/reordering_test test/trsv_test test/mtx2bin test/sptrsv_test test/ilu_test test/pcg test/ic_test


$(LIBRARY): $(OBJS)
	rm -f $(LIBRARY)
	ar qvs $(LIBRARY) $(OBJS) 

all: clean
	$(MAKE) $(LIBRARY)

test: $(out)

test/%: test/%.o $(LIBRARY)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -fPIC -c $< -o $@

clean:
	rm -f $(LIBRARY) $(OBJS) test/*.o $(out)
