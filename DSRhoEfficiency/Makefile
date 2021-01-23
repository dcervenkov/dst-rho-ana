TARGET := DSRhoEfficiency

CXX_FILES := $(wildcard src/*.cc) src/gitversion.cc
OBJ_FILES := $(addprefix obj/,$(notdir $(CXX_FILES:.cc=.o)))
DEP_FILES := $(OBJ_FILES:.o=.d)

CXX := g++
LD_FLAGS := -Wl,-rpath $(HOME)/bin/belle/packages/Meerkat/lib
CXX_FLAGS := -std=c++11 -g -Wall -Wno-non-virtual-dtor
ROOT_FLAGS := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs) -lRooFitCore -lRooFit
LD_LIBS := -lboost_system -lboost_filesystem -L$(HOME)/bin/belle/packages/Meerkat/lib -lMeerkat -L$(HOME)/programs/boost_1_74_0/stage/lib
INCS := -Isrc -I$(HOME)/bin/belle/packages/Meerkat/inc -I$(HOME)/programs/boost_1_74_0

.PHONY: all clean

all:: $(TARGET)

-include $(DEP_FILES)

$(TARGET): $(OBJ_FILES)
	$(CXX) -o $@ $^ $(LD_FLAGS) $(ROOT_LIBS) $(LD_LIBS)

obj/%.o: src/%.cc | obj
	$(CXX) $(CXX_FLAGS) $(INCS) $(ROOT_FLAGS) -c -MMD -o $@ $<

obj:
	mkdir obj

src/gitversion.cc: ../.git/HEAD ../.git/index $(filter-out src/gitversion.cc, $(CXX_FILES))
	echo "const char* gitversion = \"$(shell git rev-parse HEAD)$(shell git diff --quiet || echo '-dirty')\";" > $@

tidy: $(CXX_FILES)
clean::
	rm -f $(OBJ_FILES) $(DEP_FILES) $(TARGET) src/gitversion.cc

