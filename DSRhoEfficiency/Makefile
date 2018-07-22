TARGET := DSRhoEfficiency

CXX_FILES := $(wildcard src/*.cc)
OBJ_FILES := $(addprefix obj/,$(notdir $(CXX_FILES:.cc=.o)))
DEP_FILES := $(OBJ_FILES:.o=.d)

CXX := g++
LD_FLAGS := -Wl,-rpath $(HOME)/bin/belle/packages/Meerkat/lib
CXX_FLAGS := -std=c++11 -g -Wall -Wno-non-virtual-dtor
#BOOST_INC_DIR := /usr/local/include/boost-1_37
ROOT_FLAGS := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs) -lRooFitCore -lRooFit
LIBS := -L$(HOME)/bin/belle/packages/Meerkat/lib -lMeerkat
INCS := -Isrc -I$(HOME)/bin/belle/packages/Meerkat/inc

.PHONY: all clean

all:: $(TARGET)

-include $(DEP_FILES)

$(TARGET): $(OBJ_FILES)
	$(CXX) -o $@ $^ $(LD_FLAGS) $(ROOT_LIBS) $(LIBS)

obj/%.o: src/%.cc | obj
	$(CXX) $(CXX_FLAGS) $(INCS) $(ROOT_FLAGS) -c -MMD -o $@ $<

obj:
	mkdir obj

clean::
	rm -f $(OBJ_FILES) $(DEP_FILES) $(TARGET)

