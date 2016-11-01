TARGET := DSRhoEfficiency

CXX_FILES := $(wildcard src/*.cc)
OBJ_FILES := $(addprefix obj/,$(notdir $(CXX_FILES:.cc=.o)))
DEP_FILES := $(OBJ_FILES:.o=.d)

CXX := g++
#LD_FLAGS := ...
CXX_FLAGS := -std=c++11 -g -Wall -Wno-non-virtual-dtor
#BOOST_INC_DIR := /usr/local/include/boost-1_37
ROOT_FLAGS := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs) -lRooFitCore -lRooFit
INCS := -Isrc

.PHONY: all clean

all:: $(TARGET)

-include $(DEP_FILES)

$(TARGET): $(OBJ_FILES)
	g++ -o $@ $^ $(LD_FLAGS) $(ROOT_LIBS) 

obj/%.o: src/%.cc
	g++ $(CXX_FLAGS) $(INCS) $(ROOT_FLAGS) -c -MMD -o $@ $<

clean::
	rm -f $(OBJ_FILES) $(DEP_FILES) $(TARGET)

