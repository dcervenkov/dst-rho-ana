TARGET := DSRhoLifetime

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

BELLE_LIB_DIR ?= /belle
BELLE_INCS := -I$(BELLE_LIB_DIR)/belle/b20091203_2126/include
BELLE_LIBS := -L$(BELLE_LIB_DIR)/belle/b20091203_2126/x86_64-unknown-linux-gnu/opt/lib -ltatami -lbelleutil
CERN_LIBS := -L$(BELLE_LIB_DIR)/cern/2006/lib64 -lmathlib

.PHONY: all clean

all:: $(TARGET)

-include $(DEP_FILES)

$(TARGET): $(OBJ_FILES)
	g++ -o $@ $^ $(LD_FLAGS) $(ROOT_LIBS) $(BELLE_LIBS) $(CERN_LIBS)

obj/%.o: src/%.cc | obj
	g++ $(CXX_FLAGS) $(INCS) $(BELLE_INCS) $(ROOT_FLAGS) -c -MMD -o $@ $<

obj:
	mkdir obj

clean::
	rm -f $(OBJ_FILES) $(DEP_FILES) $(TARGET)

