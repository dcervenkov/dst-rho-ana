TARGET := DSRhoLifetime

CXX_FILES := $(wildcard src/*.cc) src/gitversion.cc
OBJ_FILES := $(addprefix obj/,$(notdir $(CXX_FILES:.cc=.o)))
DEP_FILES := $(OBJ_FILES:.o=.d)

CXX := g++
#LD_FLAGS := ...
CXX_FLAGS := -std=c++11 -g -Wall -Wno-non-virtual-dtor
ROOT_FLAGS := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs) -lRooFitCore -lRooFit
LD_LIBS := -lboost_system -lboost_filesystem
INCS := -Isrc

BELLE_LIB_DIR ?= /belle
BELLE_INCS := -isystem $(BELLE_LIB_DIR)/belle/b20091203_2126/include
BELLE_LIBS := -L$(BELLE_LIB_DIR)/belle/b20091203_2126/x86_64-unknown-linux-gnu/opt/lib -ltatami -lbelleutil
CERN_LIBS := -L$(BELLE_LIB_DIR)/cern/2006/lib64 -lmathlib

.PHONY: all clean

all:: $(TARGET)

-include $(DEP_FILES)

$(TARGET): $(OBJ_FILES)
	$(CXX) -o $@ $^ $(LD_FLAGS) $(ROOT_LIBS) $(BELLE_LIBS) $(CERN_LIBS) $(LD_LIBS)

obj/%.o: src/%.cc | obj
	$(CXX) $(CXX_FLAGS) $(INCS) $(BELLE_INCS) $(ROOT_FLAGS) -c -MMD -o $@ $<

obj:
	mkdir obj

src/gitversion.cc: ../.git/HEAD ../.git/index $(filter-out src/gitversion.cc, $(CXX_FILES))
	echo "const char* gitversion = \"$(shell git rev-parse HEAD)$(shell git diff --quiet || echo '-dirty')\";" > $@

tidy: $(CXX_FILES)

clean::
	rm -f $(OBJ_FILES) $(DEP_FILES) $(TARGET) src/gitversion.cc

