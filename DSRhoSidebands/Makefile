TARGET := DSRhoSidebands

CXX_FILES := $(wildcard src/*.cc) src/gitversion.cc
OBJ_FILES := $(addprefix obj/,$(notdir $(CXX_FILES:.cc=.o)))
DEP_FILES := $(OBJ_FILES:.o=.d)

CXX := g++
CXX_FLAGS := -std=c++11 -g -Wall -Wno-non-virtual-dtor
LD_LIBS := -lboost_system -lboost_filesystem
ROOT_FLAGS := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs) -lRooFitCore -lRooFit
INCS := -Isrc

.PHONY: all clean tests

all:: $(TARGET)

-include $(DEP_FILES)

$(TARGET): $(OBJ_FILES)
	$(CXX) -o $@ $^ $(LD_LIBS) $(ROOT_LIBS)

obj/%.o: src/%.cc | obj
	$(CXX) $(CXX_FLAGS) $(INCS) $(ROOT_FLAGS) -c -MMD -o $@ $<

obj:
	mkdir obj

src/gitversion.cc: ../.git/HEAD ../.git/index $(filter-out src/gitversion.cc, $(CXX_FILES))
	echo "const char* gitversion = \"$(shell git rev-parse HEAD)$(shell git diff --quiet || echo '-dirty')\";" > $@

tests:: $(TARGET)
	cd tests && ./run_tests.py

tidy: $(CXX_FILES)
clean::
	rm -f $(OBJ_FILES) $(DEP_FILES) $(TARGET) src/gitversion.cc
