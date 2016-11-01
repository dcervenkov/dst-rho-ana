#
# Makefile for libRooTatami
#

TARGET_LRT = libRooTatami.so
TARGET_LRT_PORTABLE = libRooTatamiPortable.so
OBJS_LRT = RooDtLifetime.o RooDtBkg.o RooDtCPSignal.o RooTatamiHelper.o RooTatami.o 
HDRS_LRT = RooDtLifetime.h RooDtBkg.h RooDtCPSignal.h 

# Compiler
CXX = g++

# Optimization and Debug flags
OPT = -O2 #-g -v

# Other compiler flags
CFLAGS   = $(OPT) -DHEP_SHORT_NAMES -DBELLE_SHORT_NAMES -DDSTXX_NOINLINE -fPIC -Wall -Wextra
CXXFLAGS = $(OPT) -DHEP_SHORT_NAMES -DBELLE_SHORT_NAMES -DDSTXX_NOINLINE -fPIC -Wall -Wextra
SOFLAGS  = $(OPT) -shared -Wl,-export-dynamic

# ROOT flags and libs
ROOT_FLAGS_INC = $(shell root-config --cflags)
ROOT_LIB = $(shell root-config --libs)
ROOFIT_LIB = -lRooFitCore -lRooFit -lHtml -lMinuit -lFoam -lGui
# Link against libMathMore if it exists
ROOFIT_LIB += $(shell if [ -e $(ROOTSYS)/lib/libMathMore.so ]; then echo -lMathMore; fi)

# Belle and cernlib specific include and libs (for libRooTatami only)
ifdef BELLE_LEVEL
INCS_TATAMI = -I$(BELLE_TOP_DIR)/include
LIBS_TATAMI = -L$(BELLE_LIB_DIR)/so -L$(BELLE_LIB_DIR) -ltatami -lbelleutil -lbelleCLHEP -L$(CERN_ROOT)/lib64/so -L$(CERN_ROOT)/lib64 -lpacklib -lmathlib
endif

# BOOST directories
BOOST_INC =
BOOST_LIB =

# All Includes
INCS = -I. $(ROOT_FLAGS_INC) $(INCS_TATAMI) $(BOOST_INC)

# All libs 
LIBS = -L. $(ROOT_LIB) $(ROOFIT_LIB) $(BOOST_LIB)

.SUFFIXES: .c .cc

# C++ source files
.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCS) -o $@ $<

# ------------------------------------------------------------------------

all::
ifdef BELLE_LEVEL
	@echo "BELLE_LEVEL "$(BELLE_LEVEL)" detected"
	$(MAKE) $(TARGET_LRT)
	$(MAKE) $(TARGET_LRT_PORTABLE)
else
	@echo "-----------------------------------------------------------------------"
	@echo "Belle software seems not to be installed on this machine. At least the"
	@echo "BELLE_LEVEL environment variable does not exist. Probably you need to"
	@echo "source the bashrc_general or cshrc_general."
	@echo "This might also be on purpose if you don't want to use BELLE specific"
	@echo "resolution functions, e.g. in case of BELLE 2."
	@echo "Skipping built of $(TARGET_LRT) and $(TARGET_LRT_PORTABLE)"
	@echo "-----------------------------------------------------------------------"
endif

ifndef BELLE_LEVEL
$(TARGET_LRT):
	@echo "-----------------------------------------------------------------------"
	@echo "Belle software seems not to be installed on this machine. At least the"
	@echo "BELLE_LEVEL environment variable does not exist. Probably you need to"
	@echo "source the bashrc_general or cshrc_general."
	@echo "This might also be on purpose if you don't want to use BELLE specific"
	@echo "resolution functions, e.g. in case of BELLE 2."	
	@echo "-----------------------------------------------------------------------"	
else
$(TARGET_LRT): $(OBJS_LRT)
	$(CXX) -o $@ $(OBJS_LRT) $(LIBS) $(LIBS_TATAMI) $(SOFLAGS) 
endif

ifndef BELLE_LEVEL
$(TARGET_LRT_PORTABLE):
	@echo "-----------------------------------------------------------------------"
	@echo "Belle software seems not to be installed on this machine. At least the"
	@echo "BELLE_LEVEL environment variable does not exist. Probably you need to"
	@echo "source the bashrc_general or cshrc_general."
	@echo "This might also be on purpose if you don't want to use BELLE specific"
	@echo "resolution functions, e.g. in case of BELLE 2."	
	@echo "-----------------------------------------------------------------------"	
else
$(TARGET_LRT_PORTABLE): $(OBJS_LRT)
	$(CXX) -o $@ $(OBJS_LRT) $(SOFLAGS) -Wl,-Bstatic $(LIBS_TATAMI) -Wl,-Bdynamic $(LIBS) #-Wl,--verbose
endif

RooTatami.o: $(HDRS_LRT)
	rm -f RooTatami.cc RooTatami.h RooTatami.o
	rootcint RooTatami.cc -c $(HDRS_LRT)
	$(CXX) -c $(CXXFLAGS) $(INCS) -o RooTatami.o RooTatami.cc

clean::
	rm -f $(OBJS_LRT)

distclean::
	$(MAKE) clean
	rm -f $(TARGET_LRT)
ifndef BELLE_LEVEL
	@echo "-----------------------------------------------------------------------"
	@echo $(TARGET_LRT_PORTABLE)" was not deleted as it probably can't be rebuilt on this system"
	@echo "-----------------------------------------------------------------------"
else
	rm -f $(TARGET_LRT_PORTABLE)
endif
	rm -f RooTatami.cc RooTatami.h RooTatami.o
