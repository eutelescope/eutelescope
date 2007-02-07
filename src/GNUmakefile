###########################################################
#
# Example GNUmakefile to build an Application with Marlin
#
# @author Frank Gaede
# @author Jan Engels
# @version $Id: GNUmakefile,v 1.1.1.1 2007-02-07 10:53:12 bulgheroni Exp $
#
###########################################################


#----------------------------------------------------------------------------
#        user variables
#----------------------------------------------------------------------------

# name of the program
PROGNAME = Eutelescope

# source files for the binary
# leave this commented out if your Application does not produce binary files
#BINSRCS = $(PROGNAME).cc test.cc

# marlin working directory
ifndef MARLINWORKDIR
 MARLINWORKDIR=$(MARLIN)
 export MARLINWORKDIR
endif

# additional user include/libraries
ifndef USERINCLUDES
 USERINCLUDES +=
endif

ifndef USERLIBS
 USERLIBS +=
endif

#----------------------------------------------------------------------------
#        do not modify anything below in order to use top level makefile
#----------------------------------------------------------------------------

# directories
INCDIR = ../include
LIBDIR = $(MARLINWORKDIR)/lib
BINDIR = $(MARLINWORKDIR)/bin
DOCDIR = $(MARLINWORKDIR)/doc/$(PROGNAME)
TMPDIR = $(MARLINWORKDIR)/tmp/$(PROGNAME)
OBJDIR = $(TMPDIR)/obj
DEPDIR = $(TMPDIR)/dep

# library file without path
LIBNAME = lib$(PROGNAME).a

# library file with path
LIB = $(LIBDIR)/$(LIBNAME)

ifdef BINSRCS
# binary files with path
BIN = $(patsubst %.cc,$(BINDIR)/%,$(BINSRCS))

# object files for the binaries without path
BINOBJS = $(patsubst %.cc,%.o,$(BINSRCS))

# object files for the binaries with path
BINOBJ = $(patsubst %,$(OBJDIR)/%,$(BINOBJS))
endif

# source files without path
SRCS = $(filter-out $(BINSRCS),$(wildcard *.cc))

# object files without path
OBJS = $(patsubst %.cc,%.o,$(SRCS))

# object files with path
OBJ = $(patsubst %,$(OBJDIR)/%,$(OBJS))

# dependency files without path
DEPS = $(patsubst %.o,%.d,$(OBJS) $(BINOBJS))

# dependency files with path
DEP = $(patsubst %,$(DEPDIR)/%,$(DEPS))

# c++ compiler
ifdef MY_CXX
  CXX = $(MY_CXX)
else
  CXX = g++
endif

# compiler flags
ifdef MARLINDEBUG
 CXXFLAGS = -g -Wall -ansi
else
 CXXFLAGS = -O3 -Wall -ansi
endif

# pre-processor flags
CPPFLAGS = -I $(INCDIR)
CPPFLAGS += $(USERINCLUDES)
CPPFLAGS += `. $(MARLINWORKDIR)/bin/marlin_includes.sh`

# libraries
#LIBS = -Wl,-whole-archive -L $(LIBDIR) -l$(PROGNAME) -Wl,-no-whole-archive
LIBS += $(USERLIBS)
LIBS += `. $(MARLINWORKDIR)/bin/marlin_libs.sh`

# linker
LD = $(CXX)

# linker flags
LDFLAGS = -Wl,--demangle

# archiever
AR = ar cru

# check if a directory exists
DIR_EXISTS = test -d

.PHONY: all lib doc clean distclean

all: lib

ifdef BINSRCS
.PHONY: bin
bin: lib $(BINOBJ) $(BIN)
endif

lib: $(LIB)

doc:
	@$(DIR_EXISTS) $(DOCDIR) || mkdir -p $(DOCDIR) ; \
	export DOX_PROGNAME=$(PROGNAME) ; \
	export DOX_DOCDIR=$(DOCDIR) ; \
	cd ../doc ; doxygen ;

clean:
	rm -f $(OBJ) $(BINOBJ) $(BIN) $(LIB)

distclean:
	rm -rf $(BIN) $(LIB) $(TMPDIR) $(DOCDIR)

# tell make to include dependencies
-include $(DEP)

# binary rule
$(BIN) : $(OBJ) $(BINOBJ)
	@$(DIR_EXISTS) $(BINDIR) || mkdir -p $(BINDIR)
	$(LD) $(OBJDIR)/$(notdir $@).o $(LDFLAGS) $(LIBS) -o $@

# objects rule
$(OBJDIR)/%.o: %.cc
	@$(DIR_EXISTS) $(OBJDIR) || mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# library rule
$(LIB) : $(OBJ)
	@$(DIR_EXISTS) $(LIBDIR) || mkdir -p $(LIBDIR)
	$(AR) $@ $?

# dependency rule
$(DEPDIR)/%.d : %.cc
	@$(DIR_EXISTS) $(DEPDIR) || mkdir -p $(DEPDIR)
	@echo "Creating dependencies for $<" ; \
	set -e; $(CXX) -MM $(CPPFLAGS) $< \
	| sed 's|\($*\)\.o[ :]*|$(OBJDIR)/\1.o $(DEPDIR)/$*.d : |g'\
	> $@; [ -s $@ ] || rm -f $@

# end
