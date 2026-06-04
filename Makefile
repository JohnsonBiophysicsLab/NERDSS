#
#  Update 2020-01-29:
#  o Uses required argument of serial, omp, or mpi.
#  o Use VPATH for finding cpp file in different directories -- this simplifies rules
#  o Abort if gsl-config isn't available
#  o Fixed (INTEL) compiler search 0=found | 1=notfound ; make conditional simple (ifeq 0|1)
#  o Also use conditional for GCC
#  o For objects, use basename to get base file name
#  o Clean up directory prefix (shorten variable names and group)
#  o Simplified obj and bin rule logic and readability.
#  o put rules in canonical order
#  o Now has PROF for profiling. (This is by default overrided with empty PROF.)
#  o Now uses INCS. CXXFLAGS is used for C++ specific options.
#  o Make executables with suffixes ( nerdss_serial | nerdss_mpi | nerdss_omp).
#  --  a bit cleaner                                            Kent milfeld@tacc.utexas.edu
#
# TODO: use function to create VPATH
# TODO: Fix MPI after learning purpose
# TODO: Make rules for *.hpp's
#
# Set terminal width to 220 to avoid viewing wrapped lines in output. A width of 200 avoids most wrapping.
#
# Update 2025-08-25:
# o dded new PHONY: "debug" and "profile". 
# o Use `make serial debug` to debug with gdb
# o use `make serial profile` to profile
#

BDIR   = bin
ODIR   = obj
SDIR   = src
EDIR   = EXEs

PROF   =

.PHONY: any debug profile clean

# ---------------- REQUIREMENTS: gsl and directories
hasGSL = $(shell type gsl-config >/dev/null 2>&1; echo $$?)
ifeq ($(hasGSL),1)
$(error " GSL must be installed, and gsl-config must be in path.")
else
$(shell mkdir -p bin)
$(shell mkdir -p obj)
endif

# ---------------- EXECUTABLE SETUP
INCLUDE_FOLDERS = boundary_conditions classes error math parser reactions system_setup trajectory_functions io

ifneq (,$(filter serial,$(MAKECMDGOALS)))
	_EXEC = nerdss
endif

ifneq (,$(filter mpi,$(MAKECMDGOALS)))
	_EXEC = nerdss_mpi
	DEFS = -Dmpi_
	INCLUDE_FOLDERS += debug io_mpi mpi
endif

ifneq (,$(filter clean,$(MAKECMDGOALS)))
	MAKECMDGOALS = dummy
endif

ifneq (,$(filter debug,$(MAKECMDGOALS)))
	ENABLE_DEBUG = true
endif

ifneq (,$(filter profile,$(MAKECMDGOALS)))
	ENABLE_PROFILING = true
endif

SRCS = $(foreach dir,$(INCLUDE_FOLDERS),$(wildcard $(SDIR)/$(dir)/*.cpp))
EXEC = $(patsubst %,$(BDIR)/%,$(_EXEC))

OS    := $(shell uname)
INTEL = $(shell type icpc  >/dev/null 2>&1; echo $$?)
GCC   = $(shell type g++   >/dev/null 2>&1; echo $$?)

INCS    = $(shell gsl-config --cflags) -Iinclude
CXXFLAGS = -std=c++0x
LIBS     = $(shell gsl-config --libs)

# ---------------- COMPILER SETUP
PROF   =

ifeq ($(GCC),0)
	CC      = g++
	ifeq (mpi,$(MAKECMDGOALS))
		CC = mpicxx
	endif
	CFLAGS  = -O3 # use -O2 if profiling is confused by optimization
endif

ifeq ($(INTEL),0)
	CC      = icpc
	ifeq (mpi,$(MAKECMDGOALS))
		CC = mpicxx
	endif
	CFLAGS  = -O3 # use -O2 if profiling is confused by optimization
endif

# ---------------- Feature toggles
# Set debug flags if ENABLE_DEBUG is true
ifdef ENABLE_DEBUG
	CFLAGS = -g -O0 -fsanitize=address -fno-omit-frame-pointer
	CXXFLAGS += -DDEBUG
endif

# Set profiling flags if ENABLE_PROFILING is true
ifdef ENABLE_PROFILING
	PROF += -pg
	CFLAGS += -DENABLE_PROFILING 
	LIBS += $(shell pkg-config --libs libprofiler)
endif

# ---------------- OBJECT FILES
OBJS = $(patsubst $(SDIR)/%.cpp,$(ODIR)/%.o,$(SRCS))

# ---------------- RULES
syntax:
	@echo "------------------------------------"
	@printf '\033[31m%s\033[0m\n' " USAGE: make serial|mpi [debug] [profile]"
	@echo "------------------------------------"
	exit 0

$(MAKECMDGOALS): $(EXEC)
	@echo "Finished making (re-)building $(MAKECMDGOALS) version, $(EXEC)."

$(EXEC): $(OBJS)
	@echo "Compiling $(EDIR)/$(@F).cpp"
	$(CC) $(CFLAGS) $(CXXFLAGS) $(INCS) $(PROF) -o $@ $(EDIR)/$(@F).cpp $(OBJS) $(LIBS) $(PLANG)
	@echo "------------"

$(ODIR)/%.o: $(SDIR)/%.cpp
	@echo "Compiling $< to $@"
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(CXXFLAGS) $(INCS) $(PROF) -c $< -o $@ $(PLANG) $(DEFS)
	@echo "------------"

clean:
	rm -rf $(ODIR) bin


# Reference: https://www.gnu.org/software/make/manual/html_node/Quick-Reference.html
#            https://www.gnu.org/software/make/
#            https://www.cmcrossroads.com/article/basics-vpath-and-vpath
#            https://www.gnu.org/software/make/manual/html_node/Implicit-Variables.html