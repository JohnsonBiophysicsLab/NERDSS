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
	ODIR := $(ODIR)/serial
endif

ifneq (,$(filter mpi,$(MAKECMDGOALS)))
	_EXEC = nerdss_mpi
	DEFS = -Dmpi_
	INCLUDE_FOLDERS += debug io_mpi mpi
	ODIR := $(ODIR)/mpi
endif

ifneq (,$(filter clean,$(MAKECMDGOALS)))
	MAKECMDGOALS = dummy
endif

# debug and profile get their own object trees for the same reason serial and
# mpi do: they change the flags without changing the sources, so make would
# otherwise reuse whatever the last build left behind.  debug is the sharper
# case -- it adds -fsanitize=address, and mixing instrumented objects with
# uninstrumented ones is its own class of run-time nonsense.
ifneq (,$(filter debug,$(MAKECMDGOALS)))
	ENABLE_DEBUG = true
	ODIR := $(ODIR)-debug
endif

ifneq (,$(filter profile,$(MAKECMDGOALS)))
	ENABLE_PROFILING = true
	ODIR := $(ODIR)-profile
endif

SRCS = $(foreach dir,$(INCLUDE_FOLDERS),$(wildcard $(SDIR)/$(dir)/*.cpp))
EXEC = $(patsubst %,$(BDIR)/%,$(_EXEC))

OS    := $(shell uname)
INTEL = $(shell type icpc  >/dev/null 2>&1; echo $$?)
GCC   = $(shell type g++   >/dev/null 2>&1; echo $$?)

INCS    = $(shell gsl-config --cflags) -Iinclude
CXXFLAGS = -std=c++0x
LIBS     = $(shell gsl-config --libs)
# Emit a .d file beside each .o listing the headers it included, so that editing
# a header triggers a rebuild of the objects that use it.  -MP adds phony targets
# for those headers so a deleted or renamed header does not wedge the build.
DEPFLAGS = -MMD -MP

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
# Each mode gets its own object tree (obj/serial, obj/mpi).  Both modes compile
# the same sources, but -Dmpi_ changes what those sources mean: it selects a
# different subvolume grid in class_SimulVolume.cpp, a different free-protein
# count in check_dissociation_implicitlipid.cpp, and -- the dangerous one --
# adds MPI_Request/MPI_Status members to MpiContext in include/split.cpp.
#
# With a single shared obj/, the second build reused the first mode's objects,
# because their .cpp files had not changed since those .o files were written.
# `make serial && make mpi` therefore linked an MPI binary in which only the 15
# MPI-only sources knew about -Dmpi_ and the other 208 did not.  The two halves
# disagreed about the layout of MpiContext -- xOffset sits at byte 56 without
# mpi_ and at byte 136 with it -- so the parallel cell assignment in
# class_SimulVolume.cpp read domain bounds from the wrong offsets.  It compiled,
# it linked, it ran, and it reported zero of every species.
OBJS = $(patsubst $(SDIR)/%.cpp,$(ODIR)/%.o,$(SRCS))

# ---------------- RULES
syntax:
	@echo "------------------------------------"
	@printf '\033[31m%s\033[0m\n' " USAGE: make serial|mpi [debug] [profile]"
	@echo "------------------------------------"
	exit 0

$(MAKECMDGOALS): $(EXEC)
	@echo "Finished making (re-)building $(MAKECMDGOALS) version, $(EXEC)."

# The executable source is a prerequisite, not just an argument on the recipe
# line.  Without it listed here, make compared bin/nerdss only against the
# objects, so editing EXEs/nerdss.cpp -- which holds the whole timestep loop --
# rebuilt nothing and left the previous binary in place.  A build would report
# success and silently produce the old program, which is the worst possible
# failure for a repository whose claims rest on measured timings.
# DEPFLAGS is applied here too, for the same reason the object rule needs it:
# the link step compiles a translation unit, and that unit includes headers.
$(EXEC): $(OBJS) $(EDIR)/$(_EXEC).cpp
	@echo "Compiling $(EDIR)/$(@F).cpp"
	$(CC) $(CFLAGS) $(CXXFLAGS) $(INCS) $(PROF) $(DEPFLAGS) -MF $(ODIR)/$(@F).d -o $@ $(EDIR)/$(@F).cpp $(OBJS) $(LIBS) $(PLANG) $(DEFS)
	@echo "------------"

$(ODIR)/%.o: $(SDIR)/%.cpp
	@echo "Compiling $< to $@"
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(CXXFLAGS) $(INCS) $(PROF) $(DEPFLAGS) -c $< -o $@ $(PLANG) $(DEFS)
	@echo "------------"

clean:
	rm -rf $(ODIR) bin

# ---------------- HEADER DEPENDENCIES
# Without these, editing a header rebuilt nothing: make only compared each .o
# against its .cpp.  Adding a member to a struct in a header therefore produced
# a binary in which some translation units used the new layout and the rest
# still used the old one -- it linked, and then misbehaved at run time.  The
# generated .d files list every header each object actually included, so a
# header edit now recompiles exactly the objects that read it.
-include $(OBJS:.o=.d)

# The same, for the executable's own translation unit.  Its .d is written into
# $(ODIR) so `make clean` removes it with everything else.
-include $(ODIR)/$(_EXEC).d


# Reference: https://www.gnu.org/software/make/manual/html_node/Quick-Reference.html
#            https://www.gnu.org/software/make/
#            https://www.cmcrossroads.com/article/basics-vpath-and-vpath
#            https://www.gnu.org/software/make/manual/html_node/Implicit-Variables.html