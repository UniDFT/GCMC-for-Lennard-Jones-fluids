COMPFLAGS     = -O # -check_bounds  -g -trapeuv

CFLAGS	      = $(COMPFLAGS)
PFLAGS	      = $(COMPFLAGS)
FFLAGS	      = $(COMPFLAGS)
CCFLAGS       = $(COMPFLAGS)
CXXFLAGS      = $(COMPFLAGS)

DEST	      = .

# Uncomment this if you want system header files to be expanded
#
# SYSHDRS       =

EXTHDRS	      =

HDRS	      =

INSTALL	      = install

LD	      = mpif90

LDFLAGS	      = $(COMPFLAGS)

LIBS	      =

LINTLIBS      =

LINTFLAGS     = -u $(CFLAGS)

MAKEFILE      = Makefile

OBJS	      = adjust.o \
		corp.o \
		coru.o \
		ener.o \
		eneri.o \
		lattice.o \
		mc_grand.o \
		mcexch.o \
		mcmove.o \
		ranf.o \
		rantest.o \
		readdat.o \
		sample.o \
		store.o \
		toterg.o

PRINT	      = pr

PRINTFLAGS    =

LP	      = lp

LPFLAGS       = 

PROGRAM       = run

SHELL	      = /bin/sh

SRCS	      = adjust.f \
		corp.f \
		coru.f \
		ener.f \
		eneri.f \
		lattice.f \
		mc_grand.f \
		mcexch.f \
		mcmove.f \
		ranf.f \
		rantest.f \
		readdat.f \
		sample.f \
		store.f \
		toterg.f

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS) $(MAKEFILE)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
		@echo "done"

clean:;		@rm -f $(OBJS) core

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags

depend:;	@mkmf -f $(MAKEFILE) ROOT=$(ROOT)

echo:;		@echo $(HDRS) $(SRCS)

-1:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@-strip $(PROGRAM)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) -f $(DEST) $(PROGRAM)); fi

lint:		$(LINTLIBS) $(SRCS) $(HDRS) $(EXTHDRS)
		lint $(LINTFLAGS) $(LINTLIBS) $(SRCS)

print:;		@$(PRINT) $(PRINTFLAGS) $(HDRS) $(SRCS) | $(LP) $(LPFLAGS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

touch:;		@touch $(SRCS) $(HDRS)
