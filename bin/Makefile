SYSOS   := $(shell uname)
CURFWD  := $(shell pwd)
DATE    := $(shell date)
WHO     := $(shell whoami)
HOST    := $(shell hostname)
SHELL   := /bin/sh
#BITS is also used in source.h to define whether code is 32 bit or 64 bit
#if code is 64 bit, must use BITS = -64
BITS = -64

# for the linux boxes in the trailer...
ifeq "$(HOST)" "popo" 
override HOST = Trailer
endif
ifeq "$(HOST)" "colima" 
override HOST = Trailer
endif
ifeq "$(HOST)" "elchichon" 
override HOST = Trailer
endif
ifeq "$(HOST)" "rainier"
override HOST = Trailer
endif
ifeq "$(HOST)" "sthelens"
override HOST = Trailer
endif
ifeq "$(HOST)" "casita" 
override HOST = Trailer
endif

include ../source.h

#Nothing serious. Just provide some build information
#should use for loop

DW2 = $(word 2, $(DATE))
DW3 = $(word 3, $(DATE))
DW4 = $(word 4, $(DATE))
DW6 = $(word 6, $(DATE))
BUILDINFO = LastBuild_$(DW2)_$(DW3)_$(DW4)_$(DW6)_By_$(WHO)


#First clear off, then specify
.SUFFIXES:
.SUFFIXES: .C .f .o

DIRNAME = $(notdir $(CURPWD))


#default config
CC          = mpicc 
CCC         = mpicxx
F77         = mpif90
OPTFLAGS    = 
CFLAGS      = $(OPTFLAGS)
FFLAGS      = $(OPTFLAGS)
INCLUDES    =
SYSLIBS	    = 
LDFLAGS	    = $(OPTFLAGS)
LINKER	    = 
CCC_SUFFIXES = 
F77_SUFFIXES = 
#config for different systems

#1. ganesha, Sun workstation
ifeq "$(SYSOS)" "SunOS"
MPI_HOME   = /eng/mae/ganesha4/mpich-1.2.1
override CCC = $(MPI_HOME)/bin/mpiCC 
override CC  = $(MPI_HOME)/bin/mpicc
override F77 = $(MPI_HOME)/bin/mpif77
ARCH        = solaris
COMM        = ch_shmem
AR          = ar cr
RANLIB      = ranlib
LIB_LIST    = 
override CFLAGS += -DFORTRANUNDERSCORE -DHAVE_STDLIB_H=1 \
                    -DHAVE_UNISTD_H=1 -DHAVE_STDARG_H=1  \
	            -DUSE_STDARG=1 -DMALLOC_RET_VOID=1 
override CCC_SUFFIXES += -mt 
override F77_SUFFIXES +=
override SYSLIBS  = -lF77 -lM77 -lsunmath -L/eng/mae/ganesha3/linlib -lblas #-xlic_lib=sunperf
override LDFLAGS += $(LIB_LIST) 
override LINKER = $(CCC)
#RcLIB = 
#RCLIB = ../cobj/libcobj_$(SYSOS)_C.a
#RFLIB = ../fobj/libfobj_$(SYSOS)_f.a
endif

#2. IBM SP
ifeq "$(SYSOS)" "AIX"
ARCH        = rs6000
COMM        = ch_eui
AR          = ar crl
RANLIB      = ranlib
override OPTFLAGS += -bmaxdata:0x80000000 -bmaxstack:0x10000000 -qarch=pwr3 -qtune=pwr3
override FFLAGS += -qextname
override CC  = mpCC
override CCC = mpCC
override F77 = mpxlf
override CFLAGS += -DFORTRANUNDERSCORE   -DHAVE_STDLIB_H=1 \
                   -DHAVE_UNISTD_H=1 -DHAVE_STDARG_H=1 \
                   -DUSE_STDARG=1 -DUSE_STDARG=1 \
                   -DMALLOC_RET_VOID=1 -DMPI_$(ARCH)
override SYSLIBS += -lc -lC  -lm -lmpi -lblas
override LINKER = $(F77)
endif

#3. joy
ifeq "$(SYSOS)" "IRIX64"
override CFLAGS +=-woff all -LANG:std
override SYSLIBS += -lmpi -lc -lC -lcomplib.sgimath -lnsl -lz -lm -lc -lC -lmpi 
override LINKER = CC
override LDFLAGS = $(BITS) -LANG:std -w 
endif

#4. crosby (SGI O2K)
ifeq "$(HOST)" "crosby.ccr.buffalo.edu"

LIBTOOL = /usr/freeware/bin/libtool 
INCLUDES += -I. 
DEFS += -DHAVE_CONFIG_H
override CFLAG +=  -woff all
override SYSLIBS += -lmpi -lc -lC -lcomplib.sgimath -lnsl -lz -lm -lc -lC -lscs -lmpi++ -lmpi  
override LINKER = CC
override LDFLAGS = $(BITS) -LANG:std -w 
endif

#5. golden
ifeq "$(HOST)" "golden"
override FFLAGS += -m 4
endif

#6. Linux boxes in the trailer
ifeq "$(HOST)" "Trailer"
LIBTOOL = /usr/bin/libtool 
INCLUDES += -I. 
DEFS += -DHAVE_CONFIG_H 

MPI_HOME   = /usr/local/mpich-1.2.4/ch_p4
override CCC = $(MPI_HOME)/bin/mpiCC -w
override CC  = $(MPI_HOME)/bin/mpicc
override F77 = $(MPI_HOME)/bin/mpif77
#override CFLAGS += -DFORTRANUNDERSCORE -DHAVE_STDLIB_H=1                     -DHAVE_UNISTD_H=1 -DHAVE_STDARG_H=1  	    -DUSE_STDARG=1 -DMALLOC_RET_VOID=1 
override SYSLIBS  = -lnsl -lz -lm -lg2c -L/usr/lib -lf77blas -latlas -lc -L$(MPI_HOME)/lib -lmpich -lpmpich

override LDFLAGS += $(LIB_LIST)
override LINKER = /usr/bin/libtool --mode=link g++  #$(CCC)
endif

#7. nash.ccr.buffalo.edu 
ifeq "$(HOST)" "nash.ccr.buffalo.edu"
hdflibdir = /util/hdf/lib
hdfincdir = /util/hdf/include 
hdfhlincdir = /util/hdf_hl/include
hdfhllibdir= /util/hdf_hl/lib
LIBTOOL = /usr/bin/libtool 
INCLUDES += -I. -I$(hdfincdir) -I$(hdfhlincdir)
DEFS += -DHAVE_CONFIG_H 

MPI_HOME   = /util/mpich-gm/gnu/ch_gm
override CCC = $(MPI_HOME)/bin/mpiCC -w
override CC  = $(MPI_HOME)/bin/mpicc
override F77 = $(MPI_HOME)/bin/mpif77
override SYSLIBS  = -lnsl -lz -lm -lg2c -L/util/ATLAS/v3.2.0/lib/Linux_PIIISSE1 -lf77blas -latlas -lc -L$(hdfhllibdir) -lhdf5_hl -L$(hdflibdir) -lhdf5
override LDFLAGS += $(LIB_LIST)
override LINKER = /usr/bin/libtool --mode=link $(CCC)
#override LINKER = /usr/bin/libtool --mode=link g++  #$(CCC)
endif
#8. joplin.ccr.buffalo.edu & joplin.ccr.buffalo.edu
ifeq "$(HOST)" "joplin.ccr.buffalo.edu"
LIBTOOL = /usr/bin/libtool 
INCLUDES += -I.
DEFS += -DHAVE_CONFIG_H 

MPI_HOME   = /util/mpich-gm/gnu/ch_gm
override CCC = $(MPI_HOME)/bin/mpiCC -w
override CC  = $(MPI_HOME)/bin/mpicc
override F77 = $(MPI_HOME)/bin/mpif77
#override CFLAGS += -DFORTRANUNDERSCORE -DHAVE_STDLIB_H=1                     -DHAVE_UNISTD_H=1 -DHAVE_STDARG_H=1  	    -DUSE_STDARG=1 -DMALLOC_RET_VOID=1 
override SYSLIBS  =  -lnsl -lz -lm -lg2c -L/util/atlas/Linux_P4SSE2_2_cn_gnu/ -lf77blas -latlas -lc   #-L$(MPI_HOME)/lib -lampe -llmpe -lfmpich -lmpe -lmpe_nompi -lmpichfsup -lpmpich++ -ltmpe -lmpich -lpmpich

override LDFLAGS += $(LIB_LIST)
override LINKER = /usr/bin/libtool --mode=link $(CCC)
#override LINKER = /usr/bin/libtool --mode=link g++  #$(CCC)
endif

#9. Linux boxes that are standalone
ifeq "$(HOST)" "eruption.home"
hdflibdir = /usr/local/5-1.4.5-linux/lib
hdfincdir = /usr/local/5-1.4.5-linux/include 
hdfhlincdir = /usr/local/hdf5_hl-ia32/include
hdfhllibdir= /usr/local/hdf5_hl-ia32/lib
LIBTOOL = /usr/bin/libtool
INCLUDES += -I. -I$(hdfincdir) -I$(hdfhlincdir)
DEFS += -DHAVE_CONFIG_H

MPI_HOME   = /usr/local/mpich-1.2.5
#MPI_HOME   = /usr/local/mpich2-0.92
override CCC = $(MPI_HOME)/bin/mpiCC -w
#override CCC = $(MPI_HOME)/bin/mpicc -w
override CC  = $(MPI_HOME)/bin/mpicc
override F77 = $(MPI_HOME)/bin/mpif77
override SYSLIBS  =  -lnsl    -L/usr/lib -lf77blas -latlas -lc -L$(MPI_HOME)/lib -lmpich  -L$(hdflibdir) -lhdf5 -lz -lm  -L/usr/lib/gcc-lib/i386-redhat-linux/3.2 -lg2c -lpmpich -L$(hdfhllibdir) -lhdf5_hl

override LDFLAGS += $(LIB_LIST)
override LINKER = /usr/bin/libtool --mode=link  $(CCC)
endif


CCFLAGS     = $(CFLAGS)

CDIRLIST = $(addprefix :$(CSRCDIR),$(CUSERDIR))
FDIRLIST = $(addprefix :$(FSRCDIR),$(FUSERDIR))

VPATH = .$(CDIRLIST):$(FDIRLIST)

COBJ  = $(CSRC:.C=_$(SYSOS)_C.o)
cOBJ  = $(cSRC:.c=_$(SYSOS)_c.o)
FOBJ  = $(FSRC:.f=_$(SYSOS)_f.o)
RCOBJ = $(addprefix $(COBJDIR), $(COBJ))
RcOBJ = $(addprefix $(cOBJDIR), $(cOBJ))
RFOBJ = $(addprefix $(FOBJDIR), $(FOBJ))

export CCC
export CC
export F77
export CCFLAGS
export CFLAGS
export FFLAGS
export SYSOS
export COBJ
export cOBJ
export FOBJ
export VPATH
export BITS
export CCC_SUFFIXES 
export CC_SUFFIXES 
export F77_SUFFIXES
export INCLUDES

hpfem: Allobj
	$(LINKER) $(LDFLAGS) \
	$(RCOBJ) \
	$(RcOBJ) \
	$(RFOBJ) \
	$(SYSLIBS) 

cleanC:   
	rm -f hpfem
	rm -rf ../cobj/*_$(SYSOS)_C.o
cleanf:
	rm -rf ../fobj/*_$(SYSOS)_f.o  
cleanc:
	rm -rf ../ccobj/*_$(SYSOS)_c.o
clean:
	rm -f hpfem
	rm -rf ../cobj/*_$(SYSOS)_C.o
	rm -rf ../fobj/*_$(SYSOS)_f.o  
	rm -rf ../ccobj/*_$(SYSOS)_c.o
	rm -rf ../csrc/*/*~
	rm -rf ../fsrc/*/*~
Allobj:	
	(cd $(FOBJDIR); $(MAKE))
	(cd $(COBJDIR); $(MAKE))
	(cd $(cOBJDIR); $(MAKE))
