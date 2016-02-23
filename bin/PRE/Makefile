SHELL = /bin/sh
SYSOS := $(shell uname)
HOST := $(shell hostname)

################################
# system specific information
################################

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

# Linux OS in the trailer
ifeq "$(HOST)" "Trailer"

CC = gcc
LIBTOOL = /usr/bin/libtool 
LIBS = -lnsl -lz -lm -lg2c -lc
INCLUDES = -I. -I../../csrc/header
DEFS = -DHAVE_CONFIG_H 
#CCLD = $(CC)
CCLD = g++
II_FILES =
endif
#Linux standalone
ifeq "$(HOST)" "localhost"

CC = gcc
LIBTOOL = /usr/bin/libtool
LIBS = -lnsl -lz -lm -lg2c -lc
INCLUDES = -I.
DEFS = -DHAVE_CONFIG_H
#CCLD = $(CC)
CCLD = g++
II_FILES =
endif

# nash.ccr.buffalo.edu
ifeq "$(HOST)" "nash.ccr.buffalo.edu"

CC = gcc
LIBTOOL = /usr/bin/libtool 
LIBS = -lnsl -lz -lm -lg2c -lc 
INCLUDES = -I. -I../../csrc/header
DEFS = -DHAVE_CONFIG_H
#CCLD = $(CC)
CCLD = g++
II_FILES =
endif

# joplin.ccr.buffalo.edu
ifeq "$(HOST)" "joplin.ccr.buffalo.edu"

CC = gcc
LIBTOOL = /usr/bin/libtool 
LIBS = -lnsl -lz -lm -lg2c -lc
INCLUDES = -I. -I../../csrc/header
DEFS = -DHAVE_CONFIG_H
#CCLD = $(CC)
CCLD = g++
II_FILES = 
endif

# crosby.ccr.buffalo.edu
ifeq "$(HOST)" "crosby.ccr.buffalo.edu"

#CC = cc -w
CC = CC -w
#LIBTOOL = /usr/freeware/bin/libtool
LIBS = -lnsl -lz -lm -lc -lftn
INCLUDES = -I. -I../../csrc/header
DEFS = -DHAVE_CONFIG_H
#CCLD = $(F77)
CCLD = CC
CPPFLAGS = -LANG:std
II_FILES = ./ii_files
endif

################################
# everything below shouldn't
# need to be changed
################################
OPT_FLAGS = -O2
#GisApi_SOURCES = create_funky.c

#CPPFLAGS =
LDFLAGS =
GisApi_OBJECTS = create_funky.o GisApi.o  GisAscFile.o  GisRasterHdr.o  GisBinFile.o  preproc.o
GisApi_LDFLAGS =
CFLAGS = -O3
COMPILE = $(CC) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
LTCOMPILE = $(CC) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS) 
LINK = $(LIBTOOL) --mode=link $(CCLD) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@

#SOURCES = $(GisApi_SOURCES)
OBJECTS = $(GisApi_OBJECTS)

.c.o:
	$(COMPILE) $(OPT_FLAGS) -c $<

.C.o: 
	$(COMPILE) $(OPT_FLAGS) -c $< 

.f.o:
	$(F77) $(OPT_FLAGS) -c $<

preproc.x: $(GisApi_OBJECTS)
	@rm -f
	$(LINK) $(GisApi_LDFLAGS) $(GisApi_OBJECTS) $(LIBS)
	@rm -rf .libs

GisAscFile.o: ../../csrc/geoflow/GisAscFile.C ../../csrc/header/GisAscFile.h
	 $(COMPILE) $(OPT_FLAGS) -c ../../csrc/geoflow/GisAscFile.C

GisRasterHdr.o: ../../csrc/geoflow/GisRasterHdr.C ../../csrc/header/GisRasterHdr.h \
		    ../../csrc/header/GisAscFile.h
	 $(COMPILE) $(OPT_FLAGS) -c ../../csrc/geoflow/GisRasterHdr.C

GisBinFile.o: ../../csrc/geoflow/GisBinFile.C ../../csrc/header/GisBinFile.h
	 $(COMPILE) $(OPT_FLAGS) -c ../../csrc/geoflow/GisBinFile.C

GisApi.o: ../../csrc/geoflow/GisApi.C ../../csrc/header/GisApi.h \
	    ../../csrc/header/GisAscFile.h ../../csrc/header/GisRasterHdr.h ../../csrc/header/GisBinFile.h
	 $(COMPILE) $(OPT_FLAGS) -c ../../csrc/geoflow/GisApi.C

create_funky.o: create_funky.C GisApi.o ../../csrc/header/GisApi.h preproc.o

preproc.o: preproc.f

clean:
	rm -rf *.o *~ preproc.x .libs pile.data frict.data funky.dat $(II_FILES)
