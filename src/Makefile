# This is a -*- Makefile -*- 
# Warning: this file contains tabs which cannot be converted to spaces

# TOP_SRCDIR must be relative
TOP_SRCDIR = ../..
SRCDIR = .

#CONFIG=debug

#CC	=	$(CROSS)gcc
CXX	=	$(CROSS)g++
#CPP	=	$(CROSS)gcc -E
LD	=	$(CXX)

# CFLAGS_OPT is normally in the config$(CONFIG).mk file,
# but here's a reasonable default value
# add -ftree-vectorizer-verbose=3 to get information about vectorized loops
#CFLAGS_OPT=-Wall -g -O2 -ftree-vectorize -msse3 -mssse3 -ffast-math
#CFLAGS_OPT=-Wall -O3 -DNDEBUG -march=core2 -ftree-vectorize -msse3 -mssse3 -ffast-math
#CFLAGS_OPT=-Wall -g -O2 -march=pentium4 -ftree-vectorize -msse3 -mssse3 -ffast-math -DDISPARITYTAGGER_DISABLE_TEXTURE_RECTANGLE -DDISPARITYTAGGER_DISABLE_PBO
#CFLAGS_OPT=-Wall -g -O2 -march=pentium4 -ftree-vectorize -msse3 -mssse3 -ffast-math -DVIDEOINPUT_WITHOUT_NVSDI
CFLAGS_OPT= -Wall  -Wmissing-declarations
#CFLAGS_OPT= -Wall -g -O2
CFLAGS_OPENMP = -fopenmp

# the flag for .cpp files
#CXXFLAGS = -std=c++11

PNG_CPPFLAGS = $(shell pkg-config libpng --cflags)
PNG_LDFLAGS = $(shell pkg-config libpng --libs-only-L)
PNG_LIBS = $(shell pkg-config libpng --libs-only-l) $(shell pkg-config libpng --libs-only-other)
JPEG_CPPFLAGS =
JPEG_LDFLAGS =
JPEG_LIBS = -ljpeg
TIFF_CPPFLAGS = $(shell pkg-config libtiff-4 --cflags)
TIFF_LDFLAGS = $(shell pkg-config libtiff-4 --libs-only-L)
TIFF_LIBS = $(shell pkg-config libtiff-4 --libs-only-l) $(shell pkg-config libtiff-4 --libs-only-other)

LIBTOOL = libtool

-include config$(CONFIG).mk

CPPFLAGS= \
	$(PNG_CPPFLAGS) \
	$(JPEG_CPPFLAGS) \
	$(TIFF_CPPFLAGS)

LIBS = \
	$(PNG_LDFLAGS) $(PNG_LIBS) \
	$(JPEG_LDFLAGS) $(JPEG_LIBS) \
	$(TIFF_LDFLAGS) $(TIFF_LIBS)


CFLAGS_EXTRA=$(CFLAGS_OPENMP)
CXXFLAGS=$(CFLAGS_OPT) $(CFLAGS_EXTRA)
CFLAGS=$(CFLAGS_OPT) $(CFLAGS_EXTRA)
LDFLAGS=$(CFLAGS_OPT) $(CFLAGS_EXTRA) $(LDFLAGS_EXTRA) $(OPENCL_LDFLAGS)

# for gcov profiling add:
#-fprofile-arcs -ftest-coverage

COMPILE.c=$(CC) -c $(CFLAGS) $(CPPFLAGS)
COMPILE.cpp=$(CXX) -c $(CXXFLAGS) $(CPPFLAGS)

PROGRAMS = brox_spatial brox_temporal

LIBRARY = 

#VPATH = $(SRCDIR):$(SRCDIR)/$(TOP_SRCDIR)/tools/log:$(SRCDIR)/$(TOP_SRCDIR)/tools/iniParser:$(SRCDIR)/$(TOP_SRCDIR)/tools/files

.PHONY: all lib

all: $(LIBRARY) $(PROGRAMS)

lib: $(LIBRARY)

brox_spatial_SOURCES_CPP = \
	brox_optic_flow_spatial.cpp \
	brox_spatial_main.cpp \
	bicubic_interpolation.cpp \
	brox_spatial_mask.cpp \
	gaussian.cpp \
	iio.cpp \
	mask.cpp \
	zoom.cpp \
	xmalloc.cpp
brox_spatial_SOURCES_C =
brox_spatial_HEADERS = \
	bicubic_interpolation.h \
	brox_optic_flow.h \
	brox_spatial_mask.h \
	gaussian.h \
	iio.h \
	mask.h \
	zoom.h \
	xmalloc.h
brox_spatial_OBJS =  $(brox_spatial_SOURCES_CPP:.cpp=.o) $(brox_spatial_SOURCES_C:.c=.o)
brox_spatial_LIBS = $(LIBS)

brox_temporal_SOURCES_CPP = \
	brox_optic_flow_temporal.cpp \
	brox_temporal_main.cpp \
	bicubic_interpolation.cpp \
	brox_temporal_mask.cpp \
	gaussian.cpp \
	iio.cpp \
	mask.cpp \
	zoom.cpp \
	xmalloc.cpp
brox_temporal_SOURCES_C =
brox_temporal_HEADERS = \
	bicubic_interpolation.h \
	brox_optic_flow.h \
	brox_temporal_mask.h \
	gaussian.h \
	iio.h \
	mask.h \
	zoom.h \
	xmalloc.h
brox_temporal_OBJS =  $(brox_temporal_SOURCES_CPP:.cpp=.o) $(brox_temporal_SOURCES_C:.c=.o)
brox_temporal_LIBS = $(LIBS)


SRCS_CPP = \
	$(brox_spatial_SOURCES_CPP) \
	$(brox_temporal_SOURCES_CPP)

SRCS_C = \
	$(brox_spatial_SOURCES_C) \
	$(brox_temporal_SOURCES_C)

HEADERS = \
	$(brox_spatial_HEADERS) \
	$(brox_temporal_HEADERS)

brox_spatial: $(brox_spatial_OBJS) $(LIBRARY)
	$(LD) -o $@ $^ $(LDFLAGS) $(brox_spatial_LIBS) $(LDADD)

brox_temporal: $(brox_temporal_OBJS) $(LIBRARY)
	$(LD) -o $@ $^ $(LDFLAGS) $(brox_temporal_LIBS) $(LDADD)


#$(LIBRARY): $(util_OBJS)
#	$(AR) rvu $@ $^
#	$(RANLIB) $@

.SUFFIXES: .c .o .cpp

## gcc-only version:
%.o : %.c
	$(COMPILE.c) -MD -o $@ $<
	@cp $*.d $*.P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	    rm -f $*.d

%.o : %.cpp
	$(COMPILE.cpp) -MD -o $@ $<
	@cp $*.d $*.P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	    rm -f $*.d

## general version:
# %.o : %.c
# 	@$(MAKEDEPEND); \
# 	    cp $*.d $*.P; \
# 	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
# 		-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
# 	    rm -f $*.d
# 	$(COMPILE.c) -o $@ $<

# %.o : %.cpp
# 	@$(MAKEDEPEND); \
# 	    cp $*.d $*.P; \
# 	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
# 		-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
# 	    rm -f $*.d
# 	$(COMPILE.cpp) -o $@ $<

.PHONY: clean distclean
clean:
	-rm -f $(PROGRAMS) $(LIBRARY) *.o  *~
	-rm -rf *.dSYM

distclean: clean
	-rm -f $(SRCS_CPP:.cpp=.P) $(SRCS_C:.c=.P)

count:
	 wc -l $(SRCS_CPP) $(HEADERS)

-include $(SRCS_CPP:.cpp=.P) $(SRCS_C:.c=.P)
