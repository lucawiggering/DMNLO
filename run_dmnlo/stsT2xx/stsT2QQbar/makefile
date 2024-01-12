FC = gfortran
FFLAGS = -ffixed-line-length-none -std=legacy
 
CC = gcc-11
CFLAGS = -O3 -fomit-frame-pointer -ffast-math -DHAVE_UNDERSCORE
 
.PHONY: all default force mostlyclean clean distclean

all: libstsT2QQbar.a 
	cp libstsT2QQbar.a ../../.
	
SRC = \
  stsT2QQbar_Main.F

LIS = $(wildcard *.F)

LIBS = ../../util.a 

$(LIBS): force
	$(MAKE) -C $*

clean:
	rm *.a

libstsT2QQbar.a: $(SRC) $(LIBS)
	$(FC) -c $(FFLAGS) -o $@ $<


