FC = gfortran
FFLAGS = -ffixed-line-length-none
 
CC = gcc
CFLAGS = -O3 -fomit-frame-pointer -ffast-math -DHAVE_UNDERSCORE
 
.PHONY: all default force mostlyclean clean distclean

all: libDD.a 
	cp libDD.a ../.
	
SRC = \
  DD_Main.F DD_Looptest.F DD_Vertex.F DD_VertexCounter.F DD_Tree.F DD_Prop.F DD_PropCounter.F DD_Box.F DD_SetCounterterms.F DD_SetCouplings.F DD_Init.F DD_SetCouplingsHO.F DD_EFTVertex.F DD_Running.F DD_EFTVertexCounter.F DD_SigmaNucleon.F DD_PrintModel.F

LIS = $(wildcard *.F)

LIBS = ../util.a 

$(LIBS): force
	$(MAKE) -C $*

clean:
	rm *.a

libDD.a: $(SRC) $(LIBS)
	$(FC) -c $(FFLAGS) -o $@ $<


