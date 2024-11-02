#
# Makefile for data analysis software v1.6
# MAC 981008
#

FC = ifc
#FFLAGS = -fpe4 -vms -check bounds -g -p
FFLAGS = -quiet #-fpe4 -vms -check bounds -g 

CC = cc
CFLAGS = -O 
#CFLAGS = -O -g3
#CFLAGS= -O -p

LIBDIR = -L/home/sfegan/Software/Numerical/CERNLIB/2003/lib \
	 -L/home/sfegan/Devel/Dublin/gdf/v083
 
Cernlib    = -lgraflib -lgrafX11 -lpacklib -lmathlib -lkernlib \
             -lgenlib -lpawlib -lgrafGKS -lbvsl -lphtools
GDFlib     = -lgdf
Whipplelib = #-lwhipple -lrecipes
Xlib       = -lX11 -lDXm -lXaw -lXm
Clib       = -lc -lm 
LIBS       = $(GDFlib) $(Cernlib) $(Whipplelib) $(Xlib) $(Clib)


all: fz2red gcpeds gn2gains gparamdat rmdbe lsdbe gparamdat_331

fz2red: fz2red.o Makefile whipClib.o 
	$(FC) fz2red.o whipClib.o $(FFLAGS) -o fz2red \
	$(LIBDIR) $(LIBS)
fz2red.o: Makefile fz2red.f

gcpeds: gcpeds.o Makefile whipClib.o
	$(FC) gcpeds.o whipClib.o $(FFLAGS) -o gcpeds \
	$(LIBDIR) $(LIBS)
gcpeds.o: Makefile gcpeds.f

gn2gains: gn2gains.o Makefile whipClib.o
	$(FC) gn2gains.o whipClib.o $(FFLAGS) -o gn2gains \
	$(LIBDIR) $(LIBS)
gn2gains.o: Makefile gn2gains.f

gparamdat: gparamdat.o Makefile whipClib.o
	$(FC) gparamdat.o whipClib.o $(FFLAGS) -o gparamdat \
	$(LIBDIR) $(LIBS)
gparamdat.o: Makefile gparamdat.f

gparamdat_331: gparamdat_331.o Makefile whipClib.o
	$(FC) gparamdat_331.o whipClib.o $(FFLAGS) -o gparamdat_331 \
	$(LIBDIR) $(LIBS)
gparamdat_331.o: Makefile gparamdat_331.f

rmdbe: rmdbe.o Makefile whipClib.o
	$(CC) rmdbe.o whipClib.o $(CFLAGS) -o rmdbe -L/lib $(Clib)
rmdbe.o: Makefile rmdbe.c

lsdbe: lsdbe.o Makefile whipClib.o
	$(CC) lsdbe.o whipClib.o $(CFLAGS) -o lsdbe -L/lib $(Clib)
lsdbe.o: Makefile rmdbe.c

whipClib.o: Makefile whipClib.c

clean:
	rm -f *~ 

# end of Makefile
