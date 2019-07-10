.SUFFIXES :
.SUFFIXES : .f .f90 .o
FC = lf95 
#FFLAGS=  -g -C  -std95 -warn all
#FFLAGS=  -g -C  -w95 
FFLAGS=  --nfix --o2
#FFLAGS = --chkglobal --f95 --g --in --maxfatals 10 --trace --trap --nfix
# --in is for implicit none 
#FFLAGS = --chkglobal --f95 -g --maxfatals 10 --trace --trap --nfix
#FFLAGS77=   -g  -w95
#FFLAGS77=   -O3 -cm -vec-report0 -w95
OBJF = header.o  he3props.o dlinpk.o dwtmeta.o
rg33opt : he4state.o  globmod.o compute_mod.o output_mod.o  input_mod.o  main.o \
	  he4props.o $(OBJF) 
	$(FC) $(FFLAGS)   he4state.o globmod.o he4props.o compute_mod.o  \
	input_mod.o  output_mod.o main.o $(OBJF) -o rg33opt
main.o $(OBJF) input_mod.o  output_mod.o compute_mod.o : globmod.o 
compute_mod.o  : he4props.o
he4props.o  :  he4state.o
main.o $(OBJF) input_mod.o  output_mod.o  :   compute_mod.o
input_mod.o : output_mod.o 
main.o : input_mod.o
dlinpk.o : dlinpk.f90
	$(FC) $(FFLAGS) -c dlinpk.f90 
.f.o :
	$(FC) $(FFLAGS) -FR -c  $<
.f90.o :
	$(FC) $(FFLAGS)  -c  $<
