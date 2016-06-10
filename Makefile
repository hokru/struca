
  PROG = ~/bin/struca
#  OBJS=modules.o main.o io.o string.o dist.o rmsd.o prtim.o math.o intcoords.o bond_matrix.o align.o  eval_opt.o molecule.o version.o

  SOURCES=$(wildcard *.f90)
  OBJS=$(SOURCES:.f90=.o)

BUILID:=$(shell date)
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
$(info Building: $(GIT_VERSION))


  FC = gfortran -static 
#  FLAGS= -O3 -ffree-line-length-none -m64 
  FLAGS= -Og -g -fbounds-check -ffree-line-length-none -m64 
  LIBS= -llapack -lblas

# targets:
.PHONY: all
.PHONY: clean

all: version $(PROG)


version:
	@touch version.f90
	@echo 'writing new version.f90'
	@echo "subroutine version" > version.f90
	@echo     " print*,  'Build info:'" >> version.f90
	@echo     " print*,  ' build date    : $(BUILID)     '"      >> version.f90
	@echo     " print*,  ' git version   : $(GIT_VERSION)'"      >> version.f90
	@echo     " print*,  ' '          "     >> version.f90
	@echo "end subroutine" >> version.f90


%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FLAGS) -c $< -o $@


$(PROG):$(OBJS) 
	$(FC) $(LINK) $(OBJS) $(LIBS) -o $(PROG)





clean:
	rm -f *.o *.mod $(PROG) 

