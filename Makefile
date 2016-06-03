
  PROG = ~/bin/struca

#  OBJS=modules.o main.o io.o  string.o eval_opt.o dist.o prtim.o
  OBJS=modules.o main.o io.o string.o dist.o rmsd.o prtim.o math.o intcoords.o bond_matrix.o align.o 

#  FC = gfortran -fopenmp
  FC = gfortran -fopenmp
  LINK =  
  FLAGS= -O3 -fbounds-check -ffree-line-length-none -m64 
#  FLAGS= -Og -g -fbounds-check -ffree-line-length-none -m64 
  LIBS= -llapack -lblas

# targets:
.PHONY: all
.PHONY: clean


%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FLAGS) -c $< -o $@

$(PROG):$(OBJS) 
	$(FC) $(LINK) $(OBJS) $(LIBS) -o $(PROG)


clean:
	rm -f *.o *.mod $(PROG) 

