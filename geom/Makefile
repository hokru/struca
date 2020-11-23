# program name
PROG     = libgeom.a

# compiler options
# FC= gfortran #-static
# FFLAGS = -Og -g -ffree-line-length-none -fbounds-check
LINKER  = ar

OBJS=main.o geometry.o string.o tools.o ebe.o

%.o: %.f
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

# link
$(PROG): $(OBJS)
	$(LINKER) rc $(PROG) $(OBJS) 
	ranlib $(PROG)

clean:
	rm -f *.o $(PROG)
