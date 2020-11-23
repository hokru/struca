
 PROG = ~/bin/struca

#  MKLROOT=/home/kruse/intel/oneapi/mkl/2021.1-beta09/
  FC = gfortran #-static 
#   FC = ifort
  FLAGS= -O -fbounds-check -ffree-line-length-none -m64 
  LIBS= -llapack -lblas 
#  LIBS = LIBS= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
#  LIBS= -L$(OPENBLAS)/lib/ -lopenblas -lpthread
#  LIBS= -L$(OPENBLAS)/ -lopenblas -lpthread

 SOURCES=\
 modules.f90\
 align.f90\
 dist.f90\
 intcoords.f90\
 main.f90\
 prtim.f90\
 rmsd.f90\
 version.f90\
 bond_matrix.f90\
 eval_opt.f90\
 io.f90\
 math.f90\
 molecule.f90\
 rdf.f90\
 string.f90\
 hbonds.f90\
 single_intcoords.f90\
 pdbread.f90

OBJS=$(SOURCES:.f90=.o)


# SIMPLE INSERT FOR GEOMERTY LIBRARY
GEOM_DIR=./geom/
LIB_GEOM=$(GEOM_DIR)/libgeom.a
TMP=main.o geometry.o tools.o ebe.o
OBJS_GEOM+=$(foreach dir, $(TMP), $(addprefix $(GEOM_DIR), $(dir)))

LIBS+=$(LIB_GEOM)

BUILID:=$(shell date)
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
$(info Building: $(GIT_VERSION))



# targets:
.PHONY: all
.PHONY: clean
.PHONY: libgeom

all: version $(PROG) $(GEOM)
geom: $(GEOM)

version:
	@touch version.f90
	@echo 'writing new version.f90'
	@echo "subroutine version" > version.f90
	@echo     " print*,  'Build info:'" >> version.f90
	@echo     " print*,  ' github        : https://github.com/hokru/struca '">> version.f90
	@echo     " print*,  ' build date    : $(BUILID)     '"      >> version.f90
	@echo     " print*,  ' git version   : $(GIT_VERSION)'"      >> version.f90
	@echo     " print*,  ' '          "     >> version.f90
	@echo "end subroutine" >> version.f90


%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FLAGS) -c $< -o $@


$(PROG):$(OBJS) libgeom
	$(FC) $(LINK) $(OBJS) $(LIBS) -o $(PROG)

libgeom:
	FC='$(FC)' FFLAGS='$(FLAGS)' make -C $(GEOM_DIR)
# 	ar rc $(LIB_GEOM) $(OBJS_GEOM) 
# 	ranlib $(PROG)
	

clean:
	rm -f *.o $(GEOM_DIR)/*.o *.mod $(PROG) $(LIB_GEOM)

