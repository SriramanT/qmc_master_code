FLAGS=$(FC_FLAGS) -I../../SRC

.SUFFIXES : .mod .o .F90

test:   ../../libdqmc.a ggeom.F90
	$(FC) $(FLAGS) -o ggeom ggeom.F90 ../../libdqmc.a $(LIB)

clean:
	rm -f ggeom
