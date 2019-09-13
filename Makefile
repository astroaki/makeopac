all:		makeopac_henyey makeopac_scatmat
makeopac_henyey:	make_ca_cs_g.o bhmie.o Makefile
		gfortran make_ca_cs_g.o bhmie.o -o makeopac_henyey
make_ca_cs_g.o:	make_ca_cs_g.f90 Makefile
		gfortran -c make_ca_cs_g.f90
makeopac_scatmat:	make_scatmat_distr.o bhmie_anggrid.o Makefile
		gfortran make_scatmat_distr.o bhmie_anggrid.o -o makeopac_scatmat
make_scatmat.o:	make_scatmat.f90 Makefile
		gfortran -c make_scatmat.f90
make_scatmat_distr.o:	make_scatmat_distr.f90 Makefile
		gfortran -c make_scatmat_distr.f90
bhmie.o:	bhmie.f Makefile
		gfortran -c bhmie.f
bhmie_anggrid.o:	bhmie_anggrid.f Makefile
		gfortran -c bhmie_anggrid.f

clean:
	@rm -f	*.o *.mod *~
	@echo OBJECT and MODULE files removed.

cleanall:
	@rm -f	*.o *.mod *~ dustkapscatmat*.inp makeopac_henyey makeopac_scatmat *.out
