# $Id: makefile,v 1.18.2.29 2008/06/11 19:44:16 jqi Exp $
# $Name: New_Input $
# $Revision: 1.18.2.29 $
#-----------BEGIN MAKEFILE-------------------------------------------------
            SHELL         = /bin/sh
            EXEC          = fvcom

#--------------------------------------------------------------------------
#  Preprocessing and Compilation Directives
#--------------------------------------------------------------------------
.SUFFIXES: .o .f90 .F .F90 

.c.o:
	$(CC) -c $(CFLAGS) $(INCS) $*.c

.F.o:
	$(CPP) $(CPPARGS) $(PETSC_FC_INCLUDES) $*.F > $*.f90
	$(FC)  -c $(FFLAGS) $(INCS) $*.f90
	\rm $*.f90

#.F.o:
#	$(CPP) $(CPPARGS) $(INCS) $*.F > $*.f90
#	$(FC)  -c $(FFLAGS) $(INCS) $*.f90
#	\rm $*.f90
.F90.o:
	$(CPP) $(CPPFLAGS) $(CPPDEFS) $(INCLDIR) $< > $*.f  
	$(FC) -c $(FFLAGS) $(FIXEDFLAGS) $(INCS) $(INCLDIR) $*.f  
	\rm $*.f90


#--------------------------------------------------------------------------
#  FVCOM Source Code.
#--------------------------------------------------------------------------

# Note: Change the order of modules at your own risk!

MODS  = mod_prec.F	mod_types.F	mod_time.F	mod_main.F	\
	mod_utils.F	mod_clock.F	eqs_of_state.F	mod_bulk.F	\
	mod_interp.F	mod_par.F	mod_spherical.F	mod_ncll.F	\
	mod_nctools.F	mod_wd.F	mod_sng.F	mod_force.F	\
	mod_obcs.F	mod_petsc.F	mod_semi_implicit.F	mod_non_hydro.F	\
	mod_input.F 	mod_set_time.F	ice_kinds_mod.F	ice_model_size.F\
	ice_domain.F	ice_constants.F	ice_fileunits.F	ice_state.F	\
	ice_work.F	ice_grid.F	ice_albedo.F	ice_calendar.F	\
	ice_flux.F	ice_flux_in.F	ice_itd.F	ice_mechred.F	\
	ice_itd_linear.F		ice_scaling.F	ice_atmo.F	\
	ice_ocean.F	ice_therm_vertical.F		ice_init.F	\
	ice_therm_itd.F	mod_ice2d.F	mod_ice.F	mod_startup.F	\
	mod_visit.F	mod_wqm.F	mod_ncdio.F	mod_setup.F	\
	particle.F	linklist.F	mod_lag.F	mod_northpole.F \
	mod_dye.F	mod_assim.F	mod_report.F	mod_probe.F	\
	mod_gotm.F	mod_balance_2d.F		mod_tridiag.F	\
	mod_sed.F	mod_enkf.F	mod_etkf.F	mod_rrk.F	\
	mod_rrkassim.F	mod_optimal_interpolation.F	mod_assim_oi.F	\
	mod_nesting.F

MAIN  =	fvcom.F		genmap.F	tge.F		longshore_flow.F\
	cell_area.F 	open_all.F	load_grid.F	allocate_all.F	\
	setup_domain.F	genmap_obc.F	genmap_lsf.F	genmap_nest.F	\
	coords_n_const.F	shape_coef_gcn.F shape_coef_gcy.F	\
	depth_grad.F	grid_metrics.F	cntrl_prmtrs.F	sinter.F	\
	internal_step.F	bcond_gcn.F	bcond_gcy.F	adjust2d3d.F	\
	brough.F	advection_edge_gcn.F	advection_edge_gcy.F	\
	ghostuv.F	advave_edge_gcn.F	advave_edge_gcy.F	\
	phy_baropg.F	baropg.F	external_step.F	extel_edge.F	\
	extuv_edge.F	depth_check.F	vertvl_edge.F	adv_uv_edge_gcn.F\
	adv_uv_edge_gcy.F	vdif_uv.F		wreal.F		\
	viscofh.F	adv_q.F		fct_q2.F	fct_q2l.F	\
	vdif_q.F	adv_t.F		adv_s.F		fct_t.F		\
	vdif_ts.F	vdif_ts_gom.F	fct_s.F		bcond_ts.F	\
	adjust_ts.F	conv_over.F	visitsim.F	startup_type.F	\
	edge_len.F	adcor.F		icing.F		rho_pmean.F	\
	calc_vort.F	namelist.F

COBJS = func_pointer.c


BIOGEN	= mod_bio_3D.F	mod_onedtide.F

 SRCS = $(MODS) $(MAIN) $(BIOGEN)

# object list is composed of free (*.F) and fixed (*.F90) source
 FVCOM_OBJS = $(COBJS:.c=.o) $(SRCS:.F=.o) $(FIXED:.F90=.o) 


#--------------------------------------------------------------------------
#  Linking Directives
#--------------------------------------------------------------------------

$(EXEC):	$(FVCOM_OBJS)
		$(FC) $(FFLAGS) $(VISOPT) $(LDFLAGS) $(FVCOM_OBJS) $(LIBS) -o $(EXEC)


libfvcom: $(FVCOM_OBJS)
	$(AR) $@.a $(FVCOM_OBJS)
	$(RANLIB) $@.a

#	$(FC) -shared -o $@.so $(FVCOM_OBJS)
#       PROBLEMS WITH SHARED OPTION, REQUIRES -fPIC in All sub libraries as well when building on x86_64


libs_ext: 
	(cd libs; 	\
	tar xzf proj.tgz;	\
	tar xzf metis.tgz; 	\
	tar xzf julian.tgz;	\
	tar xzf fproj.tgz;	\
	tar xzf netcdf.tgz)
	\
	(cd libs/netcdf/;\
	./configure CC=$(CC) CFLAGS=-O3 CXX=$(CC) CXXFLAGS=-O3 F77=$(FC) F90=$(FC) FFLAGS=-O3 --prefix=$(INSTALLDIR) --build=$(MACHTYPE); \
	make;	\
	make install)
	\
	(cd libs/proj/;	\
	./configure CC=$(CC) CFLAGS=-O3 CXX=$(CC) CXXFLAGS=-O3 F77=$(FC) FFLAGS=-O3 --prefix=$(INSTALLDIR); \
	make;				\
	make install)		
	\
	(cd libs/fproj/;	\
	./configure CPPFLAGS='$(COMPILER)' CC=$(CC) CFLAGS=-O3 CXX=$(CXX) CXXFLAGS=-O3 FC=$(FC) FFLAGS=-O3 --with-proj4=../install --prefix=$(INSTALLDIR); \
	make;				\
	make install)		
	\
	(cd libs/metis/; make; cp libmetis.a $(INSTALLDIR)/lib/.)
	(cd libs/julian/; make; cp julian.a $(INSTALLDIR)/lib/libjulian.a; cp fjulian.inc $(INSTALLDIR)/include/.; cp *.h $(INSTALLDIR)/include/.)

libs_clean:
		(cd libs; \
		rm -rf proj metis julian fproj netcdf $(INSTALLDIR); )


#--------------------------------------------------------------------------
#  Target to create dependecies.
#--------------------------------------------------------------------------

depend:
		makedepf90  $(SRCS) >> makedepends

#--------------------------------------------------------------------------
#  Target to make all.
#--------------------------------------------------------------------------

all:    
	( cd BIO_source/             ; make )
	make 
	make libfvcom
	( cd testing/nctools         ; make )
	( cd testing/parallel        ; make )
	( cd utilities/cnvrt_coords/ ; make )
	( cd utilities/gen_nml/      ; make )
	( cd utilities/interp_time/  ; make )
	( cd utilities/weather_data/ ; make )
	( cd input/obc/              ; make )
	( cd input/restart_file/     ; make )
	( cd input/rivers/           ; make )
	( cd input/sst/              ; make )
	( cd input/surface_forcing/  ; make )
	( cd input/aux_input/        ; make )
	( cd input/init_file/        ; make )

#--------------------------------------------------------------------------
#  Tar Up Code                           
#--------------------------------------------------------------------------

tarfile:
	tar cvf fvcom.tar *.F *.F90  makefile exa_run.dat makedepends RELEASE_NOTES 

#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

clean:
		/bin/rm -f *.o *.mod *.f90

clobber:	
		make clean
		/bin/rm -f fvcom libfvcom.so libfvcom.a

allclean:
		make clean
		( cd testing/nctools/        ; make clean )
		( cd testing/parallel/       ; make clean )
		( cd utilities/cnvrt_coords/ ; make clean )
		( cd utilities/gen_nml/      ; make clean )
		( cd utilities/interp_time/  ; make clean )
		( cd utilities/weather_data/ ; make clean )
		( cd input/obc/              ; make clean )
		( cd input/restart_file/     ; make clean )
		( cd input/rivers/           ; make clean )
		( cd input/sst/              ; make clean )
		( cd input/surface_forcing/  ; make clean )
		( cd input/aux_input/        ; make clean )
	        ( cd BIO_source/             ; make clean )


allclobber:
		make clean
		/bin/rm -f fvcom libfvcom.so libfvcom.a
		( cd testing/nctools/        ; make clobber )
		( cd testing/parallel/       ; make clobber )
		( cd utilities/cnvrt_coords/ ; make clobber )
		( cd utilities/gen_nml/      ; make clobber )
		( cd utilities/interp_time/  ; make clobber )
		( cd utilities/weather_data/ ; make clobber )
		( cd input/obc/              ; make clobber )
		( cd input/restart_file/     ; make clobber )
		( cd input/rivers/           ; make clobber )
		( cd input/sst/              ; make clobber )
		( cd input/surface_forcing/  ; make clobber )
		( cd input/aux_input/        ; make clobber )
	        ( cd BIO_source/             ; make clobber )

#--------------------------------------------------------------------------
#  Common rules for all Makefiles - do not edit.
#--------------------------------------------------------------------------

emptyrule::

#--------------------------------------------------------------------------
#  Empty rules for directories that do not have SUBDIRS - do not edit.
#--------------------------------------------------------------------------

install::
	@echo "install in $(CURRENT_DIR) done"

install.man::
	@echo "install.man in $(CURRENT_DIR) done"

Makefiles::

includes::
include ./makedepends
include make.inc
