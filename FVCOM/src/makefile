#/===========================================================================/
# Copyright (c) 2007, The University of Massachusetts Dartmouth 
# Produced at the School of Marine Science & Technology 
# Marine Ecosystem Dynamics Modeling group
# All rights reserved.
#
# FVCOM has been developed by the joint UMASSD-WHOI research team. For 
# details of authorship and attribution of credit please see the FVCOM
# technical manual or contact the MEDM group.
#
# 
# This file is part of FVCOM. For details, see http://fvcom.smast.umassd.edu 
# The full copyright notice is contained in the file COPYRIGHT located in the 
# root directory of the FVCOM code. This original header must be maintained
# in all distributed versions.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED.  
#
#/---------------------------------------------------------------------------/
# CVS VERSION INFORMATION
# $Id$
# $Name$
# $Revision$
#/===========================================================================/

include ./make.inc

#-----------BEGIN MAKEFILE-------------------------------------------------
            SHELL         = /bin/sh
            EXEC          = fvcom

# Default target
default: $(EXEC)

#--------------------------------------------------------------------------
#  Preprocessing and Compilation Directives
#--------------------------------------------------------------------------
.SUFFIXES: .o .f90 .F .F90 

.c.o:
	$(CC) -c $(CFLAGS) $(INCS) $*.c

.F.o:
	$(CPP) $(CPPARGS) $(if $(PETSC_FC_INCLUDES), $(PETSC_FC_INCLUDES) $(INCS)) $*.F > $*.f90
	$(FC)  -c $(FFLAGS) $(INCS) $*.f90
	$(RM) $*.f90

.F90.o:
	$(CPP) $(CPPFLAGS) $(CPPDEFS) $(INCLDIR) $< > $*.f  
	$(FC) -c $(FFLAGS) $(FIXEDFLAGS) $(INCS) $(INCLDIR) $*.f  
	$(RM) $*.f90

.f90.o:
	$(FC) -c $(FFLAGS) $(INCS) $*.f90


#--------------------------------------------------------------------------
#  FVCOM Source Code.
#--------------------------------------------------------------------------

# Note: Change the order of modules at your own risk!

MODS  = mpi_constants.f90 mod_prec.F	sinter.F	mod_types.F	mod_time.F	\
	mod_main.F	mod_spherical.F	mod_utils.F	mod_clock.F	\
	eqs_of_state.F	mod_interp.F	mod_par.F	\
	mod_par_special.F		mod_ncll.F	mod_nctools.F	\
	mod_wd.F	mod_sng.F	mod_heatflux.F	mod_solar.F	\
	mod_bulk.F	mod_input.F	mod_force.F	mod_obcs.F	mod_petsc.F	\
	mod_semi_implicit.F		mod_non_hydro.F	mod_set_time.F	\
	ice_kinds_mod.F	ice_model_size.F		ice_domain.F	\
	ice_constants.F	ice_fileunits.F	ice_state.F	ice_work.F	\
	ice_grid.F	ice_albedo.F	ice_calendar.F	ice_flux.F	\
	ice_flux_in.F	ice_itd.F	ice_mechred.F	ice_itd_linear.F\
	ice_scaling.F	ice_atmo.F	ice_ocean.F	ice_therm_vertical.F\
	ice_init.F	ice_therm_itd.F	mod_ice2d.F	mod_ice.F	\
	mod_startup.F	mod_wqm.F	mod_ncdio.F	mod_setup.F	\
	mod_newinp.F	particle.F	linklist.F	mod_lag.F	\
	mod_northpole.F mod_pwp.F	mod_dye.F	mod_optimal_interpolation.F\
	mod_report.F	mod_probe.F	mod_gotm.F	mod_balance_2d.F\
	mod_tridiag.F	mod_scal.F	mod_meanflow.F	mod_obcs2.F	\
	mod_obcs3.F	mod_sed.F	mod_enkf.F	mod_etkf.F	\
	mod_rrk.F	mod_rrkf_obs.F	mod_rrkassim.F	mod_enkf_ncd.F	\
	enkf_ncdio.F	mod_enkf_obs.F	mod_enkfassim.F mod_assim.F	\
	mod_nesting.F	mod_visit.F     mod_plbc.F mod_dam.F		\
	mod_station_timeseries.F	mod_sparse_timeseries.F 	\
	mod_boundschk.F mod_esmf_nesting.F	\
	mod_cstms_vars.F	mod_flocmod.F	mod_sed_cstms.F		\
	mod_fluid_mud.F mod_tvd.F mod_mld_rho.F mod_heatflux_sediment.F \
	mod_vector_projection.F
# lwang added mod_mld_rho.F
MODS_W = mod_main_wave.F	swmod1.F	swmod2.F	swmod3.F\
	mod_action_im.F mod_action_ex.F mod_wavesetup.F 		\
	mod_wave_current_interaction.F  mod_bbl.F

MAIN  =	fvcom.F		genmap.F	tge.F		longshore_flow.F\
	cell_area.F 	open_all.F	load_grid.F	allocate_all.F	\
	setup_domain.F	genmap_obc.F	genmap_lsf.F	print_vals.F    \
	coords_n_const.F	shape_coef_gcn.F shape_coef_gcy.F	\
	depth_grad.F	grid_metrics.F	cntrl_prmtrs.F	init_sed.F	\
	internal_step.F	bcond_gcn.F	bcond_gcy.F	adjust2d3d.F	\
	brough.F	advection_edge_gcn.F	advection_edge_gcy.F	\
	ghostuv.F	advave_edge_gcn.F	advave_edge_gcy.F	\
	phy_baropg.F	baropg.F	external_step.F	extel_edge.F	\
	extuv_edge.F	depth_check.F	vertvl_edge.F	adv_uv_edge_gcn.F\
	adv_uv_edge_gcy.F  vdif_uv.F    extelpf_edge.F	wreal.F		\
	viscofh.F	adv_q.F		fct_q2.F	fct_q2l.F	\
	vdif_q.F	adv_t.F		adv_s.F		fct_t.F		\
	vdif_ts.F	fct_s.F		bcond_ts.F	\
	adjust_ts.F	conv_over.F	visitsim.F	startup_type.F	\
	edge_len.F	adcor.F		icing.F		rho_pmean.F	\
	calc_vort.F	namelist.F	nh_set_nesting.F		\
	coare26z.F	coare40vn.F	adv_uv_edge_gcn_rk.F		\
	adv_uv_edge_gcy_rk.F		adv_t_rk.F	adv_s_rk.F	\
	adv_q_rk.F

MAIN_W = ocpmix.F	ocpcre.F	ocpids.F	\
	swanpre1.F	swanpre2.F	swanser.F	swanmain.F	\
	swancom1.F	swancom2.F	swancom3.F	swancom4.F	\
	swancom5.F	w3part.F

COBJS = $(if $(FLAG_411), func_pointer.c partition.c, func_pointer.c)

BIOGEN	= mod_bio_3D.F	mod_onedtide.F

postFVCOM = postFVCOM_ncnest_reorder.F

SRCS = $(MODS) $(MODS_W) $(MAIN) $(BIOGEN) $(MAIN_W) 

# object list is composed of free (*.F) and fixed (*.F90) source
FVCOM_OBJS = $(COBJS:.c=.o) $(SRCS:.F=.o) $(FIXED:.F90=.o) $(SRCS:.f90=.o)

# Explicit module dependencies based on MODS order and actual USE statements
# These ensure modules are compiled in the correct order

# mpi_constants.f90 - no dependencies (base module)

# mod_prec.F - no dependencies (base module)

# sinter.F - depends on mod_prec
sinter.o: mod_prec.o

# mod_types.F - depends on mod_prec
mod_types.o: mod_prec.o

# mod_time.F - depends on mod_prec
mod_time.o: mod_prec.o

# mod_main.F - depends on mod_types, mod_time, mod_prec
mod_main.o: mod_types.o mod_time.o mod_prec.o

# mod_spherical.F - depends on mod_main, mod_prec
mod_spherical.o: mod_main.o mod_prec.o

# mod_utils.F - depends on mod_main, mod_spherical, mod_prec
mod_utils.o: mod_main.o mod_spherical.o mod_prec.o

# mod_clock.F - depends on mod_utils, mod_main, mod_time, mod_prec
mod_clock.o: mod_utils.o mod_main.o mod_time.o mod_prec.o

# eqs_of_state.F - depends on mod_utils, mod_main, mod_prec
eqs_of_state.o: mod_utils.o mod_main.o mod_prec.o

# mod_interp.F - depends on mod_main, mod_spherical, mod_clock, mod_utils, mod_prec
mod_interp.o: mod_main.o mod_spherical.o mod_clock.o mod_utils.o mod_prec.o

# mod_par.F - depends on mod_main, mod_time, mod_utils, mod_types, mod_prec
mod_par.o: mod_main.o mod_time.o mod_utils.o mod_types.o mod_prec.o

# mod_par_special.F - depends on mod_main, mod_time, mod_utils, mod_types, mod_prec
mod_par_special.o: mod_main.o mod_time.o mod_utils.o mod_types.o mod_prec.o

# mod_ncll.F - depends on mod_main, mod_interp, mod_utils, mod_time, mod_prec
mod_ncll.o: mod_main.o mod_interp.o mod_utils.o mod_time.o mod_prec.o

# mod_nctools.F - depends on mod_main, mod_par_special, mod_ncll, mod_prec
mod_nctools.o: mod_main.o mod_par_special.o mod_ncll.o mod_prec.o

# mod_wd.F - depends on mod_par_special, mod_main, mod_prec
mod_wd.o: mod_par_special.o mod_main.o mod_prec.o

# mod_sng.F - depends on mod_utils, mod_prec
mod_sng.o: mod_utils.o mod_prec.o

# mod_heatflux.F - depends on mod_main, mod_utils, mod_prec
mod_heatflux.o: mod_main.o mod_utils.o mod_prec.o

# mod_solar.F - depends on mod_main, mod_utils, mod_prec
mod_solar.o: mod_main.o mod_utils.o mod_prec.o

# mod_bulk.F - depends on mod_main_wave, mod_par_special, mod_main, mod_prec
mod_bulk.o: mod_main_wave.o mod_par_special.o mod_main.o mod_prec.o

# mod_input.F - depends on mod_main, mod_solar, mod_heatflux, mod_sng, mod_utils, mod_nctools, mod_prec, mpi_constants
mod_input.o: mod_main.o mod_solar.o mod_heatflux.o mod_sng.o mod_utils.o mod_nctools.o mod_prec.o mpi_constants.o

# mod_force.F - depends on swmod1, mod_solar, mod_heatflux, mod_input, mod_par_special, mod_spherical, mod_utils, mod_ncll, mod_nctools, mod_time, mod_main, mod_interp, mod_prec
mod_force.o: swmod1.o mod_solar.o mod_heatflux.o mod_input.o mod_par_special.o mod_spherical.o mod_utils.o mod_ncll.o mod_nctools.o mod_time.o mod_main.o mod_interp.o mod_prec.o

# mod_obcs.F - depends on mod_spherical, mod_main, mod_par_special, mod_main, mod_force, mod_time, mod_utils, mod_prec
mod_obcs.o: mod_spherical.o mod_main.o mod_par_special.o mod_main.o mod_force.o mod_time.o mod_utils.o mod_prec.o

# mod_petsc.F - depends on mod_main, mod_par_special, mod_prec
mod_petsc.o: mod_main.o mod_par_special.o mod_prec.o

# cntrl_prmtrs.F - depends on mod_heatflux, mod_nesting, mod_ncdio, mod_time, mod_assim, mod_utils, mod_main, mod_prec, mpi_constants
cntrl_prmtrs.o: mod_heatflux.o mod_nesting.o mod_ncdio.o mod_time.o mod_assim.o mod_utils.o mod_main.o mod_prec.o mpi_constants.o

# mod_semi_implicit.F - depends on mod_par_special, mod_main, mod_prec
mod_semi_implicit.o: mod_par_special.o mod_main.o mod_prec.o

# mod_non_hydro.F - depends on mod_main, mod_prec
mod_non_hydro.o: mod_main.o mod_prec.o

# mod_set_time.F - depends on mod_main, mod_prec
mod_set_time.o: mod_main.o mod_prec.o

# ice_kinds_mod.F - depends on mod_prec
ice_kinds_mod.o: mod_prec.o

# ice_model_size.F - depends on mod_prec
ice_model_size.o: mod_prec.o

# ice_domain.F - depends on mod_prec
ice_domain.o: mod_prec.o

# ice_constants.F - depends on mod_prec
ice_constants.o: mod_prec.o

# ice_fileunits.F - depends on mod_prec
ice_fileunits.o: mod_prec.o

# ice_state.F - depends on mod_prec
ice_state.o: mod_prec.o

# ice_work.F - depends on mod_prec
ice_work.o: mod_prec.o

# ice_grid.F - depends on mod_prec
ice_grid.o: mod_prec.o

# ice_albedo.F - depends on mod_prec
ice_albedo.o: mod_prec.o

# ice_calendar.F - depends on mod_prec
ice_calendar.o: mod_prec.o

# ice_flux.F - depends on mod_prec
ice_flux.o: mod_prec.o

# ice_flux_in.F - depends on mod_prec
ice_flux_in.o: mod_prec.o

# ice_itd.F - depends on mod_prec
ice_itd.o: mod_prec.o

# ice_mechred.F - depends on mod_prec
ice_mechred.o: mod_prec.o

# ice_itd_linear.F - depends on mod_prec
ice_itd_linear.o: mod_prec.o

# ice_scaling.F - depends on mod_prec
ice_scaling.o: mod_prec.o

# ice_atmo.F - depends on mod_prec
ice_atmo.o: mod_prec.o

# ice_ocean.F - depends on mod_prec
ice_ocean.o: mod_prec.o

# ice_therm_vertical.F - depends on mod_prec
ice_therm_vertical.o: mod_prec.o

# ice_init.F - depends on mod_prec
ice_init.o: mod_prec.o

# ice_therm_itd.F - depends on mod_prec
ice_therm_itd.o: mod_prec.o

# mod_ice2d.F - depends on mod_prec
mod_ice2d.o: mod_prec.o

# mod_ice.F - depends on mod_prec
mod_ice.o: mod_prec.o

# mod_startup.F - depends on mod_main, mod_prec
mod_startup.o: mod_main.o mod_prec.o

# mod_wqm.F - depends on mod_main, mod_prec
mod_wqm.o: mod_main.o mod_prec.o

# mod_ncdio.F - depends on mod_main, mod_prec
mod_ncdio.o: mod_main.o mod_prec.o

# mod_setup.F - depends on mod_main, mod_prec
mod_setup.o: mod_main.o mod_prec.o

# mod_newinp.F - depends on mod_main, mod_prec
mod_newinp.o: mod_main.o mod_prec.o

# particle.F - depends on mod_main, mod_prec
particle.o: mod_main.o mod_prec.o

# linklist.F - depends on mod_main, mod_prec
linklist.o: mod_main.o mod_prec.o

# mod_lag.F - depends on mod_main, mod_prec
mod_lag.o: mod_main.o mod_prec.o

# mod_northpole.F - depends on mod_main, mod_prec
mod_northpole.o: mod_main.o mod_prec.o

# mod_pwp.F - depends on mod_main, mod_prec
mod_pwp.o: mod_main.o mod_prec.o

# mod_dye.F - depends on mod_main, mod_prec
mod_dye.o: mod_main.o mod_prec.o

# mod_optimal_interpolation.F - depends on mod_main, mod_prec
mod_optimal_interpolation.o: mod_main.o mod_prec.o

# mod_report.F - depends on mod_main, mod_prec
mod_report.o: mod_main.o mod_prec.o

# mod_probe.F - depends on mod_main, mod_prec
mod_probe.o: mod_main.o mod_prec.o

# mod_gotm.F - depends on mod_main, mod_prec
mod_gotm.o: mod_main.o mod_prec.o

# mod_balance_2d.F - depends on mod_main, mod_prec
mod_balance_2d.o: mod_main.o mod_prec.o

# mod_tridiag.F - depends on mod_main, mod_prec
mod_tridiag.o: mod_main.o mod_prec.o

# mod_scal.F - depends on mod_main, mod_prec
mod_scal.o: mod_main.o mod_prec.o

# mod_meanflow.F - depends on mod_main, mod_prec
mod_meanflow.o: mod_main.o mod_prec.o

# mod_obcs2.F - depends on mod_main, mod_prec
mod_obcs2.o: mod_main.o mod_prec.o

# mod_obcs3.F - depends on mod_main, mod_prec
mod_obcs3.o: mod_main.o mod_prec.o

# mod_sed.F - depends on mod_main, mod_prec
mod_sed.o: mod_main.o mod_prec.o

# mod_enkf.F - depends on mod_main, mod_prec
mod_enkf.o: mod_main.o mod_prec.o

# mod_etkf.F - depends on mod_main, mod_prec
mod_etkf.o: mod_main.o mod_prec.o

# mod_rrk.F - depends on mod_main, mod_prec
mod_rrk.o: mod_main.o mod_prec.o

# mod_rrkf_obs.F - depends on mod_main, mod_prec
mod_rrkf_obs.o: mod_main.o mod_prec.o

# mod_rrkassim.F - depends on mod_main, mod_prec
mod_rrkassim.o: mod_main.o mod_prec.o

# mod_enkf_ncd.F - depends on mod_main, mod_prec
mod_enkf_ncd.o: mod_main.o mod_prec.o

# enkf_ncdio.F - depends on mod_main, mod_prec
enkf_ncdio.o: mod_main.o mod_prec.o

# mod_enkf_obs.F - depends on mod_main, mod_prec
mod_enkf_obs.o: mod_main.o mod_prec.o

# mod_enkfassim.F - depends on mod_main, mod_prec
mod_enkfassim.o: mod_main.o mod_prec.o

# mod_assim.F - depends on mod_main, mod_prec
mod_assim.o: mod_main.o mod_prec.o

# mod_nesting.F - depends on mod_main, mod_prec
mod_nesting.o: mod_main.o mod_prec.o

# mod_visit.F - depends on mod_main, mod_par, mod_ncdio, mod_utils, mod_prec
mod_visit.o: mod_main.o mod_par.o mod_ncdio.o mod_utils.o mod_prec.o

# mod_plbc.F - depends on mod_main, mod_prec
mod_plbc.o: mod_main.o mod_prec.o

# mod_dam.F - depends on mod_main, mod_prec
mod_dam.o: mod_main.o mod_prec.o

# mod_station_timeseries.F - depends on mod_main, mod_prec
mod_station_timeseries.o: mod_main.o mod_prec.o

# mod_sparse_timeseries.F - depends on mod_main, mod_prec
mod_sparse_timeseries.o: mod_main.o mod_prec.o

# mod_boundschk.F - depends on mod_main, mod_prec
mod_boundschk.o: mod_main.o mod_prec.o

# mod_esmf_nesting.F - depends on mod_main, mod_prec
mod_esmf_nesting.o: mod_main.o mod_prec.o

# mod_cstms_vars.F - depends on mod_main, mod_prec
mod_cstms_vars.o: mod_main.o mod_prec.o

# mod_flocmod.F - depends on mod_main, mod_prec
mod_flocmod.o: mod_main.o mod_prec.o

# mod_sed_cstms.F - depends on mod_main, mod_prec
mod_sed_cstms.o: mod_main.o mod_prec.o

# mod_fluid_mud.F - depends on mod_main, mod_prec
mod_fluid_mud.o: mod_main.o mod_prec.o

# mod_tvd.F - depends on mod_main, mod_prec
mod_tvd.o: mod_main.o mod_prec.o

# mod_mld_rho.F - depends on mod_main, mod_prec
mod_mld_rho.o: mod_main.o mod_prec.o

# mod_heatflux_sediment.F - depends on mod_main, mod_prec
mod_heatflux_sediment.o: mod_main.o mod_prec.o

# mod_vector_projection.F - depends on mod_main, mod_prec
mod_vector_projection.o: mod_main.o mod_prec.o 


#--------------------------------------------------------------------------
#  Linking Directives
#--------------------------------------------------------------------------

$(EXEC):	$(FVCOM_OBJS)
		$(FC) $(FFLAGS) $(VISOPT) $(LDFLAGS) $(FVCOM_OBJS) $(LIBS) -o $(EXEC)
		make libfvcom


libfvcom: $(FVCOM_OBJS)
	$(AR) $@.a $(FVCOM_OBJS)
	$(RANLIB) $@.a

#	$(FC) -shared -o $@.so $(FVCOM_OBJS)
#       PROBLEMS WITH SHARED OPTION, REQUIRES -fPIC in All sub libraries as well when building on x86_64

postFVCOM:	$(postFVCOM:.F=.o)
		$(FC) $< $(FFLAGS) $(LDFLAGS) $(LIBS) -o postFVCOM_ncnest_reorder

#--------------------------------------------------------------------------
#  Target to create dependecies.
#--------------------------------------------------------------------------

depend:
		makedepf90 -d mod_prec.o $(SRCS) >> makedepends

#--------------------------------------------------------------------------
#  Target to make all.
#--------------------------------------------------------------------------

SUBDIRS = testing input utilities


all:    
	make depend
	cd BIO_source && make
	make 
	for dir in $(SUBDIRS); do (cd $$dir && make) || exit 1; done

#--------------------------------------------------------------------------
#  Tar Up Code                           
#--------------------------------------------------------------------------

tarfile:
	tar cvf fvcom.tar *.F   makefile  makedepends RELEASE_NOTES 

#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

clean:
		/bin/rm -f *.o *.mod *.f90 fvcom libfvcom.a

clobber:	
		make clean
		/bin/rm -f fvcom libfvcom.so libfvcom.a postFVCOM_ncnest_reorder

allclean:
		make clean
		for dir in $(SUBDIRS); do (cd $$dir && make clean) || exit 1; done
	        cd BIO_source && make clean


allclobber:
		make clean
		/bin/rm -f fvcom libfvcom.so libfvcom.a
		for dir in $(SUBDIRS); do (cd $$dir && make clobber) || exit 1; done
	        cd BIO_source && make clobber

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
# Temporarily disable makedepends due to circular dependencies
# include ./makedepends
include ./make.inc
