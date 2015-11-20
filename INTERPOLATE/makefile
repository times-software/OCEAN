include ../Makefile.arch

ifdef OBF_DIR

SHIRLEY_DIR = $(OBF_DIR)/SHIRLEY
DEBUG =

include $(OBF_DIR)/make.sys

FFLAGS+= -g -traceback -warn all #-check all -init=snan -e03
FFLAGS+= $(OPTIONS)

#.f90.o:
#	$(MPIF90) $(F90FLAGS) -I$(OBF_DIR)/include -c $<

MODFLAGS= $(MOD_FLAG)/$(OBF_DIR)/iotk/src $(MOD_FLAG)/$(OBF_DIR)/Modules $(MOD_FLAG)/$(OBF_DIR)/PW $(MOD_FLAG).  $(MOD_FLAG)/$(OBF_DIR)/SHIRLEY

IFLAGS = -I$(OBF_DIR)/include

LIBOBJS = $(OBF_DIR)/flib/ptools.a $(OBF_DIR)/flib/flib.a $(OBF_DIR)/clib/clib.a $(OBF_DIR)/iotk/src/libiotk.a


OCEAN_BUILDER_EXE=ocean_builder$(DEBUG).x
SHIRLEY_HAM_EXE=shirley_ham_o$(DEBUG).x
OCEAN_QDIAGP_EXE=ocean_qdiagp$(DEBUG).x

all: OCEAN_timer.o $(OCEAN_BUILDER_EXE) $(SHIRLEY_HAM_EXE) $(OCEAN_QDIAGP_EXE)

install:
	cp $(OCEAN_BUILDER_EXE) $(INSTDIR)/ 
	cp $(SHIRLEY_HAM_EXE)   $(INSTDIR)/ 
	cp $(OCEAN_QDIAGP_EXE)  $(INSTDIR)/
	cp $(OBF_DIR)/bin/pw.x $(INSTDIR)/obf_pw.x
	cp $(OBF_DIR)/bin/pp.x $(INSTDIR)/obf_pp.x
	cp $(OBF_DIR)/bin/shirley_basis.x $(INSTDIR)/


OCEAN_O2L = \
../CNBSE/zcks/zobj/getprefs.o \
../CNBSE/zcks/zobj/snatch.o \
../CNBSE/zcks/zobj/nbseprjprep.o \
../CNBSE/zcks/zobj/nbseprjnproj.o \
../CNBSE/zcks/zobj/nbseftread.o \
../CNBSE/zcks/zobj/setlm.o  \
../CNBSE/zcks/zobj/getylm.o 

#$(OBF_DIR)/OCEAN/getprefs.o \
$(OBF_DIR)/OCEAN/snatch.o \
$(OBF_DIR)/OCEAN/nbseprjprep.o \
$(OBF_DIR)/OCEAN/nbseprjnproj.o \
$(OBF_DIR)/OCEAN/nbseftread.o \
$(OBF_DIR)/OCEAN/setlm.o \
$(OBF_DIR)/OCEAN/getylm.o

NEWDIAGOBJ= \
$(OBF_DIR)/SHIRLEY/bspline90_22.o \
$(OBF_DIR)/SHIRLEY/corerepair.o \
$(OBF_DIR)/SHIRLEY/diag.o \
$(OBF_DIR)/SHIRLEY/freeunit.o \
$(OBF_DIR)/SHIRLEY/hamq_pool.o \
$(OBF_DIR)/SHIRLEY/hamq_shirley.o \
$(OBF_DIR)/SHIRLEY/kpt.o \
$(OBF_DIR)/SHIRLEY/mpio.o \
$(OBF_DIR)/SHIRLEY/plot.o \
$(OBF_DIR)/SHIRLEY/pwmat_module.o \
$(OBF_DIR)/SHIRLEY/scalapack_module.o \
$(OBF_DIR)/SHIRLEY/shirley_constants.o \
$(OBF_DIR)/SHIRLEY/shirley_input.o \
$(OBF_DIR)/SHIRLEY/start_shirley.o 

OCEAN_QDIAGP = OCEAN_qdiagp.o $(NEWDIAGOBJ) fix_fermi.o OCEAN_timer.o


SHIRLEY_HAM_O = \
shirley_ham_o.o \
OCEAN_bofr_multi.o \
OCEAN_obf2localbyk.o \
OCEAN_obf2localbyk_subs.o \
../zbridge/zobj/chkfftreal.o \
../zbridge/zobj/chkfftrecp.o \
../zbridge/zobj/facpowfind.o \
../zbridge/zobj/optim.o \
../CNBSE/zcmsnmult/zobj/jlmfft.o \
OCEAN_mod_bofx.o \
OCEAN_timer.o \
OCEAN_mod_obf2loc.o \
OCEAN_bofx.o 
#OCEAN_par_gentoreal.o \
OCEAN_par_gen_setup.o \

#$(OBF_DIR)/OCEAN/jlmfft.o \
#$(OBF_DIR)/OCEAN/OCEAN_2bphase.o\
#$(OBF_DIR)/OCEAN/OCEAN_bigbofr.o \
#$(OBF_DIR)/OCEAN/OCEAN_bradial.o \
#$(OBF_DIR)/OCEAN/OCEAN_bofj.o \
$(OBF_DIR)/OCEAN/chkfftreal.o \
$(OBF_DIR)/OCEAN/chkfftrecp.o \
$(OBF_DIR)/OCEAN/facpowfind.o \
$(OBF_DIR)/OCEAN/optim.o \
$(OBF_DIR)/OCEAN/OCEAN_bofr.o \

SHIRLEY_HAM = \
$(OBF_DIR)/SHIRLEY/append_missing_slash.o \
$(OBF_DIR)/SHIRLEY/atomic_wfc_shirley.o \
$(OBF_DIR)/SHIRLEY/B0.o \
$(OBF_DIR)/SHIRLEY/bspline90_22.o \
$(OBF_DIR)/SHIRLEY/elph_shirley.o \
$(OBF_DIR)/SHIRLEY/freeunit.o \
$(OBF_DIR)/SHIRLEY/get_dnl.o \
$(OBF_DIR)/SHIRLEY/get_proj_ham.o \
$(OBF_DIR)/SHIRLEY/gk_l2gmap.o \
$(OBF_DIR)/SHIRLEY/hamq.o \
$(OBF_DIR)/SHIRLEY/hamq_pool.o \
$(OBF_DIR)/SHIRLEY/hamq_shirley.o \
$(OBF_DIR)/SHIRLEY/init_us_2_shirley.o \
$(OBF_DIR)/SHIRLEY/interaction.o \
$(OBF_DIR)/SHIRLEY/interaction_test.o \
$(OBF_DIR)/SHIRLEY/kpt.o \
$(OBF_DIR)/SHIRLEY/mpio.o \
$(OBF_DIR)/SHIRLEY/orthoatwfc_shirley.o \
$(OBF_DIR)/SHIRLEY/qinteraction.o \
$(OBF_DIR)/SHIRLEY/plot.o \
$(OBF_DIR)/SHIRLEY/plot_basis.o \
$(OBF_DIR)/SHIRLEY/pole_strength.o \
$(OBF_DIR)/SHIRLEY/proj_shirley.o \
$(OBF_DIR)/SHIRLEY/print_clock_shirley.o \
$(OBF_DIR)/SHIRLEY/pwmat.o \
$(OBF_DIR)/SHIRLEY/read_file_shirley.o \
$(OBF_DIR)/SHIRLEY/reallocate_nlpot.o \
$(OBF_DIR)/SHIRLEY/rgrid.o \
$(OBF_DIR)/SHIRLEY/scalapack_module.o \
$(OBF_DIR)/SHIRLEY/shirley_ham_input.o \
$(OBF_DIR)/SHIRLEY/splines.o \
$(OBF_DIR)/SHIRLEY/start_shirley.o \
$(OBF_DIR)/SHIRLEY/stop_shirley.o \
$(OBF_DIR)/SHIRLEY/triples.o \
$(OBF_DIR)/SHIRLEY/trnvect_wrapper.o \
$(OBF_DIR)/SHIRLEY/update_pseudo.o

PWOBJS = \
$(OBF_DIR)/PW/a2fmod.o \
$(OBF_DIR)/PW/add_bfield.o \
$(OBF_DIR)/PW/add_efield.o \
$(OBF_DIR)/PW/add_paw_to_deeq.o \
$(OBF_DIR)/PW/add_vuspsi.o \
$(OBF_DIR)/PW/addusdens.o \
$(OBF_DIR)/PW/allocate_bp_efield.o \
$(OBF_DIR)/PW/allocate_fft.o \
$(OBF_DIR)/PW/allocate_locpot.o \
$(OBF_DIR)/PW/allocate_nlpot.o \
$(OBF_DIR)/PW/allocate_wfc.o \
$(OBF_DIR)/PW/atomic_rho.o \
$(OBF_DIR)/PW/atomic_wfc.o \
$(OBF_DIR)/PW/average_pp.o \
$(OBF_DIR)/PW/becmod.o \
$(OBF_DIR)/PW/bp_c_phase.o \
$(OBF_DIR)/PW/bp_calc_btq.o \
$(OBF_DIR)/PW/bp_qvan3.o \
$(OBF_DIR)/PW/bp_strings.o \
$(OBF_DIR)/PW/buffers.o \
$(OBF_DIR)/PW/ccgdiagg.o \
$(OBF_DIR)/PW/cdiagh.o \
$(OBF_DIR)/PW/cdiaghg.o \
$(OBF_DIR)/PW/cegterg.o \
$(OBF_DIR)/PW/clean_pw.o \
$(OBF_DIR)/PW/close_files.o \
$(OBF_DIR)/PW/compute_dip.o \
$(OBF_DIR)/PW/compute_rho.o \
$(OBF_DIR)/PW/compute_ux.o \
$(OBF_DIR)/PW/coset.o \
$(OBF_DIR)/PW/d_matrix.o \
$(OBF_DIR)/PW/data_structure.o \
$(OBF_DIR)/PW/deriv_drhoc.o \
$(OBF_DIR)/PW/divide.o \
$(OBF_DIR)/PW/divide_class.o \
$(OBF_DIR)/PW/divide_class_so.o \
$(OBF_DIR)/PW/divide_et_impera.o \
$(OBF_DIR)/PW/dqvan2.o \
$(OBF_DIR)/PW/drhoc.o \
$(OBF_DIR)/PW/dvloc_of_g.o \
$(OBF_DIR)/PW/dynamics_module.o \
$(OBF_DIR)/PW/efermig.o \
$(OBF_DIR)/PW/efermit.o \
$(OBF_DIR)/PW/eqvect.o \
$(OBF_DIR)/PW/ewald.o \
$(OBF_DIR)/PW/ewald_dipole.o \
$(OBF_DIR)/PW/exx.o \
$(OBF_DIR)/PW/find_group.o \
$(OBF_DIR)/PW/g_psi.o \
$(OBF_DIR)/PW/g_psi_mod.o \
$(OBF_DIR)/PW/g2_kin.o \
$(OBF_DIR)/PW/gen_at_dj.o \
$(OBF_DIR)/PW/gen_at_dy.o \
$(OBF_DIR)/PW/gen_us_dj.o \
$(OBF_DIR)/PW/gen_us_dy.o \
$(OBF_DIR)/PW/get_locals.o \
$(OBF_DIR)/PW/gk_sort.o \
$(OBF_DIR)/PW/gradcorr.o \
$(OBF_DIR)/PW/gweights.o \
$(OBF_DIR)/PW/h_1psi.o \
$(OBF_DIR)/PW/h_epsi_her_apply.o \
$(OBF_DIR)/PW/h_epsi_her_set.o \
$(OBF_DIR)/PW/h_psi.o \
$(OBF_DIR)/PW/h_psi_meta.o \
$(OBF_DIR)/PW/hinit0.o \
$(OBF_DIR)/PW/hinit1.o \
$(OBF_DIR)/PW/init_at_1.o \
$(OBF_DIR)/PW/init_ns.o \
$(OBF_DIR)/PW/init_run.o \
$(OBF_DIR)/PW/init_us_1.o \
$(OBF_DIR)/PW/init_us_2.o \
$(OBF_DIR)/PW/init_vloc.o \
$(OBF_DIR)/PW/interpolate.o \
$(OBF_DIR)/PW/io_rho_xml.o \
$(OBF_DIR)/PW/irrek.o \
$(OBF_DIR)/PW/iweights.o \
$(OBF_DIR)/PW/kernel_table.o \
$(OBF_DIR)/PW/kpoint_grid.o \
$(OBF_DIR)/PW/lchk_tauxk.o \
$(OBF_DIR)/PW/make_pointlists.o \
$(OBF_DIR)/PW/martyna_tuckerman.o \
$(OBF_DIR)/PW/memory_report.o \
$(OBF_DIR)/PW/multable.o \
$(OBF_DIR)/PW/n_plane_waves.o \
$(OBF_DIR)/PW/new_ns.o \
$(OBF_DIR)/PW/new_occ.o \
$(OBF_DIR)/PW/newd.o \
$(OBF_DIR)/PW/noncol.o \
$(OBF_DIR)/PW/ns_adj.o \
$(OBF_DIR)/PW/offset_atom_wfc.o \
$(OBF_DIR)/PW/openfil.o \
$(OBF_DIR)/PW/ortho_wfc.o \
$(OBF_DIR)/PW/orthoatwfc.o \
$(OBF_DIR)/PW/output_tau.o \
$(OBF_DIR)/PW/para.o \
$(OBF_DIR)/PW/paw_init.o \
$(OBF_DIR)/PW/paw_onecenter.o \
$(OBF_DIR)/PW/paw_symmetry.o \
$(OBF_DIR)/PW/potinit.o \
$(OBF_DIR)/PW/print_clock_pw.o \
$(OBF_DIR)/PW/punch.o \
$(OBF_DIR)/PW/pwcom.o \
$(OBF_DIR)/PW/pw_restart.o \
$(OBF_DIR)/PW/qvan2.o \
$(OBF_DIR)/PW/rcgdiagg.o \
$(OBF_DIR)/PW/rdiaghg.o \
$(OBF_DIR)/PW/read_pseudo.o \
$(OBF_DIR)/PW/realus.o \
$(OBF_DIR)/PW/report_mag.o \
$(OBF_DIR)/PW/rho2zeta.o \
$(OBF_DIR)/PW/rotate_wfc.o \
$(OBF_DIR)/PW/rotate_wfc_gamma.o \
$(OBF_DIR)/PW/rotate_wfc_k.o \
$(OBF_DIR)/PW/ruotaijk.o \
$(OBF_DIR)/PW/s_1psi.o \
$(OBF_DIR)/PW/s_psi.o \
$(OBF_DIR)/PW/save_in_cbands.o \
$(OBF_DIR)/PW/save_in_electrons.o \
$(OBF_DIR)/PW/save_in_ions.o \
$(OBF_DIR)/PW/scale_h.o \
$(OBF_DIR)/PW/scf_mod.o \
$(OBF_DIR)/PW/seqopn.o \
$(OBF_DIR)/PW/set_hubbard_l.o \
$(OBF_DIR)/PW/set_kplusq.o \
$(OBF_DIR)/PW/set_kup_and_kdw.o \
$(OBF_DIR)/PW/set_rhoc.o \
$(OBF_DIR)/PW/set_vrs.o \
$(OBF_DIR)/PW/setlocal.o \
$(OBF_DIR)/PW/setqf.o \
$(OBF_DIR)/PW/setup.o \
$(OBF_DIR)/PW/spinor.o \
$(OBF_DIR)/PW/sph_ind.o \
$(OBF_DIR)/PW/struct_fact.o \
$(OBF_DIR)/PW/sum_band.o \
$(OBF_DIR)/PW/sumkg.o \
$(OBF_DIR)/PW/sumkt.o \
$(OBF_DIR)/PW/summary.o \
$(OBF_DIR)/PW/symm_base.o \
$(OBF_DIR)/PW/symme.o \
$(OBF_DIR)/PW/tabd.o \
$(OBF_DIR)/PW/transform_becsum_nc.o \
$(OBF_DIR)/PW/transform_becsum_so.o \
$(OBF_DIR)/PW/trnvecc.o \
$(OBF_DIR)/PW/tweights.o \
$(OBF_DIR)/PW/update_pot.o \
$(OBF_DIR)/PW/usnldiag.o \
$(OBF_DIR)/PW/v_of_rho.o \
$(OBF_DIR)/PW/vhpsi.o \
$(OBF_DIR)/PW/vloc_of_g.o \
$(OBF_DIR)/PW/vloc_psi.o \
$(OBF_DIR)/PW/wannier_check.o \
$(OBF_DIR)/PW/wannier_clean.o \
$(OBF_DIR)/PW/wannier_init.o \
$(OBF_DIR)/PW/weights.o \
$(OBF_DIR)/PW/wfcinit.o \
$(OBF_DIR)/PW/write_ns.o \
$(OBF_DIR)/PW/wsweight.o \
$(OBF_DIR)/PW/xc_vdW_DF.o


NEWDIAGOBJ= \
$(SHIRLEY_DIR)/bspline90_22.o \
$(SHIRLEY_DIR)/corerepair.o \
$(SHIRLEY_DIR)/diag.o \
$(SHIRLEY_DIR)/freeunit.o \
$(SHIRLEY_DIR)/hamq_pool.o \
$(SHIRLEY_DIR)/hamq_shirley.o \
$(SHIRLEY_DIR)/kpt.o \
$(SHIRLEY_DIR)/mpio.o \
$(SHIRLEY_DIR)/plot.o \
$(SHIRLEY_DIR)/pwmat_module.o \
$(SHIRLEY_DIR)/scalapack_module.o \
$(SHIRLEY_DIR)/shirley_constants.o \
$(SHIRLEY_DIR)/shirley_input.o \
$(SHIRLEY_DIR)/start_shirley.o

PPOBJS = $(OBF_DIR)/PP/xsf.o $(OBF_DIR)/PP/cube.o

PW_QDIAGP = \
$(OBF_DIR)/flib/date_and_tim.o \
$(OBF_DIR)/flib/flush_unit.o \
$(OBF_DIR)/flib/inpfile.o

MODULES = \
$(OBF_DIR)/Modules/atom.o \
$(OBF_DIR)/Modules/basic_algebra_routines.o \
$(OBF_DIR)/Modules/cell_base.o \
$(OBF_DIR)/Modules/check_stop.o \
$(OBF_DIR)/Modules/clocks.o \
$(OBF_DIR)/Modules/constants.o \
$(OBF_DIR)/Modules/constraints_module.o \
$(OBF_DIR)/Modules/control_flags.o \
$(OBF_DIR)/Modules/coulomb_vcut.o \
$(OBF_DIR)/Modules/descriptors.o \
$(OBF_DIR)/Modules/dspev_drv.o \
$(OBF_DIR)/Modules/electrons_base.o \
$(OBF_DIR)/Modules/environment.o \
$(OBF_DIR)/Modules/error_handler.o \
$(OBF_DIR)/Modules/fft_base.o \
$(OBF_DIR)/Modules/fft_interfaces.o \
$(OBF_DIR)/Modules/fft_parallel.o \
$(OBF_DIR)/Modules/fft_scalar.o \
$(OBF_DIR)/Modules/fft_types.o \
$(OBF_DIR)/Modules/funct.o \
$(OBF_DIR)/Modules/griddim.o \
$(OBF_DIR)/Modules/image_io_routines.o \
$(OBF_DIR)/Modules/input_parameters.o \
$(OBF_DIR)/Modules/io_files.o \
$(OBF_DIR)/Modules/io_global.o \
$(OBF_DIR)/Modules/ions_base.o \
$(OBF_DIR)/Modules/kind.o \
$(OBF_DIR)/Modules/mm_dispersion.o \
$(OBF_DIR)/Modules/mp.o \
$(OBF_DIR)/Modules/mp_base.o \
$(OBF_DIR)/Modules/mp_global.o \
$(OBF_DIR)/Modules/mp_wave.o \
$(OBF_DIR)/Modules/parallel_include.o \
$(OBF_DIR)/Modules/parameters.o \
$(OBF_DIR)/Modules/parser.o \
$(OBF_DIR)/Modules/paw_variables.o \
$(OBF_DIR)/Modules/printout_base.o \
$(OBF_DIR)/Modules/pseudo_types.o \
$(OBF_DIR)/Modules/ptoolkit.o \
$(OBF_DIR)/Modules/radial_grids.o \
$(OBF_DIR)/Modules/random_numbers.o \
$(OBF_DIR)/Modules/read_ncpp.o \
$(OBF_DIR)/Modules/read_upf_v1.o \
$(OBF_DIR)/Modules/read_upf_v2.o \
$(OBF_DIR)/Modules/read_uspp.o \
$(OBF_DIR)/Modules/recvec.o \
$(OBF_DIR)/Modules/recvec_subs.o \
$(OBF_DIR)/Modules/splinelib.o \
$(OBF_DIR)/Modules/stick_base.o \
$(OBF_DIR)/Modules/stick_set.o \
$(OBF_DIR)/Modules/timestep.o \
$(OBF_DIR)/Modules/uspp.o \
$(OBF_DIR)/Modules/upf.o \
$(OBF_DIR)/Modules/upf_to_internal.o \
$(OBF_DIR)/Modules/version.o \
$(OBF_DIR)/Modules/wannier.o \
$(OBF_DIR)/Modules/wannier_new.o \
$(OBF_DIR)/Modules/wavefunctions.o \
$(OBF_DIR)/Modules/wrappers.o \
$(OBF_DIR)/Modules/write_upf_v2.o \
$(OBF_DIR)/Modules/ws_base.o \
$(OBF_DIR)/Modules/xml_io_base.o \
$(OBF_DIR)/Modules/zhpev_drv.o

ifdef NIST
MODULES += $(OBF_DIR)/Modules/netcdf_io.o 
endif


OCEAN_BUILDER = OCEAN_builder.o mkvipt.o intval.o getomega.o OCEAN_timer.o realspace_paw.o getprefs.o getylm.o OCEAN_build_chi.o $(NEWDIAGOBJ) 

$(OCEAN_BUILDER_EXE) :  $(OCEAN_BUILDER) $(PPOBJS) $(PW_QDIAGP) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
	$(OCEAN_BUILDER) $(PPOBJS) $(PWOBJS) $(OBF_DIR)/SHIRLEY/read_file_shirley.o \
	$(PW_QDIAGP) $(MODULES) $(LIBOBJS) $(LIBS) 

$(SHIRLEY_HAM_EXE) : $(SHIRLEY_HAM) $(SHIRLEY_HAM_O) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(OCEAN_O2L)
	$(MPIF90) $(LDFLAGS) -o $@ \
	$(SHIRLEY_HAM) $(SHIRLEY_HAM_O) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(OCEAN_O2L) $(LIBS)

$(OCEAN_QDIAGP_EXE) : $(OCEAN_QDIAGP) $(PPOBJS) $(PW_QDIAGP) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
	$(OCEAN_QDIAGP) $(PPOBJS) $(PWOBJS) $(OBF_DIR)/SHIRLEY/read_file_shirley.o $(PW_QDIAGP) $(MODULES) $(LIBOBJS) $(LIBS)

else
all :
	@echo "Your build does not have OBF support."

install:
	@echo "Your build does not have OBF support."
  
endif

clean : 
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i work.pc *__genmod.f90



