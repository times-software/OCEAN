ESPRESSO_DIR = /home/jtv1/shirley_QE4.3-gnu-atlas
SHIRLEY_DIR = $(ESPRESSO_DIR)/SHIRLEY
DEBUG =

include $(ESPRESSO_DIR)/make.sys

FFLAGS+= -Wall

#.f90.o:
#	$(MPIF90) $(F90FLAGS) -I$(ESPRESSO_DIR)/include -c $<

MODFLAGS= $(MOD_FLAG)/$(ESPRESSO_DIR)/iotk/src $(MOD_FLAG)/$(ESPRESSO_DIR)/Modules $(MOD_FLAG)/$(ESPRESSO_DIR)/PW $(MOD_FLAG).  $(MOD_FLAG)/$(ESPRESSO_DIR)/SHIRLEY

IFLAGS = -I$(ESPRESSO_DIR)/include

LIBOBJS = $(ESPRESSO_DIR)/flib/ptools.a $(ESPRESSO_DIR)/flib/flib.a $(ESPRESSO_DIR)/clib/clib.a $(ESPRESSO_DIR)/iotk/src/libiotk.a


OCEAN_BUILDER_EXE=ocean_builder$(DEBUG).x
SHIRLEY_HAM_EXE=shirley_ham_o$(DEBUG).x

all: $(OCEAN_BUILDER_EXE) $(SHIRLEY_HAM_EXE)

OCEAN_O2L = \
$(ESPRESSO_DIR)/OCEAN/getprefs.o \
$(ESPRESSO_DIR)/OCEAN/snatch.o \
$(ESPRESSO_DIR)/OCEAN/nbseprjprep.o \
$(ESPRESSO_DIR)/OCEAN/nbseprjnproj.o \
$(ESPRESSO_DIR)/OCEAN/nbseftread.o \
$(ESPRESSO_DIR)/OCEAN/setlm.o \
$(ESPRESSO_DIR)/OCEAN/getylm.o


SHIRLEY_HAM_O = \
$(ESPRESSO_DIR)/OCEAN/shirley_ham_o.o \
$(ESPRESSO_DIR)/OCEAN/OCEAN_bofx.o \
$(ESPRESSO_DIR)/OCEAN/OCEAN_bofr.o \
$(ESPRESSO_DIR)/OCEAN/OCEAN_bofr_multi.o \
$(ESPRESSO_DIR)/OCEAN/OCEAN_bigbofr.o \
$(ESPRESSO_DIR)/OCEAN/OCEAN_bradial.o \
$(ESPRESSO_DIR)/OCEAN/OCEAN_bofj.o \
OCEAN_obf2localbyk.o \
$(ESPRESSO_DIR)/OCEAN/chkfftreal.o \
$(ESPRESSO_DIR)/OCEAN/chkfftrecp.o \
$(ESPRESSO_DIR)/OCEAN/facpowfind.o \
$(ESPRESSO_DIR)/OCEAN/optim.o \
$(ESPRESSO_DIR)/OCEAN/jlmfft.o \
$(ESPRESSO_DIR)/OCEAN/OCEAN_2bphase.o\
OCEAN_timer.o

SHIRLEY_HAM = \
$(ESPRESSO_DIR)/SHIRLEY/append_missing_slash.o \
$(ESPRESSO_DIR)/SHIRLEY/atomic_wfc_shirley.o \
$(ESPRESSO_DIR)/SHIRLEY/B0.o \
$(ESPRESSO_DIR)/SHIRLEY/bspline90_22.o \
$(ESPRESSO_DIR)/SHIRLEY/elph_shirley.o \
$(ESPRESSO_DIR)/SHIRLEY/freeunit.o \
$(ESPRESSO_DIR)/SHIRLEY/get_dnl.o \
$(ESPRESSO_DIR)/SHIRLEY/get_proj_ham.o \
$(ESPRESSO_DIR)/SHIRLEY/gk_l2gmap.o \
$(ESPRESSO_DIR)/SHIRLEY/hamq.o \
$(ESPRESSO_DIR)/SHIRLEY/hamq_pool.o \
$(ESPRESSO_DIR)/SHIRLEY/hamq_shirley.o \
$(ESPRESSO_DIR)/SHIRLEY/init_us_2_shirley.o \
$(ESPRESSO_DIR)/SHIRLEY/interaction.o \
$(ESPRESSO_DIR)/SHIRLEY/interaction_test.o \
$(ESPRESSO_DIR)/SHIRLEY/kpt.o \
$(ESPRESSO_DIR)/SHIRLEY/mpio.o \
$(ESPRESSO_DIR)/SHIRLEY/orthoatwfc_shirley.o \
$(ESPRESSO_DIR)/SHIRLEY/qinteraction.o \
$(ESPRESSO_DIR)/SHIRLEY/plot.o \
$(ESPRESSO_DIR)/SHIRLEY/plot_basis.o \
$(ESPRESSO_DIR)/SHIRLEY/pole_strength.o \
$(ESPRESSO_DIR)/SHIRLEY/proj_shirley.o \
$(ESPRESSO_DIR)/SHIRLEY/print_clock_shirley.o \
$(ESPRESSO_DIR)/SHIRLEY/pwmat.o \
$(ESPRESSO_DIR)/SHIRLEY/read_file_shirley.o \
$(ESPRESSO_DIR)/SHIRLEY/reallocate_nlpot.o \
$(ESPRESSO_DIR)/SHIRLEY/rgrid.o \
$(ESPRESSO_DIR)/SHIRLEY/scalapack_module.o \
$(ESPRESSO_DIR)/SHIRLEY/shirley_ham_input.o \
$(ESPRESSO_DIR)/SHIRLEY/splines.o \
$(ESPRESSO_DIR)/SHIRLEY/start_shirley.o \
$(ESPRESSO_DIR)/SHIRLEY/stop_shirley.o \
$(ESPRESSO_DIR)/SHIRLEY/triples.o \
$(ESPRESSO_DIR)/SHIRLEY/trnvect_wrapper.o \
$(ESPRESSO_DIR)/SHIRLEY/update_pseudo.o

PWOBJS = \
$(ESPRESSO_DIR)/PW/a2fmod.o \
$(ESPRESSO_DIR)/PW/add_bfield.o \
$(ESPRESSO_DIR)/PW/add_efield.o \
$(ESPRESSO_DIR)/PW/add_paw_to_deeq.o \
$(ESPRESSO_DIR)/PW/add_vuspsi.o \
$(ESPRESSO_DIR)/PW/addusdens.o \
$(ESPRESSO_DIR)/PW/allocate_bp_efield.o \
$(ESPRESSO_DIR)/PW/allocate_fft.o \
$(ESPRESSO_DIR)/PW/allocate_locpot.o \
$(ESPRESSO_DIR)/PW/allocate_nlpot.o \
$(ESPRESSO_DIR)/PW/allocate_wfc.o \
$(ESPRESSO_DIR)/PW/atomic_rho.o \
$(ESPRESSO_DIR)/PW/atomic_wfc.o \
$(ESPRESSO_DIR)/PW/average_pp.o \
$(ESPRESSO_DIR)/PW/becmod.o \
$(ESPRESSO_DIR)/PW/bp_c_phase.o \
$(ESPRESSO_DIR)/PW/bp_calc_btq.o \
$(ESPRESSO_DIR)/PW/bp_qvan3.o \
$(ESPRESSO_DIR)/PW/bp_strings.o \
$(ESPRESSO_DIR)/PW/buffers.o \
$(ESPRESSO_DIR)/PW/ccgdiagg.o \
$(ESPRESSO_DIR)/PW/cdiagh.o \
$(ESPRESSO_DIR)/PW/cdiaghg.o \
$(ESPRESSO_DIR)/PW/cegterg.o \
$(ESPRESSO_DIR)/PW/clean_pw.o \
$(ESPRESSO_DIR)/PW/close_files.o \
$(ESPRESSO_DIR)/PW/compute_dip.o \
$(ESPRESSO_DIR)/PW/compute_rho.o \
$(ESPRESSO_DIR)/PW/compute_ux.o \
$(ESPRESSO_DIR)/PW/coset.o \
$(ESPRESSO_DIR)/PW/d_matrix.o \
$(ESPRESSO_DIR)/PW/data_structure.o \
$(ESPRESSO_DIR)/PW/deriv_drhoc.o \
$(ESPRESSO_DIR)/PW/divide.o \
$(ESPRESSO_DIR)/PW/divide_class.o \
$(ESPRESSO_DIR)/PW/divide_class_so.o \
$(ESPRESSO_DIR)/PW/divide_et_impera.o \
$(ESPRESSO_DIR)/PW/dqvan2.o \
$(ESPRESSO_DIR)/PW/drhoc.o \
$(ESPRESSO_DIR)/PW/dvloc_of_g.o \
$(ESPRESSO_DIR)/PW/dynamics_module.o \
$(ESPRESSO_DIR)/PW/efermig.o \
$(ESPRESSO_DIR)/PW/efermit.o \
$(ESPRESSO_DIR)/PW/eqvect.o \
$(ESPRESSO_DIR)/PW/ewald.o \
$(ESPRESSO_DIR)/PW/ewald_dipole.o \
$(ESPRESSO_DIR)/PW/exx.o \
$(ESPRESSO_DIR)/PW/find_group.o \
$(ESPRESSO_DIR)/PW/g_psi.o \
$(ESPRESSO_DIR)/PW/g_psi_mod.o \
$(ESPRESSO_DIR)/PW/g2_kin.o \
$(ESPRESSO_DIR)/PW/gen_at_dj.o \
$(ESPRESSO_DIR)/PW/gen_at_dy.o \
$(ESPRESSO_DIR)/PW/gen_us_dj.o \
$(ESPRESSO_DIR)/PW/gen_us_dy.o \
$(ESPRESSO_DIR)/PW/get_locals.o \
$(ESPRESSO_DIR)/PW/gk_sort.o \
$(ESPRESSO_DIR)/PW/gradcorr.o \
$(ESPRESSO_DIR)/PW/gweights.o \
$(ESPRESSO_DIR)/PW/h_1psi.o \
$(ESPRESSO_DIR)/PW/h_epsi_her_apply.o \
$(ESPRESSO_DIR)/PW/h_epsi_her_set.o \
$(ESPRESSO_DIR)/PW/h_psi.o \
$(ESPRESSO_DIR)/PW/h_psi_meta.o \
$(ESPRESSO_DIR)/PW/hinit0.o \
$(ESPRESSO_DIR)/PW/hinit1.o \
$(ESPRESSO_DIR)/PW/init_at_1.o \
$(ESPRESSO_DIR)/PW/init_ns.o \
$(ESPRESSO_DIR)/PW/init_run.o \
$(ESPRESSO_DIR)/PW/init_us_1.o \
$(ESPRESSO_DIR)/PW/init_us_2.o \
$(ESPRESSO_DIR)/PW/init_vloc.o \
$(ESPRESSO_DIR)/PW/interpolate.o \
$(ESPRESSO_DIR)/PW/io_rho_xml.o \
$(ESPRESSO_DIR)/PW/irrek.o \
$(ESPRESSO_DIR)/PW/iweights.o \
$(ESPRESSO_DIR)/PW/kernel_table.o \
$(ESPRESSO_DIR)/PW/kpoint_grid.o \
$(ESPRESSO_DIR)/PW/lchk_tauxk.o \
$(ESPRESSO_DIR)/PW/make_pointlists.o \
$(ESPRESSO_DIR)/PW/martyna_tuckerman.o \
$(ESPRESSO_DIR)/PW/memory_report.o \
$(ESPRESSO_DIR)/PW/multable.o \
$(ESPRESSO_DIR)/PW/n_plane_waves.o \
$(ESPRESSO_DIR)/PW/new_ns.o \
$(ESPRESSO_DIR)/PW/new_occ.o \
$(ESPRESSO_DIR)/PW/newd.o \
$(ESPRESSO_DIR)/PW/noncol.o \
$(ESPRESSO_DIR)/PW/ns_adj.o \
$(ESPRESSO_DIR)/PW/offset_atom_wfc.o \
$(ESPRESSO_DIR)/PW/openfil.o \
$(ESPRESSO_DIR)/PW/ortho_wfc.o \
$(ESPRESSO_DIR)/PW/orthoatwfc.o \
$(ESPRESSO_DIR)/PW/output_tau.o \
$(ESPRESSO_DIR)/PW/para.o \
$(ESPRESSO_DIR)/PW/paw_init.o \
$(ESPRESSO_DIR)/PW/paw_onecenter.o \
$(ESPRESSO_DIR)/PW/paw_symmetry.o \
$(ESPRESSO_DIR)/PW/potinit.o \
$(ESPRESSO_DIR)/PW/print_clock_pw.o \
$(ESPRESSO_DIR)/PW/punch.o \
$(ESPRESSO_DIR)/PW/pwcom.o \
$(ESPRESSO_DIR)/PW/pw_restart.o \
$(ESPRESSO_DIR)/PW/qvan2.o \
$(ESPRESSO_DIR)/PW/rcgdiagg.o \
$(ESPRESSO_DIR)/PW/rdiaghg.o \
$(ESPRESSO_DIR)/PW/read_pseudo.o \
$(ESPRESSO_DIR)/PW/realus.o \
$(ESPRESSO_DIR)/PW/report_mag.o \
$(ESPRESSO_DIR)/PW/rho2zeta.o \
$(ESPRESSO_DIR)/PW/rotate_wfc.o \
$(ESPRESSO_DIR)/PW/rotate_wfc_gamma.o \
$(ESPRESSO_DIR)/PW/rotate_wfc_k.o \
$(ESPRESSO_DIR)/PW/ruotaijk.o \
$(ESPRESSO_DIR)/PW/s_1psi.o \
$(ESPRESSO_DIR)/PW/s_psi.o \
$(ESPRESSO_DIR)/PW/save_in_cbands.o \
$(ESPRESSO_DIR)/PW/save_in_electrons.o \
$(ESPRESSO_DIR)/PW/save_in_ions.o \
$(ESPRESSO_DIR)/PW/scale_h.o \
$(ESPRESSO_DIR)/PW/scf_mod.o \
$(ESPRESSO_DIR)/PW/seqopn.o \
$(ESPRESSO_DIR)/PW/set_hubbard_l.o \
$(ESPRESSO_DIR)/PW/set_kplusq.o \
$(ESPRESSO_DIR)/PW/set_kup_and_kdw.o \
$(ESPRESSO_DIR)/PW/set_rhoc.o \
$(ESPRESSO_DIR)/PW/set_vrs.o \
$(ESPRESSO_DIR)/PW/setlocal.o \
$(ESPRESSO_DIR)/PW/setqf.o \
$(ESPRESSO_DIR)/PW/setup.o \
$(ESPRESSO_DIR)/PW/spinor.o \
$(ESPRESSO_DIR)/PW/sph_ind.o \
$(ESPRESSO_DIR)/PW/struct_fact.o \
$(ESPRESSO_DIR)/PW/sum_band.o \
$(ESPRESSO_DIR)/PW/sumkg.o \
$(ESPRESSO_DIR)/PW/sumkt.o \
$(ESPRESSO_DIR)/PW/summary.o \
$(ESPRESSO_DIR)/PW/symm_base.o \
$(ESPRESSO_DIR)/PW/symme.o \
$(ESPRESSO_DIR)/PW/tabd.o \
$(ESPRESSO_DIR)/PW/transform_becsum_nc.o \
$(ESPRESSO_DIR)/PW/transform_becsum_so.o \
$(ESPRESSO_DIR)/PW/trnvecc.o \
$(ESPRESSO_DIR)/PW/tweights.o \
$(ESPRESSO_DIR)/PW/update_pot.o \
$(ESPRESSO_DIR)/PW/usnldiag.o \
$(ESPRESSO_DIR)/PW/v_of_rho.o \
$(ESPRESSO_DIR)/PW/vhpsi.o \
$(ESPRESSO_DIR)/PW/vloc_of_g.o \
$(ESPRESSO_DIR)/PW/vloc_psi.o \
$(ESPRESSO_DIR)/PW/wannier_check.o \
$(ESPRESSO_DIR)/PW/wannier_clean.o \
$(ESPRESSO_DIR)/PW/wannier_init.o \
$(ESPRESSO_DIR)/PW/weights.o \
$(ESPRESSO_DIR)/PW/wfcinit.o \
$(ESPRESSO_DIR)/PW/write_ns.o \
$(ESPRESSO_DIR)/PW/wsweight.o \
$(ESPRESSO_DIR)/PW/xc_vdW_DF.o


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

PPOBJS = $(ESPRESSO_DIR)/PP/xsf.o $(ESPRESSO_DIR)/PP/cube.o

PW_QDIAGP = \
$(ESPRESSO_DIR)/flib/date_and_tim.o \
$(ESPRESSO_DIR)/flib/flush_unit.o \
$(ESPRESSO_DIR)/flib/inpfile.o

MODULES = \
$(ESPRESSO_DIR)/Modules/atom.o \
$(ESPRESSO_DIR)/Modules/basic_algebra_routines.o \
$(ESPRESSO_DIR)/Modules/cell_base.o \
$(ESPRESSO_DIR)/Modules/check_stop.o \
$(ESPRESSO_DIR)/Modules/clocks.o \
$(ESPRESSO_DIR)/Modules/constants.o \
$(ESPRESSO_DIR)/Modules/constraints_module.o \
$(ESPRESSO_DIR)/Modules/control_flags.o \
$(ESPRESSO_DIR)/Modules/coulomb_vcut.o \
$(ESPRESSO_DIR)/Modules/descriptors.o \
$(ESPRESSO_DIR)/Modules/dspev_drv.o \
$(ESPRESSO_DIR)/Modules/electrons_base.o \
$(ESPRESSO_DIR)/Modules/environment.o \
$(ESPRESSO_DIR)/Modules/error_handler.o \
$(ESPRESSO_DIR)/Modules/fft_base.o \
$(ESPRESSO_DIR)/Modules/fft_interfaces.o \
$(ESPRESSO_DIR)/Modules/fft_parallel.o \
$(ESPRESSO_DIR)/Modules/fft_scalar.o \
$(ESPRESSO_DIR)/Modules/fft_types.o \
$(ESPRESSO_DIR)/Modules/funct.o \
$(ESPRESSO_DIR)/Modules/griddim.o \
$(ESPRESSO_DIR)/Modules/image_io_routines.o \
$(ESPRESSO_DIR)/Modules/input_parameters.o \
$(ESPRESSO_DIR)/Modules/io_files.o \
$(ESPRESSO_DIR)/Modules/io_global.o \
$(ESPRESSO_DIR)/Modules/ions_base.o \
$(ESPRESSO_DIR)/Modules/kind.o \
$(ESPRESSO_DIR)/Modules/mm_dispersion.o \
$(ESPRESSO_DIR)/Modules/mp.o \
$(ESPRESSO_DIR)/Modules/mp_base.o \
$(ESPRESSO_DIR)/Modules/mp_global.o \
$(ESPRESSO_DIR)/Modules/mp_wave.o \
$(ESPRESSO_DIR)/Modules/parallel_include.o \
$(ESPRESSO_DIR)/Modules/parameters.o \
$(ESPRESSO_DIR)/Modules/parser.o \
$(ESPRESSO_DIR)/Modules/paw_variables.o \
$(ESPRESSO_DIR)/Modules/printout_base.o \
$(ESPRESSO_DIR)/Modules/pseudo_types.o \
$(ESPRESSO_DIR)/Modules/ptoolkit.o \
$(ESPRESSO_DIR)/Modules/radial_grids.o \
$(ESPRESSO_DIR)/Modules/random_numbers.o \
$(ESPRESSO_DIR)/Modules/read_ncpp.o \
$(ESPRESSO_DIR)/Modules/read_upf_v1.o \
$(ESPRESSO_DIR)/Modules/read_upf_v2.o \
$(ESPRESSO_DIR)/Modules/read_uspp.o \
$(ESPRESSO_DIR)/Modules/recvec.o \
$(ESPRESSO_DIR)/Modules/recvec_subs.o \
$(ESPRESSO_DIR)/Modules/splinelib.o \
$(ESPRESSO_DIR)/Modules/stick_base.o \
$(ESPRESSO_DIR)/Modules/stick_set.o \
$(ESPRESSO_DIR)/Modules/timestep.o \
$(ESPRESSO_DIR)/Modules/uspp.o \
$(ESPRESSO_DIR)/Modules/upf.o \
$(ESPRESSO_DIR)/Modules/upf_to_internal.o \
$(ESPRESSO_DIR)/Modules/version.o \
$(ESPRESSO_DIR)/Modules/wannier.o \
$(ESPRESSO_DIR)/Modules/wannier_new.o \
$(ESPRESSO_DIR)/Modules/wavefunctions.o \
$(ESPRESSO_DIR)/Modules/wrappers.o \
$(ESPRESSO_DIR)/Modules/write_upf_v2.o \
$(ESPRESSO_DIR)/Modules/ws_base.o \
$(ESPRESSO_DIR)/Modules/xml_io_base.o \
$(ESPRESSO_DIR)/Modules/zhpev_drv.o


OCEAN_BUILDER = OCEAN_builder.o mkvipt.o intval.o getomega.o OCEAN_timer.o realspace_paw.o getprefs.o getylm.o OCEAN_build_chi.o $(NEWDIAGOBJ)

$(OCEAN_BUILDER_EXE) :  $(OCEAN_BUILDER) $(PPOBJS) $(PW_QDIAGP) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
	$(OCEAN_BUILDER) $(PPOBJS) $(PWOBJS) $(ESPRESSO_DIR)/SHIRLEY/read_file_shirley.o \
	$(PW_QDIAGP) $(MODULES) $(LIBOBJS) $(LIBS) 

$(SHIRLEY_HAM_EXE) : $(SHIRLEY_HAM) $(SHIRLEY_HAM_O) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(OCEAN_O2L)
	$(MPIF90) $(LDFLAGS) -o $@ \
	$(SHIRLEY_HAM) $(SHIRLEY_HAM_O) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(OCEAN_O2L) $(LIBS)


clean : 
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i work.pc *__genmod.f90

