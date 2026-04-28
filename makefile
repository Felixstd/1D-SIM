## INSTRUCTION: make EXE='executable_filename' ice
# Makefile for my sea-ice model
FC      = gfortran
LIBS    = 
EXE     = 'zoupa'
BUILDDIR = build

# Module files must come BEFORE files that use them.
# ice.f90 (main program) goes last.
SRC = src/parameter_mod.f90 \
      src/constants_mod.f90 \
      src/option_mod.f90 \
      src/diag_stress_mod.f90 \
      src/global_var_mod.f90 \
      src/shallow_water_mod.f90 \
      src/ini_get.f90 \
      src/wind_forcing.f90 \
      src/ice_strength.f90 \
      src/viscouscoefficient.f90 \
      src/Cw_coefficient.f90 \
      src/Fu.f90 \
      src/SOR.f90 \
      src/advection.f90 \
      src/util.f90 \
      src/prep_fgmres_NK.f90 \
      src/identity.f90 \
      src/Jacobian.f90 \
      src/fgmresD.f90 \
      src/dcopy.f90 \
      src/ddot.f90 \
      src/daxpy.f90 \
      src/bvect.f90 \
      src/output_results.f90 \
      src/output_file.f90 \
      src/calc_scaling.f90 \
      src/EVP2solver.f90 \
      src/mu_phi.f90 \
      src/energy_dissipation.f90 \
      src/mechanical_energy.f90 \
      src/par_get.f90 \
      src/ice.f90

ice: $(BUILDDIR) $(SRC)
	$(FC) $(FFLAGS) -J $(BUILDDIR) -I $(BUILDDIR) -o $(EXE) $(SRC) #$(LIBS)
	rm -f *.mod
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

tidy:
	rm -f ./*.o src/*~ ./zoupa

clean:
	rm -rf $(BUILDDIR)
	rm -f ./*.o ./zoupa
