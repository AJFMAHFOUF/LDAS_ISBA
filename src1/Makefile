SOURCE =\
mod_const.f90 mod_forc.f90 mod_soil.f90 mod_surf1.f90 mod_setup.f90 \
master.f90 main_oi.f90 main_enkf.f90 main.f90 \
isba.f90 energy_budget.f90 \
clean_exit.f90 drag_coeff_z0h.f90 fluxes.f90 \
init_surf1.f90 interpol_forcing.f90 read_forcing.f90 \
rs_soil.f90 rs_veg.f90 soil_def.f90 soil_prop.f90 \
vdfppcfls.f90 water_budget.f90 thermo_functions.f90 \
print_output_oi.f90 print_output_enkf.f90 cholsl.f90 \
choldc.f90 solar_angle.f90 \
oi_coeffs_ec.f90 oi_coeffs_mf.f90 gasdev.f90 out_product.f90 \
inverse_matrix.f90 cls_interpol.f90 vdfppcfls2.f90
#FC = g95
#FFLAGS = -r8
FC = gfortran
FFLAGS = -fdefault-real-8	

OBJECTS	=	$(SOURCE:.f90=.o)

.SUFFIXES:
.SUFFIXES:	.o .f90

.f90.o :
	$(FC)  $(FFLAGS) -c $*.f90

main:	$(OBJECTS)
	$(FC) -o main-ref ${OBJECTS} 

clean:
	 \rm -f $(OBJECTS) fort.* *~

fort:
	 \rm -f  fort.* *~
