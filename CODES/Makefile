#Set compiler
FC  = ifort
LD  = ifort
NETCDF_INCDIR=/usr/local/netcdf-fortran-4.4.4/include
NETCDF_LIBDIR=-L/usr/local/netcdf-fortran-4.4.4/lib -lnetcdff -L/usr/local/hdf5/lib -L/usr/local/netcdf4/lib -lnetcdf -lhdf5_hl -lhdf5 -lcurl -I/usr/local/netcdf4/include -lnetcdf -L/usr/local/netcdf-4.7.3/lib
#-L/usr/local/netcdf-fortran-4.4.4/lib -lnetcdff -L/usr/local/hdf5/lib -L/usr/local/netcdf4/lib -L/usr/local/lib -lnetcdf -L/opt/intel/lib/intel64_lin -lhdf
FFLAGS += -I$(NETCDF_INCDIR)
#FFLAGS +=-g -check all,noarg_temp_created
LDFLAGS += $(NETCDF_LIBDIR) -lnetcdf -lnetcdff 
LDFLAGS += -fpp -qopenmp 
TARGET = 1d.out
TARGET_LIN = do_lin.out
OBJ  = param.o \
       arrays_sub.o \
       calendar_sub.o \
       mod_hist.o \
       mod_avg.o \
       mod_diags.o \
       mod_init.o \
       mod_rst.o \
       mod_atm.o \
       mod_clm.o \
       mod_adv.o \
       bulk_sub.o \
       eos_sub.o \
       solve_diag.o \
       input_files.o \
       kc_sub.o \
       kpp_sub.o \
       nnf_sub.o
OBJ_adv  = mod_adv.o
OBJ_ncdf = ncdf_write.o
OBJ_main = main.o
OBJ_lin = mod_lin.o \
	lin_main.o
MODS  = ncdf_write.mod \
        param.mod \
        bulk_sub.mod \
        eos_sub.mod \
        calendar_sub.mod \
        solve_diag.mod \
        arrays_sub.mod \
        mod_hist.mod \
        mod_avg.mod \
        mod_diags.mod \
        mod_init.mod \
        mod_rst.mod \
        mod_atm.mod \
        mod_clm.mod \
        mod_adv.mod \
        kc_sub.mod \
        nnf_sub.mod \
        kpp_sub.mod \
        input_files.mod 
all: $(TARGET)
$(TARGET): $(OBJ_ncdf) $(OBJ) $(OBJ_main)
	$(LD) -o $@ $(OBJ_ncdf) $(OBJ) $(OBJ_main) ${LDFLAGS} ${FFLAGS}
$(TARGET_LIN): $(OBJ_ncdf) $(OBJ) $(OBJ_lin)
	$(LD) -o $@ $(OBJ_ncdf) $(OBJ) $(OBJ_lin) ${LDFLAGS} ${FFLAGS}
#Objects--------------------------------------
ncdf_write.o:ncdf_write.f90
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
param.o:param.f90
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
arrays_sub.o:arrays_sub.f90 param.f90
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
mod_hist.o:mod_hist.f90 param.f90 CPPLISTS.h
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
mod_avg.o:mod_avg.f90 calendar_sub.f90 param.f90 CPPLISTS.h
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
mod_diags.o:mod_diags.f90 param.f90 CPPLISTS.h
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
mod_init.o:mod_init.f90 param.f90
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
mod_rst.o:mod_rst.f90 param.f90
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
mod_lin.o:mod_lin.f90 param.f90
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
mod_atm.o:mod_atm.f90 param.f90 CPPLISTS.h
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
mod_clm.o:mod_clm.f90 param.f90
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
mod_adv.o:mod_adv.f90 param.f90
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
bulk_sub.o:bulk_sub.f90 
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
eos_sub.o:eos_sub.f90 
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
calendar_sub.o:calendar_sub.f90
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
input_files.o:input_files.f90 param.f90
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
output_files.o:output_files.f90
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
solve_diag.o: solve_diag.f90 CPPLISTS.h
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
kc_sub.o: kc_sub.f90 eos_sub.o param.o
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
nnf_sub.o: nnf_sub.f90 eos_sub.o param.o
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
kpp_sub.o: kpp_sub.f90 eos_sub.o param.o
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
main.o:main.f90 bulk_sub.o eos_sub.o param.o CPPLISTS.h
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
lin_main.o:lin_main.f90 bulk_sub.o eos_sub.o param.o CPPLISTS.h
	$(FC) -c $< ${LDFLAGS} ${FFLAGS} ${CPPFLAGS}
clean:
	rm -f $(TARGET) $(OBJ) $(OBJ_ncdf) $(OBJ_main) $(OBJ_nudge) $(MODS)
