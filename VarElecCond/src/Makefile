DOBJS     = types.o comms.o params.o dimens.o trans_params.o control.o fft.o \
	arrays.o timing.o initial.o constants.o netcdf.o

OBJS = utils.o dbydz.o magnetic.o hydro_deriv.o

SRC = main.F90 ranset.f90

DERIVSRC = test_dbydxy.F90 ranset.f90

all: prog derivtest

prog : $(DOBJS) $(OBJS) $(SRC)
	$(FC) -o $(EXEC) $(FFLAGS) $(FPPFLAGS) $(SRC) $(DOBJS) $(OBJS) $(INCLUDE) $(FFTLIB) $(NETCDFLIB)
	mv $(EXEC) ../bin/


derivtest    : $(DOBJS) $(DERIVSRC)
	$(FC) -o $(EXECTEST) $(FFLAGS) $(FPPFLAGS) $(DERIVSRC) $(DOBJS) $(INCLUDE) $(FFTLIB) $(NETCDFLIB)
	mv $(EXECTEST) ../bin/


initial.o : initial.F90 types.o dimens.o comms.o params.o trans_params.o \
	arrays.o fft.o control.o constants.o netcdf.o
	$(FC) $(FFLAGS) $(FPPFLAGS) -c initial.F90

fft.o	: fft.F90 types.o dimens.o comms.o timing.o \
	params.o trans_params.o control.o
	$(FC) $(FFLAGS) $(FPPFLAGS) -c fft.F90

timing.o : timing.F90 types.o control.o
	$(FC) $(FFLAGS) $(FPPFLAGS) -c timing.F90

control.o : control.f90
	 $(FC) $(FFLAGS) $(FPPFLAGS) -c control.f90

dimens.o : dimens.f90
	 $(FC) $(FFLAGS) $(FPPFLAGS) -c dimens.f90

params.o : params.f90 types.o
	 $(FC) $(FFLAGS) $(FPPFLAGS) -c params.f90

trans_params.o : trans_params.f90  types.o dimens.o comms.o params.o \
	constants.o arrays.o
	 $(FC) $(FFLAGS) $(FPPFLAGS) -c trans_params.f90

comms.o	: comms.F90 types.o
	$(FC) $(FFLAGS) $(FPPFLAGS) -c comms.F90

types.o	: types.f90
	$(FC) $(FFLAGS) $(FPPFLAGS) -c types.f90

constants.o : constants.f90 types.o dimens.o params.o trans_params.o
	$(FC) $(FFLAGS) $(FPPFLAGS) -c constants.f90

arrays.o : arrays.f90 types.o dimens.o comms.o
	$(FC) $(FFLAGS) $(FPPFLAGS) -c arrays.f90

utils.o : utils.F90 types.o dimens.o params.o trans_params.o constants.o \
	arrays.o comms.o fft.o dbydz.o
	$(FC) $(FFLAGS) $(FPPFLAGS) -c utils.F90

dbydz.o : dbydz.F90 types.o dimens.o comms.o timing.o params.o trans_params.o \
	constants.o arrays.o control.o
	$(FC) $(FFLAGS) $(FPPFLAGS) -c dbydz.F90

hydro_deriv.o : hydro_deriv.f90 types.o dimens.o comms.o timing.o params.o trans_params.o \
	constants.o arrays.o control.o fft.o dbydz.o
	$(FC) $(FFLAGS) $(FPPFLAGS) -c hydro_deriv.f90

magnetic.o : magnetic.f90 types.o dimens.o comms.o timing.o params.o trans_params.o \
	constants.o arrays.o control.o fft.o dbydz.o
	$(FC) $(FFLAGS) $(FPPFLAGS) -c magnetic.f90

netcdf.o : netcdf.F90 types.o dimens.o params.o trans_params.o arrays.o comms.o
	$(FC) $(FFLAGS) $(FPPFLAGS) -c netcdf.F90 $(INCLUDE)

clean:
	rm -f *.o *.oo *~ *.d work.pc work.pcl *.mod $(EXEC) $(EXECTEST)
