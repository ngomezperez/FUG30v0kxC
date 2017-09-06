#######################################
# Upper level Makefile for parallel code
########################################

# First set executable names
PROG="p3d"
PROGTEST="derivtest"
INTERP="interpolate"
GLUE="glue"
CONV="convert"

# Now set environment variables. There are (currently) six possibilities:
# 1) Newcastle home directory (serial jobs)
#FC="gfortran"
#FFLAGS="-O3"
#FPPFLAGS=""
#FFTLIB="-L/data/npjb7/libs/fftw2/lib -lrfftw -lfftw"
#INCLUDE="-I/data/npjb7/libs/netcdf4/include"
#NETCDFLIB="-L/data/npjb7/libs/netcdf4/lib -lnetcdff"

# 2) Newcastle quantum cluster (parallel jobs)
#FC="gfortran"
#FFLAGS="-O3 -I/usr/lib/openmpi/include/"
#FPPFLAGS="-DMPI "
#FFTLIB="-L/data/npjb7/libs/fftw2/lib -lrfftw -lfftw -limf -lm -llapack"
#INCLUDE="-I/usr/lib/openmpi/include -I/data/npjb7/libs/netcdf3/include"
#NETCDFLIB="-L/data/npjb7/libs/netcdf3/lib -lnetcdf -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -ldl -lnsl -lutil -lm" 

# 3) HECToR XT6 (parallel)
#FC="ftn"
#FFLAGS="-fastsse -Mipa=fast,inline"
#FPPFLAGS="-DMPI " 
#FFTLIB="-L/opt/fftw/2.1.5/cnos/lib -ldrfftw -ldfftw"
#INCLUDE="-I/home/e92/e92/npjb7/libs/netcdf/include"
#NETCDFLIB="-L/home/e92/e92/npjb7/libs/netcdf/lib -lnetcdf"

# 4) UKMHD (serial)
#FC="ifort"
#FFLAGS="-O3"
#FPPFLAGS="" 
#FFTLIB="-L/home/mhd/pjb40/libs/fftw2/lib -lrfftw -lfftw -limf -lm -llapack"
#INCLUDE="-I/home/mhd/pjb40/libs/netcdf3/include"
#NETCDFLIB="-L/home/mhd/pjb40/libs/netcdf3/lib -lnetcdf"

# 5) UKMHD (parallel)                                                         
#FC="mpif90"                                                                   
#FFLAGS="-O3"                                                                 
#FPPFLAGS="-DMPI"                                                             
#FFTLIB="-L/home/mhd/pjb40/libs/fftw2/lib -lrfftw -lfftw -limf -lm -llapack"           
#INCLUDE="-I/home/mhd/pjb40/libs/netcdf3/include"                             
#NETCDFLIB="-L/home/mhd/pjb40/libs/netcdf3/lib -lnetcdf" 

# 6) Laptop (serial)                                                           
FC="gfortran"                                                         
FFLAGS="-O3"
FPPFLAGS=""                                                               
FFTLIB="-L/Users/ngomez/libs/fftw2/lib -lrfftw -lfftw"                         
INCLUDE="-I/Users/ngomez/libs/netcdf3/include"                                 
NETCDFLIB="-L/Users/ngomez/libs/netcdf3/lib -lnetcdf"   

# Now we specify the target subdirectories
SUBDIRS = src interpolate
SUBDIRSP3D = src
SUBDIRSUTIL = interpolate

###################################
# The Make targets are defined here
all: p3d

p3d:
	@for i in $(SUBDIRSP3D) ; do \
		cd $$i ; \
		$(MAKE) 	\
			FC=$(FC) \
			EXEC=$(PROG) \
			EXECTEST=$(PROGTEST) \
			FFLAGS=$(FFLAGS) \
			FPPFLAGS=$(FPPFLAGS) \
			FFTLIB=$(FFTLIB) \
			INCLUDE=$(INCLUDE) \
			NETCDFLIB=$(NETCDFLIB) \
			all ; cd .. ;\
	done

tools:
	@for i in $(SUBDIRSUTIL) ; do \
		cd $$i ; \
		$(MAKE) \
			FC=$(FC) \
			INTERP=$(INTERP) \
			EXECGLUE=$(GLUE) \
			EXECCONVERT=$(CONV) \
			EXECTEST=$(PROGTEST) \
			FFLAGS=$(FFLAGS) \
			FPPFLAGS=$(FPPFLAGS) \
			FFTLIB=$(FFTLIB) \
			INCLUDE=$(INCLUDE) \
			NETCDFLIB=$(NETCDFLIB) \
			all ; cd .. ;\
	done

clean:
	@for i in $(SUBDIRS); do \
		(cd $$i; \
		$(MAKE) clean ); \
	done
	rm -f *~
	rm -f ./bin/$(PROG)
	rm -f ./bin/$(PROGTEST)
