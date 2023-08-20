# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ${MOAB_DIR}/share/examples/makefile.config

default: ExchangeHalos
all: ExchangeHalos

ExchangeHalos: Driver.o ExchangeHalos.o ${MOAB_LIBDIR}/libMOAB.la
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	@echo "  [LD]   ExchangeHalos..."
	${VERBOSE}${MOAB_CXX} Driver.o ExchangeHalos.o ${MOAB_LIBS_LINK} -o ExchangeHalos
endif

run: ExchangeHalos
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	${RUNSERIAL} ./ExchangeHalos
	${RUNPARALLEL} ./ExchangeHalos
endif

clean: clobber
	rm -rf ExchangeHalos
