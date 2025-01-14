if (DEBUG)
  string(APPEND CFLAGS " -g")
endif()
string(APPEND CFLAGS " -std=c99")
string(APPEND CPPDEFS " -DFORTRANUNDERSCORE -DNO_CRAY_POINTERS -DNO_SHR_VMATH -DCPRNAG")
set(FC_AUTO_R8 "-r8")
set(FFLAGS "-wmismatch=mpi_send,mpi_recv,mpi_bcast,mpi_allreduce,mpi_reduce,mpi_isend,mpi_irecv,mpi_irsend,mpi_rsend,mpi_gatherv,mpi_gather,mpi_scatterv,mpi_allgather,mpi_alltoallv,mpi_file_read_all,mpi_file_write_all,mpibcast,mpiscatterv,mpi_alltoallw,nfmpi_get_vara_all,NFMPI_IPUT_VARA,NFMPI_GET_VAR_ALL,NFMPI_PUT_VARA,NFMPI_PUT_ATT_REAL,NFMPI_PUT_ATT_DOUBLE,NFMPI_PUT_ATT_INT,NFMPI_GET_ATT_REAL,NFMPI_GET_ATT_INT,NFMPI_GET_ATT_DOUBLE,NFMPI_PUT_VARA_DOUBLE_ALL,NFMPI_PUT_VARA_REAL_ALL,NFMPI_PUT_VARA_INT_ALL    -convert=BIG_ENDIAN")
if (NOT DEBUG)
  string(APPEND FFLAGS " -ieee=full -O2")
endif()
if (DEBUG)
  string(APPEND FFLAGS " -g -time -f2003 -ieee=stop")
endif()
if (DEBUG)
  string(APPEND FFLAGS " -C=all -g -time -f2003 -ieee=stop")
endif()
if (DEBUG AND NOT compile_threaded)
  string(APPEND FFLAGS " -gline")
endif()
if (compile_threaded)
  string(APPEND FFLAGS " -openmp")
endif()
set(FIXEDFLAGS "-fixed")
set(FREEFLAGS "-free")
set(HAS_F2008_CONTIGUOUS "FALSE")
if (compile_threaded)
  string(APPEND LDFLAGS " -openmp")
endif()
set(MPICC "mpicc")
set(MPIFC "mpif90")
set(SCC "gcc")
set(SFC "nagfor")
