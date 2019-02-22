set(T3NSLIB_SOURCE_FILES
    "Heff.c"
    "Wigner.c"
    "bookkeeper.c"
    "davidson.c"
    "hamiltonian.c"
    "hamiltonian_nn_hubbard.c"
    "hamiltonian_qc.c"
    "hamiltonian_doci.c"
    "instructions.c"
    "instructions_nn_hubbard.c"
    "instructions_qc.c"
    "instructions_doci.c"
    "io.c"
    "io_to_disk.c"
    "macros.c"
    "network.c"
    "opType.c"
    "opType_qc.c"
    "optScheme.c"
    "optimize_network.c"
    "rOperators_bUpdate.c"
    "rOperators_init.c"
    "rOperators_misc.c"
    "rOperators_pUpdate.c"
    "siteTensor_decompose.c"
    "siteTensor_decomposeNEW.c"
    "siteTensor_init.c"
    "siteTensor_misc.c"
    "sort.c"
    "sparseblocks.c"
    "symmetries.c"
    "symmetry_pg.c"
    "symmetry_su2.c"
    "symmetry_u1.c"
    "symmetry_z2.c"
    "symsecs.c"
    "tensorproducts.c"
    "wrapper_solvers.c"
    "RedDM.c"
    )

add_library(T3NS-shared SHARED ${T3NSLIB_SOURCE_FILES})
target_link_libraries(T3NS-shared ${LAPACK_LIBRARIES} ${HDF5_LIBRARIES})
set_target_properties(T3NS-shared PROPERTIES OUTPUT_NAME "T3NS" EXPORT_NAME "T3NS")

add_executable(T3NS-bin executable.c)
target_link_libraries(T3NS-bin T3NS-shared)
set_target_properties(T3NS-bin PROPERTIES OUTPUT_NAME "T3NS")

install(TARGETS T3NS-shared T3NS-bin
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
