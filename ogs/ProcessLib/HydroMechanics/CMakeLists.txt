ogs_add_process(HydroMechanics)

if(OGS_BUILD_TESTING AND OGS_USE_MFRONT AND NOT OGS_USE_MPI AND NOT OGS_USE_PETSC)
    ogs_add_executable(hm-substepping-tests
        ${PROJECT_SOURCE_DIR}/Tests/MaterialLib/MFront/TestHMSubstepping.cpp
        ${PROJECT_SOURCE_DIR}/Tests/TestTools.cpp)
    target_link_libraries(hm-substepping-tests PRIVATE
        MaterialLib_SolidModels_MFront ParameterLib MeshLib gtest)
    add_test(NAME HMSubsteppingTransactions COMMAND hm-substepping-tests)
    set_tests_properties(HMSubsteppingTransactions PROPERTIES TIMEOUT 60
        ENVIRONMENT "OMP_NUM_THREADS=1;OGS_ASM_THREADS=1")
endif()
