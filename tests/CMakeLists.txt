set(TESTDIR ${CMAKE_BINARY_DIR}/tests)

set(TESTLIST "test1" "test2" "test3" "test4" "test5")
if(PERFORMANCETEST)
    set(TEST_INIT_OPTION c)
    set(TEST_PREFIX performance)
else()
    set(TEST_INIT_OPTION r)
    set(TEST_PREFIX )
endif()

foreach(ITEM ${TESTLIST})
    configure_file(${CMAKE_SOURCE_DIR}/tests/${ITEM}.c.in ${TESTDIR}/${ITEM}.c)
    add_executable(${TEST_PREFIX}${ITEM} ${TESTDIR}/${ITEM}.c)
    target_link_libraries(${TEST_PREFIX}${ITEM} T3NS-shared)
    add_test(${TEST_PREFIX}${ITEM} ${TESTDIR}/${TEST_PREFIX}${ITEM})
    set_tests_properties(${TEST_PREFIX}${ITEM}
	PROPERTIES
	TIMEOUT 1200)
endforeach()
