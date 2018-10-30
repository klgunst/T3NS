set(TESTDIR ${CMAKE_BINARY_DIR}/tests)

set(TESTLIST "test1" "test2" "test3" "test4")

foreach(ITEM ${TESTLIST})
    configure_file(${CMAKE_SOURCE_DIR}/tests/${ITEM}.c.in ${TESTDIR}/${ITEM}.c)
    add_executable(${ITEM} ${TESTDIR}/${ITEM}.c)
    target_link_libraries(${ITEM} T3NS-shared)
    add_test(${ITEM} ${TESTDIR}/${ITEM})
    set_tests_properties(${ITEM}
	PROPERTIES
	PASS_REGULAR_EXPRESSION "Test passed")
    set_tests_properties(${ITEM}
	PROPERTIES
	FAIL_REGULAR_EXPRESSION "(Error|Test failed)")
    set_tests_properties(${ITEM}
	PROPERTIES
	TIMEOUT 1200)
endforeach()
