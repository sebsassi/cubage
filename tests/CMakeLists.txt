macro(create_test TESTNAME)
    add_executable(${TESTNAME}.test ${TESTNAME}.cpp)
    target_include_directories(${TESTNAME}.test
        PRIVATE ${PROJECT_SOURCE_DIR}/include
    )
    target_compile_options(${TESTNAME}.test
        PRIVATE
            $<$<CONFIG:Debug>:-Wall;-Wextra;-Wpedantic;-Wconversion;-Werror;-g;-O0>
    )
    target_link_libraries(${TESTNAME}.test cubage)
    add_test(NAME ${TESTNAME} COMMAND ${TESTNAME}.test)
endmacro()

create_test(test_box)
create_test(test_cubage)
create_test(test_gauss_kronrod)
create_test(test_genz_malik)