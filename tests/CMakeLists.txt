# Copyright (c) 2024 Sebastian Sassi

# Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the "Software"), to deal # in the Software without restriction, including without limitation the rights # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in 
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
# SOFTWARE.
macro(create_test TESTNAME)
    add_executable(${TESTNAME}.test ${TESTNAME}.cpp)
    target_include_directories(${TESTNAME}.test
        PRIVATE ${PROJECT_SOURCE_DIR}/include
    )
    target_compile_options(${TESTNAME}.test
        PRIVATE
            $<$<CONFIG:Debug>:-Wall;-Wextra;-Wpedantic;-Wconversion;-Werror;-g;-O0>
    )

    target_compile_features(${TESTNAME}.test PUBLIC cxx_std_20)

    target_link_libraries(${TESTNAME}.test cubage)
    add_test(NAME ${TESTNAME} COMMAND ${TESTNAME}.test)
endmacro()

create_test(test_box)
create_test(test_cubage)
create_test(test_gauss_kronrod)
create_test(test_genz_malik)