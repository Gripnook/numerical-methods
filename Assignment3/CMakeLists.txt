cmake_minimum_required (VERSION 3.1)
set (CMAKE_CXX_STANDARD 14)

project (ecse543-a3)

# Catch library
add_library (Catch INTERFACE)
target_include_directories (Catch INTERFACE lib/catch)

# Main executable
include (src/CMakeLists.txt)
set (MAIN_EXE_SOURCES src/${MAIN_ENTRY})
foreach (FILE ${MAIN_SOURCES})
	set (MAIN_EXE_SOURCES ${MAIN_EXE_SOURCES} src/${FILE})
endforeach (FILE)

add_executable (ecse543-a3 ${MAIN_EXE_SOURCES})

# Test executable
include (test/CMakeLists.txt)
set (TEST_EXE_SOURCES test/${TEST_ENTRY})
foreach (FILE ${TEST_SOURCES})
	set (TEST_EXE_SOURCES ${TEST_EXE_SOURCES} test/${FILE})
endforeach (FILE)
foreach (FILE ${MAIN_SOURCES})
	set (TEST_EXE_SOURCES ${TEST_EXE_SOURCES} src/${FILE})
endforeach (FILE)

include_directories (src)

add_executable (tests ${TEST_EXE_SOURCES})
target_link_libraries (tests Catch)

enable_testing ()
add_test (tests tests)
