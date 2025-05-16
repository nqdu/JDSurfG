#  CMakeLists.txt

cmake_minimum_required(VERSION 3.10.0)

################## compiler options ##################
# compiler
set(CMAKE_CXX_COMPILER          ${CXX})
set(CMAKE_FORTRAN_COMPILER     ${FC})

# compiler flags
set(CMAKE_CXX_FLAGS             ${CXXFLAGS}) 

# build function test
set(BUILD_TEST                  ${BUILD_TEST})

################## dependencies ##################
# preferred installation prefix of dependencies
set(EIGEN_INC                 ${EIGEN_INC})

# set default values
macro(setDefault var val)
if(NOT DEFINED ${var})
set(${var} ${val})
endif()
endmacro()
setDefault(EIGEN_INC   "~/software/eigen-3.4.0")
setDefault(BUILD_TEST   FALSE)

# project
project(JDSurfG LANGUAGES CXX Fortran)

# build type 
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/obj)

################# cmake setup #################
set(CMAKE_CXX_FLAGS  "-g  -DEIGEN_DONT_PARALLELIZE")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(CMAKE_CXX_FLAGS_RELEASE " ${CMAKE_CXX_FLAGS_RELEASE} -march=native")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native -ffixed-line-length-none")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE " ${CMAKE_CXX_FLAGS_RELEASE} -xHost -DEIGEN_STRONG_INLINE=inline")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xHost -extend-source")
else () # mac
  set(CMAKE_CXX_FLAGS_RELEASE " ${CMAKE_CXX_FLAGS_RELEASE} -march=native")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native -ffixed-line-length-none")
endif()

set(CMAKE_VERBOSE_MAKEFILE on)

################# find packages #################
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  include_directories(${OpenMP_CXX_INCLUDE_DIRS})
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# eigen
include_directories(${EIGEN_INC})

################# include paths #################
include_directories(${PROJECT_SOURCE_DIR}/src  ${PROJECT_SOURCE_DIR}/include)

# subs
add_subdirectory(${PROJECT_SOURCE_DIR}/src/shared)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/gravity)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/SWD)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/surftomo)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/JSurfGTomo)

# exe
add_executable(mkmat ${PROJECT_SOURCE_DIR}/src/gravity/main.cpp)
add_executable(syngrav ${PROJECT_SOURCE_DIR}/utils/syn_gravity.cpp)
add_executable(DSurfTomo ${PROJECT_SOURCE_DIR}/src/surftomo/main.cpp)
add_executable(JointTomo ${PROJECT_SOURCE_DIR}/src/JSurfGTomo/main.cpp)

if (${BUILD_TEST}) 
  message(STATUS "build test examples = ON ...")
  add_executable(test_fmst ${PROJECT_SOURCE_DIR}/test/fmst_deriv/main.cpp)
  target_link_libraries(test_fmst surftomo SWD share)
  add_executable(test_smooth ${PROJECT_SOURCE_DIR}/test/smooth/main.cpp)
  target_link_libraries(test_smooth share)
  add_executable(test_swd ${PROJECT_SOURCE_DIR}/test/swd/main.cpp)
  target_link_libraries(test_swd SWD share)
endif()

target_link_libraries(mkmat gravity share OpenMP::OpenMP_CXX)
target_link_libraries(syngrav gravity SWD share OpenMP::OpenMP_CXX)
target_link_libraries(DSurfTomo surftomo SWD share OpenMP::OpenMP_CXX)
target_link_libraries(JointTomo jointomo SWD gravity share OpenMP::OpenMP_CXX)