#  CMakeLists.txt

cmake_minimum_required(VERSION 3.1.0)

################## compiler options ##################
# compiler
set(CMAKE_CXX_COMPILER          ${CXX})
set(CMAKE_FORTRAN_COMPILER     ${FC})

# compiler flags
set(CMAKE_CXX_FLAGS             ${CXXFLAGS}) 

################## dependencies ##################
# preferred installation prefix of dependencies
set(EIGEN_INC                 ${EIGEN_INC})

# set default values
macro(setDefault var val)
if(NOT DEFINED ${var})
set(${var} ${val})
endif()
endmacro()
setDefault(EIGEN_INC   "/mnt/c/linuxapp/eigen-3.4.0")

# project
project(JDSurfG LANGUAGES CXX Fortran)

# build type 
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/build/obj)

################# cmake setup #################
set(CMAKE_CXX_FLAGS  "-g -DEIGEN_NO_DEBUG -DEIGEN_DONT_PARALLELIZE")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math")
  set(CMAKE_Fortran_FLAGS "-O3 -march=native -ffixed-line-length-none")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -xHost -DEIGEN_STRONG_INLINE=inline")
  set(CMAKE_Fortran_FLAGS "-O3 -xHost -extend-source")
endif()

set(CMAKE_VERBOSE_MAKEFILE off)

################# find packages #################
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
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

target_link_libraries(mkmat gravity share)
target_link_libraries(syngrav gravity SWD share)
target_link_libraries(DSurfTomo surftomo SWD share)
target_link_libraries(JointTomo jointomo SWD gravity share)