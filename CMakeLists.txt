cmake_minimum_required(VERSION 3.15)

project(test_cubed_sphere_read
        VERSION 1.0.0
        LANGUAGES C Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/Modules")

include(${PROJECT_NAME}_compiler_flags)
include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Intel)$")
  message(WARNING "${CMAKE_Fortran_COMPILER_ID} is not supported.")
endif()

option(ENABLE_NCDIAG_SERIAL "Enable Serial NetCDF diagnostics" ON)

# Dependencies
find_package(MPI REQUIRED)
find_package(ESMF REQUIRED)
find_package(NetCDF REQUIRED COMPONENTS Fortran)

# Get list of files
list(APPEND READ_CUBED_SPHERE_SRC driver.F90  module_read_netcdf.F90)

add_executable(test_read ${READ_CUBED_SPHERE_SRC})
target_link_libraries(test_read esmf NetCDF::NetCDF_Fortran /data/jedi/ufs-build/FV3/atmos_cubed_sphere/libfv3.a /home/mpotts/spack-stack/envs/ufs-jedi.gnu/install/gcc/9.3.0/fms-2022.01-37a7p6b/lib/libfms_r8.a /data/jedi/ufs-build/FV3/ccpp/libfv3ccpp.a /data/jedi/ufs-build/FV3/ccpp/physics/libccpp_physics.a /home/mpotts/spack-stack/envs/ufs-jedi.gnu/install/gcc/9.3.0/w3nco-2.4.1-4frea23/lib/libw3nco_d.a /home/mpotts/spack-stack/envs/ufs-jedi.gnu/install/gcc/9.3.0/sp-2.3.3-oj6qi5u/lib/libsp_d.a /home/mpotts/spack-stack/envs/ufs-jedi.gnu/install/gcc/9.3.0/bacio-2.4.1-jnhzejy/lib/libbacio_8.a   )
target_include_directories(test_read PRIVATE ${ESMF_F90COMPILEPATHS} /data/jedi/ufs-build/FV3/atmos_cubed_sphere/include/fv3 /data/jedi/ufs-build/FV3/ccpp/physics /home/mpotts/spack-stack/envs/ufs-jedi.gnu/install/gcc/9.3.0/fms-2022.01-37a7p6b/include_r8 )
