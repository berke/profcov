cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
enable_language(Fortran)
enable_language(C)
project(scarbo NONE)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    if(DEBUG)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all -Wall -g -Wimplicit-interface \
                         -Waliasing -Wconversion -Wintrinsic-shadow -march=native -cpp \
                         -std=f2018")
    else()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=bounds -Wall -Wimplicit-interface \
                         -Waliasing -Wconversion -Wintrinsic-shadow -march=native -O3 -cpp \
                         -std=f2018")
    endif()
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I/usr/include/x86_64-linux-gnu -Wall -g -O3")

find_package(sexptran-${CMAKE_Fortran_COMPILER_ID} 1.3.1 REQUIRED )
include_directories("${sexptran_INCLUDE_DIRS}" )

set(SRCS profcov.f90)

set(LAPACK_LIBS blas lapack)

include_directories(/usr/local/hdf5/include)
link_directories(/usr/local/hdf5/lib)
set(HDF_LIBS hdf5_fortran hdf5 z dl)

set(ALL_LIBS ${HDF_LIBS} ${LAPACK_LIBS})

add_executable(profcov profcov.f90 ${SRCS})
target_link_libraries(profcov ${ALL_LIBS})
