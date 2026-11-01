cmake_minimum_required(VERSION 3.20)

project(ADCSTest)
set(EXECUTABLE adcs-test)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# glob CMSIS-DSP library sources only (skip Examples/ComputeLibrary/Lapack, which need embedded/ARM-only headers)
file(GLOB_RECURSE SRC CONFIGURE_DEPENDS linalg/Source/*.c)
# each Source/<Group> folder has an amalgamation file (e.g. BasicMathFunctions.c) that #includes
# every arm_*.c in that folder, so compiling the individual arm_*.c files too causes duplicate symbols
list(FILTER SRC EXCLUDE REGEX "/Source/[^/]+/arm_[^/]*\\.c$")
# exclude underscore-prefixed files: these are #include-only fragments (e.g. Neon kernel templates),
# not standalone translation units, and fail to compile on their own
list(FILTER SRC EXCLUDE REGEX "/_[^/]*\\.c$")
# manual specified sources:
set(ADCS_SRC src/main.c src/test.c src/quat.c src/laextension.c src/bdot.c src/body.c src/down_quat.c src/ecef2eci.c src/filters.c src/iterate.c src/kepler.c src/magnetosphere.c src/moments2currents.c src/orbital2eci.c src/pd.c src/photodiode_determination.c src/quest.c src/sunvec.c src/torque2moments.c src/sgp4_init.c src/sgp4_propagate.c src/sgp4_step.c src/closest_approach.c)
list(APPEND SRC ${ADCS_SRC})
# extra warnings on our own code only (CMSIS is noisy). -Wshadow / -Wfloat-equal match the PVDX
# firmware build, which compiles with -Werror.
set_source_files_properties(${ADCS_SRC} PROPERTIES COMPILE_OPTIONS "-Wall;-Wextra;-Wshadow;-Wfloat-equal")
set(CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -Wpedantic -g" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)

add_executable(${EXECUTABLE} ${SRC})
target_include_directories(${EXECUTABLE} PRIVATE .)
target_include_directories(${EXECUTABLE} PRIVATE linalg
    linalg/Include
    linalg/PrivateInclude
    linalg/LinearAlgebra
    linalg/Lapack/Blas
    linalg/Lapack/F2c
    linalg/Lapack/Include
    linalg/Lapack/Install
    linalg/Lapack/Scr
)

target_link_libraries(adcs-test m)

target_compile_definitions(${EXECUTABLE} PRIVATE __GNUC_PYTHON__)

# expose headers in an include directory if desired:
# target_include_directories(${EXECUTABLE} PUBLIC include)
