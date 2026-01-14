cmake_minimum_required(VERSION 3.10)
project(PiCalculation)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- 1. Define Flags as Lists ---
set(OPT_FLAGS  -g -Ofast -ffast-math -funroll-loops)
set(CHIP_FLAGS -march=native -mtune=native)
set(PIC_FLAGS  -fPIC)

# Combine them into one list (CMake automatically handles lists)
set(MY_CXX_FLAGS ${OPT_FLAGS} ${CHIP_FLAGS} ${PIC_FLAGS})

# --- 2. Find Dependencies ---
find_package(OpenMP REQUIRED)

# --- 3. Define Executable ---
add_executable(main pi.cpp xtime.cpp matrix_vector.cpp main.cpp)

# --- 4. Apply Flags (Using the List) ---
# CMake expands the list properly here into separate arguments
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(main PRIVATE ${MY_CXX_FLAGS})
endif()

# --- 5. Link Libraries ---
# Link OpenMP and Math library (-lm) here
if(OpenMP_CXX_FOUND)
    target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX m) # 'm' is libm (-lm)
endif()
