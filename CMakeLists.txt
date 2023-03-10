cmake_minimum_required(VERSION 3.20)

project(fluid_sim)

set(CMAKE_CXX_STANDARD 20)

if (EMSCRIPTEN)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s USE_GLFW=3 -s ASSERTIONS=1 -s WASM=1 -s EXPORTED_FUNCTIONS=\"['_main', '_malloc']\" -s \"EXPORTED_RUNTIME_METHODS=['ccall']\"")
endif ()

add_subdirectory(lib/raylib-4.2.0)
add_subdirectory(lib/raylib-cpp-4.2.7)

set(LIB_INCLUDES
        lib/raygui-3.2/include
        lib/thread-pool-3.3.0/include
        )

set(SOURCE_FILES
        src/main.cpp
        src/util/fixed_loop.cpp
        )

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_link_libraries(${PROJECT_NAME} -static-libgcc -static-libstdc++)
endif ()

target_include_directories(${PROJECT_NAME} PRIVATE ${LIB_INCLUDES})

target_link_libraries(${PROJECT_NAME} raylib raylib_cpp)