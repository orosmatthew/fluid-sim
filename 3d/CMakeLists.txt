cmake_minimum_required(VERSION 3.24)

project(fluid_sim_3d)

set(CMAKE_CXX_STANDARD 20)

function(add_shaders TARGET)
    find_program(GLSLANGVALIDATOR glslangvalidator)
    foreach (SHADER ${ARGN})
        set(current-shader-path ${CMAKE_CURRENT_SOURCE_DIR}/src/shader/${SHADER})
        set(current-output-path ${CMAKE_CURRENT_SOURCE_DIR}/res/bin/shader/${SHADER}.spv)

        # Add a custom command to compile GLSL to SPIR-V using glslangvalidator.
        get_filename_component(current-output-dir ${current-output-path} DIRECTORY)
        file(MAKE_DIRECTORY ${current-output-dir})

        add_custom_command(
                OUTPUT ${current-output-path}
                COMMAND ${GLSLANGVALIDATOR} -V -o ${current-output-path} ${current-shader-path}
                DEPENDS ${current-shader-path}
                IMPLICIT_DEPENDS CXX ${current-shader-path}
                VERBATIM)

        # Make sure our build depends on this output.
        set_source_files_properties(${current-output-path} PROPERTIES GENERATED TRUE)
        target_sources(${TARGET} PRIVATE ${current-output-path})
    endforeach ()
endfunction()


find_package(Vulkan REQUIRED)

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

set(FT_DISABLE_ZLIB ON CACHE BOOL "" FORCE)
set(FT_DISABLE_BZIP2 ON CACHE BOOL "" FORCE)
set(FT_DISABLE_PNG ON CACHE BOOL "" FORCE)
set(FT_DISABLE_HARFBUZZ ON CACHE BOOL "" FORCE)
set(FT_DISABLE_BROTLI ON CACHE BOOL "" FORCE)

add_subdirectory(lib/glfw-3.3.8)
add_subdirectory(lib/freetype-2.13.0)

set(LIB_SOURCE_FILES
        lib/spirv-reflect-1.3.236/src/spirv_reflect.cpp)

set(SOURCE_FILES
        src/main.cpp
        src/mve/monitor.cpp
        src/mve/renderer.cpp
        src/mve/shader.cpp
        src/mve/vertex_data.cpp
        src/mve/window.cpp
        src/logger.cpp
        )

set(LIBS
        Vulkan::Vulkan
        glfw
        freetype)

set(LIB_INCLUDES
        lib/spdlog-1.11.0/include
        lib/vma-3.0.1/include
        lib/stb-8b5f1f3/include
        lib/tinyobjloader-45f683c/include
        lib/spirv-reflect-1.3.236/include
        lib/thread-pool-3.3.0/include)

add_executable(${PROJECT_NAME})

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_link_libraries(${PROJECT_NAME} -static-libgcc -static-libstdc++)
endif ()

target_sources(${PROJECT_NAME} PRIVATE
        ${LIB_SOURCE_FILES}
        ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} ${LIBS})

target_include_directories(${PROJECT_NAME} PRIVATE ${LIB_INCLUDES})

add_shaders(fluid_sim_3d
        cloud.frag
        cloud.vert)

