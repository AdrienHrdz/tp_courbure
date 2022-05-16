cmake_minimum_required(VERSION 3.7.0)

project(pgm LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)


find_package(Qt5 COMPONENTS Widgets OpenGL REQUIRED)
find_package(OpenGL REQUIRED)

find_package(GLEW)
if(NOT ${GLEW_FOUND} OR APPLE)

  message(WARNING "Glew not found, so manual compilation of glew added to the project")
  
  option(glew-cmake_BUILD_SHARED "Build the shared glew library" OFF)
  option(ONLY_LIBS "Do not build executables" ON)

  add_subdirectory(extern/glew)
  include_directories(${PROJECT_PATH}/extern/glew/include)
  set(USED_GLEW_LIBRARY libglew_static)
else()
  include_directories(${GLEW_INCLUDE_DIRS})
  set(USED_GLEW_LIBRARY ${GLEW_LIBRARIES})
endif()

file(
GLOB_RECURSE
source_files
"src/*.cpp"
)

file(
GLOB_RECURSE
header_files
"src/*.h"
"src/*.hpp"
)

file(
GLOB_RECURSE
shader_files
"shaders/*.frag"
"shaders/*.vert"
)



message("Source files")
foreach(line IN LISTS source_files)
    message("${line}")
endforeach() 

message("header files")
foreach(line IN LISTS header_files)
    message("${line}")
endforeach() 

message("Shader files")
foreach(line IN LISTS shader_files)
    message("${line}")
endforeach() 

if (WIN32)
  set(LINK_FLAGS_PERSO )
elseif (APPLE)
  set(LINK_FLAGS_PERSO )
else ()
  set(LINK_FLAGS_PERSO -lm -ldl -fopenmp)
  add_definitions( -Wall -Wextra -Wno-comment -Wno-unused-parameter -Wno-unused-function -Wno-unused-variable -fopenmp  -DCORRECTION)
endif (WIN32)

add_executable( 
  pgm
  ${source_files}
  ${header_files}
  ${shader_files}
)

target_link_libraries(pgm ${USED_GLEW_LIBRARY} ${OPENGL_LIBRARIES} Qt5::Widgets Qt5::OpenGL ${LINK_FLAGS_PERSO}	)
