# CMakeLists.txt for basic analysis package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(agana)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
if(EXISTS $ENV{ROOTSYS}/cmake/ROOTConfig.cmake)
    if(NOT DEFINED ENV{CMAKE_PREFIX_PATH})
        list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
    endif()
endif()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{ROOTSYS}/cmake ")

# Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS)
# specify the location of the ROOT libraries
link_libraries(${ROOT_LIBRARIES})
# Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})


find_package(Git REQUIRED)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
   option(GIT_SUBMODULE "Check submodules during build" ON)
endif()

function(get_submodule MODULE_NAME)
   if(GIT_SUBMODULE)
      if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE_NAME}/.git")
         message("${MODULE_NAME} submodule ok")
      else()
         message(STATUS "Submodule ${MODULE_NAME} update --init")
         execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init ${MODULE_NAME}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
         if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
         endif()
      endif()
    endif()
endfunction()


# You need to tell CMake where to find the ROOTANA installation. This can be done in a number of ways, this is a painful one.
if(EXISTS $ENV{ROOTANASYS}/rootana-config.cmake)
    if(NOT DEFINED $ENV{CMAKE_PREFIX_PATH})
        list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTANASYS})
    endif()
    message(STATUS $ENV{ROOTANASYS}/rootana-config.cmake " Found")
endif()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{ROOTANASYS} ")
find_package(ROOTANA)

include_directories(${ROOTANA_INCLUDE_DIRS})
# specify the location of the ROOTANA library
link_directories($ENV{ROOTANASYS}/lib)

# determine whether ROOTANA is a submodule
string(COMPARE EQUAL "$ENV{ROOTANASYS}" "${CMAKE_SOURCE_DIR}/rootana" _rootana_default )

SET(_AGANA_LIBS rootana z m dl pthread Thread MultiProc Gui XMLParser XMLIO RHTTP)

# MIDAS -- live analyzer
set(MIDASSYS "" CACHE PATH "")
if(MIDASSYS STREQUAL "")
    set(MIDASSYS "$ENV{MIDASSYS}")
endif()
if(MIDASSYS STREQUAL "")
    find_package(Midas QUIET)
else()
    set(MIDAS_FOUND TRUE)
    set(MIDAS_INCLUDE_DIRS "$ENV{MIDASSYS}/include")
    set(MIDAS_LIBRARIES "$ENV{MIDASSYS}/lib/libmidas.a" -lutil -lrt)
endif()
if(MIDAS_FOUND)
    message(DEBUG "[${PROJECT_NAME}] MIDAS_INCLUDE_DIRS = ${MIDAS_INCLUDE_DIRS}")
    message(DEBUG "[${PROJECT_NAME}] MIDAS_LIBRARIES = ${MIDAS_LIBRARIES}")
    list(APPEND _AGANA_LIBS midas)
    link_libraries(${MIDAS_LIBRARIES})
    link_directories($ENV{MIDASSYS}/lib)
endif()

file(GLOB Sources "*.cxx")

get_submodule( agcfmdb )

# determine whether agana is a submodule or standalone
if( _AGANA_STATIC )

  message(STATUS "building agana static library")
  add_library(agana STATIC ${Sources})
  if(_rootana_default)
    add_dependencies(agana rootana)
  endif()

file(GLOB HEADER_FILES include/*.h)
install(FILES ${HEADER_FILES} DESTINATION include)

else()
  message(STATUS "building agana standalone")
  set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

  add_executable(agana.exe $ENV{ROOTANASYS}/manalyzer/manalyzer_main ${Sources} )
  target_link_libraries(agana.exe ${_AGANA_LIBS})

  install(
    TARGETS agana.exe
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    )
endif()
