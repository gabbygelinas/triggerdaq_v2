# CMakeLists.txt for basic analysis package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(agana)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# ROOT
if(EXISTS $ENV{ROOTSYS}/cmake/ROOTConfig.cmake)
    if(NOT DEFINED $ENV{CMAKE_PREFIX_PATH})
        list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
    endif()
endif()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{ROOTSYS}/cmake ")
# Locate the ROOT package and defines a number of variables 
find_package(ROOT REQUIRED COMPONENTS XMLIO XMLParser RHTTP Gui)
if(ROOT_FOUND)
  add_compile_definitions(HAVE_ROOT)
  include_directories(${ROOT_INCLUDE_DIRS})
  # Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
  include(${ROOT_USE_FILE})
endif()

# ROOTANA
if(EXISTS $ENV{ROOTANASYS}/rootana-config.cmake)
    if(NOT DEFINED $ENV{CMAKE_PREFIX_PATH})
        list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTANASYS})
    endif()
    message(STATUS $ENV{ROOTANASYS}/rootana-config.cmake " Found")
endif()
# Locate the ROOTANA package and defines a number of variables 
find_package(ROOTANA REQUIRED)
if(ROOTANA_FOUND)
  message(STATUS "${ROOTANA_PACKAGE_NAME} ok $ENV{ROOTANASYS}")
  set(ROOTANA_LIBRARY_DIRS "$ENV{ROOTANASYS}/lib")
  message(STATUS "ROOTANA header files: ${ROOTANA_INCLUDE_DIRS}")
else()
  message(FATAL "ROOTANA not found")
endif(ROOTANA_FOUND)
# specify the location of the ROOTANA headers
include_directories(${ROOTANA_INCLUDE_DIRS})
# specify the location of the ROOTANA library
link_directories(${ROOTANA_LIBRARY_DIRS})

##########################################
# determine whether ROOTANA is a submodule
string(COMPARE EQUAL "$ENV{ROOTANASYS}" "${CMAKE_SOURCE_DIR}/rootana" _rootana_default )
##########################################

# MIDAS -- live analyzer
set(MIDASSYS "" CACHE PATH "")
if(MIDASSYS STREQUAL "")
    set(MIDASSYS "$ENV{MIDASSYS}")
endif()
list(APPEND CMAKE_PREFIX_PATH $ENV{MIDASSYS})
# Locate the MIDAS package
find_package(Midas QUIET)
if(MIDAS_FOUND)
  set(MIDASSYS "$ENV{MIDASSYS}")
  message(STATUS "MIDAS found ${MIDASSYS}")
  message(STATUS "[${PROJECT_NAME}] MIDAS_INCLUDE_DIRS = ${MIDAS_INCLUDE_DIRS}")
  message(STATUS "[${PROJECT_NAME}] MIDAS_LIBRARIES = ${MIDAS_LIBRARIES}")
  # specify the location of the MIDAS headers
  include_directories(${MIDAS_INCLUDE_DIRS}) 
  # specify the location of the MIDAS library
  link_directories(${MIDAS_LIBRARY_DIRS})
endif()


##########################################
find_package(Git REQUIRED)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
   option(GIT_SUBMODULE "Check submodules during build" ON)
endif()

function(get_submodule MODULE_NAME)
   if(GIT_SUBMODULE)
      if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE_NAME}/.git")
         message(STATUS "${MODULE_NAME} submodule ok")
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

# electronics database
get_submodule( agcfmdb )
##########################################

# find source for this project
file(GLOB Headers ${PROJECT_SOURCE_DIR}/*.h)
file(GLOB Sources ${PROJECT_SOURCE_DIR}/*.cxx)

##########################################
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

  add_executable(agana.exe "$ENV{ROOTANASYS}/manalyzer/manalyzer_main.cxx" "$ENV{ROOTANASYS}/manalyzer/manalyzer.cxx" ${Sources} ${Headers})
  target_link_libraries(agana.exe ${ROOT_LIBRARIES} rootana z ${MIDAS_LIBRARIES} util rt) 

  install(
    TARGETS agana.exe
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    )
endif()
