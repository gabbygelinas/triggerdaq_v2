# CMakeLists.txt for basic analysis package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(agana)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

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
  add_compile_definitions(HAVE_THTTP_SERVER)
  include_directories(${ROOT_INCLUDE_DIRS})
  # Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
  include(${ROOT_USE_FILE})
endif()
##########################################

# MIDAS
if(NO_MIDAS)
   message(STATUS "${PROJECT_NAME}: MIDAS support is disabled via NO_MIDAS")
   unset(MIDAS_FOUND)
elseif(DEFINED ENV{MIDASSYS})
   include($ENV{MIDASSYS}/lib/midas-targets.cmake)
   message(STATUS "${PROJECT_NAME}: Building with MIDAS from $ENV{MIDASSYS}/lib/midas-targets.cmake")
   set(MIDAS_FOUND TRUE)
   set(MIDASSYS "$ENV{MIDASSYS}")
   set(MIDAS_LIBRARIES midas)
   set(MIDAS_LIBRARY_DIRS "${MIDASSYS}/lib")
   message(STATUS "MIDAS found ${MIDASSYS}")
else()
   message(STATUS "${PROJECT_NAME}: Building without MIDAS")
   unset(MIDAS_FOUND)
endif()


##########################################
# determine whether agana is a submodule or standalone
if( _AGANA_STATIC )
  message(STATUS "building agana static library")
  add_library(agana STATIC ${Sources})
  
  if( BUILD_MANALYZER )
    target_include_directories(agana PUBLIC ${CMAKE_SOURCE_DIR}/manalyzer)
    target_include_directories(agana PUBLIC ${CMAKE_SOURCE_DIR}/midasio)
    target_include_directories(agana PUBLIC ${CMAKE_SOURCE_DIR}/mvodb)
  elseif( MIDAS_FOUND )
    target_include_directories(agana PUBLIC ${MIDAS_INCLUDE_DIRS})
    target_include_directories(agana PUBLIC $ENV{MIDASSYS}/manalyzer)
  else()
    message(FATAL_ERROR "Failed configuration for target libagana")
  endif()

  file(GLOB HEADER_FILES include/*.h)
  install(FILES ${HEADER_FILES} DESTINATION include)

else()
  message(STATUS "building agana standalone")
  set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

  add_executable(agana.exe ${Sources} ${Headers})

  if(MIDAS_FOUND)
    message(STATUS "MIDAS found in AGANA standalone")
   
    #set(manalyzer_objs)
    target_include_directories(agana.exe PUBLIC ${MIDAS_INCLUDE_DIRS})
    target_include_directories(agana.exe PUBLIC $ENV{MIDASSYS}/manalyzer)
    target_compile_definitions(agana.exe PUBLIC -DHAVE_MIDAS -DHAVE_TMFE)
    target_link_directories(agana.exe PUBLIC ${MIDAS_LIBRARY_DIRS})
  else()
    # get midas-analyzer submodules
    get_submodule(midasio)
    get_submodule(mvodb)
    get_submodule(mjson)
    get_submodule(mxml)
    get_submodule(manalyzer)
    set(manalyzer_objs
      ${PROJECT_SOURCE_DIR}/midasio/midasio.cxx
      ${PROJECT_SOURCE_DIR}/midasio/lz4frame.c
      ${PROJECT_SOURCE_DIR}/midasio/xxhash.c
      ${PROJECT_SOURCE_DIR}/mvodb/mvodb.cxx
      ${PROJECT_SOURCE_DIR}/mvodb/nullodb.cxx
      ${PROJECT_SOURCE_DIR}/mvodb/mjsonodb.cxx
      ${PROJECT_SOURCE_DIR}/mvodb/mxmlodb.cxx
      ${PROJECT_SOURCE_DIR}/mjson/mjson.cxx
      ${PROJECT_SOURCE_DIR}/mxml/mxml.cxx
      )
    add_library(manalyzer      STATIC ${PROJECT_SOURCE_DIR}/manalyzer/manalyzer ${manalyzer_objs})
    add_library(manalyzer_main STATIC ${PROJECT_SOURCE_DIR}/manalyzer/manalyzer_main)
    
    target_include_directories(manalyzer      PUBLIC ${PROJECT_SOURCE_DIR})
    target_include_directories(manalyzer_main PUBLIC ${PROJECT_SOURCE_DIR})
    
    target_compile_features(manalyzer      PUBLIC cxx_std_11)
    target_compile_features(manalyzer_main PUBLIC cxx_std_11)
    
    target_include_directories(manalyzer      PUBLIC ${PROJECT_SOURCE_DIR}/midasio ${PROJECT_SOURCE_DIR}/mvodb ${PROJECT_SOURCE_DIR}/mjson ${PROJECT_SOURCE_DIR}/mxml ${PROJECT_SOURCE_DIR}/manalyzer)
    target_include_directories(manalyzer_main PUBLIC ${PROJECT_SOURCE_DIR}/midasio ${PROJECT_SOURCE_DIR}/mvodb ${PROJECT_SOURCE_DIR}/manalyzer)
 

    target_include_directories(manalyzer      PUBLIC ${ROOT_INCLUDE_DIRS})
    target_include_directories(manalyzer_main PUBLIC ${ROOT_INCLUDE_DIRS})
    target_compile_definitions(manalyzer      PUBLIC ${ROOT_CCX_FLAGS} -DHAVE_ROOT -DHAVE_ROOT_HTTP -DHAVE_THTTP_SERVER)
    target_compile_definitions(manalyzer_main PUBLIC ${ROOT_CCX_FLAGS} -DHAVE_ROOT -DHAVE_ROOT_HTTP -DHAVE_THTTP_SERVER)
  endif()

   target_link_libraries(agana.exe ${ROOT_LIBRARIES} manalyzer_main manalyzer ${MIDAS_LIBRARIES})
   
  install(
    TARGETS agana.exe
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    )
endif()
