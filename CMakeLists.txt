# CMakeLists.txt for basic analysis package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(dlana)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_compile_options("-O2")
add_compile_options("-g")
add_compile_options(-Wall -Wformat=2 -Wno-format-nonliteral -Wno-strict-aliasing -Wuninitialized -Wno-unused-function)

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

function(clone_module MODULE_NAME MODULE_URL)
   if(GIT_SUBMODULE)
      if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE_NAME}")
         message(STATUS "${MODULE_NAME} submodule ok")
	 execute_process(COMMAND ${GIT_EXECUTABLE} pull
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE_NAME}
           RESULT_VARIABLE GIT_SUBMOD_UPDATE)
         if(NOT GIT_SUBMOD_UPDATE EQUAL "0")
            message(FATAL_ERROR "git clone ${MODULE_NAME} ${GIT_SUBMOD_RESULT}, please check submodules")
         endif()
      else()
         message(STATUS "Cloning Submodule ${MODULE_NAME}")
         execute_process(COMMAND ${GIT_EXECUTABLE} clone ${MODULE_URL}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
         if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git clone ${MODULE_URL} FAILED with ${GIT_SUBMOD_RESULT}, please checkout submodules")
         endif()
      endif()
    endif()
endfunction()


# electronics database
get_submodule( dlcfmdb )
##########################################

##########################################
# MIDAS and manalyzer
##########################################

if(NO_MIDAS)
   message(STATUS "${PROJECT_NAME}: MIDAS support is disabled via NO_MIDAS")
   unset(MIDAS_FOUND)
elseif(DEFINED ENV{MIDASSYS} AND (NOT ("$ENV{MIDASSYS}" STREQUAL "")))
   include($ENV{MIDASSYS}/lib/midas-targets.cmake)
   message(STATUS "${PROJECT_NAME}: Building with MIDAS from $ENV{MIDASSYS}/lib/midas-targets.cmake")
   set(MIDAS_FOUND TRUE)
   include($ENV{MIDASSYS}/lib/manalyzer-targets.cmake)
   message(STATUS "${PROJECT_NAME}: Building with manalyzer from $ENV{MIDASSYS}/lib/manalyzer-targets.cmake")
else()
   message(STATUS "${PROJECT_NAME}: Building without MIDAS")
   unset(MIDAS_FOUND)
endif()

##########################################
message(STATUS "building dlana static library")

# find sources for the dlana library
#file(GLOB Headers ${PROJECT_SOURCE_DIR}/*.h)
file(GLOB DlanaSources ${PROJECT_SOURCE_DIR}/*.cxx)
##########################################

list(FILTER DlanaSources EXCLUDE REGEX ".*_module.cxx$")
message(STATUS "from ${DlanaSources}")
add_library(dlana STATIC ${DlanaSources})
target_include_directories(dlana PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  
if( BUILD_MANALYZER )
  target_include_directories(dlana PUBLIC ${CMAKE_SOURCE_DIR}/manalyzer)
  target_include_directories(dlana PUBLIC ${CMAKE_SOURCE_DIR}/midasio)
  target_include_directories(dlana PUBLIC ${CMAKE_SOURCE_DIR}/mvodb)
elseif( MIDAS_FOUND )
  #target_include_directories(dlana PUBLIC ${MIDASSYS}/manalyzer)
  #target_include_directories(dlana PUBLIC ${MIDASSYS}/midasio)
  #target_include_directories(dlana PUBLIC ${MIDASSYS}/mvodb)
  target_link_libraries(dlana PRIVATE midas)
else()
  target_include_directories(dlana PUBLIC ${CMAKE_SOURCE_DIR}/midasio)
endif()

file(GLOB HEADER_FILES ${PROJECT_SOURCE_DIR}/*.h)
install(FILES ${HEADER_FILES} DESTINATION include)

if( _DLANA_STATIC )
else()
  message(STATUS "building dlana.exe")
  set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

  set(DlanaModules unpack_cb_module.cxx cbko_module.cxx coinc_module.cxx dltdc_module.cxx dltdc8_module.cxx dltdc4_module.cxx)
  #set(DlanaModules unpack_cb_module.cxx cbko_module.cxx dltdc_module.cxx dltdc8_module.cxx dltdc4_module.cxx)
  #set(TestunpackModules testunpack_module)

  add_executable(dlana.exe ${DlanaModules})
  #add_executable(testunpack.exe ${TestunpackModules})

  if(MIDAS_FOUND)
    message(STATUS "using manalyzer from MIDAS")
    target_link_libraries(dlana.exe PRIVATE manalyzer_main dlana manalyzer midas)
    #target_link_libraries(testunpack.exe PRIVATE manalyzer_main dlana manalyzer midas)
  else()
    message(STATUS "building our own manalyzer")
    # get midas-analyzer submodules
    #get_submodule(midasio)
    #get_submodule(mvodb)
    #get_submodule(mjson)
    #get_submodule(mxml)
    #get_submodule(manalyzer)
    clone_module(midasio https://bitbucket.org/tmidas/midasio.git)
    clone_module(mvodb https://bitbucket.org/tmidas/mvodb.git)
    clone_module(mjson https://bitbucket.org/tmidas/mjson.git)
    clone_module(mxml https://bitbucket.org/tmidas/mxml.git)
    clone_module(manalyzer https://bitbucket.org/tmidas/manalyzer.git)
    add_subdirectory(manalyzer)
    target_link_libraries(dlana.exe PRIVATE manalyzer_main dlana manalyzer ${ROOT_LIBRARIES})
    #target_link_libraries(testunpack.exe PRIVATE manalyzer_main dlana manalyzer ${ROOT_LIBRARIES})
  endif()
   
  install(
    TARGETS dlana.exe
    #TARGETS testunpack.exe #everything testunpack was commented out by Gabby
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    )
endif()
