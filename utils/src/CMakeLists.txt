set(SRCFILES BinManager.cc EigenDecomp.cc OptParser.cc ToyThrower.cc XsecDial.cc XsecNorm.cc)
set(HEADERS BinManager.hh EigenDecomp.hh OptParser.hh ToyThrower.hh XsecDial.hh XsecNorm.hh)

if(CMAKE_BUILD_TYPE MATCHES DEBUG)
  add_library(xslfutils STATIC ${SRCFILES})
else(CMAKE_BUILD_TYPE MATCHES RELEASE)
  add_library(xslfutils SHARED ${SRCFILES})
endif()

set_target_properties(xslfutils PROPERTIES VERSION
  "${xsLLhFitter_VERSION_MAJOR}.${xsLLhFitter_VERSION_MINOR}.${xsLLhFitter_VERSION_REVISION}")
set_target_properties(xslfutils PROPERTIES LINK_FLAGS ${ROOT_LD_FLAGS})
target_include_directories(xslfutils PUBLIC ${CMAKE_SOURCE_DIR}/include)

install(TARGETS xslfutils DESTINATION lib)
#Can uncomment this to install the headers... but is it really neccessary?
# install(FILES ${HEADERS} DESTINATION include)

set(MODULETargets ${MODULETargets} xslfutils PARENT_SCOPE)
