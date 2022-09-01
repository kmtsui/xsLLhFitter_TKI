#xsllhCalcXsec
add_executable(xsllhCalcXsec xsllhCalcXsec.cxx)

set_target_properties(xsllhCalcXsec PROPERTIES LINK_FLAGS ${ROOT_LD_FLAGS})

foreach(mod ${MODULES})
 target_include_directories(xsllhCalcXsec PUBLIC ${CMAKE_SOURCE_DIR}/${mod}/src)
 target_include_directories(xsllhCalcXsec PUBLIC ${CMAKE_SOURCE_DIR}/include)
endforeach()

foreach(mod ${MODULETargets})
  cmessage(STATUS "Linking with: ${mod}")
  target_link_libraries(xsllhCalcXsec ${mod})
endforeach()

#xsllhDetVariations
add_executable(xsllhDetVar xsllhDetVariations.cxx)

set_target_properties(xsllhDetVar PROPERTIES LINK_FLAGS ${ROOT_LD_FLAGS})

foreach(mod ${MODULES})
 target_include_directories(xsllhDetVar PUBLIC ${CMAKE_SOURCE_DIR}/${mod}/src)
 target_include_directories(xsllhDetVar PUBLIC ${CMAKE_SOURCE_DIR}/include)
endforeach()

foreach(mod ${MODULETargets})
  cmessage(STATUS "Linking with: ${mod}")
  target_link_libraries(xsllhDetVar ${mod})
endforeach()

#xsllhDetVariations_test
add_executable(xsllhDetVar_test xsllhDetVariations_test.cxx)

set_target_properties(xsllhDetVar_test PROPERTIES LINK_FLAGS ${ROOT_LD_FLAGS})

foreach(mod ${MODULES})
 target_include_directories(xsllhDetVar_test PUBLIC ${CMAKE_SOURCE_DIR}/${mod}/src)
 target_include_directories(xsllhDetVar_test PUBLIC ${CMAKE_SOURCE_DIR}/include)
endforeach()

foreach(mod ${MODULETargets})
  cmessage(STATUS "Linking with: ${mod}")
  target_link_libraries(xsllhDetVar_test ${mod})
endforeach()

#xsllhFit
add_executable(xsllhFit xsllhFit.cxx)

set_target_properties(xsllhFit PROPERTIES LINK_FLAGS ${ROOT_LD_FLAGS})

foreach(mod ${MODULES})
 target_include_directories(xsllhFit PUBLIC ${CMAKE_SOURCE_DIR}/${mod}/src)
 target_include_directories(xsllhFit PUBLIC ${CMAKE_SOURCE_DIR}/include)
endforeach()

foreach(mod ${MODULETargets})
  cmessage(STATUS "Linking with: ${mod}")
  target_link_libraries(xsllhFit ${mod})
endforeach()

#xsllhTreeConvert
add_executable(xsllhTreeConvert xsllhTreeConvert.cxx)

set_target_properties(xsllhTreeConvert PROPERTIES LINK_FLAGS ${ROOT_LD_FLAGS})

foreach(mod ${MODULES})
 target_include_directories(xsllhTreeConvert PUBLIC ${CMAKE_SOURCE_DIR}/${mod}/src)
 target_include_directories(xsllhTreeConvert PUBLIC ${CMAKE_SOURCE_DIR}/include)
endforeach()

foreach(mod ${MODULETargets})
  cmessage(STATUS "Linking with: ${mod}")
  target_link_libraries(xsllhTreeConvert ${mod})
endforeach()

#xsllhUnitTests
add_executable(xsllhUnitTest xsllhUnitTests.cxx)

set_target_properties(xsllhUnitTest PROPERTIES LINK_FLAGS ${ROOT_LD_FLAGS})

foreach(mod ${MODULES})
 target_include_directories(xsllhUnitTest PUBLIC ${CMAKE_SOURCE_DIR}/${mod}/src)
 target_include_directories(xsllhUnitTest PUBLIC ${CMAKE_SOURCE_DIR}/include)
endforeach()

foreach(mod ${MODULETargets})
  cmessage(STATUS "Linking with: ${mod}")
  target_link_libraries(xsllhUnitTest ${mod})
endforeach()

#xsllhXseCov
add_executable(xsllhXsecCov xsllhXsecCov.cxx)

set_target_properties(xsllhXsecCov PROPERTIES LINK_FLAGS ${ROOT_LD_FLAGS})

foreach(mod ${MODULES})
 target_include_directories(xsllhXsecCov PUBLIC ${CMAKE_SOURCE_DIR}/${mod}/src)
 target_include_directories(xsllhXsecCov PUBLIC ${CMAKE_SOURCE_DIR}/include)
endforeach()

foreach(mod ${MODULETargets})
  cmessage(STATUS "Linking with: ${mod}")
  target_link_libraries(xsllhXsecCov ${mod})
endforeach()

#xsllhOneSigma
add_executable(xsllhOneSigma xsllhOneSigma.cxx)

set_target_properties(xsllhOneSigma PROPERTIES LINK_FLAGS ${ROOT_LD_FLAGS})

foreach(mod ${MODULES})
 target_include_directories(xsllhOneSigma PUBLIC ${CMAKE_SOURCE_DIR}/${mod}/src)
 target_include_directories(xsllhOneSigma PUBLIC ${CMAKE_SOURCE_DIR}/include)
endforeach()

foreach(mod ${MODULETargets})
  cmessage(STATUS "Linking with: ${mod}")
  target_link_libraries(xsllhOneSigma ${mod})
endforeach()

#install targets
install(TARGETS xsllhCalcXsec xsllhFit xsllhDetVar xsllhDetVar_test xsllhTreeConvert xsllhUnitTest xsllhXsecCov xsllhOneSigma DESTINATION bin)
