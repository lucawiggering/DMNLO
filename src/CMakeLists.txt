# -- DM@NLO library ---
add_library(libdmnlo STATIC DMNLO.cc params.cc)
SET_TARGET_PROPERTIES(libdmnlo PROPERTIES PREFIX "")
TARGET_LINK_LIBRARIES(libdmnlo util util2 ChiChi2QQ NeuQ2qx NeuQ2qx_M2Amp2V NeuQ2qx_M2AmpBox2V QQ2xx staustau2QQ stst2QQ stsT2xx DD ${LOOPTOOLS_LIBRARIES} ${SLHA_LIBRARIES})

# -- DM@NLO executable ---
add_executable(dmnlo main.cc)
TARGET_LINK_LIBRARIES(dmnlo libdmnlo)

install(TARGETS dmnlo RUNTIME DESTINATION .)
install(TARGETS libdmnlo ARCHIVE DESTINATION .)