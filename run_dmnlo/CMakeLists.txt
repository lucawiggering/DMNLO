set(CMAKE_Fortran_FLAGS  "-ffixed-line-length-none -std=legacy")
include_directories(${CMAKE_SOURCE_DIR}/external/lib/SLHALib-2.2/src)
include_directories("util/cuba/common")


# --- util ---
file(GLOB_RECURSE loop_SRC "util/loop/*.F")
file(GLOB_RECURSE mint_SRC "util/mint/*.F")
file(GLOB_RECURSE misc_SRC "util/misc/*.F")
file(GLOB_RECURSE mssm_SRC "util/mssm/*.F")
file(GLOB_RECURSE resum_SRC "util/resum/*.F")
file(GLOB_RECURSE init_SRC "util/*.F")

add_library(util OBJECT "util/cuba/vegas/Vegas.c" "util/cuba/vegas/vegas-f.c" ${loop_SRC} ${mint_SRC} ${misc_SRC} ${mssm_SRC} ${resum_SRC} ${init_SRC})

# --- util2 ---
file(GLOB_RECURSE util2_SRC "util2/*.F")
add_library(util2 OBJECT ${util2_SRC})

# --- ChiChi2QQ ---
add_library(ChiChi2QQ OBJECT "ChiChi2QQ/ChiChi2QQ_Main.F")

# --- NeuQ2qx ---
add_library(NeuQ2qx_M2Amp2V OBJECT "NeuQ2qx/NeuQ2qx_M2Amp2V.F")
add_library(NeuQ2qx_M2AmpBox2V OBJECT "NeuQ2qx/NeuQ2qx_M2AmpBox2V.F")
add_library(NeuQ2qx OBJECT "NeuQ2qx/NeuQ2qx_Main.F")

# --- QQ2xx ---
add_library(QQ2xx OBJECT "QQ2xx/QQ2xx_Main.F")

# --- staustau2QQ ---
add_library(staustau2QQ OBJECT "staustau2QQ/staustau2QQ_Main.F")

# --- stst2QQ ---
add_library(stst2QQ OBJECT "stst2QQ/stst2QQ_Main.F")

# --- stsT2xx ---
add_library(stsT2xx OBJECT "stsT2xx/stsT2xx_Main.F")

# --- stsT2xx ---
add_library(DD OBJECT "DD/DD_Main.F")
