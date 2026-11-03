# Makefile fragment for compiling the ADCS library into another Make-based firmware build (PVDXos).
# `include` it, then add the variables below to your own source/include lists:
#   ADCS_SRCS         - .c files to compile: the ADCS sources plus the CMSIS-DSP sources they use
#   ADCS_INCLUDE_DIRS - directories to pass with -I
#   ADCS_CFLAGS       - flags the sources need
# All paths are absolute, so it doesn't matter where the including Makefile runs from.
# The host test build (CMakeLists.txt) doesn't use this file.

ADCS_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))

ifeq ($(wildcard $(ADCS_DIR)/linalg/Include/arm_math.h),)
    $(error CMSIS-DSP not found in $(ADCS_DIR)/linalg. Run 'git submodule update --init --recursive')
endif

# main.c and test.c are the host test runner, not part of the library
ADCS_LIB_SRCS := $(filter-out $(ADCS_DIR)/src/main.c $(ADCS_DIR)/src/test.c,$(wildcard $(ADCS_DIR)/src/*.c))

# CMSIS-DSP: cmsis_dsp/ has one wrapper per Source/<Group>/<Group>.c amalgamation (the same files
# CMakeLists.txt compiles, minus F16, which the Cortex-M4 has no FPU for). The wrappers silence
# -Wfloat-equal, which CMSIS trips. linalg/ itself can't be globbed: it also holds every arm_*.c the
# amalgamations #include (duplicate symbols), plus examples and tests. Unused functions are dropped
# by -ffunction-sections and --gc-sections.
ADCS_CMSIS_SRCS := $(wildcard $(ADCS_DIR)/cmsis_dsp/*.c)

ADCS_SRCS := $(ADCS_LIB_SRCS) $(ADCS_CMSIS_SRCS)

# The sources include "include/foo.h" (relative to ADCS_DIR) and "arm_math.h" / "Include/dsp/foo.h"
# (relative to linalg/Include and linalg)
ADCS_INCLUDE_DIRS := $(ADCS_DIR) \
                     $(ADCS_DIR)/linalg \
                     $(ADCS_DIR)/linalg/Include \
                     $(ADCS_DIR)/linalg/PrivateInclude

ADCS_CFLAGS :=
