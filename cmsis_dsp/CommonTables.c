// CMSIS-DSP's CommonTables group, compiled into PVDXos by adcs.mk. CMSIS compares floats with ==, and
// PVDXos builds with -Werror -Wfloat-equal, so that warning is off for this file only.
#pragma GCC diagnostic ignored "-Wfloat-equal"

// These tables total ~600 KB (FFT twiddles etc.). PVDXos doesn't build with -fdata-sections, so they'd
// all share one .rodata section, and using any one of them (e.g. sinTable_f32 for arm_sin_f32) would
// keep all of them and overflow RAM. CMSIS tags every table with ARM_DSP_TABLE_ATTRIBUTE, so give each
// one its own .rodata.* section; --gc-sections then drops the unused ones. (ELF only: this file isn't
// part of the host build, but editors may still parse it for a Mach-O target.)
#ifdef __ELF__
#define ADCS_STR_(x) #x
#define ADCS_STR(x) ADCS_STR_(x)
#define ARM_DSP_TABLE_ATTRIBUTE __attribute__((section(ADCS_STR(.rodata.cmsis_dsp_table.__COUNTER__))))
#endif

#include "Source/CommonTables/CommonTables.c"
