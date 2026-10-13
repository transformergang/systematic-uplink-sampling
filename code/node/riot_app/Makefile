# RIOT application for BOARD=iotlab-m3 (Cortex-M3) or BOARD=samr21-xpro (Cortex-M0+).
# Build:  make RIOTBASE=/path/to/RIOT BOARD=iotlab-m3
# The portable sources (../alg) and the test vectors (../vectors) are copied next to main.c at build time.
APPLICATION = madow_node_bench
BOARD ?= iotlab-m3
RIOTBASE ?= $(HOME)/RIOT
$(shell cp -f $(CURDIR)/../alg/*.c $(CURDIR)/../alg/*.h $(CURDIR)/../vectors/vectors.h $(CURDIR)/../vectors/vectors_madow.h $(CURDIR)/)
SRC = main.c indep_node.c excl_node.c baselines.c extp_node.c madow_node.c
ifneq (,$(wildcard $(RIOTBASE)/sys/include/ztimer.h))
  USEMODULE += ztimer_msec
  CFLAGS += -DPAM_USE_ZTIMER=1
else
  USEMODULE += xtimer
endif
DEVHELP ?= 1
CFLAGS += -O2 -DDEVELHELP
QUIET ?= 1
include $(RIOTBASE)/Makefile.include
