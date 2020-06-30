# see Config file for compile-time settings
#include Config

include ../Make_include

CFLAGS += -std=c99 -I../include

MAKEFLAGS += --no-print-directory

# macOS compiler options (uncomment on macOS) ---------------------------------

# SOFLAGS += -undefined dynamic_lookup

# default targets
BUILD_CFP = 0
BUILD_ZFORP = 0
BUILD_UTILITIES = 0
BUILD_EXAMPLES = 0
BUILD_TESTING = 0
BUILD_SHARED_LIBS = 0

LIBRARY = static
LIBZFP = libzfp.a

# compiler options ------------------------------------------------------------

# default: build all targets enabled in Config
all:
	@echo $(LIBRARY)
	@cd src; $(MAKE) $(LIBRARY)
ifneq ($(BUILD_CFP),0)
	@cd cfp/src; $(MAKE) clean $(LIBRARY)
endif
ifneq ($(BUILD_ZFORP),0)
	@cd fortran; $(MAKE) clean $(LIBRARY)
endif
ifneq ($(BUILD_UTILITIES),0)
	@cd utils; $(MAKE) clean all
endif
ifneq ($(BUILD_TESTING),0)
	@cd tests; $(MAKE) clean all
endif
ifneq ($(BUILD_EXAMPLES),0)
	@cd examples; $(MAKE) clean all
endif

install: all
	-cp -rp include/* $I
	-cp -rp lib/* $L

# run basic regression tests
test:
	@cd tests; $(MAKE) test

# clean all
clean:
	@cd src; $(MAKE) clean
ifneq ($(BUILD_CFP),0)
	@cd cfp/src; $(MAKE) clean
endif
ifneq ($(BUILD_ZFORP),0)
	@cd fortran; $(MAKE) clean
endif
ifneq ($(BUILD_UTILITIES),0)
	@cd utils; $(MAKE) clean
endif
ifneq ($(BUILD_TESTING),0)
	@cd tests; $(MAKE) clean
endif
ifneq ($(BUILD_EXAMPLES),0)
	@cd examples; $(MAKE) clean
endif
