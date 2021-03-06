# current version of eiscor
LIBNAME := eiscor
MAJOR := 0
MINOR := 2
PATCH := 0
VERSION := $(MAJOR).$(MINOR).$(PATCH)

# change library extension based on OS
ifeq ($(findstring Windows_NT,$(OS)),Windows_NT)
	SLIB := dll
else ifeq ($(findstring darwin,$(OSTYPE)),darwin)
	SLIB := dylib
else
	SLIB := so
endif
