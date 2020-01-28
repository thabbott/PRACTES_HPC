SHELL=/bin/bash

#=====================
# User-defined options
#=====================

# Order of polynomial reconstruction within each element
# 1 = linear, 2 = quadratic, 3 = cubic, etc
DG_N=1

#=========================
# End user-defined options
#=========================

# C compiler
CC=mpicc

# dependencies
DEPS=$(shell ls *.h)

CPPFLAGS=-DDG_N=${DG_N} -O3 --std=gnu99
LDFLAGS=-lpng -lm

loon: loon.c ${DEPS}
	${CC} -o $@ ${CPPFLAGS} loon.c ${LDFLAGS}

vortices: vortices.c ${DEPS}
	${CC} -o $@ ${CPPFLAGS} vortices.c ${LDFLAGS}

create_input: create_input.c ${DEPS}
	${CC} -o $@ ${CPPFLAGS} create_input.c ${LDFLAGS}
	./create_input

clean:
	-rm create_input
	-rm loon
	-rm vortices
	-rm output/*.png
