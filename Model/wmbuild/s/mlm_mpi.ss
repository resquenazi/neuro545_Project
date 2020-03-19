#!/bin/bash
mpicc -c ./libc/$1.c
ar rv ./lib/lib$1.a $1.o
rm -f ./$1.o
