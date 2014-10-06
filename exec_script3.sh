#!/bin/bash

[[ -z $COMPILER ]] && COMPILER="gcc"

if [[ $COMPILER="gcc" ]]; then
   make -f Makefile build3;
   for i in {1000..25000..2500}
   do
   	./dtrmv 3 $i
   done
fi
