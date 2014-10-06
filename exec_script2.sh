#!/bin/bash

[[ -z $COMPILER ]] && COMPILER="gcc"

if [[ $COMPILER="gcc" ]]; then
   make -f Makefile build2;
   for i in {1000..25000..2500}
   do
      ./dtrmv 2 $i
   done
fi
