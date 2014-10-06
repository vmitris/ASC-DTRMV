#!/bin/bash

[[ -z $COMPILER ]] && COMPILER="gcc"

if [[ $COMPILER="gcc" ]]; then
	make -f Makefile build1;
	for i in {1000..25000..2500}
	do
   		./dtrmv 1 $i
   	done
fi
