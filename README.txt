Mitris Victor
334 CC ASC - dtrmv

Implementation:
	function dtrmv - the one from fortram
	function dtrmv2 - optimized function. 

Results for opteron, nehalem or quad:
	For the one created by me - mana.out
	FOr the blas one - blas.out
	For the optimized one - optimi.out

Execution:
	To execute run make.
	Make executes run.sh which should not be interrupted because it verifies the numeber of raws returned by qstat for grafic construction.
	run.sh executes exec_script.sh for opteron, exec_script1.sh for nehalem and exec_Script2.sh for quad.
	The exec_script.sh executes the build1 targets to copile for opteron, build 2 to compile for nehalem and build3 to compile for quad.
	After that the program is executed using as parameter a variable length from 1000 to 25000 with 2500 incrmentation.
	Scriptul run.sh apeleaza exec_script.sh pentru opteron, exec_script1.sh pentru nehalem si exec_Script2.sh pentru quad.
	The file.gp contains the commands for graphics.	









	
