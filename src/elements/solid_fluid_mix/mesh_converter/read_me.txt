Mesh created by cubit mesh generator has its own connectivity assumption which is not the same as Tahoe connectivity assumption. A C++ code has been written(genesis_txt_to_geom.cxx) to create “.geom” input file for Tahoe with consistent connectivity assumption. The executable version of this file which is named ”converter” is used to do this conversion. “converter”, receives “input.txt” file(text file) and generates “output.geom” file.
In order to use this program, after creating mesh with cubit do the following steps:
1- export that mesh as a genesis(input.g) file.
2- create text file associated with genesis file using the following command:
>exotxt input.g input.txt
3- run “converter” program:
>converter
Note: “converter” and “input.txt” should be in the same directory.
4- If you want, you can rename “output.geom”. Then copy this geometry file to the geometry directory.

Note1: reading and writing “sidesets” are not implemented to the “converter” program. This part is ignored in conversion.

Note2: you can look at the genesis geometry file(input.g) in “paraview” to assign correct boundary  conditions to node sets.(in “output.geom”, connectivity has been changed. Node numbers and node sets are similar to the genesis file)


