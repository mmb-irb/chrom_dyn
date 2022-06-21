# Chromatin Dynamics
## Molecular Modelling and Bioinformatics Group, IRB Barcelona

### General description

Chromatin Dynamics allows to obtain 3D structures of chromatin fibers at the base pair 
level considering nucleosome positions. It requires as input the sequence of linker DNA
and the position of the nucleosomes in the sequence. The conformational sampling is done
by a Monte Carlo (MC) algorithm, considering flexible linkers and rigid nucleosomes.
Linker DNA is represented at the base pair level by a pseudo-harmonic potential expressed
in helical parameters (rise, slide, shift, twist, roll, tilt). Debye-Huckel electrostatics
and excluded volume potentials are added to avoid overlaps.

### Precompiled version for Linux

Tested on different Ubuntu and Fedora terminals.<br/>
No installation required (see below for compiling/installing from source files)

Execution:

> sh run.sh [nucleosome positions] [linker sequence] [# of MC iterations] [output folder]

- The output folder should not exist and an absolute path should be given.
- See DESCRIPTION below for a brief explanation of required files.

Nucleosome positions file example: nucl_pos.dat<br/>
Linker DNA sequence file example: link_seq.dat

Example:

> sh run.sh nucl_pos.dat link_seq.dat 100 ~/chrom_dyn/test

*Outputs generted inside the output folder:*

"output_pos":<br/>
&emsp; xyz files containing the 3D structures obtained for the fibers (linker DNA + nucleosomes).<br/>
&emsp; They are numbered consecutively as cartnucl_000000_000000.xyz, cartnucl_000000_000001.xyz, ...).<br/>
&emsp; Each nucleosome center is represented as a N atom and each base pair center as a C atom.<br/>
&emsp; The structures can be visualized using VMD software. To make the nucleosome core visible, do<br/>
&emsp; &emsp; Representations -> Selected Atoms: name N; Drawing method:VDW; Sphere Scale: 14.<br/>

"output_nucl":<br/>
&emsp; xyz files containing 3D nucleosome configurations for each structure.<br/>
&emsp; Numbered consecutively as nucl_000000_000001.xyz, nucl_000000_000002.xyz, ...).<br/>
&emsp; Each nucleosome center is represented as a N atom.<br/>

The expected time of execution is approximately 10 seconds per MC iteration on a standard
desktop computer.

------------------------------------------------------------------------------------------
### Web version

Chromatin dynamics can be executed on a web server: http://vre.multiscalegenomics.eu/launch/

Select "Chromatin Dynamics" > "CREATE 3D REPRESENTATION OF CHROMATIN FIBER FROM SEQUENCE"

DNA sequence: [linker sequence file]<br/>
Position of nucleosomes: [nucleosome positions file]<br/>
Operations: "Create Trajectory"<br/>
Number of structures: number of structures to generate.<br/>

Help: http://vre.multiscalegenomics.eu/tools/chromatindyn/help/help.php

------------------------------------------------------------------------------------------
### Source files

The src/ directory contains the code for the two instances of the program: obtaining the
starting 3D structure from the linker sequence and the nucleosome positions (ini_str
directory) & simulating the different fiber conformations (sim directory). The main code
is in the file "main_tetra_di.c" of each directory. Type make to compile the program and
execute as indicated in the README file.

The program compiles in less than a minute on a standard desktop computer.

