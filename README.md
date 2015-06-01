# vinalooper

Python script for docking library of ligands from single SDF file to protein PDB file.
Requires the following be installed: AutoDockTools, AutoDock Vina, OpenBabel, OpenBabel Python API

Under Debian:
sudo apt-get install autodocktools
sudo apt-get install autodock-vina
sudo apt-get install openbabel
sudo apt-get install python-openbabel

OpenBabel used for extracting individual ligand PDB files from SDF file.
AutoDockTools used for preparing ligand PDB files to PDBQT format and preparing receptor PDB file to PDBQT
AutoDock Vina used for docking.

Usage:
python vinalooper.py arguments
	-i	input SDF file
	-r 	receptor PDB file
	-x	x dimension (angstroms)
	-y	y dimension (angstroms)
	-z	z dimension (angstroms)
	-xc	x offset (angstroms)
	-yc	y offset (angstroms)
	-zc	z offset (angstroms)
	-o	output directory
	-n	"any notes to save"
	
Eg: python vinalooper.py -i ligands.sdf -r receptor.pdb -x 20 -y 20 -z 20 -xc 4.5 -yc -5.6 -zc 14.2 -o myoutput -n "some notes"
OR pass no arguments and interactively define parameters. 	
	
