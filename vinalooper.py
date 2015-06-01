#!/usr/bin/env python
__author__ = 'Stefan Hinote'
import sys, argparse, os, pybel
#errno, openbabel
#from termcolor import colored
from AutoDockTools.MoleculePreparation import AD4LigandPreparation
from MolKit import Read
import subprocess, time

def dockligand(ligand,receptor,x,y,z,xc,yc,zc,o):
	command = ['vina','--ligand',o+'/ligandPDBQTs/'+ligand,'--receptor',receptor+'qt','--center_x',xc,'--center_y',yc,'--center_z',zc,'--size_x',x,'--size_y',y,'--size_z',z,'--log',o+'/logs/log_'+ligand,'--out',o+'/docked/docked_'+ligand]
	subprocess.call(command,stdout=subprocess.PIPE)

def prepareligand(dir,ligand):
    ligand_filename = dir+'ligandPDBs/'+ligand
    verbose = None
    add_bonds = False
    repairs = ""
    charges_to_add = 'gasteiger'
    preserve_charge_types=''
    cleanup  = "nphs_lps"
    allowed_bonds = "backbone"
    root = 'auto'
    outputfilename = dir+'ligandPDBQTs/'+ligand+'qt'
    check_for_fragments = False
    bonds_to_inactivate = ""
    inactivate_all_torsions = False
    attach_nonbonded_fragments = False
    attach_singletons = False
    mode = 'automatic'
    dict = None

    mols = Read(ligand_filename)
    mol = mols[0]
    coord_dict = {}
    for a in mol.allAtoms: coord_dict[a] = a.coords

    mol.buildBondsByDistance()
    if charges_to_add is not None:
        preserved = {}
        preserved_types = preserve_charge_types.split(',') 
        for t in preserved_types:
            if not len(t): continue
            ats = mol.allAtoms.get(lambda x: x.autodock_element==t)
            for a in ats:
                if a.chargeSet is not None:
                    preserved[a] = [a.chargeSet, a.charge]

    LPO = AD4LigandPreparation(mol, mode, repairs, charges_to_add, 
                            cleanup, allowed_bonds, root, 
                            outputfilename=outputfilename,
                            dict=dict, check_for_fragments=check_for_fragments,
                            bonds_to_inactivate=bonds_to_inactivate, 
                            inactivate_all_torsions=inactivate_all_torsions,
                            attach_nonbonded_fragments=attach_nonbonded_fragments,
                            attach_singletons=attach_singletons)

    if charges_to_add is not None:
        #restore any previous charges
        for atom, chargeList in preserved.items():
            atom._charges[chargeList[0]] = chargeList[1]
            atom.chargeSet = chargeList[0] 
    bad_list = []
    for a in mol.allAtoms:
        if a in coord_dict.keys() and a.coords!=coord_dict[a]: 
            bad_list.append(a)
    if len(bad_list):
        print len(bad_list), ' atom coordinates changed!'    
        for a in bad_list:
            print a.name, ":", coord_dict[a], ' -> ', a.coords
    else:
        if verbose: print "No change in atomic coordinates"
    if mol.returnCode!=0: 
        sys.stderr.write(mol.returnMsg+"\n")
	
def dock(input,receptor,x,y,z,xc,yc,zc,o,notes):
	try:					
        	os.makedirs(o)
		print "\033[94m [+]\033[0m output directory created"
	except OSError as exception:
		print "\033[94m [+]\033[0m output directory already exists, exiting..."
		sys.exit()

	os.makedirs(o+"/ligandPDBs")
	os.makedirs(o+"/ligandPDBQTs")
	os.makedirs(o+"/docked")
	os.makedirs(o+"/logs")
	
	#extract ligands from SDF to PDBs
	for mol in pybel.readfile("sdf", input):
        	#print len(mol.atoms)
        	output = pybel.Outputfile("pdb", o+"/ligandPDBs/"+mol.title+".pdb",overwrite=True)
        	output.write(mol)
        	output.close()
	print "\033[94m [+]\033[0m extracted ligands from SDF"
#rename PDB files, autodetect if from ZINC and if so then rename with zinc ids.
	
	import MolKit.molecule
	#import MolKit.protein
	from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation

	receptor_filename =  receptor
	mols = Read(receptor_filename)
	mol = mols[0]
	preserved = {}
	mol.buildBondsByDistance()
	RPO = AD4ReceptorPreparation(mol, 'automatic', '', 'gasteiger', 
                        'nphs_lps_waters_nonstdres', outputfilename=None,
                        preserved=preserved, 
                        delete_single_nonstd_residues=None,
                        dict=None)    

	for atom, chargeList in preserved.items():
		atom._charges[chargeList[0]] = chargeList[1]
		atom.chargeSet = chargeList[0]

	#loop for calling prepareligand
	liganddir = o+'/ligandPDBs'
	for root, dirs, files in os.walk(liganddir):
    		for name in files:
        		#print(os.path.join(root, name))
			prepareligand(o+'/',name)

	#loop for calling dockligand
	ligandPDBQTs = o+'/ligandPDBQTs'
	startingi = 1
	numligands = str(len(os.listdir(ligandPDBQTs)))
        for root, dirs, files in os.walk(ligandPDBQTs):
                for name in files:
                        dockligand(name,receptor,x,y,z,xc,yc,zc,o)
			print "\t\033[94m [+]\033[0m "+str(startingi)+"/"+numligands+" docked"
			startingi += 1

	#scan log files, sort docking scores, write summary file
	log = o+'/logs/'
	start = '   1';
	end = '0.000'
	ligands = []
	scores = []

	for root, dirs, files in os.walk(log):
		for name in files:
			data = open(log+name,'r')
			ligand = name.split('log_')[1].split('.')[0]
			scoring = ''.join(data.readlines())
			scoring = scoring.split(start)[1].split(end)[0]
			scoring = float(''.join(scoring.split()))
			ligands.append(ligand)
			scores.append(scoring)
			data.close()

	scores, ligands = zip(*sorted(zip(scores, ligands)))
	summary = open(o+'/summary','w')
	top = '[VinaLooper] Virtual screening with Autodock VINA\nDATE: '+time.strftime("%m/%d/%y")+'\nRECEPTOR: '+receptor+"\nX: "+x+" Y: "+y+" Z: "+z+"\nOffset: X: "+xc+" Y: "+yc+" Z: "+zc+"\nexhaustiveness: 8\nNOTES: "+notes+"\n"
	summary.write(top+"\n\n")
	for i in range(0,len(scores)):
		summary.write(ligands[i]+"\t"+str(scores[i])+"\n")	
	summary.close()
        print "\033[94m [+]\033[0m docking output: "+o+"/"
        print "\033[94m [+]\033[0m docking summary: "+o+"/summary"
	print "           VinaLooper Homepage: http://bloc.pw/vinalooper"
	print " #################################################################"
	print " # If you used AutoDock Vina in your work, please cite:          #"
	print " #                                                               #"
	print " # O. Trott, A. J. Olson,                                        #"
	print " # AutoDock Vina: improving the speed and accuracy of docking    #"
	print " # with a new scoring function, efficient optimization and       #"
	print " # multithreading, Journal of Computational Chemistry 31 (2010)  #"
	print " # 455-461                                                       #"
	print " #                                                               #"
	print " # DOI 10.1002/jcc.21334                                         #"
	print " #                                                               #"
	print " # Please see http://vina.scripps.edu for more information.      #"
	print " #################################################################"
	#End of docking
#later add support to upload results to MySQL DB.

#main banner
print " \033[92mVinaLooper\033[0m \033[1;33m////\033[0m virtual screening \033[1;33m/\033[0m http://bloc.pw/vinalooper"
print "          \033[1;33m /\033[0mrequirements:"                  
print "          \033[1;33m///\033[0mAutodock Vina http://vina.scripps.edu"
print "         \033[1;33m////\033[0mAutodock Tools http://autodock.scripps.edu"
print "        \033[1;33m/////\033[0mOpenBabel http://openbabel.org"

class menuparser(argparse.ArgumentParser):
	def error(self, message):
		print "\033[92m [+]\033[0m usage: "
		print "\t-i	input SDF file"
		print "\t-r 	receptor PDB file"
		print "\t-x	x dimension (angstroms)"
		print "\t-y	y dimension (angstroms)"
		print "\t-z	z dimension (angstroms)"
		print "\t-xc	x offset (angstroms)"
		print "\t-yc	y offset (angstroms)"
		print "\t-zc	z offset (angstroms)"
		print "\t-o	output directory"
		print '\t-n	"any notes to save"'
		print "\033[94m [+]\033[0m crtl+Z to kill, or continue in interactive mode"
		i = raw_input(' Input ligand SDF file > ')
		r = raw_input(' Receptor PDB file > ')
		x = raw_input(' X dimension in angstroms > ')
		y = raw_input(' Y dimension in angstroms > ')
		z = raw_input(' Z dimension in angstroms > ')
		xc = raw_input(' X offset > ')
		yc = raw_input(' Y offset > ')
		zc = raw_input(' Z offset > ')
		o = raw_input(' Write output to directory > ')
		n = raw_input(' Notes > ')
		#message = '%s' % i
		#sys.stderr.write('%s' % message)
		dock(i,r,x,y,z,xc,yc,zc,o,n)
		sys.exit()		
parser = menuparser()
parser.add_argument('-i', help='Input ligands SDF format',required=True)
parser.add_argument('-r', help='Receptor PDB file',required=True)
parser.add_argument('-x', help='X dimension in angstroms',required=True)
parser.add_argument('-y', help='Y dimension in angstroms',required=True)
parser.add_argument('-z', help='Z dimension in angstroms',required=True)
parser.add_argument('-xc', help='X offset',required=True)
parser.add_argument('-yc', help='Y offset',required=True)
parser.add_argument('-zc', help='Z offset',required=True)
parser.add_argument('-o', help='Output folder', required=True)
parser.add_argument('-n', help='Notes', required=True)
args = parser.parse_args()


#print colored("input file: %s" % input_file, 'red')
dock(args.i,args.r,args.x,args.y,args.z,args.xc,args.yc,args.zc,args.o,args.n)
