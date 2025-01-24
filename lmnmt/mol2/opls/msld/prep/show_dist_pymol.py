from pymol import cmd, util
import sys

# input variables
resid1 = 'XRAY'
resid2 = 'REVR'

remove_solv = 'True'
#=========================
fx = open(f'noe_{resid1.lower()}.str', 'r')
anchors1 = []
for lines in fx.readlines():
    line = lines.split()
    try:
        if line[0] == 'assign':
            prot_atom = f'{line[3]}, {line[4]}, {line[5]}'
            lig_atom = f'{line[9]}, {line[10]}, {line[11]}'
            #print(prot_atom, lig_atom)
            anchors = prot_atom, lig_atom
            anchors1.append(anchors)
    except IndexError: None

ff = open(f'noe_{resid2.lower()}.str', 'r')
anchors2 = []
for lines in ff.readlines():
    line = lines.split()
    try:
        if line[0] == 'assign':
            prot_atom = f'{line[3]}, {line[4]}, {line[5]}'
            lig_atom = f'{line[9]}, {line[10]}, {line[11]}'
            #print(prot_atom, lig_atom)
            anchors = prot_atom, lig_atom
            anchors2.append(anchors)
    except IndexError: None
#========================

cmd.load('minimized.pdb')
pdb = 'minimized'
lig_sele = f'resname {resid1} or resname {resid2}'
#prot_sele = 'resid 330 or resid 376 or resid 396 or resid 90 or resid 345'

util.performance(0)
cmd.set('specular', 0)
#cmd.set('two_sided_lighting', 1)
cmd.extract (f'{resid1}', f'minimized and resname {resid1}')
cmd.extract (f'{resid2}', f'minimized and resname {resid2}')

cmd.set_color("c0", [0.2, 1.0, 0.2]); cmd.color("c0", "elem C")
#cmd.set_color("c0", [0.8, 0.8, 0.8]); cmd.color("c0", "elem C")
cmd.set_color("c1", [1.0, 1.0, 0.0]); cmd.color("c1", f"resname {resid1} and elem C") # Yellow for carbon
cmd.set_color("c2", [0.9, 0.9, 0.9]); cmd.color("c2", f"resname {resid2} and elem C") # Grey for carbon
##cmd.set_color("c3", [1.0, 0.0, 0.0]); cmd.color("c3", "elem O") # Red for oxygen
##cmd.set_color("c4", [0.0, 0.0, 1.0]); cmd.color("c4", "elem N") # Blue for nitrogen
##cmd.set_color("c5", [0.9, 0.84, 0.0]); cmd.color("c5", "elem S") # Dark yellow for sulfur


#cmd.show('spheres', f'{lig_sele}') # or {prot_sele}')
#cmd.hide('spheres', 'name H*')
#cmd.set('sphere_scale', 0.2, f'{lig_sele}') # or {prot_sele}')
##cmd.ray ('ray_trace_mode', 1)
##cmd.ray(1200, 900)

cmd.set('cartoon_color', 'grey80')
##cmd.show('sticks', prot_sele)
cmd.bg_color('white')
if remove_solv == 'True':
    cmd.remove('resname cof or resname pot or resname sod or resname cla or resname tip3')
else:
    cmd.extract('junk', 'resname tip3')
    cmd.remove('resname cof or resname pot or resname sod or resname cla')
#=================================================================================

# input variable
pose = anchors1

for i in range(len(pose)):
    prot = pose[i][0].split(',')
    lig = pose[i][-1].split(',')
    print(prot, lig)
    #cmd.hide('sticks', f"{lig[0]}")
    #cmd.hide('cartoon')
    #cmd.show('lines', f"{lig[0]}")
    cmd.show('lines', f'segid PROT and resid {prot[1]}')
    #cmd.show('cartoon', f'segid PROT and resid {prot[1]}')
    cmd.show('spheres', f"segid PROT and resid {prot[1]} and name {prot[-1]}")
    #cmd.label(selection=f"segid PROT and resid {prot[1]} and name O", expression='resn, resi')
    #cmd.set('label_size', 24)
    #cmd.label(selection=f"segid PROT and resid {prot[1]}", expression='resi')
    cmd.show('spheres', f"{lig[0]} and name {lig[-1]}")
    #cmd.distance(f'd{i+1}', f"segid PROT and resid {prot[1]} and name {prot[-1]}", f"{lig[0]} and name {lig[-1]}")

cmd.hide('sticks', 'name H*')
cmd.hide('lines', 'name H*')
cmd.hide('spheres', 'name H*')
cmd.set('sphere_scale', 0.2)

print('FINISHED')

#=================================================================================
## Setting: cartoon_transparency set to 0.60000.
## set ray_trace_mode, 1
## png xray_mdr.png, dpi=600, ray=1

