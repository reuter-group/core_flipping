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

# input variable
pose = anchors2

cmd.hide('lines', f"segid PROT and resid {prot[1]}")
cmd.hide('spheres', f"segid PROT and resid {prot[1]} and name {prot[-1]}")

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

cmd.hide('lines', 'name H*')
cmd.hide('spheres', 'name H*')
cmd.set('sphere_scale', 0.2)

print('FINISHED')

#=================================================================================
## set ray_trace_mode, 1
## hide cartoon
## png xray_mdr.png, dpi=600, ray=1

