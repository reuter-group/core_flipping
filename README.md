# core_flipping
The directories are organized as:

```
HNE directory:
    molX:  MSLD prep, varibales.inp analysis and nbond.str files for compound X
           compound paramter files are inside the prep directory
    msld_files: CHARMM input scripts for flattening and production
                bash scripts for running MSLD flattening and production
    toppar: CHARMM36m force field files
```
```
lmnmt:
    molX: compounds 1-3 follow the same structure
        charmm-cgenff:
                      msld: MSLD prep directory, varibales.inp analysis and nbond.str files for compound X
                      osp : prep directory, CHARMM and bash run scripts, python analysis file
                      simulation: protein, ligand psf, crd and force field stream files
        opls:
            msld:
            osp :
            simulation:

    msld_files: CHARMM input scripts for flattening and production
                bash scripts for running MSLD flattening and production

    multiple_distance_restraints: python script to generate protein-ligand multiple distance restraints

    simularion_files: CHARMM scripts for running small MD simulation for a protein-ligand complex in a cubic water box
``` 
