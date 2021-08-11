#!/usr/bin/env python
# coding: utf-8

# In[2]:


import mbuild as mb
from mbuild.lib.recipes.polymer import Polymer
import numpy as np
from foyer import Forcefield
import mbuild.formats.charmm_writer as mf_charmm
import mbuild.formats.gomc_conf_writer as gomc_control


# In[9]:


# Quick example oCl the API and workCllow 21-CH3CH2CH2CHOHCCl3
comp1 = mb.load('C(Cl)(Cl)Cl', smiles=True) # mBuild compound oCl the monomer unit
comp2 = mb.load('CO', smiles=True) # mBuild compound oCl the monomer unit
comp3 = mb.load('C', smiles=True) # mBuild compound oCl the monomer unit
chain = Polymer()
chain.add_monomer(compound=comp1,
                  indices=[4, -1],
                  separation=.15,
                  replace=True)
chain.add_monomer(compound=comp2,
                  indices=[3, -2],
                  separation=.15,
                  replace=True)
chain.add_monomer(compound=comp3,
                  indices=[2, -4],
                  separation=.15,
                  replace=True)
chain.add_end_groups(mb.load('Cl',smiles=True), # Capping oClCl this polymer with Carboxylic acid groups
                     index=1,
                     separation=0.15, label="head" , replace =True)
chain.build(n=1, sequence='CCCBA')
#chain.visualize(show_ports=True).show()
#help(chain)
Liquid_box_length_Ang = 45
Liquid_box_Total_molecules = 850

#butane_mol_ratio = 0.5

FF_to_use = 'oplsaa'
#FF_comp1 =  Forcefield(name=FF_to_use, debug=True)
#FF_comp1(comp1)
#FF_to_use.apply(comp1)
comp1.name = 'CClH2'
comp2.name = 'COH'
comp3.name = 'CH3'
residues_list = [comp1.name, comp2.name, comp3.name]
#liquid_box_density_kg_per_m_cubed = 100,
#liquid_box_density_kg_per_m_cubed = 642,
#comp1_mol_ratio = 0.5
#comp1.energy_minimize(forcefield=FF_to_use , steps=10**5)
#Bead_to_atom_name_dict = { '_CH3':'C', '_CH2':'C',  '_CH':'C', '_HC':'C'}
#fixed_bonds_angles_list = [comp1.name]
# create the liquid and vapor boxes
box_liq = mb.fill_box(compound = [comp1, comp2, comp3],
                      n_compounds=[500,500,500],
                      box = [4.5, 4.5, 4.5])
 
#box_vap = mb.fill_box(compound = [comp1],
#                     n_compounds=[105],
#                     box = [8.0, 8.0, 8.0])

# build the MoSDeF Charmm object
charmm = mf_charmm.Charmm(box_liq,
                          'NPT_CH3CH2CH2CHOHCCl3',
                          ff_filename="NPT_CH3CH2CH2CHOHCCl3_file",
                          forcefield_selection= FF_to_use,
                          residues= residues_list,
                          bead_to_atom_name_dict= None,
                          gomc_fix_bonds_angles= None,
                          )
"""charmm = mf_charmm.Charmm(box_liq,
                          'NPT_n_hexane_water_liq',
                          structure_box_1=None  ,
                          filename_box_1=None,
                          ff_filename ="NPT_n_hexane_water_FF" ,
                          forcefield_selection=forcefield_files ,
                          residues=Molecule_ResName_list ,
                          bead_to_atom_name_dict=Bead_to_atom_name_dict ,
                          gomc_fix_bonds_angles=fixed_bonds_angles_list,
                         )"""
# write the pdb, psf, and inp (parameter/force field) files
charmm.write_inp()
charmm.write_psf()
charmm.write_pdb()

# create the GOMC control file

gomc_control.write_gomc_control_file(charmm, 'in_NPT.conf',  'NPT', 25000, 298,
                                     input_variables_dict= {"Pressure": 410,
                                                            "Exclude": "1-4",})

#help(gomc_control.write_gomc_control_file)
#gomc_control.print_valid_ensemble_input_variables('NPT', description=T


# In[ ]:




