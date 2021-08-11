#!/usr/bin/env python
# coding: utf-8

# In[2]:


import mbuild as mb
from mbuild.lib.recipes.polymer import Polymer
import numpy as np
from foyer import Forcefield
import mbuild.formats.charmm_writer as mf_charmm
import mbuild.formats.gomc_conf_writer as gomc_control


# In[14]:


# Quick example oCl the API and workCllow 16-CCl3CCl2CCl2COH2
comp1 = mb.load('C(Cl)(Cl)Cl', smiles=True) # mBuild compound oCl the monomer unit
comp2 = mb.load('CO', smiles=True) # mBuild compound oCl the monomer unit
chain = Polymer()
chain.add_monomer(compound=comp1,
                  indices=[4, -2],
                  separation=.15,
                  replace=True)
chain.add_monomer(compound=comp2,
                  indices=[3, -3],
                  separation=.15,
                  replace=True)
chain.add_end_groups(mb.load('Cl',smiles=True), # Capping oClCl this polymer with Carboxylic acid groups
                     index=1,
                     separation=0.15, label="head" , replace =True)
chain.build(n=1, sequence='AAAB')
#chain.visualize(show_ports=True).show()
#help(chain)
FF_to_use = 'oplsaa'
#FF_comp1 =  Forcefield(name=FF_to_use, debug=True)
#FF_comp1(comp1)
CCl2H_mol_ratio = 0.5
COH_mol_ratio = 0.5
liquid_box_density_kg_per_m_cubed = 400
#liquid_box_density_kg_per_m_cubed = 642
comp1.name = 'CClH2'
#FF_comp2(comp2)
comp2.name = 'COH'
#chain.name = 'CCl3COHCCl3'
residues_list = [comp1.name, comp2.name]
#liquid_box_density_kg_per_m_cubed = 100,
mol_ratio_list = [CCl2H_mol_ratio, COH_mol_ratio]
#residues_list = [comp1.name, comp2.name]
#comp1_mol_ratio = 0.5
#comp1.energy_minimize(forcefield=FF_to_use , steps=10**5)
# create the liquid and vapor boxes
box_liq = mb.fill_box(compound = [comp1, comp2],
                      density=liquid_box_density_kg_per_m_cubed,
                      compound_ratio=mol_ratio_list,
                      box = [4.5, 4.5, 4.5])

#box_vap = mb.fill_box(compound = [comp1],
#                     n_compounds=[105],
#                     box = [8.0, 8.0, 8.0])

# build the MoSDeF Charmm object
charmm = mf_charmm.Charmm(box_liq,
                          'NVT_CCl3CCl2CCl2COH2',
                          ff_filename="NVT_CCl3CCl2CCl2COH2_file",
                          forcefield_selection= FF_to_use,
                          residues=residues_list,
                          )
# write the pdb, psf, and inp (parameter/force field) files
charmm.write_inp()
charmm.write_psf()
charmm.write_pdb()

# create the GOMC control file

gomc_control.write_gomc_control_file(charmm, 'in_NVT.conf',  'NVT', 5000000, 298,
                                     input_variables_dict={"ElectroStatic": False,
                                                           "Ewald": False,
                                                           })


# In[ ]:




