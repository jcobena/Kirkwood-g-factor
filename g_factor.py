'''
Created on Nov 1, 2015

@author: Jose_Cobena
'''
import sys, csv
import copy
import numpy as np
import pandas as pd
import itertools
import math as mt
import collections
from itertools import count, zip_longest, islice
#from numpy.core.setup_common import fname2def
#from pandas.core.frame import DataFrame

def main():
    time_snapshot = 33000000
    fixed_dipole = 2.32
    init_r = 1
    dr = 1
    max_dr = 20
    N = 851
    
    
    
#########################################################################################################
                
    with open('kirkwood-data-{:08}.txt'.format(time_snapshot), 'r') as file_t0:
        
        test_x = positions_list(file_t0)
        structures_molecule = recreating_molecule(test_x, fixed_dipole)
        
        #print(structures_molecule)
        
        test2 = get_kirkwood_f(structures_molecule, init_r, dr, max_dr, fixed_dipole, N)
        
        #print(test2)
        
def get_kirkwood_f(structures_molecule, init_r, dr, max_dr, fixed_dipole, N):
    r_for_rdfs = []
    copied_list = copy.deepcopy(structures_molecule)
    #results = []
    for x in np.arange(init_r, max_dr, dr):
        r_for_rdfs.append(x)
    #print(r_for_rdfs)
    
    for rs in r_for_rdfs:
        momentum_i = 0
        summed_resultant_momenta = 0
        for molecule_vector_i in structures_molecule:
            neigh_list_centered = []
            
            #center hs from molecule i:
            hi1_centered_x = molecule_vector_i[2][0]- molecule_vector_i[1][0]
            hi1_centered_y = molecule_vector_i[2][1]- molecule_vector_i[1][1]
            hi1_centered_z = molecule_vector_i[2][2]- molecule_vector_i[1][2]
            hi1c = [hi1_centered_x, hi1_centered_y, hi1_centered_z]
            
            hi2_centered_x = molecule_vector_i[3][0]- molecule_vector_i[1][0]
            hi2_centered_y = molecule_vector_i[3][1]- molecule_vector_i[1][1]
            hi2_centered_z = molecule_vector_i[3][2]- molecule_vector_i[1][2]
            hi2c = [hi2_centered_x, hi2_centered_y, hi2_centered_z]
            
            resultante_ix = hi1_centered_x + hi2_centered_x
            resultante_iy = hi1_centered_y + hi2_centered_y
            resultante_iz = hi1_centered_z + hi2_centered_z
            
            alfa_ii = resultante_ix/ mt.sqrt(((resultante_ix*resultante_ix) + (resultante_iy*resultante_iy) + (resultante_iz*resultante_iz)))
            beta_ii = resultante_iy/ mt.sqrt(((resultante_ix*resultante_ix) + (resultante_iy*resultante_iy) + (resultante_iz*resultante_iz)))
            gama_ii = resultante_iz/ mt.sqrt(((resultante_ix*resultante_ix) + (resultante_iy*resultante_iy) + (resultante_iz*resultante_iz)))
            
            x_final_i = fixed_dipole*alfa_ii
            y_final_i = fixed_dipole*beta_ii
            z_final_i = fixed_dipole*gama_ii
            #print(x_final_i, y_final_i, z_final_i)
            
            #print(x_final_i, y_final_i, z_final_i)
            
            
            for molecule_vector_j in copied_list:
                
                if molecule_vector_i[0] != molecule_vector_j[0]:
                
                #find the oxygens atoms near the radius:
                    r_test = (molecule_vector_j[1][0] - molecule_vector_i[1][0])*(molecule_vector_j[1][0] - molecule_vector_i[1][0]) + \
                             (molecule_vector_j[1][1] - molecule_vector_i[1][1])*(molecule_vector_j[1][1] - molecule_vector_i[1][1]) + \
                             (molecule_vector_j[1][2] - molecule_vector_i[1][2])*(molecule_vector_j[1][2] - molecule_vector_i[1][2])
                             
                    
                    
        
                    if r_test < (rs*rs):
                        #h1_and_h2 = [molecule_vector_j[2],molecule_vector_j[3]]
                        #neigh_list_centered.append(h1_and_h2)
                        #move the Hs to the origin:
                        h1_centered_x = molecule_vector_j[2][0]- molecule_vector_j[1][0]
                        h1_centered_y = molecule_vector_j[2][1]- molecule_vector_j[1][1]
                        h1_centered_z = molecule_vector_j[2][2]- molecule_vector_j[1][2]
                        h1c = [h1_centered_x, h1_centered_y, h1_centered_z]
                        
                        h2_centered_x = molecule_vector_j[3][0]- molecule_vector_j[1][0]
                        h2_centered_y = molecule_vector_j[3][1]- molecule_vector_j[1][1]
                        h2_centered_z = molecule_vector_j[3][2]- molecule_vector_j[1][2]
                        h2c = [h2_centered_x, h2_centered_y, h2_centered_z]
                        
                        hs = [h1c, h2c]
                        neigh_list_centered.append(hs)
                    #print(neigh_list_centered) #checked
                            
                    #getting the vector resultante:
                    #all_resultant_in_i = []
                        for vectors in neigh_list_centered:
                            all_resultant_in_i = []
                            resultante_x = vectors[0][0] + vectors[1][0]
                            resultante_y = vectors[0][1] + vectors[1][1]
                            resultante_z = vectors[0][2] + vectors[1][2]
                            
                            resultante_v = [resultante_x, resultante_y, resultante_z]
                            all_resultant_in_i.append(resultante_v)
                            
                        #getting directional cosines:
                        #dir_cosines = []
                        for coords in all_resultant_in_i:
                            dir_cosines = []
                            alfa_i = coords[0]/ mt.sqrt(((coords[0]*coords[0]) + (coords[1]*coords[1]) + (coords[2]*coords[2])))
                            beta_i = coords[1]/ mt.sqrt(((coords[0]*coords[0]) + (coords[1]*coords[1]) + (coords[2]*coords[2])))
                            gama_i = coords[2]/ mt.sqrt(((coords[0]*coords[0]) + (coords[1]*coords[1]) + (coords[2]*coords[2])))
                            dir_cosines_j = [alfa_i, beta_i, gama_i]
                            dir_cosines.append(dir_cosines_j)
                            
                        #getting the vectors:
                        for angles_i in dir_cosines:
                            final_coordinates = []
                            x_final = fixed_dipole*angles_i[0]
                            y_final = fixed_dipole*angles_i[1]
                            z_final = fixed_dipole*angles_i[2]
                            final_coords = [x_final, y_final, z_final]
                            final_coordinates.append(final_coords)
                        #print(final_coordinates) #checked
                

            #summing up the dipole momentum around molecule i:
                        zipped_list = zip(*final_coordinates)
                        
                        summed_vector_dipole_momenta = [sum(col) for col in zip(*final_coordinates)]
                        
                        dot_product_x = x_final_i*summed_vector_dipole_momenta[0]
                        dot_product_y = y_final_i*summed_vector_dipole_momenta[1]
                        dot_product_z = z_final_i*summed_vector_dipole_momenta[2]
                        t_dot = dot_product_x + dot_product_y + dot_product_z
                    #t_dot =  dot_product_z
                    
                    
                    #resultant_momenta = mt.sqrt(dot_product_x*dot_product_x + dot_product_y*dot_product_y + dot_product_z*dot_product_z)
                    #resultant_momenta = mt.sqrt(dot_product_x + dot_product_y + dot_product_z)
                    #print(resultant_momenta)
                    
                    #summed_resultant_momenta += resultant_momenta
                        summed_resultant_momenta += t_dot
            
            #print(summed_resultant_momenta)
        ave_resultant_momenta = summed_resultant_momenta/N
        #print(ave_resultant_momenta)
        
        ave_resultant_momenta_norm = (ave_resultant_momenta/(fixed_dipole*fixed_dipole))+1.0 #fixed_dipole  #1.84
        
        print(ave_resultant_momenta_norm)
    

        
def recreating_molecule(test_x, fixed_dipole):
  
    final_coordinates = []
    oxys = islice(test_x, 0, None, 3)
    h1s = islice(test_x, 1, None, 3)
    h2s = islice(test_x, 2, None, 3)
    
    for oxy, h1, h2 in zip(oxys, h1s, h2s):
        molecules_structure = [oxy[0],oxy[1],h1[1],h2[1]]
        final_coordinates.append(molecules_structure)
        
    return final_coordinates
        
        
def positions_list(file_tx):
    dict_atoms = {}
    coords = []
    for _skip in range(9):
        next(file_tx)
    for rows in file_tx:
        values = rows.split()
        keys = int(values[0])
        coords = [float(values[1]), float(values[2]), float(values[3])]
        dict_atoms[keys] = coords   
        
    positions = sorted(dict_atoms.items())    
    return positions


if __name__ == "__main__": main()
