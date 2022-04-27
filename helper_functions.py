#!/usr/bin/env python
# coding: utf-8

"""
This file is part of BrainMolecularAtlas.

This module provide helper functions 

Copyright (c) 2021-2022 Blue Brain Project/EPFL 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import pandas as pd
import numpy as np

import timeit
from collections import Counter

def variability_score(levels_list):
    """
    return: variability score as a measure for the visualisation 
    """
    var_score = np.abs((1/len(levels_list))*((np.std(levels_list))/(np.mean(levels_list)))) # this one is not for the log scale data
    return(var_score)


def get_samples(arr,sample_size,num_samples):
    """
    return: samples: np.ndarray (num_samples, sample_size)
    """
    samples = np.zeros((num_samples + 1, sample_size), np.float64)
    
    for elem in range(0, num_samples):
        sample = np.random.choice(arr,sample_size,replace=False)
        samples[elem] = sample
    return samples


def split_df_column_list_to_multiple_rows(df,column_of_lists,output_type=float):
    ''' 
    return: df where list elements in specified column are separated into new rows
    '''
    new_rows = []
    
    def split_list(row):
        current_value = row[column_of_lists]
        
        if isinstance(current_value, list):
            
            if current_value == []:
                new_row = row.to_dict()
                new_row[column_of_lists] = None
                new_rows.append(new_row)
                
            else:    
                for s in current_value:
                    new_row = row.to_dict()
                    new_row[column_of_lists] = s
                    new_rows.append(new_row)

            
        else:
            new_row = row.to_dict()
            new_row[column_of_lists] = current_value
            new_rows.append(new_row)
            
    
    df.apply(split_list, axis=1)
    
    df_no_lists = pd.DataFrame(new_rows)
    return df_no_lists
