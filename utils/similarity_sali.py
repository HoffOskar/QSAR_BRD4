### Imports

import pandas as pd 
import numpy as np
import math

from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs

from sklearn.neighbors import NearestNeighbors


### Module class

class MorganSimilarity:

    """
    Class to perform similarity operations using Morgan fingerprints. 
    
    Attributes:
        df (pd.DataFrame): The dataset containing molecular structures as mol objects.
        act_col (str): Column name for activity values.
        mol_col (str): Column name for molecules.
        radius (int): Radius for Morgan fingerprint generation.
        nBits (int): Bit length of Morgan fingerprints.
        useChirality (bool): Whether to consider chirality for Morgan fingerprint generation.
        fp (list): List of Morgan fingerprints as Bit Vectors. Can be formatted as a NumPy array or DataFrame.
        bit_info (list): List of dictionaries containing bit information for inspection.
        SALI_dict (dict): Dictionary containing NumPy arrays of the same shape describing pairwise 
            - distance (float)
            - activity values (float; e.g. pIC50) - name depends on the nature of the original dataset
            - delta activity values (float; e.g. ΔpIC50) - name depends on the nature of the original dataset
            - SALI (float): Similarity Activity Landscape Index
            - Positional indices of the k-nearest neighbors (int) 
            - Positional indices of the molecules as integers (int)
        info (dict): Dictionary containing class metadata.
        dist_metric (str): The distance metric used for the k-nearest neighbors.
        pair_analysis (pd.DataFrame): DataFrame containing molecule pairs with SALI values above the chosen threshold.
    """
    

    def __init__(self, df: pd.DataFrame, act_col='pIC50', mol_col='Mol', radius=2, nBits=1024, useChirality=True):
        
        self.df = df.copy()
        self.act_col = act_col
        self.mol_col = mol_col
        self.radius = radius
        self.nBits = nBits
        self.useChirality = useChirality
        self.info = {}


    def __get_fp__(self, mol):
        """
        Create a tuple of Morgan fingerprints and bit information for the mol object.
        """

        bit_info = {}
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol,
            radius=self.radius,
            nBits=self.nBits,
            useChirality=self.useChirality,
            bitInfo=bit_info
        )
        return fp, bit_info


    def get_morgan_fingerprint(self):
        """
        Create Morgan fingerprints for the classes DataFrame containing a Mol column.

        Returns:
            A list of Morgan fingerprints as Bit Vectors.
        
        Side effects:
                self.bit_info: A list of dictionaries containing bit information for inspection.
                self.fp: A list of Morgan fingerprints as Bit Vectors.
                self.info is updated. 
        """

        df = self.df

        ## Process each molecule in the DataFrame in a vectorized manner
        fps_bit_info = df['Mol'].apply(self.__get_fp__)

        ### Unzip the results into separate lists
        fingerprints, bit_info_list = zip(*fps_bit_info)
        
        ### Update the class attributes
        self.bit_info = list(bit_info_list)
        self.fp = list(fingerprints)
        self.info['fp_format'] = 'List of BitVectors'
        self.info['radius'] = self.radius
        self.info['nBits'] = self.nBits
        self.info['useChirality'] = self.useChirality

        return list(fingerprints)
    

    def format_fp(self, format='boolean_np_array'):
        """
        Reformat the Morgan fingerprints as a boolean NumPy array, DataFrame, or NumPy array. 

        Parameters:
        format (str): The format to return the fingerprints in:
                  - 'boolean_np_array' 
                  - 'DataFrame'
                  - 'np_array'

        Returns:
            Morgan fingerprints in the specified format.

        Side effects:
            self.fp: The fingerprints are updated in the specified format.
            self.info['fp_format'] is updated. 

        """

        if format == 'boolean_np_array':
            fp = np.array(self.fp).astype(bool)
        
        elif format == 'DataFrame':
            fp = pd.DataFrame(np.array(self.fp), index=self.df.index)

        elif format == 'np_array':
            fp = np.array(self.fp)            
        
        else:
            raise ValueError("Invalid format. Please use 'boolean_np_array' (default), 'np_array' or 'DataFrame'.")
        
        ### Update the class attributes
        self.fp = fp
        self.info['fp_format'] = format

        return fp
    

    def get_SALI_matrix(self, k=10, dist_metric='jaccard', stat_summary=False):
        """
        Compute the Similarity Activity Landscape Indexes (SALI) for the k-nearest neighbors.

        Parameters:
            k (int): The number of nearest neighbors to consider.
            dist_metric (str): The distance metric to use for the k-nearest neighbors.
        
        Returns:
            A NumPy array containing the SALI values for each molecule's k-nearest neighbors.
        
        Side effects:
            self.SALI_dict: A dictionary containing NumPy arrays of identical shape: 
                - distance (float)
                - activity values (float; e.g. pIC50) - name depends on the nature of the original dataset
                - delta activity values (float; e.g. ΔpIC50) - name depends on the nature of the original dataset
                - SALI (float): Similarity Activity Landscape Index
                - Positional indices of the k-nearest neighbors (int) 
                - Positional indices of the molecules as integers (int)
        """
        
        ### Instantiate the KNN model
        neigh = NearestNeighbors(n_neighbors=k+1, metric=dist_metric)  

        ### Calculate the distances and indices of the k-nearest neighbors
        neigh = neigh.fit(self.fp)
        distances_np, neigb_idx = neigh.kneighbors(self.fp, return_distance=True)

        ### Removing the first column (distance to each molecule itself)
        distances_np = distances_np[:, 1:]
        neigb_idx = neigb_idx[:, 1:]

        ### Convert pIC50 column to NumPy array for fast indexing
        pIC50_np = self.df[self.act_col].values

        ### Compute the absolute differences for each molecule's k-nearest neighbors
        delta_pIC50_np = np.abs(pIC50_np.reshape(-1,1) - pIC50_np[neigb_idx])

        ### Compute the SALI matrix
        SALI_np = delta_pIC50_np / distances_np

        ### Distribution analysis
        if stat_summary:
        
            SALI_df = pd.Series(SALI_np.flatten())
            print(f'{len(SALI_df):.2e} pairs analyzed (including duplicates).')
            print(f'{math.comb(len(self.df),2):.2e} pairs are possible.')

            print('Quantiles (SALI):')
            print(SALI_df.quantile([0, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999]))
                
        ### Summarizing all NumPy arrays in one dictionary 
        SALI_dict = {
            'distance': distances_np,
            self.act_col: pIC50_np,
            f'delta_{self.act_col}': delta_pIC50_np,
            'SALI': SALI_np,
            'neigb_idx': neigb_idx,
            'row_idx': np.array([list(range(neigb_idx.shape[0])) for i in range(neigb_idx.shape[1])]).T
        }
        
        ### Update the class attributes
        self.SALI_dict = SALI_dict
        self.dist_metric = dist_metric
        self.info['k'] = k
        self.info['dist_metric'] = dist_metric

        return SALI_np


    def get_sali_pairs(self, SALI_threshold, min_activity_diff=None):
        """
        List molecule pairs with SALI values higher than the threshold. 

        Parameters:
            SALI_threshold (int or float): The threshold for the SALI value.
            min_activity_diff (None, int, or float): A minimum difference in activity can be defined.

        Returns:
            A DataFrame containing the molecule pairs with high SALI values.

        Side effects:
            self.pair_analysis (pd.DataFrame): A DataFrame containing the molecule pairs with high SALI values.
        """
        
        ### Unpack the SALI dictionary
        distances_np = self.SALI_dict['distance']
        pIC50_np = self.SALI_dict[f'{self.act_col}']
        delta_pIC50_np = self.SALI_dict[f'delta_{self.act_col}']
        SALI_np = self.SALI_dict['SALI']
        neigb_idx = self.SALI_dict['neigb_idx']
        row_idx = self.SALI_dict['row_idx']

        ### Create a map for boolean indexing of the NumPy arrays
        SALI_mask = SALI_np > SALI_threshold

        ### Formatting the columns for the first molecule of each pair
        row_df = self.df.reset_index().iloc[row_idx[SALI_mask]]
        row_df.rename(columns={name: name + '_1' for name in row_df.columns}, inplace=True)
        row_df.reset_index(inplace=True, drop=True)

        ### Formatting the columns for the second molecule
        neigb_df = self.df.reset_index().iloc[neigb_idx[SALI_mask]]
        neigb_df.rename(columns={name: name + '_2' for name in neigb_df.columns}, inplace=True)
        neigb_df.reset_index(inplace=True, drop=True)

        ### Combining all columns
        inspect_df = pd.concat([row_df, neigb_df], axis=1)

        ### Adding the comparison metrics 
        inspect_df[self.dist_metric] = distances_np[SALI_mask]
        inspect_df[f'Δ{self.act_col}'] = delta_pIC50_np[SALI_mask]
        inspect_df['SALI'] = SALI_np[SALI_mask]

        ### Create a new Series where each element is a tuple of the sorted IDs
        sorted_pairs = inspect_df.apply(lambda row: tuple(sorted((row['ID_1'], row['ID_2']))), axis=1)

        ### Generate a boolean mask that is True for any duplicated pair
        cross_dupl_mask = sorted_pairs.duplicated(keep='first')
        
        ### Remove duplicates
        inspect_df = inspect_df.loc[cross_dupl_mask]
        
        ### Filter by minimum activity difference
        if min_activity_diff:
            inspect_df = inspect_df.loc[inspect_df[f'Δ{self.act_col}'] > min_activity_diff]                                                                            

        ### Order by SALI values and reset the index
        inspect_df = inspect_df.sort_values('SALI', ascending=False)
        inspect_df.reset_index(drop=True, inplace=True)

        ### Update the class attribute
        self.pair_analysis = inspect_df
        self.info['SALI_threshold'] = SALI_threshold
        
        return inspect_df


    def count_cliffs(self, SALI_threshold, min_activity_diff=None):
        """
        Count the number of activity cliffs for each molecule. 

        Parameters:
            SALI_threshold (int or float): The threshold for the SALI value.
            min_activity_diff (None, int, or float): A minimum difference in activity can be defined.

        Returns:
            A DataFrame containing the molecules sorted by the number of cliffs.
        """
        
        ### Unpack the SALI dictionary
        distances_np = self.SALI_dict['distance']
        pIC50_np = self.SALI_dict[self.act_col]
        delta_pIC50_np = self.SALI_dict[f'delta_{self.act_col}']
        SALI_np = self.SALI_dict['SALI']
        neigb_idx = self.SALI_dict['neigb_idx']
        row_idx = self.SALI_dict['row_idx']

        ### Set the minimum activity difference
        if min_activity_diff is None:
            min_activity_diff = 0
        
        ### Create a boolean masks for indexing
        SALI_mask = SALI_np > SALI_threshold
        delta_pIC50_mask = delta_pIC50_np > min_activity_diff
        combined_mask = np.logical_and(SALI_mask, delta_pIC50_mask)

        ### Count the number of cliffs for each molecule
        cliff_count = np.sum(combined_mask, axis=1)

        ### Create a DataFrame for the results
        inspect_df = self.df.copy()
        inspect_df['CliffCount'] = cliff_count

        return inspect_df.sort_values('CliffCount', ascending=False)