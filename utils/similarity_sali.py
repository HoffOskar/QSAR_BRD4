### Imports

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pd_mol import standardize_mol_col
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina
from sklearn.neighbors import NearestNeighbors

####################################################################################################
### Utility functions: Bemis-Murcko Scaffold Analysis


def get_bemis_murcko(df, mol_col="Mol", return_smiles=True, return_cluster=False):
    """
    Generate standardized generic Bemis-Murcko Fragments.

    Parameters:
    ----------
    df : pandas.DataFrame
        DataFrame containing a column of RDKit Mol objects.
    mol_col : str
        Column name containing RDKit Mol objects.
    return_smiles : bool
        If True, an additional column 'BM_smiles' with SMILES strings is added to the DataFrame.
    return_cluster : bool
        If True, an additional column 'BM_cluster' with cluster IDs is added to the DataFrame.

    Returns:
    --------
    df : pandas.DataFrame
        Updated DataFrame.

    Notes:
    ------
    - Generic Bemis-Murcko: Side chains, stereochemistry, and bond orders are removed. Every atom is set to carbon.
    - Standardized Bemis-Murcko: Standardization results in identical mol objects for compounds with the same scaffold.
    """

    ### Generate generic Bemis-Murcko scaffold
    print("MurckoScaffold.MakeScaffoldGeneric")
    BM_scaf = df[mol_col].progress_apply(MurckoScaffold.MakeScaffoldGeneric)

    ### Generate Bemis-Murcko scaffold
    #   Apparently this sequence adds robustness
    print("MurckoScaffold.GetScaffoldForMol")
    df["BM_mol"] = BM_scaf.progress_apply(MurckoScaffold.GetScaffoldForMol)

    ### Standardize Bemis-Murcko scaffold
    df["BM_mol"] = standardize_mol_col(
        df["BM_mol"],
        largest_frag=False,
        remove_charge=False,
    )

    ### Generate scaffold SMILES
    if return_smiles or return_cluster:
        ### Generate SMILES
        print("Chem.MolToSmiles")
        df["BM_smiles"] = df["BM_mol"].progress_apply(
            lambda mol: Chem.MolToSmiles(mol, canonical=True)
        )

    ### Generate cluster IDs
    if return_cluster:
        ### Get unique scaffolds
        unique_scaffolds = df["BM_smiles"].unique()

        ### Dictionary to map each scaffold to a cluster ID
        scaffold_to_cluster = {
            scaffold: i for i, scaffold in enumerate(unique_scaffolds)
        }

        ### Map each scaffold to its cluster ID
        df["BM_cluster"] = df["BM_smiles"].map(scaffold_to_cluster)

    ### Drop SMILES column if not needed
    if not return_smiles:
        df.drop(columns=["BM_smiles"], inplace=True)

    return df


#############################################################################################
### Utility Functions: Butina Clustering


def butina_clustering(
    df,
    distance_matrix,
    threshold=0.25,
    analyze=True,
    inspect_n_members=5,
):
    """
    Butina clustering using rdkit.ML.Cluster.Butina.

    Parameters:
    -----------
    df : pd.DataFrame
        The dataset containing molecular structures as mol objects.

    distance_matrix : np.array
        The distance matrix to cluster. Needs to be a symmetrical matrix with distances
        (e.g. Jaccard, Euclidean) not similarities (e.g. Tanimoto).
        Can be calculated using MorganSimilarity.get_tanimoto_similarity_matrix().

    threshold : float
        The threshold for clustering. Default is 0.25.

    analyze : bool
        Whether to analyze the clusters by calling analyze_clusters(). Default is True.

    inspect_n_members : int
        The number of members to inspect in each cluster. Default is 5.

    Returns:
    --------
    df : pd.DataFrame
        The input DataFrame with additional columns (one compound per row):
        - Cluster_ID: The cluster ID assigned to each molecule. Cluster are inherently
            sorted by cluster size.
        - is_centroid: Boolean indicating if the molecule is a centroid of a cluster.
        - is_singleton: Boolean indicating if the molecule is a singleton.

    cluster_df : pd.DataFrame
        A DataFrame containing the clusters with the following columns:
        - Size: The size of the cluster.
        - Centroid (Mol objects): The centroid of the cluster.
        - Member_1, Member_2, ..., Member_n (Mol objects): The members of the cluster.
            The number of members is defined by the inspect_n_members parameter.

    """
    ### Copy the DataFrame
    df = df.copy()

    ### Create additional columns
    df["Cluster_ID"] = -1
    df["is_centroid"] = False
    df["is_singleton"] = False

    ### Extract upper triangle without diagonal as vector
    distances = distance_matrix[np.triu_indices(len(distance_matrix), k=1)]

    ### Clustering as a tuple of tuples (the first element is the centroid)
    clusters = Butina.ClusterData(
        data=distances, nPts=len(df), distThresh=threshold, isDistData=True
    )

    ### Instantiate DataFrame for clusters
    member_cols = [f"Member_{i + 1}" for i in range(inspect_n_members)]
    columns = ["Size", "Centroid"] + member_cols
    cluster_df = pd.DataFrame(index=range(len(clusters)), columns=columns)

    ### Loop through the clusters
    for i, cluster in enumerate(clusters):
        ### Get indexes
        members = list(cluster)
        centroid_idx = members[0]
        member_idxs = members[1:]

        ### ### df ### ###

        ### Assign cluster ID to all cluster members
        df.loc[df.index[members], "Cluster_ID"] = i

        ### Assign centroid tags (first cluster member)
        df.loc[df.index[centroid_idx], "is_centroid"] = True

        ### Check if the cluster is a singleton
        if len(members) == 1:
            df.loc[df.index[centroid_idx], "is_singleton"] = True

        ### ### cluster_df ### ###

        ### Cluster size
        cluster_df.at[i, "Size"] = len(members)

        ### Populate centroid
        cluster_df.at[i, "Centroid"] = df.loc[df.index[centroid_idx], "Mol"]

        ### Select random members for inspection
        if len(members) > inspect_n_members:
            selected = np.random.choice(member_idxs, inspect_n_members, replace=False)

        ### None padding for smaller clusters
        else:
            selected = member_idxs + [None] * (inspect_n_members - len(member_idxs))

        ### Populate members
        cluster_df.loc[i, member_cols] = [
            df.loc[df.index[idx], "Mol"] if idx is not None else None
            for idx in selected
        ]

    ### Analyse the clusters
    if analyze:
        analyze_clusters(df, cluster_df["Size"])

    return df, cluster_df


def analyze_clusters(
    df,
    clust_count,
    plot_histogram=True,
    count_cutoff_1=5,
    count_cutoff_2=20,
):
    """
    Helper function to analyze Butina clusters.

    """
    ### Copy DataFrame
    df = df.copy()

    ### Cluster Statistics
    # clust_count = df['Cluster_ID'].value_counts().sort_index()

    ### Print statistics
    print(f"{len(clust_count)} clusters.")
    print(f"{df['is_singleton'].sum()} sigletons.")
    print(f"{clust_count[0]} compounds in the largest cluster.")
    print(
        f"Number of clusters with more than {count_cutoff_1} compound: {clust_count[clust_count > count_cutoff_1].count()}"
    )
    print(
        f"Number of clusters with more than {count_cutoff_2} compound: {clust_count[clust_count > count_cutoff_2].count()}"
    )

    ### Countplot
    if plot_histogram:
        plt.figure(figsize=(6, 3))
        plt.hist(
            clust_count.values,
            bins=range(1, clust_count.values.max() + 2),
            align="left",
        )
        plt.xlabel("Cluster Size")
        plt.ylabel("Number of Clusters")
        plt.title("Distribution of Cluster Sizes")
        plt.yscale("log")
        plt.tight_layout()
        plt.show()


#############################################################################################
### Utility Class: MorganSimilarity


class MorganSimilarity:
    """
    Class to perform similarity operations using Morgan fingerprints.

    Attributes:
    -----------
        df : pd.DataFrame
            The dataset containing molecular structures as mol objects.

        act_col : str
            Column name for activity values.

        mol_col : str
                Column name for molecules.

        radius : int
                        Radius for Morgan fingerprint generation.

        nBits : int
                        Bit length of Morgan fingerprints.

        useChirality : bool
                        Whether to consider chirality for Morgan fingerprint generation.

        fp : list
                        List of Morgan fingerprints as Bit Vectors. Can be formatted as a NumPy array or DataFrame.

        bit_info : list
                        List of dictionaries containing bit information for inspection.

        SALI_dict : dict
                        Dictionary containing NumPy arrays of the same shape describing pairwise
            - distance (float)
            - activity values (float; e.g. pIC50) - name depends on the nature of the original dataset
            - delta activity values (float; e.g. ΔpIC50) - name depends on the nature of the original dataset
            - SALI (float): Similarity Activity Landscape Index
            - Positional indices of the k-nearest neighbors (int)
            - Positional indices of the molecules as integers (int)

        info : dict
                        Dictionary containing class metadata.

        dist_metric : str
                        The distance metric used for the k-nearest neighbors.

        pair_analysis : pd.DataFrame
                        DataFrame containing molecule pairs with SALI values above the chosen threshold.
    """

    def __init__(
        self,
        df: pd.DataFrame,
        act_col="pIC50",
        mol_col="Mol",
        radius=2,
        nBits=1024,
        useChirality=True,
    ):
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
            bitInfo=bit_info,
        )
        return fp, bit_info

    def get_morgan_fingerprint(self):
        """
        Create Morgan fingerprints for the classes DataFrame containing a Mol column.

        Returns:
        --------
            A list of Morgan fingerprints as Bit Vectors.

        Side effects:
        -------------
        self.bit_info : list
            A list of dictionaries containing bit information for inspection.

        self.fp : list
            A list of Morgan fingerprints as Bit Vectors.

        self.info : dict
            is updated.
        """

        df = self.df

        ## Process each molecule in the DataFrame in a vectorized manner
        fps_bit_info = df["Mol"].apply(self.__get_fp__)

        ### Unzip the results into separate lists
        fingerprints, bit_info_list = zip(*fps_bit_info)

        ### Update the class attributes
        self.bit_info = list(bit_info_list)
        self.fp = list(fingerprints)
        self.info["fp_format"] = "List of BitVectors"
        self.info["radius"] = self.radius
        self.info["nBits"] = self.nBits
        self.info["useChirality"] = self.useChirality

        return list(fingerprints)

    def get_tanimoto_similarity_matrix(self, return_jaccard_distances=False):
        """
        Compute the full Tanimoto similarity matrix for the Morgan fingerprints.

        Parameters:
        -----------
        return_jaccard_distance : bool
            Whether to return the Jaccard distance matrix instead of the similarity matrix.

        Returns:
        --------
        np.array
            The Tanimoto similarity matrix or Jaccard distance matrix.

        Raises:
        -------
        ValueError
            If the Morgan fingerprints are not in a list format.

        """
        ### Check fingerprint format
        if type(self.fp) is not list:
            raise ValueError(
                "Morgan fingerprints must be in a list format. Please run get_morgan_fingerprint() first."
            )

        ### Number of compounds
        n_fps = len(self.fp)

        ### Instantiate the Tanimoto similarity matrix
        similarity_matrix = np.zeros((n_fps, n_fps))

        ### Compute the Tanimoto similarity matrix
        for i in range(n_fps):
            similarity_matrix[i, :] = DataStructs.BulkTanimotoSimilarity(
                self.fp[i], self.fp
            )

        ### Return the Jaccard distance matrix
        if return_jaccard_distances:
            return 1 - similarity_matrix

        ### Return the Tanimoto similarity matrix
        else:
            return similarity_matrix

    def format_fp(self, format="boolean_np_array"):
        """
        Reformat the Morgan fingerprints as a boolean NumPy array, DataFrame, or NumPy array.

        Parameters:
        -----------
        format : str
                        The format to return the fingerprints in:
            - 'boolean_np_array'
            - 'DataFrame'
            - 'np_array'

        Returns:
        --------
        np.array or pd.DataFrame
            The formatted Morgan fingerprints in the specified format.


        Side effects:
        -------------
        self.fp : np.array or pd.DataFrame
            The fingerprints are updated in the specified format.

        self.info['fp_format'] : str
            is updated.

        """

        ### Change the format of the Morgan fingerprints
        if format == "boolean_np_array":
            fp = np.array(self.fp).astype(bool)

        elif format == "DataFrame":
            fp = pd.DataFrame(np.array(self.fp), index=self.df.index)

        elif format == "np_array":
            fp = np.array(self.fp)

        else:
            raise ValueError(
                "Invalid format. Please use 'boolean_np_array' (default), 'np_array' or 'DataFrame'."
            )

        ### Update the class attributes
        self.fp = fp
        self.info["fp_format"] = format

        return fp

    def get_SALI_matrix(self, k=10, dist_metric="jaccard", stat_summary=False):
        """
        Compute the Similarity Activity Landscape Indexes (SALI) for the k-nearest neighbors.

        Parameters:
        -----------
        k : int
                        The number of nearest neighbors to consider.

        dist_metric : str
                        The distance metric to use for the k-nearest neighbors.

        Returns:
        --------
        np.array
            SALI values for each molecule's k-nearest neighbors.

        Side effects:
        -------------
        self.SALI_dict : dict
            Containing NumPy arrays of identical shape:
            - distance (float)
            - activity values (float; e.g. pIC50) - name depends on the nature of the original dataset
            - delta activity values (float; e.g. ΔpIC50) - name depends on the nature of the original dataset
            - SALI (float): Similarity Activity Landscape Index
            - Positional indices of the k-nearest neighbors (int)
            - Positional indices of the molecules as integers (int)
        """

        ### Instantiate the KNN model
        neigh = NearestNeighbors(n_neighbors=k + 1, metric=dist_metric)

        ### Calculate the distances and indices of the k-nearest neighbors
        neigh = neigh.fit(self.fp)
        distances_np, neigb_idx = neigh.kneighbors(self.fp, return_distance=True)

        ### Removing the first column (distance to each molecule itself)
        distances_np = distances_np[:, 1:]
        neigb_idx = neigb_idx[:, 1:]

        ### Convert pIC50 column to NumPy array for fast indexing
        pIC50_np = self.df[self.act_col].values

        ### Compute the absolute differences for each molecule's k-nearest neighbors
        delta_pIC50_np = np.abs(pIC50_np.reshape(-1, 1) - pIC50_np[neigb_idx])

        ### Compute the SALI matrix
        SALI_np = delta_pIC50_np / distances_np

        ### Distribution analysis
        if stat_summary:
            SALI_df = pd.Series(SALI_np.flatten())
            print(f"{len(SALI_df):.2e} pairs analyzed (including duplicates).")
            print(f"{math.comb(len(self.df), 2):.2e} pairs are possible.")

            print("Quantiles (SALI):")
            print(SALI_df.quantile([0, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999]))

        ### Summarizing all NumPy arrays in one dictionary
        SALI_dict = {
            "distance": distances_np,
            self.act_col: pIC50_np,
            f"delta_{self.act_col}": delta_pIC50_np,
            "SALI": SALI_np,
            "neigb_idx": neigb_idx,
            "row_idx": np.array(
                [list(range(neigb_idx.shape[0])) for i in range(neigb_idx.shape[1])]
            ).T,
        }

        ### Update the class attributes
        self.SALI_dict = SALI_dict
        self.dist_metric = dist_metric
        self.info["k"] = k
        self.info["dist_metric"] = dist_metric

        return SALI_np

    def get_sali_pairs(self, SALI_threshold, min_activity_diff=None):
        """
        List molecule pairs with SALI values higher than the threshold.

        Parameters:
        -----------
        SALI_threshold : int or float
            The threshold for the SALI value.

        min_activity_diff : None, int, or float
            A minimum difference in activity can be defined.

        Returns:
        --------
        pd.DataFrame
            A DataFrame containing the molecule pairs with high SALI values.

        Side effects:
        -------------
        self.pair_analysis : pd.DataFrame
            A DataFrame containing the molecule pairs with high SALI values.
        """

        ### Unpack the SALI dictionary
        distances_np = self.SALI_dict["distance"]
        pIC50_np = self.SALI_dict[f"{self.act_col}"]
        delta_pIC50_np = self.SALI_dict[f"delta_{self.act_col}"]
        SALI_np = self.SALI_dict["SALI"]
        neigb_idx = self.SALI_dict["neigb_idx"]
        row_idx = self.SALI_dict["row_idx"]

        ### Create a map for boolean indexing of the NumPy arrays
        SALI_mask = SALI_np > SALI_threshold

        ### Formatting the columns for the first molecule of each pair
        row_df = self.df.reset_index().iloc[row_idx[SALI_mask]]
        row_df.rename(
            columns={name: name + "_1" for name in row_df.columns}, inplace=True
        )
        row_df.reset_index(inplace=True, drop=True)

        ### Formatting the columns for the second molecule
        neigb_df = self.df.reset_index().iloc[neigb_idx[SALI_mask]]
        neigb_df.rename(
            columns={name: name + "_2" for name in neigb_df.columns}, inplace=True
        )
        neigb_df.reset_index(inplace=True, drop=True)

        ### Combining all columns
        inspect_df = pd.concat([row_df, neigb_df], axis=1)

        ### Adding the comparison metrics
        inspect_df[self.dist_metric] = distances_np[SALI_mask]
        inspect_df[f"Δ{self.act_col}"] = delta_pIC50_np[SALI_mask]
        inspect_df["SALI"] = SALI_np[SALI_mask]

        ### Create a new Series where each element is a tuple of the sorted IDs
        sorted_pairs = inspect_df.apply(
            lambda row: tuple(sorted((row["ID_1"], row["ID_2"]))), axis=1
        )

        ### Generate a boolean mask that is True for any duplicated pair
        cross_dupl_mask = sorted_pairs.duplicated(keep="first")

        ### Remove duplicates
        inspect_df = inspect_df.loc[cross_dupl_mask]

        ### Filter by minimum activity difference
        if min_activity_diff:
            inspect_df = inspect_df.loc[
                inspect_df[f"Δ{self.act_col}"] > min_activity_diff
            ]

        ### Order by SALI values and reset the index
        inspect_df = inspect_df.sort_values("SALI", ascending=False)
        inspect_df.reset_index(drop=True, inplace=True)

        ### Update the class attribute
        self.pair_analysis = inspect_df
        self.info["SALI_threshold"] = SALI_threshold

        return inspect_df

    def count_cliffs(self, SALI_threshold, min_activity_diff=None):
        """
        Count the number of activity cliffs for each molecule.

        Parameters:
        -----------
        SALI_threshold : int or float
                        The threshold for the SALI value.

        min_activity_diff : None, int, or float
                        A minimum difference in activity can be defined.

        Returns:
        --------
        pd.DataFrame
            A DataFrame containing the molecules sorted by the number of cliffs.
        """

        ### Unpack the SALI dictionary
        distances_np = self.SALI_dict["distance"]
        pIC50_np = self.SALI_dict[self.act_col]
        delta_pIC50_np = self.SALI_dict[f"delta_{self.act_col}"]
        SALI_np = self.SALI_dict["SALI"]
        neigb_idx = self.SALI_dict["neigb_idx"]
        row_idx = self.SALI_dict["row_idx"]

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
        inspect_df["CliffCount"] = cliff_count

        return inspect_df.sort_values("CliffCount", ascending=False)
