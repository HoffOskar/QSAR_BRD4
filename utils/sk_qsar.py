import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib_venn import venn3
from sklearn.manifold import TSNE
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm

###########################################
### Visualization functions


def plot_tsne(
    X_0,
    X_1,
    label_0="Train",
    label_1="Test",
    show=True,
    SEED=0,
    learning_rate="auto",
    perplexity=30,
):
    """
    Plot the t-SNE of the train and test sets.

    Parameters:
    ----------
    X_0: np.array (2D)
        Training data.
    X_1: np.array (2D)
        Test data.
    label_0: str
        Label for the training set.
    label_1: str
        Label for the test set.
    show: bool
        Show the plot.
    SEED: int
        Random seed.
    learning_rate: float or 'auto'
        Learning rate for the t-SNE algorithm.
    perplexity: int
        Perplexity for the t-SNE algorithm
    """

    ### Concatenate the train and test sets
    X_combined = np.concatenate((X_0, X_1), axis=0)

    ### Apply TSNE to the combined dataset
    tsne = TSNE(
        n_components=2,
        random_state=SEED,
        learning_rate=learning_rate,
        perplexity=perplexity,
    )
    X_tsne = tsne.fit_transform(X_combined)

    ### Split the transformed data back into train and test sets
    X_tsne_0 = X_tsne[: len(X_0)]
    X_tsne_1 = X_tsne[len(X_0) :]

    ### Plot the TSNE results
    sns.scatterplot(x=X_tsne_0[:, 0], y=X_tsne_0[:, 1], label=label_0, alpha=0.9)
    sns.scatterplot(x=X_tsne_1[:, 0], y=X_tsne_1[:, 1], label=label_1, alpha=0.9)
    plt.legend()
    plt.title(f"t-SNE plot of {label_0} and {label_1} data")
    plt.xlabel("TSNE Component 1")
    plt.ylabel("TSNE Component 2")

    if show:
        plt.show()


def plot_venn3(df, threshold_dict, labels):
    """
    Function to plot a Venn diagram based on numeric or categorical filters.

    Parameters:
    ----------
    df : pd.DataFrame
        DataFrame containing the columns to filter.
    threshold_dict : dict
        Dictionary where keys are column names and values are `[lower, upper]` thresholds.
    labels : list
        List containing three labels corresponding to `threshold_dict.keys()`.

    Returns:
    -------
    full_mask : np.array
        Boolean mask of compounds that meet all three criteria.
    """

    ### Validate that dictionary keys match DataFrame columns
    dict_keys = list(threshold_dict.keys())
    for key in dict_keys:
        if key not in df.columns:
            raise ValueError(f"Error: Column '{key}' not found in DataFrame!")

    ### Extract column names dynamically
    col_A, col_B, col_C = dict_keys  # Extract the three column names

    ### Create Boolean masks dynamically
    A_mask = (df[col_A] >= threshold_dict[col_A][0]) & (
        df[col_A] <= threshold_dict[col_A][1]
    )
    B_mask = (df[col_B] >= threshold_dict[col_B][0]) & (
        df[col_B] <= threshold_dict[col_B][1]
    )
    C_mask = (df[col_C] >= threshold_dict[col_C][0]) & (
        df[col_C] <= threshold_dict[col_C][1]
    )

    ### Define sets for Venn diagram (counts of each condition)
    A_only = (A_mask & ~B_mask & ~C_mask).sum()
    B_only = (~A_mask & B_mask & ~C_mask).sum()
    C_only = (~A_mask & ~B_mask & C_mask).sum()

    A_B = (A_mask & B_mask & ~C_mask).sum()
    A_C = (A_mask & ~B_mask & C_mask).sum()
    B_C = (~A_mask & B_mask & C_mask).sum()

    A_B_C = (A_mask & B_mask & C_mask).sum()

    ### Calculate the total number of compounds in the Venn diagram
    total_cmpds = A_only + B_only + C_only + A_B + A_C + B_C + A_B_C

    ### Aggregate masks
    masks = [A_mask, B_mask, C_mask]

    ### Use Seaborn's color palette
    colors = sns.color_palette("deep")[:3]

    ### Create Venn diagram
    plt.figure(figsize=(6, 6))
    venn = venn3(
        subsets=(A_only, B_only, A_B, C_only, A_C, B_C, A_B_C),
        set_labels=labels,
        set_colors=colors,
    )

    ### Change label colors to match fill colors
    for i, label in enumerate(venn.set_labels):
        label.set_color(colors[i])  # Match label text to set color
        label.set_fontsize(12)  # Adjust label size

    ### Add threshold values dynamically to the title
    title = (
        f"Venn Diagram of {total_cmpds} ({100 * total_cmpds / len(df):.0f}%) cmpds:\n"
    )
    for i in range(3):
        title += f"{labels[i]}: {threshold_dict[dict_keys[i]][0]} - {threshold_dict[dict_keys[i]][1]} ({100 * masks[i].sum() / len(df):.0f}%)\n"

    plt.title(title, fontsize=14)
    plt.show()

    ### Return full mask
    return np.array(A_mask & B_mask & C_mask)


############################################
### Classes


class RfClfAppDom:
    """
    Class to determine the applicability domain of a Random Forest Classifier based on tree agreement.

    Attributes:
    ----------
    rf_model : RandomForestClassifier
        Random Forest Classifier model.
    df_test : pd.DataFrame
        Test set with unique .index and RDKit.Mol column
    threshold : float or None
        Threshold for the applicability domain
    mol_col : str
        Column name for the RDKit.Mol objects in df_test
    pca_model : PCA or None
        PCA model used for feature transformation
    neigb_sim_np : np.array (2D) or None
        pairwise agreement between train and test compounds
    neigb_dist_np : np.array (1D) or None
        maximum agreement between test compounds to the closest train compound
    X_train : np.array (2D) or None
        Training data
    """

    def __init__(
        self, rf_model, df_test, threshold=None, mol_col="Mol", pca_model=None
    ):
        self.rf_model = rf_model
        self.n_estimators = rf_model.n_estimators
        self.estimators_np = np.array(rf_model.estimators_)
        self.df_test = df_test
        self.threshold = threshold
        self.mol_col = mol_col
        self.pca_model = pca_model
        self.neigb_sim_np = None
        self.neigb_dist_np = None
        self.X_train = None

    def fit(self, X_train, X_test):
        """
        Get the maximum agreement between the training set to the respective test compound.

        Parameters:
        ----------
        X_train: np.array (2D)
            Training data.
        X_test: np.array (1D)
            Test data.

        Side-effects - updates the class attributes:
        -------
        neigb_sim_np: np.array (2D)
            pairwise agreement between train and test compounds
        neigb_dist_np: np.array (1D)
            maximum agreement between test compounds to the closest train compound
        X_train: np.array (2D)
            Training data
        """
        ### Vectorized predictions of the individual trees - shape: (n_samples, n_trees)
        Train_array = np.stack(
            [tree.predict(X_train) for tree in self.estimators_np], axis=1
        )
        Test_array = np.stack(
            [tree.predict(X_test) for tree in self.estimators_np], axis=1
        )

        ### Compute pairwise agreement between Train_array and Test_array
        neigb_sim_np = (
            np.sum(Train_array[None, :, :] == Test_array[:, None, :], axis=2)
            / self.n_estimators
        )

        ### Compute the maximum agreement
        neigb_dist_np = np.max(neigb_sim_np, axis=1)

        ### Update class attributes
        self.neigb_sim_np = neigb_sim_np
        self.neigb_dist_np = neigb_dist_np
        self.X_train = X_train

    def def_threshold(self, threshold):
        """
        (Re)define the threshold for the applicability domain.
        """
        self.threshold = threshold

    def predict(self, threshold=None, print_stats=True):
        """
        Determine if the test compounds are within the applicability domain.

        Parameters:
        ----------
        threshold: float or None
            Threshold for the applicability domain. If None, the class attribute is used.
        print_stats: bool
            Print statistics about the applicability domain.

        Returns:
        -------
        AD_df: pd.DataFrame
            DataFrame of the test compounds with the AD status, the maximum agreement value, index and Mol objects.
        """

        ### Update threshold if provided
        if threshold is not None:
            self.threshold = threshold

        ### Check if the threshold is defined
        if self.threshold is None:
            print("No threshold saved or provided.")
            return

        ### Get indexes and mol objects
        AD_df = self.df_test[[self.mol_col]].copy()

        ### Boolean array for the applicability domain
        AD_df["AD_status"] = self.neigb_dist_np > self.threshold
        AD_df["AD_status"] = AD_df["AD_status"].astype(bool)

        ### Add the maximum agreement value
        AD_df["AD_value"] = self.neigb_dist_np
        AD_df["AD_value"] = AD_df["AD_value"].astype(float)

        ### Print statistics
        if print_stats:
            print("Threshold:", self.threshold)
            print(
                f"{len(AD_df) - AD_df['AD_status'].sum()} ({100 * (1 - AD_df['AD_status'].mean())} %) test cmpds outside the applicability domain"
            )
            print(AD_df["AD_value"].describe())

        return AD_df


class KNeighborAppDom:
    """
    Class to determine the applicability domain based on the k-NN distances.

    Attributes:
    ----------
    k : int or None
        Number of nearest neighbors to use.
    X_train_df : pd.DataFrame or None
        Training data.
    RefVal : float or None
        Reference value for the applicability domain.
    AD_thresholds : np.array or None
        AD thresholds for the training data.
    neigh : NearestNeighbors or None
        NearestNeighbors object for the training data.
    """

    def __init__(self, k=None):
        self.k = k
        self.X_train_df = None
        self.RefVal = None
        self.AD_thresholds = None
        self.neigh = None

    def explore_k(self, X_train_df, k_values, plot_ref_vals=False):
        """
        Plot k against
            - the count of training compounds with AD radius = 0 (up)
            - the count of training compounds outside AD ('quasi leave one out').

        Parameters:
        ----------
        X_train : pd.DataFrame
            Training data.
        k_values : list
            List of k values to evaluate.
        plot_ref_vals : bool
            Plot the Ref_Val against k.
        """

        ### train_1: full distances matrix from helper function
        Train_neigb_distances = self._compute_distance_matrix(X_train_df)

        ### Container
        no_ADs = []
        out_of_ADs = []
        Ref_Vals = []

        for k in tqdm(k_values):
            ### Termine the AD thresholds and global Ref_Val for the training data
            AD_thresholds, Ref_Val = self._compute_AD_thresholds(
                Train_neigb_distances, k
            )

            ### Count the number of training points with AD_threshold = 0
            no_ADs.append(np.sum(AD_thresholds == 0))

            ### Count the number of training points within AD for each training point
            AD_count = np.sum(Train_neigb_distances <= AD_thresholds, axis=1)

            ### Count the number of training points one (themselves) or less neighbors within the AD_thresholds
            out_of_ADs.append(np.sum(AD_count <= 1))

            ### Save the Ref_Val
            Ref_Vals.append(Ref_Val)

        ### Invert the values of out_of_ADs for plotting
        inverted_out_of_ADs = [-val for val in out_of_ADs]

        ### Count plot
        plt.figure(figsize=(10, 6))

        ### Plot no_ADs going up
        sns.barplot(x=k_values, y=no_ADs, label="Cmpds with AD radius = 0")

        ### Plot out_of_ADs going down
        sns.barplot(
            x=k_values,
            y=inverted_out_of_ADs,
            color="orange",
            label='Out of AD ("quasi leave-one out")',
        )

        plt.xlabel("k")
        plt.ylabel("Count")
        plt.legend()
        plt.grid(axis="y")
        plt.show()

        ### Plot Ref_Vals
        if plot_ref_vals:
            plt.figure(figsize=(10, 6))
            plt.plot(k_values, Ref_Vals, marker="o", color="red")
            plt.xlabel("k")
            plt.ylabel("Ref_Val")
            plt.grid()
            plt.show()

    def fit(self, X_train_df, k):
        """
        Train the applicability domain model using the chosen k.

        Parameters:
        ----------
        X_train : pd.DataFrame
            Training data.
        k : int
            Number of nearest neighbors to use.
        """

        ### train_1: full distances matrix from helper function
        Train_neigb_distances = self._compute_distance_matrix(X_train_df)

        ### train_2, 3 and 4 from helper function
        AD_thresholds, RefVal = self._compute_AD_thresholds(Train_neigb_distances, k)

        ### Save object attributes
        self.k = k
        self.X_train_df = X_train_df
        self.RefVal = RefVal
        self.AD_thresholds = AD_thresholds

    def predict(self, X_test_df, return_neighbors=False):
        """
        Classify new compounds to check if they are within the applicability domain.

        Parameters:
        ----------
        X_test_df : pd.DataFrame
            Test data to classify.

        Returns:
        -------
        AD_df : pd.DataFrame
            DataFrame with the AD status and the number of training compounds within the AD.
        """

        ### Check if the model is trained
        if not hasattr(self, "neigh"):
            raise ValueError("Model is not trained yet. Call `train()` first.")

        ### query 1: Calculate a distances matrix between test (rows) and training (columns) cmpds
        Test_neigb_distances, Test_neigb_idx = self.neigh.kneighbors(
            X_test_df, return_distance=True
        )
        print(
            "Shape of distance matrix:",
            Test_neigb_distances.shape,
            "(test x training cmpds)",
        )

        ### query 2: Count the number of training points within the AD_thresholds
        within_AD_mask = Test_neigb_distances <= self.AD_thresholds
        counts = np.sum(within_AD_mask, axis=1)

        ### query 3: Calculate the AD for each test point
        AD = counts != 0

        ### Summarize the results in a DataFrame
        AD_df = pd.DataFrame(
            {"AD_status": AD, "AD_density": counts}, index=X_test_df.index
        )

        ### Return neighbors, if requested
        if return_neighbors:
            neighbors_within_AD = [
                self.X_train_df.index[Test_neigb_idx[i][within_AD_mask[i]]].tolist()
                for i in range(Test_neigb_idx.shape[0])
            ]
            AD_df["train_neigbs_idx"] = neighbors_within_AD

        print(AD.sum(), f"({round(AD.mean() * 100)}%) compounds are within AD")

        return AD_df

    def _compute_distance_matrix(self, X_train_df):
        """
        Helper function to compute and store the k-NN distance matrix for the training data.

        Parameters:
        ----------
        X_train_df : pd.DataFrame
            Training data.

        Returns:
        -------
        Train_neigb_distances : np.array
            The full distance matrix for the training data.
        """

        ### train_1:
        ### Calculating the full distances matrix for the training set
        neigh = NearestNeighbors(n_neighbors=X_train_df.shape[0])
        neigh.fit(X_train_df)
        Train_neigb_distances, _ = neigh.kneighbors(X_train_df, return_distance=True)

        ### Save the NearestNeighbors object
        self.neigh = neigh

        return Train_neigb_distances

    def _compute_AD_thresholds(self, Train_neigb_distances, k):
        """
        Internal helper function to compute AD metrics using a precomputed distance matrix.

        Parameters:
        ----------
        Train_neigb_distances : np.array
            Precomputed k-NN distance matrix.
        k : int
            Number of nearest neighbors.

        Returns:
        -------
        AD thresholds : np.array
            AD thresholds for the training data.
        RefVal : float
            Reference value for the applicability domain
        """

        ### Slicing for the k nearest neighbors, excluding the data point itself
        Train_k_neigb_distances = Train_neigb_distances[:, 1 : k + 1]

        ### train_2:
        ### Calculate the mean distances to k neighbors for each training cmpd
        mean_Train_k_neigb_distances = Train_k_neigb_distances.mean(axis=1)

        ### train_3:
        ### Tukey’s outlier method for the distances to k neighbors in the training set
        Q3 = np.quantile(mean_Train_k_neigb_distances, 0.75)
        Q1 = np.quantile(mean_Train_k_neigb_distances, 0.25)
        RefVal = Q3 + 1.5 * (Q3 - Q1)

        ### train_4:
        ### For each training cmpd subset to all distances among the training set smaller/equal to RefVal
        subsetted_distances = [row[row <= RefVal] for row in Train_neigb_distances]

        ### For each training cmpd calculate the average distance to all neighbors within RefVal
        AD_thresholds = np.array([array.mean() for array in subsetted_distances])

        return AD_thresholds, RefVal
