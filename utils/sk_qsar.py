import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib_venn import venn3
from sklearn.manifold import TSNE

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
    Class to determine the applicability domain of a Random Forest Classifier.

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
