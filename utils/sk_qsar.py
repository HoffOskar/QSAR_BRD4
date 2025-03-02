import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.manifold import TSNE

###########################################
### Visualization functions


def plot_tsne(
    X_train, X_test, y_train, y_test, show=True, SEED=0, learning_rate="auto"
):
    """ """

    ### Concatenate the train and test sets
    X_combined = np.concatenate((X_train, X_test), axis=0)
    y_combined = np.concatenate((y_train, y_test), axis=0)

    ### Apply TSNE to the combined dataset
    tsne = TSNE(n_components=2, random_state=SEED, learning_rate=learning_rate)
    X_tsne = tsne.fit_transform(X_combined)

    ### Split the transformed data back into train and test sets
    X_tsne_train = X_tsne[: len(X_train)]
    X_tsne_test = X_tsne[len(X_train) :]

    ### Plot the TSNE results
    sns.scatterplot(
        x=X_tsne_train[:, 0], y=X_tsne_train[:, 1], label="Train", alpha=0.9
    )
    sns.scatterplot(x=X_tsne_test[:, 0], y=X_tsne_test[:, 1], label="Test", alpha=0.9)
    plt.legend()
    plt.title("t-SNE plot of Train and Test sets")
    plt.xlabel("TSNE Component 1")
    plt.ylabel("TSNE Component 2")

    if show:
        plt.show()


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

    def __init__(self, rf_model, df_test, threshold=None, mol_col="Mol", pca_model=None):
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

    def predict(self, threshold=None):
        """
        Determine if the test compounds are within the applicability domain.

        Parameters:
        ----------
        threshold: float or None
            Threshold for the applicability domain. If None, the class attribute is used.

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

        ### Add the maximum agreement value
        AD_df["AD_value"] = self.neigb_dist_np

        ### Print statistics
        print("Threshold:", self.threshold)
        print(
            f"{len(AD_df) - AD_df['AD_status'].sum()} ({100 * (1 - AD_df['AD_status'].mean())} %) test cmpds outside the applicability domain"
        )
        print(AD_df["AD_value"].describe())
        return AD_df
