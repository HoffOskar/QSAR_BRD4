### Imports
import copy
import math

import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools, rdChemReactions, rdCoordGen
from rdkit.Chem.Draw import rdDepictor
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from tqdm import tqdm

### Settings
### tqdm for pandas
tqdm.pandas()

####################################################################################################
### Visualization


def render_mol(df):
    """
    Renders the "Mol" column of a DataFrame for visual inspection of chemical structures.
    """
    PandasTools.ChangeMoleculeRendering(df)
    return df


def to_clipboard_smiles(df):
    """
    Loads a DataFrame with a "Mol" column to_clipboard
    """
    # Create a copy of the input DataFrame to avoid changes to the original object
    df_copy = copy.deepcopy(df)

    # Generate SMILES for each molecule in the 'Mol' column
    df_copy["Smiles"] = df_copy["Mol"].apply(lambda mol: Chem.MolToSmiles(mol))

    # Drop the 'Mol' column
    df_copy.drop(columns=["Mol"], inplace=True)

    # Copy the modified DataFrame to clipboard
    df_copy.to_clipboard(index=True)

    return


def grid_view(mol_col):
    """
    Accepts a pandas.Series of Mol objects ("Mol" column) and displays all structures in a grid
    """
    return Draw.MolsToGridImage(
        mol_col,
        molsPerRow=5,  # Adjust this number to control how many molecules appear per row
        subImgSize=(400, 400),  # Optional: adjust the size of each sub-image
        useSVG=True,  # Optional: use SVG format for better quality in Jupyter notebooks
    )


####################################################################################################
### Calculation


def convert_ic50nM_to_pic50(IC50_value):
    """
    Calculate pIC50 from IC50 values in nM.
    """
    pIC50_value = 9 - math.log10(IC50_value)
    return pIC50_value


####################################################################################################
### Standardization


def safe_apply(mol_func):
    """
    Wrap an RDKit Mol function to skip non-Mol entries and return None.
    """
    return lambda mol: mol_func(mol) if isinstance(mol, Chem.Mol) else None


def standardize_mol_col(
    mol_col,
    cleanup=False,
    sanitize=True,
    normalize=True,
    largest_frag=True,
    remove_charge=True,
    tautomerize=False,
    remove_stereo=False,
    cal_2D_coord_default=True,
    cal_2D_coord_deep=False,
):
    """
    Standardizes a pandas Series of RDKit Mol objects by applying a sequence of standardization operations.

    Parameters:
    -----------
    mol_col : pandas.Series
        A Series or DataFrame column containing RDKit Mol objects or None.

    cleanup : bool
        Cleanup structure.

    sanitize : bool
        Sanitize structure. Faster alternative to cleanup. Overruled if cleanup is set to True.

    normalize : bool
        Canonicalize functional group representation.

    largest_frag : bool
        Keep only the largest fragment.

    remove_charge : bool
        Neutralize charged molecules using `rdMolStandardize.Uncharger().uncharge`.

    tautomerize : bool
        Convert molecules to their most populated tautomer.
        (This can be slow, so use only when necessary.)

    remove_stereo : bool
        Removes stereochemical information.

    cal_2D_coord_default : bool
        Computes 2D coordinates for molecules using `rdDepictor.Compute2DCoords`.
        (Important for combinatorial libraries.)

    cal_2D_coord_deep : bool
        Generate 2D coordinate, which can better handle specific structural features (e.g., straightening alkynes to 180°).
        Overrides `cal_2D_coord_default` if set to True.

    Returns:
    --------
    pandas.Series
        A Series or DataFrame column of RDKit Mol objects after applying the selected standardization steps.

    Notes:
    ------
    - None values in the input column are preserved.
    """

    ### Work with a copy of the input column to avoid changes to the original object
    mol_col = copy.deepcopy(mol_col)

    ### Clean up visualization
    if cleanup:
        print(f"cleanup={cleanup} (rdMolStandardize.Cleanup):")
        mol_col = mol_col.progress_apply(safe_apply(rdMolStandardize.Cleanup))
        print("sanitize=not required (Chem.SanitizeMol):")

    elif sanitize:
        print(f"sanitize={sanitize} (Chem.SanitizeMol):")
        mol_col = mol_col.progress_apply(
            safe_apply(lambda mol: Chem.SanitizeMol(mol) or mol)
        )

    else:
        print(f"cleanup={cleanup} (rdMolStandardize.Cleanup)")
        print(f"sanitize={sanitize} (Chem.SanitizeMol)")

    ### Normalize (not sure if necessary)
    print(f"normalize={normalize} (rdMolStandardize.Normalize):")
    if normalize:
        mol_col = mol_col.progress_apply(safe_apply(rdMolStandardize.Normalize))

    ### Keep the largest fragment only
    print(f"largest_frag={largest_frag} (rdMolStandardize.FragmentParent):")
    if largest_frag:
        mol_col = mol_col.progress_apply(safe_apply(rdMolStandardize.FragmentParent))

    ### Remove charges
    print(f"remove_charge={remove_charge} (rdMolStandardize.Uncharger().uncharge):")
    if remove_charge:
        mol_col = mol_col.progress_apply(
            safe_apply(rdMolStandardize.Uncharger().uncharge)
        )

    #### Generate the most populated tautomer # slow only use when necessary
    print(f"tautomerize={tautomerize} (rdMolStandardize.CanonicalTautomer)")
    if tautomerize:
        mol_col = mol_col.progress_apply(safe_apply(rdMolStandardize.CanonicalTautomer))

    ### Remove stereoinformation
    print(f"remove_stereo={remove_stereo} (rdMolStandardize.StereoParent):")
    if remove_stereo:
        mol_col = mol_col.progress_apply(safe_apply(rdMolStandardize.StereoParent))

    ### 2D-coordinates are only calculated once
    if cal_2D_coord_deep:
        print("cal_2D_coord_default=not required (rdDepictor.Compute2DCoords):")
        print(f"cal_2D_coord_deep={cal_2D_coord_deep} (rdCoordGen.AddCoords):")
        # alternative to fix 2D-coordinates # straightens alkynes to 180°
        mol_col.progress_apply(safe_apply(lambda mol: rdCoordGen.AddCoords(mol)))

    elif cal_2D_coord_default:
        print(
            f"cal_2D_coord_default={cal_2D_coord_default} (rdDepictor.Compute2DCoords):"
        )
        # default to fix 2D-coordinates # important for combinatorial libraries
        mol_col.progress_apply(safe_apply(lambda mol: rdDepictor.Compute2DCoords(mol)))
        print(f"cal_2D_coord_deep={cal_2D_coord_deep} (rdCoordGen.AddCoords):")
    else:
        print(
            f"cal_2D_coord_default={cal_2D_coord_default} (rdDepictor.Compute2DCoords):"
        )
        print(f"cal_2D_coord_deep={cal_2D_coord_deep} (rdCoordGen.AddCoords):")

    return mol_col


def canonical_smiles(mol_col):
    """
    Generate canonical SMILES for a pandas Series of RDKit Mol objects.
    """
    smiles_col = mol_col.apply(Chem.MolToSmiles)
    # acids_df['smiles'] = acids_df['smiles'].apply(lambda x: rdMolStandardize.StandardizeSmiles(x))

    return smiles_col.apply(rdMolStandardize.StandardizeSmiles)


### Generate standardized SMILES from Mol


####################################################################################################
### Substructure Matching


### Function to count how often a SMARTS pattern matches for a given mol object
def _count_smarts_matches(mol, smarts):
    """
    Return number of substructure matches for a given SMARTS pattern in a molecule.

    Parameters:
    -----------
    mol : RDKit Mol object
        Molecule to search for substructure matches

    smarts : str
        SMARTS pattern to match

    Returns:
    --------
    int
        Number of substructure matches
    """
    ### Convert SMARTS to a molecule pattern
    pattern = Chem.MolFromSmarts(smarts)

    ### If mol or pattern is None, return NaN
    if mol is None or pattern is None:
        return np.nan

    ### Find all matches of the SMARTS pattern in the molecule
    matches = mol.GetSubstructMatches(pattern)
    return len(matches)


###
def count_smarts_matches(df, smarts, mol_col="Mol"):
    """
    Count the number of SMARTS matches for each molecule in a DataFrame

    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing a column of RDKit Mol objects

    smarts : str
        SMARTS pattern to match

    mol_col : str
        Column name containing RDKit Mol objects

    Returns:
    --------
    pandas.Series
        A Series containing the number of substructure matches for each molecule
    """
    return df[mol_col].apply(lambda x: _count_smarts_matches(x, smarts))


####################################################################################################
### Virtual Reactions


def _rxn_2_sms(rxn_smarts, mol_1, mol_2):
    """
    Perform a virtual reaction on two molecules and return the products.

    Parameters:
    -----------

    mol_1 : RDKit Mol object
        First starting material

    mol_2 : RDKit Mol object
        Second starting material.

    rxn_smarts : str
        SMARTS rxn pattern.

    Returns:
    --------
    RDKit Mol object
        Reaction product.
        If multiple rxns are possible, only one product is returned.
        If no rxn is possible, None is returned.
    """
    ### Define the transformation
    rxn = rdChemReactions.ReactionFromSmarts(rxn_smarts)

    ### Starting materials
    reacts = (mol_1, mol_2)

    ### Perform the virtual transformation
    products = rxn.RunReactants(reacts)  # provide starting materials as a tuple

    ### Return all products of the first reaction
    try:
        return products[0]

    except:
        return (None, None)


def rxn_2_sms(
    df,
    rxn_smarts,
    mol_col_1="Mol",
    mol_col_2=None,
    reactant=Chem.MolFromSmiles("O"),
    drop_by_product=True,
    standardize_product=True,
):
    """
    React two molecule columns in a horizontal fashion create two product columns.

    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing the molecules to react.
        If only one column is provided, the second reactant is water.

    rxn_smarts : str
        SMARTS pattern for the reaction.
        Two reactants and two products are required.
        Mass balance for heavy atoms must be maintained.
        The first product is the main product.

    mol_col_1 : str
        Column name containing the main reactant (default is "Mol").

    mol_col_2 : str or None
        Column name containing the second reactant.
        Water is used as 2nd reactant if None is provided.
    
    reactant : RDKit Mol object
        The second reactant to be used if mol_col_2 is None.

    drop_by_product : bool
        If True, the by-product column is dropped.

    standardize_product : bool
        If True, the product column is standardized using `standardize_mol_col`.

    Returns:
    --------
    pandas.DataFrame
        DataFrame containing the original columns and the products of the reaction.
        If the starting material does not fit the substructure the products are None.
    """

    ### Apply the virtual reaction to each row of the DataFrame
    rxn_products = df.copy().progress_apply(
        lambda row: _rxn_2_sms(
            rxn_smarts,
            row[mol_col_1],
            row[mol_col_2] if mol_col_2 else reactant,  # water by default
        ),
        axis=1,
    )

    ### Unpack the Series of tuples into two separate columns
    df["Product"] = rxn_products.apply(lambda x: x[0] if x else None)


    ### Drop by-product column
    if not drop_by_product:
        df["by_product"] = rxn_products.apply(lambda x: x[1] if x else None)

    ### Standardize the product column
    if standardize_product:
        df["Product"] = standardize_mol_col(
            df["Product"],
            cleanup=False,
            sanitize=True,
            normalize=True,
            largest_frag=False,
            remove_charge=False,
            tautomerize=False,
            remove_stereo=False,
            cal_2D_coord_default=True,
            cal_2D_coord_deep=False,
        )

    return df


####################################################################################################
### Scaffold Analysis


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
        df["BM_smiles"] = df["BM_mol"].progress_apply(Chem.MolToSmiles)

        ### Standardize SMILES
        print("rdMolStandardize.StandardizeSmiles")
        df["BM_smiles"] = df["BM_smiles"].progress_apply(
            Chem.MolStandardize.rdMolStandardize.StandardizeSmiles
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
