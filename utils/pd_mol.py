### Imports
import copy
import math

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools, rdChemReactions, rdCoordGen
from rdkit.Chem.Draw import rdDepictor
from rdkit.Chem.MolStandardize import rdMolStandardize

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


def standardize_mol_col(
    mol_col,
    largest_frag=True,
    remove_charge=True,
    tautomerize=False,
    normalize=True,
    remove_stereo=False,
    cal_2D_coord_default=True,
    cal_2D_coord_deep=False,
    cleanup=True,
):
    """
    Standardizes a pandas Series of RDKit Mol objects by applying a sequence of standardization operations.

    Parameters:
    -----------
    mol_col : pandas.Series
        A Series or DataFrame column containing RDKit Mol objects.

    largest_frag : bool
        Keep only the largest fragment.

    remove_charge : bool
        Neutralize charged molecules using `rdMolStandardize.Uncharger().uncharge`.

    tautomerize : bool
        Convert molecules to their most populated tautomer.
        (This can be slow, so use only when necessary.)

    normalize : bool
        Apply molecular normalization.

    remove_stereo : bool
        Removes stereochemical information.

    cal_2D_coord_default : bool
        Computes 2D coordinates for molecules using `rdDepictor.Compute2DCoords`.
        (Important for combinatorial libraries.)

    cal_2D_coord_deep : bool
        Generate 2D coordinate, which can better handle specific structural features (e.g., straightening alkynes to 180°).
        Overrides `cal_2D_coord_default` if set to True.

    cleanup : bool
        Cleanup structure.

    Returns:
    --------
    pandas.Series
        A Series or DataFrame column of RDKit Mol objects after applying the selected standardization steps.
    """

    # largest fragment only
    print(f"largest_frag={largest_frag} (rdMolStandardize.FragmentParent):")
    if largest_frag:
        mol_col = mol_col.progress_apply(rdMolStandardize.FragmentParent)

    # remove charges
    print(f"remove_charge={remove_charge} (rdMolStandardize.Uncharger().uncharge):")
    if remove_charge:
        mol_col = mol_col.progress_apply(rdMolStandardize.Uncharger().uncharge)

    # generate most populated tautomer # slow only use when necessary
    print(f"tautomerize={tautomerize} (rdMolStandardize.CanonicalTautomer)")
    if tautomerize:
        mol_col = mol_col.progress_apply(rdMolStandardize.CanonicalTautomer)

    # normalize (not sure if necessary)
    print(f"normalize={normalize} (rdMolStandardize.Normalize):")
    if normalize:
        mol_col = mol_col.progress_apply(rdMolStandardize.Normalize)

    # remove stereoinformation
    print(f"remove_stereo={remove_stereo} (rdMolStandardize.Normalize):")
    if remove_stereo:
        mol_col = mol_col.progress_apply(rdMolStandardize.StereoParent)

    ### 2D-coordinates are only calculated once
    if cal_2D_coord_deep:
        print("cal_2D_coord_default=not required (rdDepictor.Compute2DCoords):")
        print(f"cal_2D_coord_deep={cal_2D_coord_deep} (rdCoordGen.AddCoords):")
        # alternative to fix 2D-coordinates # straightens alkynes to 180°
        mol_col.progress_apply(rdCoordGen.AddCoords)

    elif cal_2D_coord_default:
        print(
            f"cal_2D_coord_default={cal_2D_coord_default} (rdDepictor.Compute2DCoords):"
        )
        # default to fix 2D-coordinates # important for combinatorial libraries
        mol_col.progress_apply(rdDepictor.Compute2DCoords)
        print(f"cal_2D_coord_deep={cal_2D_coord_deep} (rdCoordGen.AddCoords):")
    else:
        print(
            f"cal_2D_coord_default={cal_2D_coord_default} (rdDepictor.Compute2DCoords):"
        )
        print(f"cal_2D_coord_deep={cal_2D_coord_deep} (rdCoordGen.AddCoords):")

    # clean up visualization
    print(f"cleanup={cleanup} (rdMolStandardize.Cleanup):")
    if cleanup:
        mol_col = mol_col.progress_apply(rdMolStandardize.Cleanup)

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
        return None


def rxn_2_sms(df, rxn_smarts, mol_col_1="Mol", mol_col_2=None):
    """
    React two molecule columns in a horizontal fashion create products columns.

    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing the molecules to react

    rxn_smarts : str
        SMARTS pattern for the reaction

    mol_col_1 : str
        Column name containing the first molecule

    mol_col_2 : str or None
        Column name containing the second molecule. Water is used as 2nd reactant by default.

    Returns:
    --------
    pandas.DataFrame
        DataFrame containing the original columns and the products of the reaction
    """

    ### Apply the virtual reaction to each row of the DataFrame
    rxn_products = df.progress_apply(
        lambda row: _rxn_2_sms(
            rxn_smarts,
            row[mol_col_1],
            row[mol_col_2]
            if mol_col_2
            else Chem.MolFromSmiles("O"),  # water by default
        ),
        axis=1,
    )

    ### Convert the Series of tuples into a DataFrame
    products_df = pd.DataFrame(rxn_products.tolist(), index=df.index)

    ### Add new columns to the original DataFrame
    for col in products_df.columns:
        df[f"Product_{col}"] = products_df[col]

    return df
