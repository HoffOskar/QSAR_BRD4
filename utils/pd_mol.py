### Imports
import copy
import math

from rdkit import Chem
from rdkit.Chem import Draw, PandasTools, rdCoordGen
from rdkit.Chem.Draw import rdDepictor
from rdkit.Chem.MolStandardize import rdMolStandardize


def render_mol(df):
    """
    Renders the "Mol" column of a DataFrame for visual inspection of chemical structures.
    """
    PandasTools.ChangeMoleculeRendering(df)
    return df


def to_clipboard_smiles(df):
    """Loads a DataFrame with a "Mol" column to_clipboard"""
    # Create a copy of the input DataFrame to avoid changes to the original object
    df_copy = copy.deepcopy(df)

    # Generate SMILES for each molecule in the 'Mol' column
    df_copy["Smiles"] = df_copy["Mol"].apply(lambda mol: Chem.MolToSmiles(mol))

    # Drop the 'Mol' column
    df_copy.drop(columns=["Mol"], inplace=True)

    # Copy the modified DataFrame to clipboard
    df_copy.to_clipboard(index=True)

    return


def convert_ic50nM_to_pic50(IC50_value):
    pIC50_value = 9 - math.log10(IC50_value)
    return pIC50_value


### Visualize the mol columns set as Grid
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


### Functions
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
    Standardize a pandas.Series of Mol objects ("Mol" column).
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
