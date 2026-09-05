import os
import numpy as np
import pandas as pd
import MDAnalysis as mda

from MDAnalysis.lib.distances import capped_distance


# =========================================================
# USER SETTINGS
# =========================================================

# Topology / initial structure file
topology_file = "basal_arr3_0.pdb"

# Molecular dynamics trajectory
# Examples: trajectory.trr, trajectory.xtc, trajectory.pdb
trajectory_file = "basal_arr3.pdb"

# Output directory
output_directory = "outputs"

# Protein region to analyze
protein_selection = "protein and resid 1-409"

# Martini residue name used for PIP2
lipid_selection = "resname POP6"

# Contact cutoff
# 0.6 nm = 6 Angstrom
cutoff_nm = 0.6
cutoff_A = cutoff_nm * 10.0


# =========================================================
# CREATE OUTPUT DIRECTORY
# =========================================================

os.makedirs(
    output_directory,
    exist_ok=True
)


# =========================================================
# OUTPUT FILES
# =========================================================

pdb_out = os.path.join(
    output_directory,
    "arr3_pip2_contact_heatmap.pdb"
)

csv_out = os.path.join(
    output_directory,
    "arr3_pip2_contact_values.csv"
)


# =========================================================
# CHECK INPUT FILES
# =========================================================

if not os.path.exists(topology_file):

    raise FileNotFoundError(
        f"Topology file not found:\n{topology_file}"
    )


if not os.path.exists(trajectory_file):

    raise FileNotFoundError(
        f"Trajectory file not found:\n{trajectory_file}"
    )


# =========================================================
# LOAD SYSTEM
# =========================================================

u = mda.Universe(
    topology_file,
    trajectory_file
)


# =========================================================
# ATOM SELECTIONS
# =========================================================

protein = u.select_atoms(
    protein_selection
)

pip2 = u.select_atoms(
    lipid_selection
)


# =========================================================
# SYSTEM INFORMATION
# =========================================================

print("\n" + "=" * 70)
print("SYSTEM INFORMATION")
print("=" * 70)

print(
    "Topology file        :",
    topology_file
)

print(
    "Trajectory file      :",
    trajectory_file
)

print(
    "Total number of atoms:",
    len(u.atoms)
)

print(
    "Protein atoms        :",
    len(protein)
)

print(
    "Protein residues     :",
    len(protein.residues)
)

print(
    "PIP2/POP6 atoms      :",
    len(pip2)
)

print(
    "Number of frames     :",
    len(u.trajectory)
)


# =========================================================
# VALIDATE SELECTIONS
# =========================================================

if len(protein) == 0:

    raise ValueError(
        "Protein selection returned no atoms. "
        "Check protein_selection."
    )


if len(pip2) == 0:

    raise ValueError(
        "PIP2 selection returned no atoms. "
        "Check lipid_selection."
    )


# =========================================================
# INITIALIZE RESIDUE CONTACT STORAGE
# =========================================================

residue_contact_sum = {}

residue_frame_count = {}


for residue in protein.residues:

    key = (
        residue.resid,
        residue.resname
    )

    residue_contact_sum[key] = 0.0

    residue_frame_count[key] = 0


print(
    "\nNumber of protein residues analyzed:",
    len(residue_contact_sum)
)


# =========================================================
# CONTACT ANALYSIS
# =========================================================

n_frames = len(u.trajectory)


for frame_index, ts in enumerate(u.trajectory):

    # PIP2 coordinates for the current frame
    pip2_positions = pip2.positions

    # Periodic simulation box
    box = ts.dimensions


    for residue in protein.residues:

        key = (
            residue.resid,
            residue.resname
        )


        # -------------------------------------------------
        # Identify residue-PIP2 atom pairs within 6 Angstrom
        # -------------------------------------------------

        pairs = capped_distance(

            residue.atoms.positions,

            pip2_positions,

            max_cutoff=cutoff_A,

            box=box,

            return_distances=False
        )


        # -------------------------------------------------
        # Number of atom-atom contacts in this frame
        # -------------------------------------------------

        contact_count = len(pairs)


        residue_contact_sum[
            key
        ] += contact_count


        residue_frame_count[
            key
        ] += 1


    # -----------------------------------------------------
    # Progress information
    # -----------------------------------------------------

    if (
        frame_index == 0
        or (frame_index + 1) % 100 == 0
        or (frame_index + 1) == n_frames
    ):

        print(
            f"Frame "
            f"{frame_index + 1}"
            f"/{n_frames} processed"
        )


# =========================================================
# CALCULATE MEAN CONTACT VALUES
# =========================================================

residue_mean_contact = {}


for key in residue_contact_sum:

    if residue_frame_count[key] > 0:

        residue_mean_contact[key] = (

            residue_contact_sum[key]

            /

            residue_frame_count[key]
        )

    else:

        residue_mean_contact[key] = 0.0


# =========================================================
# MIN-MAX NORMALIZATION TO 0-1
# =========================================================

values = np.array(
    list(
        residue_mean_contact.values()
    ),
    dtype=float
)


vmin = values.min()

vmax = values.max()


residue_norm_score = {}


for key, value in residue_mean_contact.items():

    if vmax == vmin:

        residue_norm_score[key] = 0.0

    else:

        residue_norm_score[key] = (

            (value - vmin)

            /

            (vmax - vmin)
        )


print("\n" + "=" * 70)
print("NORMALIZATION")
print("=" * 70)


print(
    "Minimum mean contact :",
    vmin
)

print(
    "Maximum mean contact :",
    vmax
)

print(
    "Normalized range     :",
    min(residue_norm_score.values()),
    "-",
    max(residue_norm_score.values())
)


# =========================================================
# CREATE CSV OUTPUT
# =========================================================

rows = []


for key in residue_mean_contact:

    resid, resname = key


    rows.append({

        "Resid":
            resid,

        "Resname":
            resname,

        "Mean_atom_contacts_with_PIP2":
            residue_mean_contact[key],

        "Normalized_contact_0_1":
            residue_norm_score[key],

        "Frame_count":
            residue_frame_count[key]
    })


df = pd.DataFrame(
    rows
)


df = df.sort_values(
    "Resid"
)


df.to_csv(
    csv_out,
    index=False
)


print("\nCSV output:")
print(csv_out)


# =========================================================
# WRITE NORMALIZED VALUES TO PDB B-FACTOR COLUMN
# =========================================================

with open(
    topology_file,
    "r"
) as fin, open(
    pdb_out,
    "w"
) as fout:


    for line in fin:


        if line.startswith(
            ("ATOM  ", "HETATM")
        ):


            # ---------------------------------------------
            # Read residue information from PDB
            # ---------------------------------------------

            resname = (
                line[17:20].strip()
            )


            try:

                resid = int(
                    line[22:26]
                )

            except ValueError:

                fout.write(
                    line
                )

                continue


            key = (
                resid,
                resname
            )


            # ---------------------------------------------
            # Get normalized PIP2 contact score
            #
            # Non-protein / non-analyzed residues = 0
            # ---------------------------------------------

            score = (
                residue_norm_score.get(
                    key,
                    0.0
                )
            )


            line_no_newline = (
                line.rstrip("\n")
            )


            if len(line_no_newline) < 66:

                line_no_newline = (
                    line_no_newline.ljust(
                        66
                    )
                )


            # ---------------------------------------------
            # Write normalized score into
            # PDB B-factor columns 61-66
            # ---------------------------------------------

            new_line = (

                line_no_newline[:60]

                + f"{score:6.2f}"

                + line_no_newline[66:]

                + "\n"
            )


            fout.write(
                new_line
            )


        else:

            fout.write(
                line
            )


print("\nPDB output:")
print(pdb_out)


# =========================================================
# PRINT TOP 20 PIP2-CONTACTING RESIDUES
# =========================================================

print("\n" + "=" * 70)

print(
    "TOP 20 RESIDUES WITH THE "
    "HIGHEST PIP2 CONTACT VALUES"
)

print("=" * 70)


top20 = (

    df

    .sort_values(
        "Mean_atom_contacts_with_PIP2",
        ascending=False
    )

    .head(20)
)


print(

    top20[
        [
            "Resid",
            "Resname",
            "Mean_atom_contacts_with_PIP2",
            "Normalized_contact_0_1"
        ]
    ]

    .to_string(
        index=False
    )
)


# =========================================================
# VMD VISUALIZATION INFORMATION
# =========================================================

print("\n" + "=" * 70)
print("VMD VISUALIZATION")
print("=" * 70)


print(
    "\nLoad the generated PDB file:"
)

print(
    pdb_out
)


print(
    "\nGraphics -> Representations"
)

print(
    "Selected Atoms  : protein"
)

print(
    "Drawing Method  : NewCartoon"
)

print(
    "Coloring Method : Beta"
)


print(
    "\nNormalized PIP2 contact values "
    "range from 0 to 1."
)

print(
    "0 = lowest relative PIP2 contact value"
)

print(
    "1 = highest relative PIP2 contact value"
)


# =========================================================
# FINISHED
# =========================================================

print("\n" + "=" * 70)
print("ANALYSIS COMPLETED")
print("=" * 70)
