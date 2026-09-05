# β-Arrestin–PIP2 Contact Analysis

## Overview

This repository contains Python scripts for residue-level analysis of β-arrestin–PIP2 contacts from **coarse-grained molecular dynamics (CG-MD) trajectories**.

Separate scripts are provided for:

- **β-arr1 (ARR2):** residues 1–418
- **β-arr2 (ARR3):** residues 1–409

The analysis calculates the mean number of PIP2 contacts for each protein residue across the trajectory and generates both numerical contact data and a PDB file for visualization in VMD.

## Method

For each trajectory frame, the coarse-grained beads of each β-arrestin residue are compared with all PIP2 beads.

A contact is defined when a protein bead and a PIP2 bead are within:

**0.6 nm (6 Å)**

The number of bead–bead pairs satisfying this distance criterion is counted for each residue in every frame. These contact counts are then averaged over all analyzed trajectory frames to obtain the mean PIP2 contact value for each residue.

Because the analysis is performed on **coarse-grained trajectories**, the calculated contacts represent **bead–bead contacts rather than atom–atom contacts**.

In the systems used here, PIP2 is represented by the residue name `POP6`:

```python
lipid_selection = "resname POP6"
```

If a different residue name is used for PIP2 in another coarse-grained system, this selection should be modified accordingly.

## Normalization

The mean PIP2 contact values are normalized between 0 and 1 using min-max normalization:

```text
Normalized contact(i) = (x(i) - x(min)) / (x(max) - x(min))
```

where `x(i)` is the mean PIP2 contact value of residue `i`, and `x(min)` and `x(max)` are the minimum and maximum mean contact values among the analyzed protein residues.

Therefore:

- **0.0** = lowest relative PIP2 contact
- **1.0** = highest relative PIP2 contact

The normalized values represent relative PIP2 contact levels within each analyzed system.

## Input

The scripts require a structure/topology file and a coarse-grained molecular dynamics trajectory.

For example:

```python
topology_file = "structure.pdb"
trajectory_file = "trajectory.pdb"
```

Trajectory formats supported by MDAnalysis, such as multi-model PDB, TRR, or XTC, can be used provided that the topology and trajectory are compatible.

## Output

Each script generates two main output files.

### CSV Contact Data

The CSV file contains the residue-level contact results:

- `Resid` — residue number
- `Resname` — residue name
- `Mean_bead_contacts_with_PIP2` — mean number of coarse-grained bead–bead contacts with PIP2 across the trajectory
- `Normalized_contact_0_1` — min-max normalized contact score
- `Frame_count` — number of analyzed trajectory frames

### PDB Contact Map

A PDB file is also generated for visualization of the residue-level PIP2 contact distribution on the β-arrestin structure.

The normalized contact score of each residue is written into the **B-factor (Beta) field** of the corresponding PDB records:

```text
B-factor field = normalized PIP2 contact score
```

All beads belonging to the same protein residue receive the same normalized value.

The values stored in this field range from **0.00 to 1.00** and are used only for visualization. They do **not** represent experimental crystallographic B-factors.

## VMD Visualization

The generated PDB contact map can be loaded directly into VMD.

In:

```text
Graphics → Representations
```

use:

```text
Selected Atoms:  protein
Drawing Method:  NewCartoon
Coloring Method: Beta
```

VMD will color the β-arrestin structure according to the normalized PIP2 contact score stored in the B-factor field.

Residues with values closer to **1.0** have higher relative PIP2 contact levels, whereas residues with values closer to **0.0** have lower relative PIP2 contact levels.

This provides a three-dimensional residue-level map of PIP2-contacting regions on β-arrestin.

> **Note:** The calculated values represent relative coarse-grained PIP2 bead contacts. They should not be interpreted as binding affinities, interaction energies, free energies, atomistic contact numbers, or experimental B-factors.

## Requirements

- Python
- MDAnalysis
- NumPy
- pandas

Install the required packages using:

```bash
pip install MDAnalysis numpy pandas
```
