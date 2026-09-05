# beta-arrestin-pip2-contact-analysis
# β-Arrestin–PIP2 Contact Analysis

Python scripts for residue-level β-arrestin–PIP2 contact analysis using **coarse-grained molecular dynamics trajectories**.

The repository contains separate analysis scripts for β-arr1 (ARR2) and β-arr2 (ARR3).

## Scripts

- `arr2_residue_pip2_contact_analysis.py`
  - β-arr1
  - Residues 1–418

- `arr3_residue_pip2_contact_analysis.py`
  - β-arr2
  - Residues 1–409

## Analysis Overview

These scripts calculate residue-level contacts between β-arrestin and PIP2 from coarse-grained molecular dynamics simulations.

For each trajectory frame, every protein residue is compared with all PIP2 beads.

A contact is defined when a protein bead and a PIP2 bead are within:

`0.6 nm (6 Å)`

The number of protein–PIP2 bead pairs within the cutoff is counted for each residue in each frame.

The mean contact value for residue `i` is calculated as:

Mean contact(i) = total number of residue–PIP2 contacts / number of trajectory frames

The resulting mean contact values are then normalized between 0 and 1 using min-max normalization:

Normalized contact(i) = (x_i - x_min) / (x_max - x_min)

where:

- `x_i` is the mean PIP2 contact value of residue `i`
- `x_min` is the minimum mean contact value among all analyzed residues
- `x_max` is the maximum mean contact value among all analyzed residues

## Coarse-Grained Representation

The analysis is designed for **coarse-grained molecular dynamics trajectories**.

Therefore, the calculated contacts correspond to interactions between coarse-grained beads rather than individual atoms.

In the systems used for this analysis, PIP2 is represented with the residue name:

`POP6`

The corresponding MDAnalysis selection is:

```python
lipid_selection = "resname POP6"
