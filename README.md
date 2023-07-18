# IonPermeation

IonPermeation is a Python tool written for MDAnalysis that analyzes the ion permeation through ion channels.

It has the following functionalities:

* Computes the cumulative number of permeation events for a trajectory. 
* Calculates and plots the z-position of the ions and the water molecules that pass through the ion channel as a function of time.
* Calculates the average binding time of each passing ion to the channel constriction points.


# Usage

```
u = mda.Universe("test.pdb", "test.xtc")
ion = u.select_atoms("resname SOD")
in_sel = u.select_atoms("protein and resid 652 653 654 and name CA")
out_sel = u.select_atoms("protein and resid 690 694 and name CA")


ip = IonPermeation(ion, in_sel, out_sel, verbose=True).run()
```
