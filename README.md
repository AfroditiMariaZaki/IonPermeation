# IonPermeation

IonPermeation is a Python tool written for MDAnalysis that analyzes the ion permeation through ion channels.

It has the following functionalities:

* Computes the cumulative number of permeation events for a trajectory and saves data in a .csv file.
* Plots cumulative number of permeation events as a function of simulation time.
* Computes the ion conductance of the channel.

# ToDo :  
* Calculates and plots the z-position of the ions and the water molecules that pass through the ion channel as a function of time.
* Calculates the average binding time of each passing ion to the channel constriction points.


# Usage

```
u = mda.Universe("test.pdb", "test.xtc")
ions = "resname SOD"
in_sel = "protein and resid 271 272 652 653 654 and name CA"
out_sel = "protein and resid 308 312 690 694 and name CA"


IP = IonPermeation(ions, in_sel, out_sel, verbose=True)
IP.run()
IP.save()
IP.plot()
IP.conductance()
```
