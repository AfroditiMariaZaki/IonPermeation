"""
Ion Permeation Analysis
=======================================================

:Author: Afroditi Maria Zaki
:Affiliation: SBCB, University Of Oxford
:Year: 2023

This module contains the :class:'IonPermeation' that analyses the ion permeation through an ion channel in a Universe.

Input
-----

Required:
    - *universe*: an MDAnalysis Universe object
    - *atomgroup*: group of ions that permeate the ion channel
    - *in_sel*: selection of atome/residues that form the entrance of the ion channel
    - *out_sel*: selection of atome/residues that form the exit of the ion channel

Options:


Output
------
    df : A pandas DataFrame with two columns:
        - *frame*: Every frame of the trajectory
        - *events*: Cumulative number of permeation events at every frame
=======================================================

How it works
------------
    -For every frame:
        * The position of every ion across the membrane plane is read and if found to be within the cylinder defined by the pore axis and a cutoff radius
        and between the limits set by the residues that form the entry and exit points of the ion channel, it is stored in an AtomGroup (cyzone).
        * The position of every ion in the cyzone AtomGroup across the pore axis that was computed in the previous frame is read and if found to be lower than the exit point,
        the ion is considered to have permeated through the ion channel and the events counter (self.events) increases by +1.
        * If the ion index number (self.ion_id) matches the ion index number that was saved at the previous frame, the permeation event is considered to have already been counted.
          This happens because sometimes the ion will be engaged by the gate-forming residues, moving in and out of the exit point before completely passing through. 
          Counting this event multiple times will overestimate the number of permeation events and therefore the ion conductance.
        * The cyzone AtomGroup is updated to store the ions that are found in the ion channel.


        
Example use of :class: 'IonPermeation'        
--------------------------------------

import MDAnalysis as mda

u = mda.Universe(grofile, trjfile)

IP = IonPermeation(
    universe=u, 
    atomgroup='resname SOD', 
    in_sel='protein and resid 271 272 652 653 654 and name CA', 
    out_sel='protein and resid 308 312 690 694 and name CA'
    )


To compute ion permeation events::

IP.run()

To save data in a .csv file::

IP.save()

To plot data::

IP.plot()

To compute ion conductance::

IP.conductance()


Classes
-------

.. autoclass:: IonPermeation


Functions
---------

.. autofunction:: in_channel_ions

"""


import numpy as np
import MDAnalysis as mda
import pandas as pd
import matplotlib.pyplot as plt
from MDAnalysis.analysis.base import AnalysisBase


def in_channel_ions(atomgroup, in_sel, out_sel, channel_com, zmin=0, cutoff=5.0):
    """
    Creates an atomgroup with all the ions that are found in the ion conduction pore

    Parameters
    ----------
    atomgroup : str
        Selection string for the permeating ions
    in_sel : str
        Selection string for the residues or atoms that define the channel entry point
    out_sel : str
        Selection string for the residues or atoms that define the channel exit point
    channel_com : str
        Selection string for the atoms that form the lowest end of the ion conduction pore.
        Any permeating ion found to have crossed this point will be considered to have permeated
        through the channel and will increase the count of the permeation events by one.
    zmin : float, optional
        The minimum point along the z-axis of the ion channel and the cylinder that encloses the permeating ion. 
        If the channel_com contains the atoms of the exit point of the ion channel, thn zmin should be equal to zero.
    cutoff : float, optional
        Define the radius of the cylinder that encloses the ion permeation channel

    Returns
    -------
    cyzone : AtomGroup 
        Atomgroup that contains all the ions that are found in the ion conduction pore.
        This function is called at every frame, to update the selection of the permeating ions.

    Note
    ----
    The code assumes that the ion conduction pore is parallel to the z-axis. 
    It is also written for monoatomic ions. If your system contains polyatomic ions, please make sure to select only one atom type in your atomgroup selection.
    """

    in_channel = u.select_atoms(in_sel)
    out_channel = u.select_atoms(out_sel)
    zmax = in_channel.center_of_geometry()[2]-out_channel.center_of_geometry()[2]
    cyzone = u.select_atoms(f"{atomgroup} and cyzone {cutoff} {zmax} {zmin} {channel_com}")

    return cyzone


class IonPermeation(AnalysisBase):
    """
    Compute the cumulative number of ion permeation events on a given trajectory

    Parameters
    ----------
    universe : Universe
        MDAnalysis universe object
    atomgroup : str
        Selection string for the permeating ions
    in_sel : str
        Selection string for the residues or atoms that define the channel entry point
    out_sel : 
        Selection string for the residues or atoms that define the channel exit point
    verbose : bool, optional
        Turn on more logging and debugging.

    """

    def __init__(self, universe, atomgroup, in_sel, out_sel, verbose=True):
        self.u = universe
        self._trajectory = self.u.trajectory
        self.atomgroup = atomgroup
        self.in_sel = self.u.select_atoms(in_sel)
        self.out_sel = self.u.select_atoms(out_sel)


    def _prepare(self):
        """
        Prepare array to store data and initialize variables
        """
        self.cyl_zone = in_channel_ions(ions, in_sel, out_sel, channel_com)
        self.results = np.zeros((self._trajectory.n_frames, 2), dtype=np.float32)
        print(self.results.shape)
        self.events = 0
        self.ion_id = 0


    def _single_frame(self):
        """
        This function is called for every frame of the trajectory
        """
        exit_point = np.min(self.out_sel.positions[:,2])
        if len(self.cyl_zone) == 0:
            self.results[self._frame_index, 0] = self._ts.frame
            self.results[self._frame_index, 1] = self.events
        else:
            for atom in self.cyl_zone:
                if atom.position[2]<exit_point and atom.index == self.ion_id:
                    # Event has already been counted
                    self.results[self._frame_index, 0] = self._ts.frame
                    self.results[self._frame_index, 1] = self.events
                elif atom.position[2]<exit_point and atom.index != self.ion_id: 
                    self.events+=1
                    self.results[self._frame_index, 0] = self._ts.frame
                    self.results[self._frame_index, 1] = self.events
                    self.ion_id = atom.index
                else:
                    self.results[self._frame_index, 0] = self._ts.frame
                    self.results[self._frame_index, 1] = self.events
        self.cyl_zone = in_channel_ions(ions, in_sel, out_sel, channel_com)

        print(self._ts.frame, self.events)


    def conductance(self, voltage=0.15, valence=1):
        """
        Ion channel conductance calculation

        Formula
        -------
        gamma = 1 / R
        gamma = I / V
        gamma = N * valence * q / t * V

        Parameters
        ----------
        
        """
        q = 1.6e-19
        charge = valence*q
        N = self.events
        print(self._trajectory.totaltime)
        time = self._trajectory.totaltime * 1e-12
        current = (N*charge)/time
        conductance = current/voltage
        conductance_in_pS = round(conductance*1e+12, 2)
        print("Ion channel conductance is", conductance_in_pS, "pS")

        return conductance_in_pS


    def _conclude(self):
        """
        Finish up by transforming the
        results into a DataFrame.
        """
        self.df = pd.DataFrame(self.results)


    def plot(self):
        """
        Visualize the results
        """
        fig, ax = plt.subplots()
        ax.plot(self.results[:,0], self.results[:,1], color='red', label="Permeation events")
        ax.set_xlabel("Frame")
        ax.set_ylabel("Cumulative permeation events")
        ax.legend(loc=2)
        plt.show()


    def save(self, path=None, filename='permeation_events'):
        """
        Saves data in .csv file
        

        Parameters
        ----------
        path : str, optional
            Path to save the output .csv file
        filename : str, optional
            Name of the output .csv file
        """
        import os
        if path is None:
            path = os.getcwd()
        file = self.df.to_csv(f"{path}/{filename}.csv", header=["Frame", "Events"], index=False)
        return file



################################################################################

topfile = "data/test.pdb"
trjfile = "data/test.xtc"

u = mda.Universe(topfile, trjfile)

ions = "resname SOD"
channel_com = "protein and resid 312 and name OH"

in_sel = "protein and resid 271 272 652 653 654 and name CA"
out_sel = "protein and resid 690 694 308 312 and name CA"

IP = IonPermeation(u, ions, in_sel, out_sel)
IP.run()
IP.save()
IP.plot()
IP.conductance(voltage=0.75)
