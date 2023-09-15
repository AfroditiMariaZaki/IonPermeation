"""
=======================================================

:Author: Afroditi Maria Zaki
:Affiliation: SBCB, University Of Oxford
:Year: 2023

This module computes the number of ion permeation events through an ion channel in a Universe.

Input
-----

Required:
    - *Universe*: an MDAnalysis object

    Options:
    # ToDo

Output
------

    - *frame*: Frame at which an ion permeation event occured
    - *ion ID*: atom ID of the ion that crossed the channel
    - *count*: Total number of ion permeation events
=======================================================

How it works
------------
    -For every frame:
        * The position of ions across the membrane plane is read and if found to be within the cylinder defined by the pore axis and a cutoff radius,
        it is stored in a list (cylinder_ions).
        * The position of the ions in the previous list across the pore axis is read and if found to be between the upper and lower barriers set by
        the selectivity filter and the gate residues, the ion is considered to be in the channel and is stored in a list (channel_ions).
        * The current z-coordinates of the atoms in the existing (before updating) channel_ions list is compared to heir z-coordinates 
        from the previous step and if lower than selecivity filter barrier, the count of the crossing events goes up by 1. 

Example use of :class: 'IonPermeation'        
--------------------------------------

import MDAnalysis as mda
import IonPermeation as ip

u = mda.Universe(grofile, trjfile)

PermeationEvents = ip(
    universe=u, 
    ion_sel='SOD', 
    in_sel='resid 123', 
    out_sel='resid 456', 
    cutoff = 5.0)

PermeationEvents.run()

"""


import numpy as np

import logging
import warnings

import MDAnalysis as mda
import pandas as pd
import matplotlib.pyplot as plt
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.base import Results


def ion_positions(ag):
    # ag = self.atomgroup
    # print(ag.atoms.resids)
    z_pos = ag.atoms.positions[:,2]
    # sq_z_pos = z_pos
    return z_pos

def channel_limits(in_sel, out_sel):
    """
    Define entry and exit point of ion channel
    """
    z_in = in_sel.center_of_mass()[2]
    z_out = out_sel.center_of_mass()[2]

    return z_in, z_out

def pore_center(in_sel, out_sel):
    """
    Define center of pore axis
    """
    v1 = in_sel.center_of_mass()
    v2 = out_sel.center_of_mass()
    v3 = v2 - v1

    return v3

def in_cylinder(atomgroup, zmax=10, zmin=-10, cutoff=5.0):
    com = "protein and resid 690"
    cyzone = u.select_atoms(f"cyzone {cutoff} {zmax} {zmin} ({com})")
    print(cyzone.atoms)
    # in_cyl = []
    # rcp = np.array([out_sel.center_of_mass()[0], out_sel.center_of_mass()[1]])
    # for atom in atomgroup:
    #     r = np.array([atom.position[0], atom.position[1]])
    #     dist = np.linalg.norm(rcp-r)
    #     if dist <= cutoff:
    #         in_cyl.append((atom.index, atom.position[2]))
    # print(in_cyl)

# def in_channel(atomgroup, in_sel, out_sel):
#     in_chan = []

class IonPermeation(AnalysisBase):
    def __init__(self, atomgroup, in_sel, out_sel, zmin, zmax, com, cutoff, verbose=True):
        super(IonPermeation, self).__init__(atomgroup.universe.trajectory, verbose=verbose)
        # self.u = universe
        # self._trajectory = self.u.trajectory
        self.atomgroup = atomgroup
        self.in_sel = in_sel
        self.out_sel = out_sel
        self.zmin = zmin
        self.zmax = zmax
        self.com = com
        self.cutoff = cutoff

        # self.results = Results()
        # self._trajectory = self.u.trajectory
        # self._parameter = cutoff
        # self._ag = ion_sel
        # self.results = np.zeros((self._trajectory.n_frames, 2), dtype=np.float32)


    def _prepare(self):
        self.results = np.zeros((self.atomgroup.universe.trajectory.n_frames, len(self.atomgroup)+1), dtype=np.float32)
        self.channel = []
        # print(self.results)

    # def __init__(self, universe,
    #              ion_sel=None,
    #              in_sel=None,
    #              out_sel=None,
    #              cutoff=5.0):
    #     super(IonPermeation, self).__init__(ion_sel.universe.trajectory)
    #     self.u = universe
    #     self._trajectory = self.u.trajectory
    #     self.ion_sel = ion_sel.strip()
    #     self.in_sel = in_sel.strip()
    #     self.out_sel = out_sel.strip()
    #     self.cutoff = cutoff
        # self.results = Results()

    # def ion_positions(self):
    #     ions = self.atomgroup
    #     print(ions.residues)
    #     print(ions.atoms.positions[:,2])
        # print(np.array([ions.position[:,0]]))



    # def _single_frame(self):
    #     # frame = self._ts.frame
    #     ion_pos = ion_positions(self)
    #     # i = self._trajectory.ts.frame
    #     # self.results[i] = i**2
    #     # print("Single frame: ", self.results[i])
    #     # save it into self.results
    #     self.results[self._frame_index, 1] = ion_pos
    #     # the current timestep of the trajectory is self._ts
    #     self.results[self._frame_index, 0] = self._ts.frame

    def _single_frame(self):
        """
        This function is called for every frame that we choose
        in run().
        """
        # call our earlier function
        
        # z_pos = ion_positions(self)
        # res = ion_positions(self.atomgroup)
        # axis = pore_axis(self.in_sel, self.out_sel)
        pore_ions = in_cylinder(self.atomgroup, self.zmax, self.zmin, self.com, self.cutoff)
        print(pore_ions)
        for atom in pore_ions:
            if atom.index in self.channel:
                print(atom.index, True)
            els:
            print(atom.index, False)
        # save it into self.results
        # self.results[self._frame_index, 1:] = res
        # # the current timestep of the trajectory is self._ts
        # self.results[self._frame_index, 0] = self._ts.frame
        # the actual trajectory is at self._trajectory
        # self.results[self._frame_index, 1] = self._trajectory.time




    # def _conclude(self):
    #     """
    #     Finish up by transforming our
    #     results into a DataFrame.
    #     """
    #     # by now self.result is fully populated
    #     # self.average = np.mean(self.results[:, 1], axis=0)
    #     # columns = ['Frame', 'Ion Position']
    #     self.df = pd.DataFrame(self.results)
    #     print(self.df)

    def plot(self):

        fig, ax = plt.subplots()
        ax.plot(self.results[:,0], self.results[:,2], color='red', label="z-position")
        ax.set_xlabel("x axis")
        ax.set_ylabel("y axis")
        ax.legend(loc=2)
        plt.show()


    #     # self.results[self._frame_index, :] = self._trajectory.ts.frame
    #     # print(self.results[self._frame_index, :])
    #     print(self._frame_index)

#     def  _conclude(self):
#         self.results = np.asarray(self.results)
        #  self.results /= np.sum(self.results)

#     def plot(self, ax=None):
#         """Visualize the results
        
#         Parameters
#         ----------
#         ax : matplotlib.Axes, optional
#           if provided, the graph is plotted on this axis

#         Returns
#         -------
#         ax: the axis that the graph has been plotted on
#         """
#         import matplotlib as plt
#         if ax is None:
#             fig, ax = plt.subplots()
#         ax.plot(self.results
#                 'ro',
#                 label='Permeation Events')
#         ax.set_xlabel('Simulation Time [ns]')
#         ax.set_ylabel('N$_{events}$')
#         ax.legend(loc='2')

#         return ax



#     # def in_cylinder(self):
#     #     out_atoms = self.u.select_atoms(self.out_sel)
#     #     out_com = np.array([out_atoms.center_of_mass()[0], out_atoms.center_of_mass()[1]])
#     #     print(out_com)

#     # def count_by_time(self):
#     #     """Counts the number of permeation events per timestep.
        
#     #     Returns
#     #     -------
#     #     events : numpy.ndarray
#     #         Contains the total number of ion permeation events computed at each timestep.
#     #         Can be used along with :attr:`IonPermeation.times` to plot
#     #         the number of hydrogen bonds over time.
#     #     """

#     # def _conclude(self):
#     #     self.results.events = np.asarray(self.results.events).T
# ###############################################################################

path = "/Users/afroditimariazaki/Documents/GITHUB/IonPermeation"

u = mda.Universe(f"{path}/test.pdb", f"{path}/test.xtc")
ion = u.select_atoms("resname SOD")
in_sel = u.select_atoms("protein and resid 652 653 654 and name CA")
out_sel = u.select_atoms("protein and resid 690 694 and name CA")
# com = "protein and resid 690"

# cylayer innerRadius externalRadius zMax zMin selection
# cutoff = 5
# zmin = -10
# zmax= 10

# cyzone = u.select_atoms(f"cyzone {cutoff} {zmax} {zmin} (protein and resid 690 694 and name CA)")
# print(cyzone.atoms)

# for atom in cyzone:
#     print(atom.index)
# ip = IonPermeation(ion, in_sel, out_sel, -10, 10, com, 5, verbose=True).run()
in_cylinder(u)