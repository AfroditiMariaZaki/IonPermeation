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

# import logging
# import warnings

import MDAnalysis as mda
import pandas as pd
import matplotlib.pyplot as plt
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.base import Results


# def ion_positions(ag):
#     # ag = self.atomgroup
#     # print(ag.atoms.resids)
#     z_pos = ag.atoms.positions[:,2]
#     # sq_z_pos = z_pos
#     return z_pos

def in_cylinder(atomgroup, channel_upper, channel_lower, channel_com, zmin=0, cutoff=5.0):
    # cyl_atoms = []
    upper_sel = u.select_atoms(channel_upper)
    lower_sel = u.select_atoms(channel_lower)
    zmax = upper_sel.center_of_geometry()[2]-lower_sel.center_of_geometry()[2]
    cyzone = u.select_atoms(f"{atomgroup} and cyzone {cutoff} {zmax} {zmin} {channel_com}")
    # print(cyzone.names, cyzone.indices)

    return cyzone
    # cyl_atoms.append((cyzone.indices))


class IonPermeation(AnalysisBase):
    def __init__(self, universe, atomgroup, 
                 channel_upper, channel_lower, channel_com, channel_gate,
                 verbose=True):
        # super(IonPermeation, self).__init__(atomgroup.universe.trajectory, verbose=verbose)
        super(IonPermeation, self).__init__(universe.trajectory, verbose=verbose)

        self.u = universe
        # self._trajectory = self.u.trajectory
        self.atomgroup = self.u.select_atoms(atomgroup)
        self.channel_upper = channel_upper
        self.channel_lower = channel_lower
        self.channel_com = channel_com
        self.channel_gate = channel_gate
        # self.results = Results()
        # self._trajectory = self.u.trajectory
        # self._parameter = cutoff
        # self._ag = ion_sel
        # self.results = np.zeros((self._trajectory.n_frames, 2), dtype=np.float32)


    def _prepare(self):
        # self.results = np.zeros((self.atomgroup.universe.trajectory.n_frames, len(self.atomgroup)+1), dtype=np.float32)
        self.results = np.zeros((self.u.trajectory.n_frames, 2), dtype=np.float32)

        cyl_atoms = []
        events = 0
        # cum_events = []
        cyl_zone = in_cylinder(ions, channel_upper, channel_lower, channel_com)
        print("In channel atoms: ", cyl_zone)

        print(self.results)

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
        gate_sel = u.select_atoms(channel_gate)
        exit_point = gate_sel.center_of_geometry()[2]
        # print(zmin)
        for atom in cyl_zone:
            # print(atom.index, atom.position[2])
            if atom.position[2]<exit_point:
                events+=1
                print(ts.time/1000, int(events), exit_point, atom.position[2], atom.index)
        self.results[self._frame_index, 0] = self._ts.frame
        self.results[self._frame_index, 1] = events
                # cum_events.append((ts.time/1000, int(events)))
        cyl_zone = in_cylinder(ions, channel_upper, channel_lower, channel_com)

        
        # z_pos = ion_positions(self)
        # res = ion_positions(self.atomgroup)
        # save it into self.results
        # self.results[self._frame_index, 1:] = res
        # the current timestep of the trajectory is self._ts
  
        # the actual trajectory is at self._trajectory
        # self.results[self._frame_index, 1] = self._trajectory.time

    def _conclude(self):
        """
        Finish up by transforming our
        results into a DataFrame.
        """
        # by now self.result is fully populated
        # self.average = np.mean(self.results[:, 1], axis=0)
        # columns = ['Frame', 'Ion Position']
        self.df = pd.DataFrame(self.results)
        print(self.df)

    def plot(self):
        """
        Plot the results.
        """
        fig, ax = plt.subplots()
        ax.plot(self.results[:,0], self.results[:,1], color='red', label="z-position")
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

# path = "/Users/afroditimariazaki/Documents/WORK/TPC2/CONDUCTANCE/Truncated/NA/AllRestraints/V500/Repeat2"
path = "/Users/afroditimariazaki/Documents/GITHUB/IonPermeation"

# u = mda.Universe(f"{path}/conf.pdb", f"{path}/traj_fit.xtc")
u = mda.Universe(f"{path}/test.pdb", f"{path}/test.xtc")

# ion = u.select_atoms("resname SOD")
# print(len(ion))
# # print(ion_sel.positions)

ions = "resname SOD"
channel_upper = "protein and resid 654 and name CA"
channel_lower = "protein and resid 312 and name OH"
channel_com = "protein and resid 312 and name OH"
channel_gate = "protein and resid 312 and name OH"
# start, stop, step = 0, u.trajectory.n_frames, 2


ip = IonPermeation(u, ions, channel_upper, channel_lower, channel_com, channel_gate, verbose=True).run()
ip.plot()
# print(ip.ion_positions())
# print(ip._prepare())
# print(ip._conclude())
# print(ip._single_frame())


# print(IonPermeation.results)
# print(IonPermeation(universe=u, atomgroup="resname SOD")._single_frame())
# IonPermeation.
# ip.run(start=0, stop=101, step=10)
# print(ip.results)



# # permeation = IonPermeation(universe=u,
# #                            ion_sel="resname SOD",
# #                            in_sel="protein and resid 652-654",
# #                            out_sel="protein and resid 308 312 690 694 and name CA")

# # permeation.ion_positions()
# # permeation.in_cylinder()