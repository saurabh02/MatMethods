# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the ferroelectric workflow
"""

import os

from fireworks import Workflow, LaunchPad, FireTaskBase, FWAction, Firework
from fireworks.utilities.fw_utilities import explicit_serialize

from matmethods.vasp.fireworks.core import StaticFW, NonSCFFW
from matmethods.vasp.firetasks.glue_tasks import CopyVaspOutputs

from pymatgen.io.vasp.outputs import Vasprun


@explicit_serialize
class BandGapCheck(FireTaskBase):

    _fw_name = "BandGapCheck"

    def run_task(self, fw_spec):
        v = Vasprun(filename='vasprun.xml')
        (band_gap, cbm, vbm, is_direct) = v.eigenvalue_band_properties
        if band_gap > 0.5:
            return FWAction(stored_data={'band_gap': band_gap})
        else:
            return FWAction(stored_data={'band_gap': band_gap}, defuse_workflow=True)


class BandGapCheckFW(Firework):
    def __init__(self, structure, name="check_bandgap", parents=None, **kwargs):
        """
        Standard static calculation Firework for dielectric constants
        using DFPT.

        Args:
            structure (Structure): Input structure.
            name (str): Name for the Firework.
            vasp_cmd (str): Command to run vasp.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []

        t.append(CopyVaspOutputs(calc_loc=True))
        t.append(BandGapCheck())
        super(BandGapCheckFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


@explicit_serialize
class Interpolate(FireTaskBase):

    _fw_name = "Interpolate"

    def run_task(self, fw_spec):


class InterpolateFW(Firework):
    def __init__(self, structure, name="interpolate_poscar", parents=None, **kwargs):
        """
        Standard static calculation Firework for dielectric constants
        using DFPT.

        Args:
            structure (Structure): Input structure.
            name (str): Name for the Firework.
            vasp_cmd (str): Command to run vasp.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []

        interpolated_structures = structure.interpolate(end_structure=end_structure, nimages=10,
                                                        interpolate_lattices=True)
        for idx, interpolated_struct in enumerate(interpolated_structures):
            interpolated_poscar = interpolated_struct.to(fmt="poscar", filename="POSCAR_" + idx)

        t.append(CopyVaspOutputs(calc_loc=True))
        super(InterpolateFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


def get_wf_ferroelectric(structure, vasp_cmd=None, db_file=None):
    """
    Returns a workflow to calculate dipole moment.

    Firework 1 : write vasp input set for static calculation,
                 run vasp,
                 pass run location,
                 database insertion.

    Args:
        structure (Structure): input structure to be optimized and run
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.

    Returns:
        Workflow
    """
    fws = []

    fws.append(StaticFW(structure=structure, vasp_cmd=vasp_cmd, db_file=db_file))
    fws.append(NonSCFFW(structure=structure, vasp_cmd=vasp_cmd, db_file=db_file, parents=fws[0]))
    fws.append(BandGapCheckFW(structure=structure, parents=fws[1]))
    fws.append(InterpolateFW(structure=structure, parents=fws[2]))

    wfname = "{}:{}".format(structure.composition.reduced_formula, "dipole_moment")

    return Workflow(fws, name=wfname)


if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest
    from matmethods.vasp.workflows.presets.core import wf_ferroelectric

    structure = PymatgenTest.get_structure("Si")
    # wf = get_wf_ferroelectric(structure)
    my_wf = wf_ferroelectric(structure)
    # lp = LaunchPad()
    lp = LaunchPad.from_file(os.path.join("/global/homes/s/sbajaj/mm_installdir/config", "my_launchpad.yaml"))
    lp.reset('', require_password=False)
    lp.add_wf(my_wf)
