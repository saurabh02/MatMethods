# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the ferroelectric workflow
"""

import os

from fireworks import Workflow, LaunchPad

from matmethods.vasp.fireworks.core import StaticFW, NonSCFFW


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
