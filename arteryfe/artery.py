import sys
import numpy as np

from dolfinx import mesh, fem, io
import ufl

class Artery(object):
    """
    Represents an artery whose flow rate and area are calculated using the
    1D system of blood flow equations in conservation form.

    Arguments
    ---------

    root_vessel : boolean
        True if the artery is a root vessel in the artery network (has no
        parent)
    end_vessel : boolean
        True if the artery is a terminal vessel in the artery network (has no
        daughters)
    rc : float
        Characteristic radius (length)
    qc : float
        Characteristic flow
    Ru : float
        Upstream radius
    Rd : float
        Downstream radius
    L : float
        Vessel length
    k1 : float
    """
    # Rest of the class implementation
