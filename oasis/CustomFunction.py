#///////////////////////////////////////////////////////////////
#// Artery.py
#// Copyright (C) 2018 Mehdi Najafi (mnuoft@gmail.com)
#//
#// Distribution of this library is not allowed in any form.
#///////////////////////////////////////////////////////////////

__author__ = "Mehdi Najafi <mnuoft@gmail.com>"
__date__ = "2018-07-12"
__copyright__ = "Copyright (C) 2018 " + __author__
__license__  = "Private"

#///////////////////////////////////////////////////////////////

from dolfin import UserExpression, Expression, Mesh, MeshFunction, SubsetIterator, ds, MPI, assemble, Constant, sqrt,FacetNormal, as_vector, SpatialCoordinate, Point
#from dolfin import *
import sys, importlib.util, importlib.machinery
import Womersley

__all__ = ["make_custom_function_bcs", "read_function_textfile","CustomFunction"]

def error(text):
    print(text)

class CustomFunction(UserExpression):
    def __init__(self, args, degree, file_name, use_flat_profile):
        super().__init__(degree = degree)
        #print ('*** init called:', args, 'end.***')
        # Spatial args
        self.radius = args["radius"]
        self.center = Point(args["center"])
        self.normal = Point(args["normal"])
        self.normal_component = args["normal_component"]

        # Temporal args
        self.period = args["period"]

        # Physical args
        self.nu = args["nu"]

        # Internal state
        self.scale_value = args["q_mean"]

        # zeroth order Fourier coefficient
        args["Q"] = [(1-0*1j)]
        self.expr = Womersley.WomersleyComponent1(args, degree, use_flat_profile)


        # load python module
        try:
            module_name = file_name.split('/')[-1]
            spec = importlib.util.spec_from_loader(module_name,importlib.machinery.SourceFileLoader(module_name,file_name))
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            self.function_module = module
        except ImportError as error:
            error ('<!> Cannot import the flow_rate function file:'+file_name)

        self.set_t(0)

    def set_t(self, t):
        self.t = float(t) % self.period
        self.expr.set_t(t)

    def eval(self, value, x):
        self.expr.eval(value,x)
        # Scale by negative normal direction and scale_value
        value[0] *= self.function_module.flow_rate(self.t, self.period)

    def value_shape(self):
            return ()


def make_custom_function_bcs(period, q_mean, wave_form_filename, mesh, nu, area, center, radius, normal,
                       degree, scale_to=None, coeffstype="Q", flat_profile_at_intlet_bc=False, **NS_namespace):
    """Generate a list of expressions for the components of a Womersley profile."""
    # Create Expressions for each direction
    wexpressions = []
    for ncomp in normal:
        args = {
            "radius": radius,
            "center": center,
            "normal": normal,
            "normal_component": ncomp,
            "period": period,
            "nu": nu,
            "q_mean":q_mean }
        expr = CustomFunction(args, degree, wave_form_filename, flat_profile_at_intlet_bc)
        wexpressions.append(expr)

    return wexpressions

