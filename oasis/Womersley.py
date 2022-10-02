from dolfin import UserExpression, Expression, Mesh, MeshFunction, SubsetIterator, ds, MPI, assemble, Constant, sqrt,FacetNormal, as_vector, SpatialCoordinate, Point
#from dolfin import *
import numpy as np
from scipy.special import jn

__all__ = ["compute_radius", "compute_boundary_geometry_acrn", "compute_area", "make_womersley_bcs_2", "fourier_coefficients","WomersleyComponent1"]


def error(text):
    print(text)

def x_to_r2(x, c, n):
    """Compute r**2 from a coordinate x, center point c, and normal vector n.

        r is defined as the distance from c to x', where x' is
        the projection of x onto the plane defined by c and n.
        """
    # Steps:
    # rv = x - c
    # rvn = rv . n
    # rp = rv - (rv . n) n
    # r2 = ||rp||**2
    #print (x,c)
    rv = Point(x[0],x[1],x[2]) - c
    rvn = rv.dot(n)
    rp = rv - Point(rvn*n[0], rvn*n[1], rvn*n[2])
    r2 = rp.dot(rp)

    return r2

def compute_radius(mesh, facet_domains, ind, center):
    d = len(center)
    it = SubsetIterator(facet_domains, ind)
    geom = mesh.geometry()
    #maxr2 = -1.0
    maxr2 = 0
    for i, facet in enumerate(it):
        ent = facet.entities(0)
        for v in ent:
            p = geom.point(v)
            r2 = sum((p[j] - center[j])**2 for j in range(d))
            maxr2 = max(maxr2, r2)
    r = MPI.max(MPI.comm_world, sqrt(maxr2))
    return r


def compute_boundary_geometry_acrn(mesh, dsi, normals):

    d = mesh.geometry().dim()
    x = SpatialCoordinate(mesh)

    # Compute area of boundary tesselation by integrating 1.0 over all facets
    A = assemble(1.0*dsi)
    #assert A > 0.0, "Expecting positive area, probably mismatch between mesh and markers!"
    if A == 0:
        return None

    # Compute barycenter by integrating x components over all facets
    c = [assemble(x[i]*dsi) / A for i in range(d)]

    # Compute average normal (assuming boundary is actually flat)
    n = normals#FacetNormal(mesh)
    ni = np.array([assemble(n[i]*dsi) for i in range(d)])
    n_len = np.sqrt(sum([ni[i]**2 for i in range(d)])) # Should always be 1!?
    normal = ni/n_len

    # Compute radius by taking max radius of boundary points
    # (assuming boundary points are on exact geometry)
    #r = compute_radius(mesh, facet_domains, ind, c)
    # This old estimate is a few % lower because of boundary discretization errors
    r = np.sqrt(A / np.pi)

    return A, c, r, normal


def compute_area(mesh, ind, facet_domains):
    # Some convenient variables
    assert facet_domains is not None
    dsi = ds(ind, domain=mesh, subdomain_data=facet_domains)

    # Compute area of boundary tesselation by integrating 1.0 over all facets
    A = assemble(Constant(1.0, name="one")*dsi)
    assert A > 0.0, "Expecting positive area, probably mismatch between mesh and markers!"
    return A


def fourier_coefficients(x, y, T, N):
    '''From x-array and y-spline and period T, calculate N complex Fourier coefficients.'''
    omega = 2*np.pi/T
    ck = []
    ck.append(1/T*simps(y(x), x))
    for n in range(1,N):
        c = 1/T*simps(y(x)*np.exp(-1j*n*omega*x), x)

        # Clamp almost zero real and imag components to zero
        if 1:
            cr = c.real
            ci = c.imag
            if abs(cr) < 1e-14: cr = 0.0
            if abs(ci) < 1e-14: ci = 0.0
            c = cr + ci*1j

        ck.append(2*c)
    return ck

class WomersleyComponent1(UserExpression):
    def __init__(self, args, degree, use_flat_profile):
        super().__init__(degree = degree)
        #print ('*** init called:', args, 'end.***')
        # Spatial args
        self.radius = args["radius"]
        self.center = Point(args["center"])
        self.normal = Point(args["normal"])
        self.normal_component = args["normal_component"]

        self.use_flat_profile = use_flat_profile
        self.y_mean = np.sqrt(2.0) / 2.0

        # Temporal args
        self.period = args["period"]
        if "Q" in args:
            assert "V" not in args, "Cannot provide both Q and V!"
            self.Qn = args["Q"]
            self.N = len(args["Q"])
        elif "V" in args:
            self.Vn = args["V"]
            self.N = len(self.Vn)
        else:
            error("Invalid transient data type, missing argument 'Q' or 'V'.")

        # Physical args
        self.nu = args["nu"]

        # Internal state
        self.scale_value = 1.0

        # Precomputation
        self._precompute_bessel_functions()
        self._all_r_dependent_coeffs = {}

        self.set_t(0)

    def _precompute_bessel_functions(self):
        '''Calculate the Bessel functions of the Womersley profile'''
        self.omega = 2 * np.pi / self.period
        self.ns = np.arange(1, self.N)

        # Allocate for 0...N-1
        alpha = np.zeros(self.N, dtype=np.complex)
        self.beta = np.zeros(self.N, dtype=np.complex)
        self.jn0_betas = np.zeros(self.N, dtype=np.complex)
        self.jn1_betas = np.zeros(self.N, dtype=np.complex)

        # Compute vectorized for 1...N-1 (keeping element 0 in arrays to make indexing work out later)
        alpha[1:] = self.radius * np.sqrt(self.ns * (self.omega / self.nu))
        self.beta[1:] = alpha[1:] * np.sqrt(1j**3)
        self.jn0_betas[1:] = jn(0, self.beta[1:])
        self.jn1_betas[1:] = jn(1, self.beta[1:])

    def _precompute_r_dependent_coeffs(self, y):
        pir2 = np.pi * self.radius**2
        # Compute intermediate terms for womersley function
        r_dependent_coeffs = np.zeros(self.N, dtype=np.complex)
        if hasattr(self, 'Vn'):
            #r_dependent_coeffs[0] = (self.Vn[0]/2.0) * (1 - y**2)
            r_dependent_coeffs[0] = self.Vn[0] * (1 - y**2)
            for n in self.ns:
                r_dependent_coeffs[n] = self.Vn[n] * (self.jn0_betas[n] - jn(0,
                                        self.beta[n]*y)) / (self.jn0_betas[n] - 1.0)
        elif hasattr(self, 'Qn'):
            r_dependent_coeffs[0] = (2*self.Qn[0]/pir2) * (1 - y**2)
            for n in self.ns:
                bn = self.beta[n]
                j0bn = self.jn0_betas[n]
                j1bn = self.jn1_betas[n]
                r_dependent_coeffs[n] = (self.Qn[n] / pir2) * (j0bn - jn(0, bn*y)) / (j0bn - (2.0/bn)*j1bn)
        else:
            error("Missing Vn or Qn!")
        return r_dependent_coeffs

    def _get_r_dependent_coeffs(self, y):
        "Look for cached womersley coeffs."
        key = y
        r_dependent_coeffs = self._all_r_dependent_coeffs.get(key)
        if r_dependent_coeffs is None:
            # Cache miss! Compute coeffs for this coordinate the first time.
            r_dependent_coeffs = self._precompute_r_dependent_coeffs(y)
            self._all_r_dependent_coeffs[key] = r_dependent_coeffs
        return r_dependent_coeffs

    def set_t(self, t):
        self.t = float(t) % self.period
        self.expnt = np.exp((self.omega * self.t * 1j) * self.ns)
        #print ('Womersley:set_t', t)

    def eval(self, value, x):
        if self.use_flat_profile:
            y = self.y_mean
        else:
            # Compute or get cached complex coefficients that only depend on r
            y = np.sqrt(x_to_r2(x, self.center, self.normal)) / self.radius
        coeffs = self._get_r_dependent_coeffs(y)

        # Multiply complex coefficients for x with complex exponential functions in time
        wom = (coeffs[0] + np.dot(coeffs[1:], self.expnt)).real

        # Scale by negative normal direction and scale_value
        value[0] = -self.normal_component * self.scale_value * wom

    def value_shape(self):
            return ()


def make_womersley_bcs(t, Q, mesh, nu, area, center, radius, normal,
                       v_degree, scale_to=None, coeffstype="Q",
                       N=1001, num_fourier_coefficients=21, **NS_namespace):
    """Generate a list of expressions for the components of a Womersley profile."""
    # Compute transient profile as interpolation of given coefficients
    period = max(t)
    transient_profile = UnivariateSpline(t, Q, s=0, k=1)

    # Compute fourier coefficients of transient profile
    timedisc = np.linspace(0, period, N)

    #from matplotlib.pyplot import plot, savefig
    #plot(timedisc, transient_profile(timedisc))
    #savefig("series.png")
    #sys.exit(0)
    Cn = fourier_coefficients(timedisc, transient_profile, period, num_fourier_coefficients)

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
            }
        args[coeffstype] = Cn
        expr = WomersleyComponent1('1', degree=v_degree)
        expr.init(args)
        wexpressions.append(expr)
        #wexpressions.append(WomersleyComponent1('1', degree=v_degree))
        #wexpressions[-1].init(args)

    return wexpressions

def make_womersley_bcs_2(period, q_mean, wave_form_filename, mesh, nu, area, center, radius, normal,
                       degree, flat_profile_at_intlet_bc=False, scale_to=None, coeffstype="Q", **NS_namespace):
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
            "nu": nu }

        # read the Fourier coefficients
        ck = []
        for line in open( wave_form_filename, 'r').readlines():
            if line.strip()[0] in ['#','!','/']:
                continue
            abn = line.split()
            ck.append( q_mean*(float(abn[0])-float(abn[1])*1j) )

        args[coeffstype] = ck
        expr = WomersleyComponent1(args, degree, flat_profile_at_intlet_bc)
        # print out the wave-form
        # val = [float()]*3
        # f = open('wom.txt', 'w')
        # for t in range(951):
        #     expr.set_t(t/1000.0)
        #     expr.eval(val, Point(center[0],center[1],center[2]))
        #     f.write(str(t) + str(' \t ') + str(val[0]) +'\n' )
        # f.close()
        # expr.set_t(0)
        wexpressions.append(expr)

    return wexpressions

