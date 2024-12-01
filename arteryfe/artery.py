import sys
import numpy as np
from math import pi

from dolfinx import mesh, fem, io
#from dolfinx.cpp.fem import FiniteElement_float64 as FiniteElement
from dolfinx.fem import FunctionSpace, Expression, Function, dirichletbc
from dolfinx.mesh import create_interval
from dolfinx.fem.petsc import NonlinearProblem
import ufl
from ufl import split, dx, sqrt, derivative, grad, TestFunction, FiniteElement
from petsc4py import PETSc
from mpi4py import MPI

DOLFINX_EPS = 3.0e-16
def near(x, value, tol=1e-14):
    return np.isclose(x, value, atol=tol)

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
     	First constant from the relation Eh/r0
    k2 : float
        Second constant from the relation Eh/r0
    k3 : float
        Third constant from the relation Eh/R0
    rho : float
        Density of blood
    Re : float
        Reynolds' number
    nu : float
        Viscosity of blood
    p0 : float
        Diastolic pressure
    """
    def __init__(self, root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3,
                 rho, Re, nu, p0):
        self.root_vessel = root_vessel
        self.end_vessel = end_vessel
        self.rc = rc
        self.qc = qc
        self.Ru = Ru
        self.Rd = Rd
        self.L = L
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.rho = rho
        self.Re = Re
        self.nu = nu
        self.p0 = p0

    def define_geometry(self, Nx, Nt, T, N_cycles):
        """
        Initializes the artery geometry by creating the spatial refinement,
        temporal refinement, and dolfinx objects.

        Arguments
        ---------
        Nx : int
            Number of spatial points per artery
        Nt : int
            Number of time steps per cardiac cycle
        T : float
            Duration of one cardiac cycle
        N_cycles: int
            Number of cardiac cycles in the simulation
        """
        self.Nx = Nx
        self.Nt = Nt
        self.T = T
        self.N_cycles = N_cycles

        self.dx = self.L / self.Nx
        self.dt = self.T / self.Nt

        # Step for boundary condition computations, starting at normal size
        self.dex = self.dx
        self.db = np.sqrt(self.nu * self.T / (2 * np.pi))

        # Mesh and Function Spaces
        self.mesh = create_interval(MPI.COMM_WORLD, self.Nx, [0, self.L])
        element = ufl.FiniteElement("CG", self.mesh.ufl_cell(), 1)
        self.V = FunctionSpace(self.mesh, element)
        self.V2 = FunctionSpace(self.mesh, element*element)
        
        # Initial vessel-radius and derived quantities
        x = ufl.SpatialCoordinate(self.mesh)
        
        # Use ufl expressions directly
        r0_expr = self.Ru * (self.Rd / self.Ru) ** (x[0] / self.L)
        self.r0 = fem.Function(self.V)
        self.r0.interpolate(lambda x: self.Ru * (self.Rd / self.Ru) ** (x[0] / self.L))

        self.A0 = fem.Function(self.V)
        self.A0.interpolate(lambda x: pi * (self.Ru * (self.Rd / self.Ru) ** (x[0] / self.L))**2)

        f_expr = (4.0 / 3.0) * (self.k1 * ufl.exp(self.k2 * r0_expr) + self.k3)
        self.f = fem.Function(self.V)
        self.f.interpolate(lambda x: (4.0 / 3.0) * (self.k1 * np.exp(self.k2 * (self.Ru * (self.Rd / self.Ru) ** (x[0] / self.L))) + self.k3))

        dfdr_expr = (4.0 / 3.0) * self.k1 * self.k2 * ufl.exp(self.k2 * r0_expr)
        self.dfdr = fem.Function(self.V)
        self.dfdr.interpolate(lambda x: (4.0 / 3.0) * self.k1 * self.k2 * np.exp(self.k2 * (self.Ru * (self.Rd / self.Ru) ** (x[0] / self.L))))

        drdx_expr = (np.log(self.Rd / self.Ru) / self.L) * r0_expr
        self.drdx = fem.Function(self.V)
        self.drdx.interpolate(lambda x: (np.log(self.Rd / self.Ru) / self.L) * (self.Ru * (self.Rd / self.Ru) ** (x[0] / self.L)))



    def define_solution(self, q0, theta=0.5, bc_tol=1.e-14):
        """
        Defines FEniCS Function objects, boundary conditions and variational
        form of the problem.

        Arguments
        ---------
        q0 : float
            Initial flow rate in the root vessel
        theta : float
            Weighting parameter for the Crank-Nicolson method, in the interval
            [0, 1]
        bc_tol : float
            Inlet and outlet boundary thickness (tolerance)
        """
        # Initial flow value
        self.q0 = q0

        # Crank-Nicolson parameter
        self.theta = theta

        # Trial function
        self.U = Function(self.V2)
        A, q = split(self.U)

        # Test functions
        v1, v2 = split(TestFunction(self.V2))

        # Current solution, initialised
        self.Un = Function(self.V2)
        A0_init = Function(self.V)
        A0_init.interpolate(lambda x: np.pi * (self.Ru * (self.Rd / self.Ru) ** (x[0] / self.L))**2)

        self.Un.interpolate(lambda x: np.vstack((A0_init.x.array, np.full_like(x[0], q0))).T)

        # Current pressure, initialised
        self.pn = Function(self.V)
        self.pn.interpolate(lambda x: np.full_like(x[0], self.p0))

        # Boundary conditions
        def inlet_bdry(x):
            return np.isclose(x[0], 0, atol=bc_tol)

        def outlet_bdry(x):
            return np.isclose(x[0], self.L, atol=bc_tol)

        if self.root_vessel:
            self._q_in = Function(self.V)
            self._q_in.interpolate(lambda x: np.full_like(x[0], self.q0))
            bc_inlet = fem.dirichletbc(self._q_in, self.V2.sub(1))
        else:
            self._U_in = Function(self.V2)
            self._U_in.interpolate(lambda x: np.vstack((A0_init.x.array[0], q0)).T)
            bc_inlet = dirichletbc(self._U_in, self.V2)

        if self.end_vessel:
            self._A_out = Function(self.V)
            self._A_out.interpolate(lambda x: np.full_like(x[0], A0_init.x.array[-1]))
            bc_outlet = dirichletbc(self._A_out, self.V2.sub(0))
        else:
            self._U_out = Function(self.V2)
            self._U_out.interpolate(lambda x: np.vstack((A0_init.x.array[-1], q0)).T)
            bc_outlet = dirichletbc(self._U_out, self.V2)

        self.bcs = [bc_inlet, bc_outlet]

        # Define variational form
        U_v_dx = A * v1 * ufl.dx + q * v2 * ufl.dx
        Un_v_dx = self.Un[0] * v1 * ufl.dx + self.Un[1] * v2 * ufl.dx

        F2_v2_ds = (q**2 / (A + DOLFINX_EPS) 
                    + self.f * ufl.sqrt(self.A0 * (A + DOLFINX_EPS))) * v2 * ufl.ds
        F2_dv2_dx = (q**2 / (A + DOLFINX_EPS) 
                     + self.f * ufl.sqrt(self.A0 * (A + DOLFINX_EPS))) * ufl.grad(v2)[0] * ufl.dx
        dF_v_dx = ufl.grad(q)[0] * v1 * ufl.dx + F2_v2_ds - F2_dv2_dx

        Fn_v_ds = (self.Un[1]**2 / self.Un[0] 
                   + self.f * ufl.sqrt(self.A0 * self.Un[0])) * v2 * ufl.ds
        Fn_dv_dx = (self.Un[1]**2 / self.Un[0] 
                    + self.f * ufl.sqrt(self.A0 * self.Un[0])) * ufl.grad(v2)[0] * ufl.dx
        dFn_v_dx = ufl.grad(self.Un[1])[0] * v1 * ufl.dx + Fn_v_ds - Fn_dv_dx

        S_v_dx = (-2 * ufl.sqrt(pi) / self.db / self.Re * q / ufl.sqrt(A + DOLFINX_EPS) * v2 * ufl.dx 
                  + (2 * ufl.sqrt(A + DOLFINX_EPS) * (ufl.sqrt(pi) * self.f + ufl.sqrt(self.A0) * self.dfdr) 
                     - A * self.dfdr) * self.drdx * v2 * ufl.dx)

        Sn_v_dx = (-2 * ufl.sqrt(pi) / self.db / self.Re * self.Un[1] / ufl.sqrt(self.Un[0]) * v2 * ufl.dx 
                   + (2 * ufl.sqrt(self.Un[0]) * (ufl.sqrt(pi) * self.f + ufl.sqrt(self.A0) * self.dfdr) 
                
                      - self.Un[0] * self.dfdr) * self.drdx * v2 * ufl.dx)
        # Variational form
        self.variational_form = U_v_dx 
        - Un_v_dx 
        + self.dt * self.theta * dF_v_dx 
        + self.dt * (1 - self.theta) * dFn_v_dx 
        - self.dt * self.theta * S_v_dx 
        - self.dt * (1 - self.theta) * Sn_v_dx


    def solve(self):
        """
        Calls FEniCS's solve() function for the variational form.
        """
        F = self.variational_form
        J = derivative(F, self.U)
        prob = NonlinearProblem(F, self.U, self.bcs, J=J)

        # Create the PETSc SNES (Scalable Nonlinear Equations Solvers) solver
        solver = PETSc.SNES().create(MPI.COMM_WORLD)

        # Set the nonlinear problem for the solver
        solver.setFunction(prob.F, prob.b)
        solver.setJacobian(prob.J, prob.A)

        # Set solver options (similar to setting parameters in DOLFIN)
        solver.getKSP().setType("preonly")       # Set KSP type (e.g., 'preonly' for direct solvers)
        solver.getKSP().getPC().setType("lu")    # Use LU decomposition
        solver.getKSP().getPC().setFactorSolverType("mumps")  # Use MUMPS as the LU solver
        solver.setTolerances(max_it=1000)        # Set maximum number of iterations

        # Solve the nonlinear problem
        solver.solve(None, self.U.vector)



    def update_solution(self):
        """
        Stores current solution U as previous solution Un.
        """
        self.Un.assign(self.U)


    def update_pressure(self):
        """
        Calculates pressure.
        """
        self.pn.assign(Expression('p0 + f*(1-sqrt(A0/A))', degree=2, p0=self.p0,
                                  f=self.f, A0=self.A0, A=self.Un.split(True)[0]
                                  ))


    def compute_pressure(self, f, A0, A):
        """
        Computes pressure.

        Arguments
        ---------
        f : numpy.array
            Elasticity relation
        A0: numpy.array
            Area at rest
	    A : numpy.array
            Area

        Returns
        -------
        return : float
            Pressure at the outlet
        """
        return self.p0 + f*(1-np.sqrt(A0/A))


    def compute_outlet_pressure(self, A):
        """
        Computes pressure at the outlet.

        Arguments
        ---------
        A : float
            Area value at the outlet

        Returns
        -------
        return : float
            Pressure at the outlet
        """
        return self.p0 + self.f(self.L)*(1-np.sqrt(self.A0(self.L)/A))


    def CFL_term(self, x, A, q):
        """
        Computes the CFL number.

        Arguments
        ---------
        x : float
            Point at which the condition is to be checked
        A : float
            Area value at x
        q : float
            Flow rate value at x

        Returns
        -------
        return : float
            CFL number
        """
        return 1/np.abs(q/A+np.sqrt(self.f(x)/2/self.rho\
                                   *np.sqrt(self.A0(x)/A)))


    def check_CFL(self, x, A, q):
        """
        Checks the CFL condition.

        Arguments
        ---------
        x : float
            Point at which the condition is to be checked
        A : float
            Area value at x
        q : float
            Flow rate value at x

        Returns
        -------
        return : boolean
            True if the CFL condition is fulfilled
        """
        return self.dt/self.dex < self.CFL_term(x, A, q)


    def adjust_dex(self, x, A, q, margin=0.05):
        """
        Adjusts spatial step at a bifurcation to respect the CFL condition.

        Arguments
        ---------
        x : float
            Point at which the condition is to be checked
        A : float
            Area value at x
        q : float
            Flow rate value at x
        margin : float
            Margin of CFL number
        """
        M = self.CFL_term(x, A, q)
        self.dex = (1+margin)*self.dt/M


    @property
    def q_in(self):
        """
        Inlet flow rate (only for a root artery)
        """
        return self._q_in.value

    @q_in.setter
    def q_in(self, value):
        self._q_in.value = value


    @property
    def U_in(self):
        """
        Inlet boundary conditions (only for non-root arteries (daughters))
        """
        return np.array([self._U_in.A, self._U_in.q])

    @U_in.setter
    def U_in(self, U):
        self._U_in.A = U[0]
        self._U_in.q = U[1]


    @property
    def A_out(self):
        """
        Outlet area (only in use for end arteries)
        """
        return self._A_out.value

    @A_out.setter
    def A_out(self, value):
        self._A_out.value = value


    @property
    def U_out(self):
        """
        Outlet boundary conditions (only for parent arteries)
        """
        return np.array([self._U_out.A, self._U_out.q])

    @U_out.setter
    def U_out(self, U):
        self._U_out.A = U[0]
        self._U_out.q = U[1]
        
