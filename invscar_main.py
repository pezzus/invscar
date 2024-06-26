from firedrake import *
import numpy as np
from scipy.optimize import minimize

__all__ = ["invscar"]

def invscar(**params):
    """Solve the inverse scar problem.

    Parameters
    ----------

        geo         : Geometry type ['2d' (default), '3d', mesh]
        mu          : Stiffness [float] (default: 1.0)
        af          : Fiber direction [vector of floats] (default: [1,0])
        alpha_gt    : Ground truth contractility [Expression or Function]
                      Default is a circular inclusion of radius 0.2
                      at the center of the domain
        noise_level : Noise for corrupting data [float] (default: 1e-2)
        J_fide      : Type of observation ['full' (default), 'bcs']
        J_regu      : Type of regularization ['L2', 'H1' (default), 'TV']
        lmbda       : Regularization weight [float] (default: 4e-8)
        lower_bnd   : Lower bound [float] (default: 0)
        upper_bnd   : Upper bound [float] (default: inf)
        bfgs_disp   : Verbose BFGS output (default: False)
        ofile_name  : Output file name, pvd format (default: None)
        show_plot   : Show 2D plot (default: False)
    
    Return
    ------
    Structure `res` with
        `res.u`     : final displacement
        `res.alpha` : final contractility
        `res.res`   : BFGS results (iterations count, ...)
        `res.u_gt`  : ground truth displacement
        `res.alpha_gt`: ground truth contractility
        `res.J_fid` : data fidelity term
        `res.J_reg` : regularization term
    """

    #### Geometry
    geom = params.get('geom', '2d')

    if geom == '2d':
        Nx,Ny =  40, 40
        Lx,Ly = 1.0,1.0
        mm = RectangleMesh(Nx, Ny, Lx, Ly, quadrilateral=True)
    elif geom == '3d':
        Nx,Ny,Nz =  10, 10, 10
        Lx,Ly,Lz = 1.0,1.0,1.0
        mm = BoxMesh(Nx,Ny,Nz,Lx,Ly,Lz,hexahedral=False)
    else:
        mm = geom

    dim = mm.geometric_dimension()

    #### Parameters
    mu = params.get('mu', Constant(1.0,name="mu"))
    af = params.get('af', Constant(np.eye(1,dim)[0],name="fibers"))

    cnt = Constant(0.5*np.ones(dim))
    rad = Constant(0.2)
    alpha_gt = lambda x: conditional(sqrt((x-cnt)**2) < rad, 0, 1)
    alpha_gt = params.get('alpha_gt', alpha_gt)
    
    noise_level = params.get('noise_level', 1e-2)
    J_fide      = params.get('J_fide',      'full')
    J_regu      = params.get('J_regu',      'H1')
    lmbda       = params.get('lmbda', Constant(4e-8,name="lmbda"))

    lower_bnd   = params.get('lower_bnd', 0.0)
    upper_bnd   = params.get('upper_bnd', np.inf)
    
    bfgs_disp   = params.get('bfgs_disp',   False)
    ofile_name  = params.get('ofile_name',  None)
    show_plot   = params.get('show_plot',   False)

    # function space
    V = VectorFunctionSpace(mm, "P", 1)
    Q = FunctionSpace(mm, "P", 1)

    # boundary conditions
    #
    bcs = [DirichletBC(V.sub(0),Constant(0.0),1),
            DirichletBC(V.sub(1),Constant(0.0),3)]
    if dim == 3:
        bcs += [DirichletBC(V.sub(2),Constant(0.0),5)]

    #### Forward problem
    u     = Function(V, name="displacement")
    p     = Function(V, name="adjoint")
    alpha = Function(Q, name="alpha")

    I = Identity(dim)
    F = I + grad(u)
    W = mu/2 * (inner(F,F) - 2*ln(det(F))) * dx + alpha/2*inner(F*af,F*af) * dx

    G = derivative(W,u)
    fwd_prob   = NonlinearVariationalProblem(G, u, bcs,
                    form_compiler_parameters={'quadrature_degree': 2})
    fwd_solver = NonlinearVariationalSolver(fwd_prob)
                    #solver_parameters={'snes_monitor':None,'snes_rtol':1e-10})

    #### generate ground truth
    x = SpatialCoordinate(mm)
    alpha.interpolate(alpha_gt(x))
    fwd_solver.solve()

    ud = Function(V, name="displ_data")
    ud.assign(u)
    
    u_true = Function(V, name="displ_true")
    a_true = Function(Q, name="alpha_true")
    u_true.assign(u)
    a_true.assign(alpha)

    # add white noise
    u_noise = u.copy(deepcopy=True)
    sigma_u = np.std(ud.dat.data_ro)
    noise = noise_level * sigma_u
    rng = np.random.default_rng()

    ud.dat.data[:] += noise * (rng.normal(size=ud.dat.data_ro.size).reshape((-1,dim)))
    for bc in bcs: bc.apply(ud)

    #### Adjoint problem
    v = TestFunction(V)
    dG = adjoint(derivative(G,u))
    La = -dot(u-ud,v) * {'full': dx, 'bcs': ds}[J_fide]

    adj_prob   = LinearVariationalProblem(dG, La, p, bcs)
    adj_solver = LinearVariationalSolver(adj_prob)

    #### Optimality condition
    R_L2 = lambda g: 0.5*g**2*dx
    R_H1 = lambda g: 0.5*inner(grad(g),grad(g))*dx
    R_TV = lambda g: sqrt(1e-2 + inner(grad(g),grad(g)))*dx

    J_ful = lambda d: 0.5 * dot(d,d) * dx
    J_bcs = lambda d: 0.5 * dot(d,d) * ds
    J_f = {'full': J_ful, 'bcs': J_bcs}[J_fide]
    J_R = {'L2': R_L2, 'H1': R_H1, 'TV': R_TV}[J_regu]

    J = J_f(u-ud) + lmbda * J_R(alpha)

    dJ = lmbda * derivative(J_R(alpha),alpha) + derivative(action(G, p), alpha)

    #### Optimization
    alpha.interpolate(Constant(1.0))
    x0 = alpha.dat.data_ro.copy()

    lb = np.empty_like(x0)
    lb.fill(lower_bnd)
    ub = np.empty_like(x0)
    ub.fill(upper_bnd)
    bnds = np.array([lb, ub]).T

    def Jfun(x):
        alpha.dat.data[:] = x
        fwd_solver.solve()
        return assemble( J )
    
    def dJfun(x,no_forward=True):
        alpha.dat.data[:] = x

        # BFGS calls J before dJ, so we have the updated forward run already
        # we can force the recomputation here
        if not no_forward:
            fwd_solver.solve()

        adj_solver.solve()
        grd = assemble(dJ).dat.data_ro
        return grd
    
    if ofile_name is not None:
        ofile = VTKFile(ofile_name)
        ofile.write(alpha,u)

        def callback_fun(intermediate_result):
            ofile.write(alpha,u)
    else:
        callback_fun = None

    res = minimize(Jfun, x0, jac=dJfun, tol=1e-10, bounds=bnds,
                   method='L-BFGS-B',
                   callback=callback_fun,
                   options={'disp': bfgs_disp})


    class InvScarResult: pass
    result = InvScarResult()
    result.u     = u
    result.alpha = alpha
    result.res   = res
    result.u_gt      = u_true
    result.alpha_gt  = a_true
    result.J_fid     = assemble(J_f(u-ud))
    result.J_reg     = assemble(J_R(alpha))

    #print(res)

    if dim == 2 and show_plot:
        import matplotlib.pyplot as plt
        tricontourf(alpha)
        plt.show()

    return result

if __name__ == "__main__":

    # quick test
    invscar(J_fide = 'full', ofile_name='out_example.pvd')