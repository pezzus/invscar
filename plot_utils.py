from firedrake import *

import numpy as np
import matplotlib.pyplot as plt

__all__ = ["do_plot_2d","do_plot_3d"]

def do_plot_2d(res):
    fig, axs = plt.subplots(1,2,figsize=plt.figaspect(0.4))

    cols = tricontourf(res.alpha,20,axes=axs[0],cmap='inferno')
    tricontour(res.alpha,[0.5],axes=axs[0],colors=['w'])
    cb = fig.colorbar(cols)
    cb.ax.get_yaxis().labelpad = 15

    axs[0].set_title('Reconstruction')
    axs[0].set_axis_off()
    axs[0].set_aspect('equal')

    # deform the mesh
    mm = res.u.function_space().mesh()
    coord_orig = mm.coordinates.dat.data.copy()
    mm.coordinates.dat.data[:] += res.u.dat.data_ro

    cols = tricontourf(res.u, 20, axes=axs[1],cmap='inferno')
    triplot(mm, axes=axs[1], interior_kw=dict(linewidth=0.5,edgecolor='w'), boundary_kw=dict(linewidth=0))
    axs[1].set_title('Deformed mesh')
    axs[1].set_axis_off()
    axs[1].set_aspect('equal')

    fig.tight_layout()

    mm.coordinates.dat.data[:] = coord_orig

def do_plot_3d(res,backend='html'):
    import pyvista as pv
    from vtkmodules.vtkCommonDataModel import VTK_TETRA
    from firedrake.output import paraview_reordering

    # HTML visualization
    u = res.u
    mm = u.function_space().mesh()
    pts  = mm.coordinates.dat.data_ro + u.dat.data_ro
    #reorder = output.paraview_reordering.vtk_lagrange_hex_reorder(u.ufl_element())
    reorder = paraview_reordering.vtk_lagrange_tet_reorder(u.ufl_element())
    elm  = mm.coordinates.cell_node_map().values[:,reorder]

    mesh_pv = pv.UnstructuredGrid({VTK_TETRA: elm}, pts)
    #mesh_pv.point_data['u'] = u.dat.data_ro
    mesh_pv.point_data['alpha'] = res.alpha.dat.data_ro
    contours = mesh_pv.contour([0.5])

    plotter = pv.Plotter(notebook=True,window_size=[600, 600])
    plotter.add_mesh(mesh_pv, show_scalar_bar=False, show_edges=False, opacity=0.5, color='white')
    plotter.add_mesh(contours, show_scalar_bar=False, color='red')

    plotter.show(jupyter_backend=backend)