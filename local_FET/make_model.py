import numpy as np
from utilities import *


def functional_expansion(x_min, x_max, order):

    flux_tally = openmc.Tally()
    flux_tally.scores = ['flux']
    expand_filter = openmc.SpatialLegendreFilter(order, 'x', x_min, x_max)
    flux_tally.filters.append(expand_filter)

    return openmc.Tallies([flux_tally])


def make_model():
    args = argument_parser()

    bc_map = {
        "left": ['vacuum', 'transmission', 'reflective', 'vacuum', 'reflective', 'reflective'],
        "middle": ['transmission', 'transmission', 'reflective', 'vacuum', 'reflective', 'reflective'],
        "right": ['transmission', 'vacuum', 'reflective', 'vacuum', 'reflective', 'reflective'],
    }

    x_pos = np.linspace(args.x_min, args.x_max, args.Nx + 1)

    fuel = make_materials(material_dict['UO2'], percent_type='ao')
    b4c = make_materials(material_dict['Boron Carbide'], percent_type='ao')
    materials = [fuel, b4c]

    cells = []
    x_mid = (args.x_max + args.x_min) / 2
    boron_length = 5
    for i in range(args.Nx):

        material = fuel if (x_pos[i] < x_mid - boron_length or x_pos[i] > x_mid) else b4c
        boundary_conditions = bc_map["left" if i == 0 else "right" if i == args.Nx - 1 else "middle"]
        region = make_box(x_dim=[x_pos[i], x_pos[i + 1]], y_dim=[args.y_min, args.y_max], z_dim=[args.z_min, args.z_max], boundary_conditions=boundary_conditions)
        cells.append(openmc.Cell(region=region, fill=material))

    slab_model = openmc.Model()
    slab_model.geometry = openmc.Geometry(openmc.Universe(cells=cells))
    slab_model.materials = openmc.Materials(materials)
    slab_model.tallies = functional_expansion(args.x_min, args.x_max, order=10)
    slab_model.settings = simulation_settings(args, space_dist=openmc.stats.Box(lower_left=(args.x_min, args.y_min, args.z_min), upper_right=(args.x_max, args.y_max, args.z_max)))

    return slab_model


if __name__ == "__main__":
    model = make_model().export_to_model_xml()
    openmc.run(openmc_exec="/home/ebny-walid-ahammed/github/cardinal/build/openmc/bin/openmc")

