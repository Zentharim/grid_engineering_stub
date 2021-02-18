import argparse
from matplotlib import pyplot as plt
import gmsh
import os
from clean_mesh.cleaner import Cleaner
from file_handler.geo_handler import GeoWriter
from file_handler.shp_handler import ShpHandler
from exceptions.exceptions import SmoothingError

# TODO: Bugfix "taglio" poligono sbagliato (esempio: taranto qui sotto)
# TODO: controllo qualitá griglia
# TODO: documentazione con sphinx
# TODO: gestire più coste contemporaneamente (mesh type 4, 5), simile al 3


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def create_parser():
    l_parser = argparse.ArgumentParser(prefix_chars='@')
    # python bpy_test.py @i "./shp/qatar-line.shp" @s 56.441,25.533 56.525,24.8
    # or for test with bounding box:
    # python bpy_test.py @i "./shp/qatar-line.shp" @s 56.441,25.4 56.525,25.000 @b 56.332 25.4 56.525 25.000
    # same with global file:
    # python bpy_test.py @i "./shp/GSHHS_h_L1.shp" @s 56.441,25.4 56.525,25.000 @b 56.332 25.4 56.525 25.000
    # Malta (type 2):
    # python bpy_test.py @b 13.9572 36.1917 15.3044 36.1906 15.3112 35.6713 13.9476 35.6133
    # @s 13.9572,36.1917 15.3044,36.1906 15.3112,35.6713 13.9476,35.6133 @d 0.1
    # Taranto type 3
    # python bpy_test.py @b 16.5583 39.3164 17.9920 40.3917 16.7835 40.8940 16.0474 39.7845

    l_parser.add_argument("@i", "@@input_shapefile", help="Shapefile to use as input",
                          default="./shp/GSHHS_h_L1.shp")
    l_parser.add_argument("@b", "@@bbox", nargs="+",
                          help="Bounding box for the shapefile. e.g. lon_up_sx lat_up_sx lon_down_dx lat_down_dx")
    l_parser.add_argument("@s", "@@sea_points", nargs="+",
                          help="List of sea points where to close the dominion. e.g. lon1,lat1 lon2,lat2 ...",
                          default=[])
                          # default=["11.992,41.872", "12.080,41.680"])
    l_parser.add_argument("@t", "@@threshold",
                          help="Threshold for the minimum distance to merge shapefile parts and remove double points",
                          default="0.001")
    l_parser.add_argument("@o", "@@output_directory", help="Directory where to save the output",
                          default="./output")
    l_parser.add_argument("@d", "@@smoothing_degree", help="Degrees for the smoothing",
                          default="0.01")
    l_parser.add_argument("@a", "@@attractors_islands", type=str2bool, nargs='?', const=True, default=True,
                          help="All islands found will be used as attractors")
    l_parser.add_argument("@c", "@@distance_coeff", type=float, nargs='?', default=1,
                          help="Coefficient")
    l_args = l_parser.parse_args()
    return l_args


def prompt_for_boundary_data(p_args):
    change_bbox = input("\nWould you like to change bounding polygon? Y or N ")
    change_sea_points = input("\nWould you like to change sea points? Y or N ")
    change_smoothing_degree = input("\nWould you like to change smoothing degree? Y or N ")
    sea_points = []
    if change_bbox.lower() == "y":
        bbox = []
        another_point = "Y"
        while another_point.lower() == "y":
            bbox.append(input("\nInsert longitude of bounding polygon point "))
            bbox.append(input("\nInsert latitude of bounding polygon point "))
            another_point = input("\nDo you need another bounding polygon point? Y or N ")
        setattr(args, "bbox", bbox)
    if change_sea_points.lower() == "y":
        another_point = "Y"
        while another_point.lower() == "y":
            point = "{},{}".format(
                input("\nInsert longitude of sea point "),
                input("\nInsert latitude of sea point ")
            )
            sea_points.append(point)
            another_point = input("\nDo you need another sea point? Y or N ")
        setattr(args, "sea_points", sea_points)
    if change_smoothing_degree.lower() == "y":
        smoothing_degree = input("\nInsert smoothing degree ")
        setattr(args, "smoothing_degree", smoothing_degree)
    return p_args


def prompt_for_mesh_data(p_debug=True):
    if p_debug:
        dist_max = 0.1
        dist_min = 0.05
        lc_max = 0.03
        lc_min = 0.005
    else:
        dist_max = input("\nInsert L_open: ")
        dist_min = input("\nInsert L_coast: ")
        lc_max = input("\nInsert Delta_open: ")
        lc_min = input("\nInsert Delta_coast: ")
    return {
        "DistMax": float(dist_max),
        "DistMin": float(dist_min),
        "LcMax": float(lc_max),
        "LcMin": float(lc_min)
    }


def plot_shapefile(p_dataframe):
    plt.figure()
    for shape in [list(part.coords) for part in p_dataframe.geometry]:
        x = [i[0] for i in shape]
        y = [i[1] for i in shape]
        plt.plot(x, y)
    plt.show()


def generate_mesh(file=None):
    # os.system("gmsh -2 {} -o ./msh/temp.msh".format(file))
    os.system("gmsh-mac/Gmsh.app/Contents/MacOS/gmsh -2 {} -o ./msh/temp.msh".format(file))
    # os.system("gmsh-win\\gmsh.exe -2 {} -o ./msh/temp.msh".format(file))


def read_data_from_msh(file):
    l_vertices = []
    l_triangles = []
    l_quads = []
    state = 0
    with open(file, "r") as a_file:
        for line in a_file:
            if state == 0:
                if "$Nodes" in line:
                    state = 1
                    continue
            if state == 1:
                state = 2
                continue
            if state == 2:
                parts = line.split()
                if len(parts) == 1:
                    # skipping not data lines
                    continue
                if len(parts) == 4:
                    # vertices
                    l_vertices.append([parts[1], parts[2], parts[3]])
                    continue
                elif parts[1:5] == ["2", "2", "0", "1"]:
                    # triangles
                    l_triangles.append(parts[5:])
                elif parts[1] == "3" and len(parts) > 4:
                    # quads
                    l_quads.append(parts[-4, :])

    return l_vertices, l_triangles, l_quads


def msh_to_ww3(file, p_clean_vertices, p_clean_triangles):
    l_ww3_file = ".".join(file.split(".")[:-1])+"_ww3.msh"
    with open(l_ww3_file, "w") as file:
        file.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n")
        file.write("{}\n".format(len(p_clean_vertices)))
        for index, vertex in enumerate(p_clean_vertices):
            file.write("{} {} {} 0\n".format(index+1, vertex[0], vertex[1]))
        file.write("$EndNodes\n$Elements\n")
        file.write("{}\n".format(len(p_clean_triangles)))
        for index, triangle in enumerate(p_clean_triangles):
            file.write("{} 2 2 0 1 {} {} {}\n".format(index+1, triangle[0], triangle[1], triangle[2]))
        file.write("$EndElements\n")
    return l_ww3_file


def ww3_to_grd(file):
    buffer = ""
    state = 0
    nodes = []
    elements = []
    with open(file, "r") as a_file:
        for line in a_file:
            if state == 0:
                buffer += line
                if "$Nodes" in line:
                    state = 1
                    continue
            if state == 1:
                if "$Elements" in line:
                    state = 2
                    del nodes[0]
                    del nodes[-1]
                    continue
                nodes.append(line.split())
            if state == 2:
                if "$EndElements" in line:
                    del elements[0]
                    continue
                parts = line.split()
                elements.append(line.split())
    l_grd_file = "_".join(file.split("_")[:-1]) + ".grd"
    with open(l_grd_file, "w") as file:
        for node in nodes:
            file.write("1 {} 0 {} {}\n".format(node[0], node[1], node[2]))
        file.write("\n")
        for element in elements:
            file.write("2 {} 0 3 {} {} {}\n".format(element[0], element[5], element[6], element[7]))
        file.write("\n")
    return l_grd_file


if __name__ == "__main__":
    debug = True
    args = create_parser()
    input_shapefile = args.input_shapefile
    dominion_file = args.output_directory + "/result_shapefile.shp"
    mesh_type = None

    # Dominion creation
    is_fine = "N"
    while is_fine.lower() != "y":
        try:
            shp_handler = ShpHandler(input_shapefile, args.bbox)
            shp_handler.prepare_to_close()
            dominion_data_frame, extra_points = shp_handler.close_dominion(args.threshold, args.sea_points,
                                                                           float(args.smoothing_degree))
            mesh_type = shp_handler.mesh_type
            plot_shapefile(dominion_data_frame)
            is_fine = "y" if debug else input("\nIs this fine? Y or N ")
        except SmoothingError as e:
            print(e)
        if is_fine.lower() != "y":
            args = prompt_for_boundary_data(args)

    # Mesh creation
    is_fine_mesh = "N"
    mesh_params = prompt_for_mesh_data(debug)
    mesh_params["extra_points"] = extra_points
    mesh_params["islands_attractors"] = args.attractors_islands
    mesh_params["distance_coeff"] = args.distance_coeff
    mesh_params["mesh_type"] = mesh_type
    writer = GeoWriter(dominion_data_frame, mesh_params)
    while is_fine_mesh.lower() != "y":
        result_code = writer.shapefile_to_geo()
        if result_code:
            is_fine_mesh = "n"
        else:
            generate_mesh("./msh/temp.geo_unrolled")

            gmsh.initialize()
            gmsh.open("./msh/temp.msh")
            gmsh.fltk.run()
            gmsh.finalize()
            is_fine_mesh = "y" if debug else input("\nIs this fine? Y or N ")
        if is_fine_mesh.lower() != "y":
            mesh_params = prompt_for_mesh_data(debug)
            mesh_params["extra_points"] = extra_points
            mesh_params["islands_attractors"] = args.attractors_islands
            mesh_params["distance_coeff"] = args.distance_coeff
            mesh_params["mesh_type"] = mesh_type
            writer.clean(mesh_params)

    try:
        dominion_data_frame.to_file(dominion_file)
    except ValueError:
        print("Check your bbox, it seems the coastline is not there.")

    msh_result = args.output_directory+"/result_mesh.msh"
    os.replace("./msh/temp.msh", msh_result)

    # Quality checks
    vertices, triangles, quads = read_data_from_msh(msh_result)
    cleaner = Cleaner(vertices, triangles, quads)
    cleaner.clean()
    clean_vertices, clean_triangles = cleaner.get_data()

    ww3_file = msh_to_ww3(msh_result, clean_vertices, clean_triangles)
    grd_file = ww3_to_grd(ww3_file)
