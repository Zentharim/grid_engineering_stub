import geopandas
from geopy.distance import geodesic as dist
from shapely.geometry import LineString
import argparse
from matplotlib import pyplot as plt
import gmsh
import os


class SmoothingError(Exception):
    def __init__(self, *args, **kwargs):
        default_message = "The smoothing degree is too high. Retry with a smaller one."

        if args or kwargs:
            super().__init__(*args, **kwargs)
        else:
            super().__init__(default_message)


def create_parser():
    l_parser = argparse.ArgumentParser()
    # defaults are for Fiumicino for testing reasons. For Qatar:
    # python bpy_test.py -i "./shp/qatar-line.shp" -s 56.441,25.533 56.525,24.8
    # or for test with bounding box:
    # python bpy_test.py -i "./shp/qatar-line.shp" -s 56.441,25.4 56.525,25.000 -b 56.332 25.4 56.525 25.000

    l_parser.add_argument("-i", "--input_shapefile", help="Shapefile to use as input",
                          default="./shp/fiumicino-line.shp")
    l_parser.add_argument("-b", "--bbox", nargs="+",
                          help="Bounding box for the shapefile. e.g. lon_up_sx lat_up_sx lon_down_dx lat_down_dx")
    l_parser.add_argument("-s", "--sea_points", nargs="+",
                          help="List of sea points where to close the dominion. e.g. lon1,lat1 lon2,lat2 ...",
                          default=["11.992,41.872", "12.080,41.680"])
    l_parser.add_argument("-t", "--threshold",
                          help="Threshold for the minimum distance to merge shapefile parts and remove double points",
                          default="0.001")
    l_parser.add_argument("-o", "--output_directory", help="Directory where to save the output",
                          default="./output")
    l_parser.add_argument("-d", "--smoothing_degree", help="Degrees for the smoothing",
                          default="0.01")
    l_args = l_parser.parse_args()
    return l_args


def smooth_corners(p_line, p_degree):
    offset = p_line.parallel_offset(p_degree, "right", join_style=1)
    offset2 = offset.parallel_offset(p_degree, "right", join_style=1)
    return list(offset2.coords)[1:-1]


def remove_doubles(p_list, p_threshold=0.001):
    counter = 0
    while counter < len(p_list)-1:
        if dist(p_list[counter], p_list[counter+1]).km < p_threshold:
            del p_list[counter+1]
        counter += 1
    return p_list


def close_dominion(p_input_shapefile, p_bbox, p_threshold, p_sea_points, p_smoothing_degree):
    file = p_input_shapefile

    if p_bbox:
        data_frame = geopandas.read_file(file, bbox=[float(coord) for coord in p_bbox])
    else:
        data_frame = geopandas.read_file(file)

    # Merging parts of shapefile
    matrix = [list(part.coords) for part in data_frame.geometry]
    matrix = merge_parts(matrix, float(p_threshold))

    # Creating dataframe
    linestrings = [LineString(line) for line in matrix]
    new_data_frame = geopandas.GeoDataFrame({"geometry": linestrings})

    # Creation of dominion
    sea_points = []

    for string in p_sea_points:
        a, b = string.split(",")
        coords = (float(a), float(b))
        sea_points.append(coords)

    # ok for 1,2, and 3 point of dominion bbox (no model over islands, we have a coastline)
    index_of_coastline = [idx for idx, part in enumerate(new_data_frame.geometry) if not part.is_ring][0]
    matrix = [list(part.coords) for part in new_data_frame.geometry]

    coastline = matrix[index_of_coastline]
    start_coast = coastline[0]
    end_coast = coastline[-1]

    dominion_closure = [end_coast]
    distances_from_end = [dist(end_coast, sea_point).km for sea_point in sea_points]
    copy_of_distances = distances_from_end.copy()
    while copy_of_distances:
        dominion_closure.append(sea_points[distances_from_end.index(min(copy_of_distances))])
        copy_of_distances.remove(min(copy_of_distances))
    dominion_closure.append(start_coast)
    dominion_closure = smooth_corners(LineString(dominion_closure), p_smoothing_degree)
    if not dominion_closure:
        raise SmoothingError
    coastline[-1:-1] = dominion_closure
    coastline[-1] = coastline[0]
    for index, line in enumerate(matrix):
        matrix[index] = remove_doubles(matrix[index])
    # matrix[index_of_coastline] = remove_doubles(coastline)
    linestrings = [LineString(line) for line in matrix]
    types = ["island"]*len(matrix)
    types[index_of_coastline] = "coastline"
    return geopandas.GeoDataFrame({"geometry": linestrings, "types": types}), len(dominion_closure)


def prompt_for_boundary_data(p_args):
    change_bbox = input("\nWould you like to change bounding box? Y or N ")
    change_sea_points = input("\nWould you like to change sea points? Y or N ")
    change_smoothing_degree = input("\nWould you like to change smoothing degree? Y or N ")
    sea_points = []
    if change_bbox.lower() == "y":
        bbox = [
            input("\nInsert longitude of NE point "),
            input("\nInsert latitude of NE point "),
            input("\nInsert longitude of SW point "),
            input("\nInsert latitude of SW point ")
        ]
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


def prompt_for_mesh_data():
    dist_max = input("\nInsert DistMax: ")
    dist_min = input("\nInsert DistMin: ")
    lc_max = input("\nInsert LcMax: ")
    lc_min = input("\nInsert LcMin: ")
    # dist_max = 0.1
    # dist_min = 0.05
    # lc_max = 0.03
    # lc_min = 0.005
    return {
        "DistMax": float(dist_max),
        "DistMin": float(dist_min),
        "LcMax": float(lc_max),
        "LcMin": float(lc_min)
    }


def arrange_parts(p_part_one, p_part_two, p_situation):
    results = [
        p_part_two[::-1] + p_part_one,
        p_part_two + p_part_one,
        p_part_one + p_part_two,
        p_part_one + p_part_two[::-1]
    ]
    return results[p_situation]


# FUNCTION TO MERGE PARTS OF SHAPEFILE
def merge_parts(p_matrix, p_threshold=0.001):
    len_matrix = len(p_matrix)
    index = 0
    while index < len(p_matrix):
        other_index = index + 1
        while other_index < len(p_matrix):
            part = p_matrix[index]
            other_part = p_matrix[other_index]

            distances = [
                dist(part[0], other_part[0]).km,
                dist(part[0], other_part[-1]).km,
                dist(part[-1], other_part[0]).km,
                dist(part[-1], other_part[-1]).km
            ]

            if min(distances) < p_threshold:
                p_matrix[index] = arrange_parts(part, other_part, distances.index(min(distances)))
            else:
                other_index += 1
                continue

            try:
                del p_matrix[other_index]
            except IndexError:
                pass
        index += 1
    if len_matrix == len(p_matrix):
        return p_matrix
    return merge_parts(p_matrix, p_threshold)


def plot_shapefile(p_dataframe):
    plt.figure()
    for shape in [list(part.coords) for part in p_dataframe.geometry]:
        x = [i[0] for i in shape]
        y = [i[1] for i in shape]
        plt.plot(x, y)
    plt.show()


def shapefile_to_geo(p_dataframe, p_mesh_params):

    gmsh.initialize()
    gmsh.model.add("my_mesh")
    coastline_index = list(p_dataframe.types).index("coastline")
    coastline = []
    counter = 0
    curve_loops = []
    for shape in p_dataframe.geometry:
        shape_points = []
        for point in list(shape.coords):
            shape_points.append(gmsh.model.geo.addPoint(point[0], point[1], 0))
        shape_lines = []
        for index in range(len(shape_points)):
            if index == len(shape_points)-1:
                shape_lines.append(gmsh.model.geo.addLine(shape_points[index], shape_points[0]))
            else:
                shape_lines.append(gmsh.model.geo.addLine(shape_points[index], shape_points[index+1]))
        curve_loops.append(gmsh.model.geo.addCurveLoop(shape_lines))
        if counter == coastline_index:
            coastline = shape_points
        counter += 1
    gmsh.model.geo.addPlaneSurface(curve_loops)
    gmsh.model.geo.synchronize()

    coastline_len = len(coastline)

    gmsh.model.mesh.field.add("Attractor", 1)
    gmsh.model.mesh.field.setNumbers(1, "NodesList", coastline[coastline_len//100:
                                                               coastline_len - p_mesh_params["extra_points"] -
                                                               coastline_len//100])
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "IField", 1)
    gmsh.model.mesh.field.setNumber(2, "DistMax", p_mesh_params["DistMax"])
    gmsh.model.mesh.field.setNumber(2, "DistMin", p_mesh_params["DistMin"])
    gmsh.model.mesh.field.setNumber(2, "LcMax", p_mesh_params["LcMax"])
    gmsh.model.mesh.field.setNumber(2, "LcMin", p_mesh_params["LcMin"])
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.write("./msh/temp.geo_unrolled")
    gmsh.model.mesh.generate(2)
    gmsh.write("./msh/temp.msh")
    gmsh.finalize()
    return


if __name__ == "__main__":
    args = create_parser()
    input_shapefile = args.input_shapefile
    dominion_file = args.output_directory + "/result_shapefile.shp"

    # Dominion creation
    is_fine = "N"
    while is_fine.lower() != "y":
        try:
            dominion_data_frame, extra_points = close_dominion(input_shapefile, args.bbox, args.threshold,
                                                               args.sea_points, float(args.smoothing_degree))
            plot_shapefile(dominion_data_frame)
            is_fine = input("\nIs this fine? Y or N ")
            # is_fine = "y"
        except SmoothingError:
            pass
        if is_fine.lower() != "y":
            args = prompt_for_boundary_data(args)

    # Mesh creation
    is_fine_mesh = "N"
    mesh_params = prompt_for_mesh_data()
    mesh_params["extra_points"] = extra_points
    while is_fine_mesh.lower() != "y":
        shapefile_to_geo(dominion_data_frame, mesh_params)

        gmsh.initialize()
        gmsh.open("./msh/temp.msh")
        gmsh.fltk.run()
        gmsh.finalize()
        # is_fine_mesh = "y"
        is_fine_mesh = input("\nIs this fine? Y or N ")
        if is_fine_mesh.lower() != "y":
            mesh_params = prompt_for_mesh_data()
            mesh_params["extra_points"] = extra_points

    try:
        dominion_data_frame.to_file(dominion_file)
    except ValueError:
        print("Check your bbox, it seems the coastline is not there.")

    os.replace("./msh/temp.msh", args.output_directory+"/result_mesh.msh")
