import bpy
import BlenderGIS.operators.io_import_shp as io_import_shp
import BlenderGIS.operators.io_export_shp as io_export_shp
import shapefile
import geopandas
from geopy.distance import geodesic as dist
from shapely.geometry import LineString
import argparse


def arrange_parts(i_part_one, i_part_two, situation):
    results = [
        i_part_two[::-1] + i_part_one,
        i_part_two + i_part_one,
        i_part_one + i_part_two,
        i_part_one + i_part_two[::-1]
    ]
    return results[situation]


# FUNCTION TO MERGE PARTS OF SHAPEFILE
def merge_parts(io_matrix, i_threshold=0.001):
    print("iteration")
    len_matrix = len(io_matrix)
    index = 0
    while index < len(io_matrix):
        other_index = index + 1
        while other_index < len(io_matrix):
            part = io_matrix[index]
            other_part = io_matrix[other_index]

            distances = [
                dist(part[0], other_part[0]).km,
                dist(part[0], other_part[-1]).km,
                dist(part[-1], other_part[0]).km,
                dist(part[-1], other_part[-1]).km
            ]

            if min(distances) < i_threshold:
                io_matrix[index] = arrange_parts(part, other_part, distances.index(min(distances)))
            else:
                other_index += 1
                continue

            try:
                del io_matrix[other_index]
            except IndexError:
                pass
        index += 1
    if len_matrix == len(io_matrix):
        return io_matrix
    return merge_parts(io_matrix, i_threshold)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # defaults are for Fiumicino for testing reasons. For Qatar:
    # python bpy_test.py -i "./shp/qatar-line.shp" -s 56.441,25.533 56.525,24.8
    # or for test with bounding box:
    # python bpy_test.py -i "./shp/qatar-line.shp" -s 56.441,25.4 56.525,25.000 -b 56.332 25.4 56.525 25.000

    parser.add_argument("-i", "--input_shapefile", help="Shapefile to use as input", default="./shp/fiumicino-line.shp")
    parser.add_argument("-b", "--bbox", nargs="+",
                        help="Bounding box for the shapefile. e.g. lon_up_sx lat_up_sx lon_down_dx lat_down_dx")
    parser.add_argument("-s", "--sea_points", nargs="+",
                        help="List of sea points where to close the dominion. e.g. lon1,lat1 lon2,lat2 ...",
                        default=["11.992,41.872", "12.080,41.680"])
    parser.add_argument("-t", "--threshold",
                        help="Threshold for the minimum distance to merge shapefile parts and remove double points",
                        default="0.001")
    parser.add_argument("-o", "--output_shapefile", help="Shapefile to use as output",
                        default="./shp/test_dominion.shp")
    args = parser.parse_args()
    file = args.input_shapefile
    print(args.bbox)
    if args.bbox:
        data_frame = geopandas.read_file(file, bbox=[float(coord) for coord in args.bbox])
    else:
        data_frame = geopandas.read_file(file)

    # Merging parts of shapefile
    matrix = [list(part.coords) for part in data_frame.geometry]
    matrix = merge_parts(matrix, float(args.threshold))

    # Creating dataframe
    linestrings = [LineString(line) for line in matrix]
    new_data_frame = geopandas.GeoDataFrame({"geometry": linestrings})

    # Creation of dominion
    sea_points = []

    for string in args.sea_points:
        a, b = string.split(",")
        coords = (float(a), float(b))
        sea_points.append(coords)

    print(sea_points)

    # ok for 1,2, and 3 point of dominion bbox (no model over islands, we have a coastline)
    index_of_coastline = [idx for idx, part in enumerate(new_data_frame.geometry) if not part.is_ring][0]
    matrix = [list(part.coords) for part in new_data_frame.geometry]

    coastline = matrix[index_of_coastline]
    start_coast = coastline[0]
    end_coast = coastline[-1]

    distances_from_start = [dist(start_coast, sea_point).km for sea_point in sea_points]
    print(distances_from_start)
    coastline.insert(0, sea_points[distances_from_start.index(min(distances_from_start))])

    distances_from_end = [dist(end_coast, sea_point).km for sea_point in sea_points]
    print(distances_from_end)
    copy_of_distances = distances_from_end.copy()
    while copy_of_distances:
        coastline.append(sea_points[distances_from_end.index(min(copy_of_distances))])
        copy_of_distances.remove(min(copy_of_distances))

    matrix[index_of_coastline] = coastline

    linestrings = [LineString(line) for line in matrix]
    dominion_data_frame = geopandas.GeoDataFrame({"geometry": linestrings})

    try:
        dominion_file = "./shp/test_dominion.shp"
        dominion_data_frame.to_file(dominion_file)
    except ValueError:
        print("Check your bbox, it seems the coastline is not there.")

# file = "./shp/fiumicino-line.shp"

# PYSHP
# sf = shapefile.Reader(file)
# shapes = sf.shapes()
# records = sf.records()
# for record, shape in zip(records, shapes):
#     print("record {} n of points {} type {}".format(record, len(shape.points), shape.shapeType))

# GEOPANDAS
# data_frame = geopandas.read_file(file)

# FILTRO SU BBOX (punto in alto a sinistra e punto in basso a destra)

# bbox = (56.378, 25.17114, 56.380, 25.4035)
# data_frame_bbox = geopandas.read_file(file, bbox=bbox)

# TEST
# points_sf = sum([len(shape.points) for shape in sf.shapes()])
# print(points_sf)
# points_df = sum([len(part.coords) for part in data_frame.geometry])
# print(points_df)


# TEST
# matrix = [list(part.coords) for part in data_frame.geometry]
# matrix = merge_parts(matrix)
# print(data_frame)
#
# # Creating dataframe
# linestrings = [LineString(line) for line in matrix]
# new_data_frame = geopandas.GeoDataFrame({"geometry": linestrings})
# print(new_data_frame)
# try:
#     new_data_frame.to_file("./shp/test.shp")
# except ValueError:
#     print("Check your bbox, it seems the coastline is not there.")
#
# # Creation of dominion
# print(list(new_data_frame.geometry[0].coords)[0])
# print(list(new_data_frame.geometry[0].coords)[-1])
# sea_points = [(11.992, 41.872), (12.080, 41.680)]
#
# # ok for 1,2, and 3 point of dominion bbox (no model over islands, we have a coastline)
# index_of_coastline = [idx for idx, part in enumerate(new_data_frame.geometry) if not part.is_ring][0]
# matrix = [list(part.coords) for part in new_data_frame.geometry]
#
# coastline = matrix[index_of_coastline]
# start_coast = coastline[0]
# end_coast = coastline[-1]
#
# distances_from_start = [dist(start_coast, sea_point).km for sea_point in sea_points]
# coastline.insert(0, sea_points[distances_from_start.index(min(distances_from_start))])
#
# distances_from_end = [dist(end_coast, sea_point).km for sea_point in sea_points]
# copy_of_distances = distances_from_end.copy()
# while copy_of_distances:
#     coastline.append(sea_points[distances_from_end.index(min(copy_of_distances))])
#     copy_of_distances.remove(min(copy_of_distances))
#
# matrix[index_of_coastline] = coastline
#
# linestrings = [LineString(line) for line in matrix]
# dominion_data_frame = geopandas.GeoDataFrame({"geometry": linestrings})
# print(dominion_data_frame)
# for part in dominion_data_frame.geometry:
#     print(part.is_ring)
#     print(part.coords[0])
#     print(part.coords[-1])
# try:
#     dominion_file = "./shp/test_dominion.shp"
#     dominion_data_frame.to_file(dominion_file)
# except ValueError:
#     print("Check your bbox, it seems the coastline is not there.")

# Subdivide dominion on Blender
# dominion_smooth_file = "./shp/test_dominion_smooth.shp"
# # PULIZIA BLENDER
# print(list(bpy.data.objects))
# print(bpy.data.objects)
# for element in bpy.data.objects:
#     bpy.data.objects.remove(element)
# bpy.data.meshes.remove(bpy.data.meshes["Cube"])
# print(bpy.data.objects)

# BLENDERGIS
# io_import_shp.register()
# bpy.ops.importgis.shapefile("EXEC_DEFAULT", filepath=dominion_file)
#
# bpy.ops.object.mode_set(mode='OBJECT')
# obj = bpy.context.active_object
# print("EHI")
# print(obj)
# bpy.ops.object.mode_set(mode='EDIT')
# bpy.ops.mesh.select_mode(type="VERT")
# bpy.ops.mesh.select_all(action='DESELECT')
# bpy.ops.object.mode_set(mode='OBJECT')
# obj.data.vertices[-1].select = True
# obj.data.vertices[-2].select = True
# bpy.ops.object.mode_set(mode='EDIT')
# bpy.ops.mesh.subdivide(number_cuts=160)
# bpy.ops.mesh.select_all(action='SELECT')
# bpy.ops.object.mode_set(mode='OBJECT')
#
# io_export_shp.register()
# bpy.ops.exportgis.shapefile("EXEC_DEFAULT", filepath=dominion_smooth_file)
# bpy.ops.exportgis.shapefile()
