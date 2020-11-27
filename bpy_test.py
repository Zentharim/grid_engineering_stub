# import bpy
# import BlenderGIS.operators.io_import_shp as io_import_shp
import shapefile
import geopandas
from math import hypot
from geopy.distance import geodesic as dist
from shapely.geometry import LineString

# PULIZIA BLENDER
# print(list(bpy.data.objects))
# print(bpy.data.objects)
# for element in bpy.data.objects:
#     bpy.data.objects.remove(element)
# bpy.data.meshes.remove(bpy.data.meshes["Cube"])
# print(bpy.data.objects)

# BLENDERGIS
# io_import_shp.register()
# bpy.ops.importgis.importgis.shapefile_props_dialog("EXEC_DEFAULT")
file = "./shp/fiumicino-line.shp"
# bpy.ops.importgis.shapefile("EXEC_DEFAULT", filepath=file)
# print("matteo edges {}".format(len(list(bpy.data.objects['export-line'].data.vertices))))
# print(type(bpy.data.objects['export-line'].data.edges[0]))
# print(bpy.data.meshes['export-line'].vertices[0].co.xy.x)

# PYSHP
sf = shapefile.Reader(file)
# shapes = sf.shapes()
# records = sf.records()
# for record, shape in zip(records, shapes):
#     print("record {} n of points {} type {}".format(record, len(shape.points), shape.shapeType))

# GEOPANDAS
data_frame = geopandas.read_file(file)

# FILTRO SU BBOX (punto in alto a sinistra e punto in basso a destra)

bbox = (56.378, 25.17114, 56.380, 25.4035)
data_frame_bbox = geopandas.read_file(file, bbox=bbox)

# TEST
points_sf = sum([len(shape.points) for shape in sf.shapes()])
# print(points_sf)
points_df = sum([len(part.coords) for part in data_frame.geometry])
print(points_df)


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


# TEST
matrix = [list(part.coords) for part in data_frame.geometry]
matrix = merge_parts(matrix)
print(data_frame)

# Creating dataframe
linestrings = [LineString(line) for line in matrix]
new_data_frame = geopandas.GeoDataFrame({"geometry": linestrings})
print(new_data_frame)
try:
    new_data_frame.to_file("./shp/test.shp")
except ValueError:
    print("Check your bbox, it seems the coastline is not there.")

# Creation of dominion
print(list(new_data_frame.geometry[0].coords)[0])
print(list(new_data_frame.geometry[0].coords)[-1])
sea_points = [(11.992, 41.872), (12.080, 41.680)]

# ok for 1,2, and 3 point of dominion bbox (no model over islands, we have a coastline)
index_of_coastline = [idx for idx, part in enumerate(new_data_frame.geometry) if not part.is_ring][0]
matrix = [list(part.coords) for part in new_data_frame.geometry]

coastline = matrix[index_of_coastline]
start_coast = coastline[0]
end_coast = coastline[-1]

distances_from_start = [dist(start_coast, sea_point).km for sea_point in sea_points]
coastline.insert(0, sea_points[distances_from_start.index(min(distances_from_start))])

distances_from_end = [dist(end_coast, sea_point).km for sea_point in sea_points]
copy_of_distances = distances_from_end.copy()
while copy_of_distances:
    coastline.append(sea_points[distances_from_end.index(min(copy_of_distances))])
    copy_of_distances.remove(min(copy_of_distances))

matrix[index_of_coastline] = coastline

linestrings = [LineString(line) for line in matrix]
dominion_data_frame = geopandas.GeoDataFrame({"geometry": linestrings})
print(dominion_data_frame)
for part in dominion_data_frame.geometry:
    print(part.is_ring)
    print(part.coords[0])
    print(part.coords[-1])
try:
    dominion_data_frame.to_file("./shp/test_dominion.shp")
except ValueError:
    print("Check your bbox, it seems the coastline is not there.")