import geopandas
from shapely.geometry import LineString, Polygon
from geopy.distance import geodesic as dist
from exceptions.exceptions import SmoothingError
from matplotlib import pyplot as plt


class ShpHandler:

    def __init__(self, p_shapefile, p_bbox):
        self.mesh_type = None
        self.bbox = p_bbox
        if self.bbox:
            self.bbox = [float(coord) for coord in self.bbox]
            if len(self.bbox) == 4:
                self.data_frame = geopandas.read_file(p_shapefile, bbox=self.bbox)
            else:
                self.data_frame = geopandas.read_file(p_shapefile)
                self.__polygons_to_linestring__()
        else:
            self.data_frame = geopandas.read_file(p_shapefile)
        self.matrix = []

    def prepare_to_close(self):
        if isinstance(self.data_frame.geometry[0], Polygon):
            self.__polygons_to_linestring__()

    def __polygons_to_linestring__(self):
        # TODO: fix issue with intersected different from original polygon (boundary slightly different).
        #  For the moment hot fix take lat of first point of coast
        # Create a polygon from the bounding box
        if len(self.bbox) == 4:
            lons_poly = [self.bbox[0], self.bbox[2], self.bbox[2], self.bbox[0]]
            lats_poly = [self.bbox[1], self.bbox[1], self.bbox[3], self.bbox[3]]
        else:
            lons_poly = [coord for index, coord in enumerate(self.bbox) if index % 2 == 0]
            lats_poly = [coord for index, coord in enumerate(self.bbox) if index % 2 == 1]
        polygon_geom = Polygon(zip(lons_poly, lats_poly))
        polygon = geopandas.GeoDataFrame(index=[0], geometry=[polygon_geom])
        # Intersect bounding box polygon with dataframe polygon and extract the boundary
        intersected_dataframe = geopandas.overlay(self.data_frame, polygon).boundary
        polygon_boundary = polygon.boundary
        polygon_boundary = geopandas.GeoDataFrame(index=[0], geometry=polygon_boundary.geometry)
        # The difference between intersected polygon boundary and bounding box polygon boundary gives us the coastline
        intersected_dataframe = geopandas.GeoDataFrame(geometry=intersected_dataframe.geometry)
        self.data_frame = geopandas.overlay(intersected_dataframe, polygon_boundary, "difference")
        return

    def close_dominion(self, p_threshold, p_sea_points, p_smoothing_degree):
        # Merging parts of shapefile
        for part in self.data_frame.geometry:
            try:
                for line in part:
                    self.matrix.append(list(line.coords))
            except TypeError:
                self.matrix.append(list(part.coords))
        self.__merge_parts__(None, float(p_threshold))

        # Creating dataframe
        linestrings = [LineString(line) for line in self.matrix]
        new_data_frame = geopandas.GeoDataFrame({"geometry": linestrings})

        # Creation of dominion
        sea_points = []

        for string in p_sea_points:
            a, b = string.split(",")
            coords = (float(a), float(b))
            sea_points.append(coords)

        # ok for 1,2, and 3 point of dominion bbox (no model over islands, we have a coastline)
        index_of_coastlines = None
        try:
            index_of_coastlines = [idx for idx, part in enumerate(new_data_frame.geometry) if not part.is_ring]
            # This assignment is needed to trigger the exception
            trigger = index_of_coastlines[0]
            if len(index_of_coastlines) == 1:
                self.mesh_type = 1 if sea_points else 3
            else:
                self.mesh_type = 4
        except IndexError:
            self.mesh_type = 2
        self.matrix = [list(part.coords) for part in new_data_frame.geometry]
        coastline = None

        if self.mesh_type == 1:
            del self.matrix[index_of_coastlines[0]][-1]
            coastline = self.matrix[index_of_coastlines[0]]
            self.__orient_coastline__(coastline, sea_points)
        elif self.mesh_type == 3:
            del self.matrix[index_of_coastlines[0]][-1]
            coastline = self.matrix[index_of_coastlines[0]]
        elif self.mesh_type == 4:
            for i in index_of_coastlines:
                del self.matrix[index_of_coastlines[i]][-1]
                del self.matrix[index_of_coastlines[i]][0]
            coastline = [self.matrix[index_of_coastlines[i]] for i in index_of_coastlines]
        dominion_closure = self.__create_closure__(coastline, sea_points, p_smoothing_degree)
        if self.mesh_type == 1:
            coastline[-1:-1] = dominion_closure
            coastline[-1] = coastline[0]
        elif self.mesh_type == 2:
            self.matrix.append(dominion_closure)
        elif self.mesh_type == 3:
            coastline.append(dominion_closure)
        elif self.mesh_type == 4:
            self.matrix.append(dominion_closure)
        for index, line in enumerate(self.matrix):
            self.matrix[index] = self.remove_doubles(self.matrix[index])
        linestrings = [LineString(line) for line in self.matrix]
        types = ["island"] * len(self.matrix)
        if self.mesh_type == 1 or self.mesh_type == 3:
            types[index_of_coastlines[0]] = "coastline"
        elif self.mesh_type == 2:
            types[-1] = "dominion_closure"
        elif self.mesh_type == 4:
            types[-1] = "dominion_closure"
            for i in index_of_coastlines:
                types[index_of_coastlines[i]] = "coastline"
        return geopandas.GeoDataFrame({"geometry": linestrings, "types": types}), len(dominion_closure)

    def __merge_parts__(self, p_matrix=None, p_threshold=0.001):
        matrix = self.matrix if p_matrix is None else p_matrix
        len_matrix = len(matrix)
        index = 0
        while index < len(matrix):
            other_index = index + 1
            while other_index < len(matrix):
                part = matrix[index]
                other_part = matrix[other_index]

                distances = [
                    dist(part[0], other_part[0]).km,
                    dist(part[0], other_part[-1]).km,
                    dist(part[-1], other_part[0]).km,
                    dist(part[-1], other_part[-1]).km
                ]

                if min(distances) < p_threshold:
                    matrix[index] = self.__arrange_parts__(part, other_part, distances.index(min(distances)))
                else:
                    other_index += 1
                    continue

                try:
                    del matrix[other_index]
                except IndexError:
                    pass
            index += 1
        if len_matrix == len(matrix):
            return matrix
        self.__merge_parts__(matrix, p_threshold)
        return matrix

    @staticmethod
    def __arrange_parts__(p_part_one, p_part_two, p_situation):
        results = [
            p_part_two[::-1] + p_part_one,
            p_part_two + p_part_one,
            p_part_one + p_part_two,
            p_part_one + p_part_two[::-1]
        ]
        return results[p_situation]

    @staticmethod
    def __orient_coastline__(p_coastline, p_sea_points):
        mean_coastline = [
            (p_coastline[0][0] + p_coastline[-1][0]) / 2,
            (p_coastline[0][1] + p_coastline[-1][1]) / 2
        ]
        try:
            mean_sea_points = [
                sum(point[0] for point in p_sea_points) / len(p_sea_points),
                sum(point[1] for point in p_sea_points) / len(p_sea_points)
            ]
        except TypeError:
            mean_sea_points = p_sea_points

        lon_diff = abs(p_coastline[0][0] - p_coastline[-1][0])
        lat_diff = abs(p_coastline[0][1] - p_coastline[-1][1])
        if lon_diff <= lat_diff:
            if mean_sea_points[0] > mean_coastline[0]:
                # Sea is east with respect to the coastline
                if p_coastline[0][1] > p_coastline[-1][1]:
                    p_coastline.reverse()
            else:
                # Sea is west with respect to the coastline
                if p_coastline[0][1] < p_coastline[-1][1]:
                    p_coastline.reverse()
        else:
            if p_coastline[0][0] > p_coastline[-1][0]:
                p_coastline.reverse()
        return

    @staticmethod
    def __order_by_distance_(p_coastline):
        coast = p_coastline[0]
        min_distances = []
        for other_coast in p_coastline:
            min_distances.append(min([
                dist(coast[0], other_coast[0]).km,
                dist(coast[0], other_coast[-1]).km,
                dist(coast[-1], other_coast[0]).km,
                dist(coast[-1], other_coast[-1]).km
            ]))
        zipped = zip(min_distances, p_coastline)
        sorted_zip = sorted(zipped)
        p_coastline = [coast for _, coast in sorted_zip]
        return p_coastline

    def __create_closure__(self, p_coastline, p_sea_points, p_smoothing_degree):
        dominion_closure = []
        if self.mesh_type == 1:
            start_coast = p_coastline[0]
            end_coast = p_coastline[-1]

            dominion_closure = [end_coast]
            distances_from_end = [dist(end_coast, sea_point).km for sea_point in p_sea_points]
            copy_of_distances = distances_from_end.copy()
            while copy_of_distances:
                dominion_closure.append(p_sea_points[distances_from_end.index(min(copy_of_distances))])
                copy_of_distances.remove(min(copy_of_distances))
            dominion_closure.append(start_coast)
            dominion_closure = self.__smooth_corners__(LineString(dominion_closure), p_smoothing_degree)
        elif self.mesh_type == 2:
            p_sea_points.append(p_sea_points[0])
            dominion_closure = p_sea_points
            dominion_closure = self.__smooth_corners__(LineString(dominion_closure), p_smoothing_degree)
        elif self.mesh_type == 3:
            start_coast = p_coastline[0]
            dominion_closure = start_coast
        elif self.mesh_type == 4:
            p_coastline = self.__order_by_distance_(p_coastline)
            dominion_closure = self.__merge_parts__(p_coastline, 10000)[0]
            dominion_closure.append(dominion_closure[0])
        if not dominion_closure:
            raise SmoothingError
        return dominion_closure

    def __smooth_corners__(self, p_line, p_degree):
        if self.mesh_type == 1:
            offset = p_line.parallel_offset(p_degree, "right", join_style=1)

            offset2 = offset.parallel_offset(p_degree, "right", join_style=1)
            return list(offset2.coords)[1:-1]
        elif self.mesh_type == 2:
            closure = p_line.buffer(p_degree, join_style=1, mitre_limit=1).exterior
            return list(closure.coords)
        else:
            return []

    @staticmethod
    def remove_doubles(p_list, p_threshold=0.001):
        counter = 0
        while counter < len(p_list) - 1:
            if dist(p_list[counter], p_list[counter + 1]).km < p_threshold:
                del p_list[counter + 1]
            counter += 1
        return p_list
