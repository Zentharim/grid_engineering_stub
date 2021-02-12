from geopy.distance import geodesic as dist
from numpy import argmin


class GeoWriter:

    def __init__(self, p_dataframe, p_mesh_params, file="./msh/temp.geo_unrolled"):
        self.fp = open(file, "w")
        self.file = file
        self.coastline_index = list(p_dataframe.types).index("coastline")
        self.counter = 0
        self.point_counter = 1
        self.start_point = 1
        self.end_point = None
        self.line_counter = 1
        self.line_loop_counter = 1
        self.curve_loops = []
        self.points = []
        self.mesh_params = p_mesh_params
        self.dataframe = p_dataframe
        self.first_attractor = None
        self.last_attractor = None
        self.end_coastline = None
        self.start_coastline = None
        self.other_attractors = []
        self.imposed_attractors = []

    def __points_to_geo__(self, list_of_coords):
        for point in list_of_coords:
            self.fp.write("Point({}) = {{{}, {}, 0, 1.0}};\n".format(self.point_counter, point[0], point[1]))
            self.points.append(point)
            self.point_counter += 1
        self.end_point = self.point_counter - 1
        return

    def __lines_to_geo__(self):
        shape_lines = []
        self.fp.write("Line({}) = {{{}:{}}};\n".format(self.line_counter, self.start_point, self.end_point))
        shape_lines.append(self.line_counter)
        self.line_counter += 1
        self.fp.write("Line({}) = {{{},{}}};\n".format(self.line_counter, self.end_point, self.start_point))
        shape_lines.append(self.line_counter)
        self.line_counter += 1
        return shape_lines

    def __line_loops_to_geo__(self, p_shape_lines):
        loop_buffer = "Line Loop({}) = {{".format(self.line_loop_counter)
        self.curve_loops.append(self.line_loop_counter)
        self.line_loop_counter += 1
        for element in p_shape_lines:
            loop_buffer += "{}, ".format(element)
        loop_buffer = loop_buffer[:-2] + "};\n"
        self.fp.write(loop_buffer)
        return

    def __planes_to_geo__(self):
        plane_buffer = "Plane Surface(1) = {"
        for element in self.curve_loops:
            plane_buffer += "{}, ".format(element)
        plane_buffer = plane_buffer[:-2] + "};\n"
        self.fp.write(plane_buffer)
        return

    def __dicotomic_search_attractors__(self, interval, distance, point):
        distances = [abs(dist(point, interval_point).km - distance) for interval_point in interval]
        if min(distances) <= 0.001:
            return interval[argmin(distances)]
        mid_point = [(interval[0][0]+interval[1][0])/2, (interval[0][1]+interval[1][1])/2]
        if argmin(distances) == 0:
            return self.__dicotomic_search_attractors__([interval[0], mid_point], distance, point)
        else:
            return self.__dicotomic_search_attractors__([mid_point, interval[1]], distance, point)

    def __find_attractors__(self):
        attractors_min_distance_km = self.mesh_params["LcMax"] * 100 * self.mesh_params["distance_coeff"]
        coastline_points = self.points[self.start_coastline - 1:self.end_coastline - self.mesh_params["extra_points"]]
        # If there are islands add islands as attractors (if required)
        if self.mesh_params["islands_attractors"] and \
                not (len(self.points) == self.end_coastline and self.start_coastline == 1):
            if self.start_coastline == 1:
                self.other_attractors = [(self.end_coastline + 1, len(self.points))]
            elif self.end_coastline == len(self.points):
                self.other_attractors = [(1, self.start_coastline - 1)]
            else:
                self.other_attractors = [(1, self.start_coastline - 1), (self.end_coastline + 1, len(self.points))]
        for index, point in enumerate(coastline_points):
            if self.first_attractor is None and dist(coastline_points[0], point).km >= attractors_min_distance_km:
                interval_first = [coastline_points[index-1], coastline_points[index]]
                self.first_attractor = self.start_coastline + index
            if self.first_attractor is not None and dist(coastline_points[-1], point).km <= attractors_min_distance_km:
                interval_second = [coastline_points[index-1], coastline_points[index]]
                self.last_attractor = self.start_coastline + index - 1
                break

        # Create and impose attractors at "safety distance" from boundaries ON THE COAST.
        # To be revised when working on multiple coastlines
        try:
            first_point = self.__dicotomic_search_attractors__(interval_first, attractors_min_distance_km, coastline_points[0])
            second_point = self.__dicotomic_search_attractors__(interval_second, attractors_min_distance_km, coastline_points[-1])
            self.__points_to_geo__([first_point, second_point])
            self.imposed_attractors = [self.point_counter - 1, self.point_counter - 2]
        except UnboundLocalError as e:
            print("Your d_bound is off, I cannot find any legal attractor. "
                  "Please check c_bound and delta_open parameters")
            return 1
        return 0

    def __attractors_to_geo__(self):
        attractors_buffer = "Field[1].NodesList = {"
        if self.first_attractor is not None and self.last_attractor is not None:
            attractors_buffer += "{}:{}, ".format(self.first_attractor, self.last_attractor)
        for couple in self.other_attractors:
            attractors_buffer += "{}:{}, ".format(couple[0], couple[1])
        for point in self.imposed_attractors:
            attractors_buffer += "{}, ".format(point)
        attractors_buffer = attractors_buffer[:-2] + "};\n"
        return attractors_buffer

    def __mesh_params_to_geo__(self):
        self.fp.write("Field[1] = Attractor;\n")
        self.fp.write("forceParametrizablePatches = 1;\n")
        self.fp.write(self.__attractors_to_geo__())
        self.fp.write("Field[2] = Threshold;\n")
        self.fp.write("Field[2].IField = 1;\n")
        self.fp.write("Field[2].DistMax = {};\n".format(self.mesh_params["DistMax"]))
        self.fp.write("Field[2].DistMin = {};\n".format(self.mesh_params["DistMin"]))
        self.fp.write("Field[2].LcMax = {};\n".format(self.mesh_params["LcMax"]))
        self.fp.write("Field[2].LcMin = {};\n".format(self.mesh_params["LcMin"]))
        self.fp.write("Background Field = 2;\n")
        self.fp.write("Mesh.Algorithm = 6;\n")
        self.fp.write("Mesh.CharacteristicLengthExtendFromBoundary = 0;\n")
        return

    def shapefile_to_geo(self):
        for shape in self.dataframe.geometry:
            self.__points_to_geo__(list(shape.coords)[:-1])

            shape_lines = self.__lines_to_geo__()

            self.__line_loops_to_geo__(shape_lines)
            if self.counter == self.coastline_index:
                self.end_coastline = self.end_point
                self.start_coastline = self.start_point
            self.start_point = self.point_counter
            self.counter += 1

        self.__planes_to_geo__()

        result_code = self.__find_attractors__()
        if result_code:
            return result_code
        self.__mesh_params_to_geo__()
        self.fp.close()
        return result_code

    def clean(self, p_mesh_params):
        self.fp = open(self.file, 'w')
        self.counter = 0
        self.point_counter = 1
        self.start_point = 1
        self.end_point = None
        self.line_counter = 1
        self.line_loop_counter = 1
        self.curve_loops = []
        self.points = []
        self.mesh_params = p_mesh_params
        self.first_attractor = None
        self.last_attractor = None
        self.end_coastline = None
        self.start_coastline = None
        self.other_attractors = []
        self.imposed_attractors = []
