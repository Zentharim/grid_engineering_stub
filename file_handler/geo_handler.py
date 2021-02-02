from geopy.distance import geodesic as dist


class GeoWriter:

    def __init__(self, p_dataframe, p_mesh_params, file="./msh/temp.geo_unrolled"):
        self.fp = open(file, "w")
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

    def __points_to_geo__(self, p_shape):
        for point in list(p_shape.coords)[:-1]:
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

    def __find_attractors__(self):
        if not (len(self.points) == self.end_coastline and self.start_coastline == 1):
            if self.start_coastline == 1:
                self.other_attractors = [(self.end_coastline + 1, len(self.points))]
            elif self.end_coastline == len(self.points):
                self.other_attractors = [(1, self.start_coastline - 1)]
            else:
                self.other_attractors = [(1, self.start_coastline - 1), (self.end_coastline + 1, len(self.points))]
        for index, point in enumerate(self.points[self.start_coastline - 1:self.end_coastline - self.mesh_params["extra_points"]]):
            if self.first_attractor is None and \
                    dist(self.points[self.start_coastline - 1], point).km >= self.mesh_params["LcMax"] * 100 * \
                    self.mesh_params["distance_coeff"]:
                self.first_attractor = self.start_coastline - 1 + index
            if self.first_attractor is not None and self.last_attractor is None and \
                    dist(self.points[self.end_coastline - self.mesh_params["extra_points"] - 2], point).km <= \
                    self.mesh_params["LcMax"] * 100 * self.mesh_params["distance_coeff"]:
                self.last_attractor = self.start_coastline - 1 + index
                return

    def __attractors_to_geo__(self):
        attractors_buffer = "Field[1].NodesList = {"
        attractors_buffer += "{}:{}, ".format(self.first_attractor, self.last_attractor)
        for couple in self.other_attractors:
            attractors_buffer += "{}:{}, ".format(couple[0], couple[1])
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
        self.fp.write("Mesh.Algorithm = 5;\n")
        self.fp.write("Mesh.CharacteristicLengthExtendFromBoundary = 0;\n")
        return

    def shapefile_to_geo(self):
        for shape in self.dataframe.geometry:
            self.__points_to_geo__(shape)

            shape_lines = self.__lines_to_geo__()

            self.__line_loops_to_geo__(shape_lines)
            if self.counter == self.coastline_index:
                self.end_coastline = self.end_point
                self.start_coastline = self.start_point
            self.start_point = self.point_counter
            self.counter += 1

        self.__planes_to_geo__()

        self.__find_attractors__()
        self.__mesh_params_to_geo__()
        self.fp.close()
        return
