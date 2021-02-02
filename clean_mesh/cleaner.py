import numpy as np
import math


class Cleaner:

    def __init__(self, vertices, triangles, quads):
        self.vertices = vertices
        self.triangles = triangles
        self.quads = quads

    def __clean_quads__(self):
        for quad in self.quads:
            self.triangles.append(quad[0, 2])
            self.triangles.append(quad[2:])
        self.quads = []
        return

    def __clean_alone_points__(self):
        used_vertices = set(vertex for vertex_group in self.triangles for vertex in vertex_group)
        for index, vertex in enumerate(self.vertices):
            if index+1 not in used_vertices:
                self.vertices.remove(vertex)
        # [self.vertices.remove(vertex) for vertex in self.vertices if vertex not in used_vertices]
        return

    @staticmethod
    def __normalize__(x):
        return x / np.linalg.norm(x)

    # Something is seriously wrong
    def __find_bad_angles__(self, max_angle=160):
        threshold = math.cos(max_angle * math.pi / 180)
        return [triangle for triangle in self.triangles if
                self.__normalize__(self.vertices[triangle[1] - 1] - self.vertices[triangle[0] - 1]) *
                self.__normalize__(self.vertices[triangle[2] - 1] - self.vertices[triangle[0] - 1]) < threshold or
                self.__normalize__(self.vertices[triangle[2] - 1] - self.vertices[triangle[1] - 1]) *
                self.__normalize__(self.vertices[triangle[0] - 1] - self.vertices[triangle[1] - 1]) < threshold or
                self.__normalize__(self.vertices[triangle[0] - 1] - self.vertices[triangle[2] - 1]) *
                self.__normalize__(self.vertices[triangle[1] - 1] - self.vertices[triangle[2] - 1]) < threshold]
    
    def clean(self):
        self.__clean_quads__()
        self.__clean_alone_points__()
        return

    def get_data(self):
        return self.vertices, self.triangles
