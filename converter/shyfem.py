import os
import numpy as np
from scipy.spatial.qhull import Delaunay


# LOAD SHYFEM GRD AS CLASS

class ShyfemGrid:
    def __init__(self, grd_path):
        self.path = grd_path
        self.get_nodes_elements()
        self.grid = ""
        self.nodes = []
        self.elements = []
        self.triangs_array = []
        self.triangs = []
        self.centroids = []
        self.msk = None
        self.nodesInbox = []
        self.elemsInbox = []
        self.centroidsInbox = []
        self.latitude = 0
        self.longitude = 0
        self.depth = 0
        self.tri = []

    @staticmethod
    def until_emtpy_line(data):
        while True:
            line = next(data)
            if line == "" or (line[:-1].strip() == ""):
                return
            yield line

    def get_nodes_elements(self):
        self.grid = open(self.path)
        if os.path.exists(self.path + ".npz"):
            ds = np.load(self.path + ".npz")
            self.nodes, self.elements = ds["nodes"], ds["elements"]
        else:
            self.nodes = np.loadtxt(self.until_emtpy_line(self.grid))  # load first grid block
            self.elements = np.loadtxt(self.until_emtpy_line(self.grid))
            np.savez(self.path, nodes=self.nodes, elements=self.elements)

    def in_box(self):
        id_nodes = self.nodes[:, 1].astype(int)
        verts_ids = self.elements[:, 4:7].astype(int)
        node_sorter = np.argsort(id_nodes)
        sorted_idx = node_sorter[
            id_nodes.searchsorted(verts_ids, sorter=node_sorter)]  # associate original position to the named index

        self.triangs_array = self.nodes[sorted_idx]
        self.triangs = sorted_idx

        self.centroids = np.average(self.triangs_array[:, :, 3:], axis=1)

        self.msk = (
                (self.centroids[:, 0] >= -180) &
                (self.centroids[:, 0] <= 180) &
                (self.centroids[:, 1] >= -90) &
                (self.centroids[:, 1] <= 90))

        self.nodesInbox = self.triangs_array[self.msk]
        self.elemsInbox = self.elements[self.msk]
        self.centroidsInbox = self.centroids[self.msk]
        self.latitude = self.centroidsInbox[:, 1]
        self.longitude = self.centroidsInbox[:, 0]
        try:
            self.depth = self.elemsInbox[:, 7]
        except IndexError:
            self.depth = self.elemsInbox[:, -1] * 0

    def triangulate_elems(self):
        self.in_box()
        tri = Delaunay(np.vstack((self.longitude, self.latitude)).T)
        self.tri = tri.simplices.copy()

    def tri_from_grd_to_numpy(self):
        self.in_box()
        lon = self.nodes[:, -2]
        lat = self.nodes[:, -1]
        tri = self.triangs
        d = np.zeros_like(lat)
        return np.array((lon, lat, d)).T, tri

    def box_to_numpy(self, outfile):
        self.triangulate_elems()
        np.savez(outfile, v=np.array((self.longitude, self.latitude, self.depth)).T, t=self.tri)
        return np.array((self.longitude, self.latitude, self.depth)).T, self.tri