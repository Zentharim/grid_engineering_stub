import os
import numpy as np
from converter.shyfem import ShyfemGrid
import bpy
import xarray as xr


# REQUIRES SHYFEM GRD CLASS
class BlenderGrid:
    def __init__(self, name, grd_path):
        self.name = name
        self.outPath = outPath = os.path.dirname(grd_path)
        print(outPath)

        os.makedirs(outPath, exist_ok=True)

        self.grid = ShyfemGrid(grd_path)
        v, t = self.grid.tri_from_grd_to_numpy()
        self.ob = self.load_obj(self.name, v, t)

    def write_grd(self, outname):
        bpy.ops.object.mode_set(mode='OBJECT')
        m = bpy.context.object.data
        mesh = np.array([i.co for i in m.vertices])
        x = mesh[:, 0]
        y = mesh[:, 1]
        ran = [i + 1 for i, j in enumerate(x)]
        em = np.zeros_like(y)
        nodes = np.array((em + 1, ran, em, x, y)).T
        tri = np.array([list(np.array(t.vertices) + 1) for t in m.polygons])
        base = np.zeros(tri.shape[0])
        ranE = [i + 1 for i, j in enumerate(base)]
        elems = np.array((base + 2, ranE, base, base + 3, tri[:, 0], tri[:, 1], tri[:, 2])).T

        with open('%s.grd' % outname, 'w')as out:
            np.savetxt(out, nodes, fmt='%i %i %i %f %f')  ### grid.nodes /grid.elems
            np.savetxt(out, elems, fmt='%i %i %i %i %i %i %i')

    def mk_obj(self, name):
        me = bpy.data.meshes.new(name)
        ob = bpy.data.objects.new(name, me)
        scn = bpy.context.scene
        scn.objects.link(ob)
        scn.objects.active = ob
        ob.select = True
        return ob

    def load_obj(self, name, v, t):
        ob = self.mk_obj(name)
        ob.data.from_pydata(v.tolist(), [], t.tolist())
        return ob

    def selection_to_mask(self):
        obj = bpy.context.object
        mesh = obj.data
        obj.update_from_editmode()  # Loads edit-mode data into object data
        msk = [1 if e.select else 0 for e in mesh.vertices]
        # np.save(os.path.join(self.outPath, '{}_mask'.format(self.name)), msk)
        xr.Dataset({'mask': (['node'], msk)},
                   coords={'node': np.arange(len(msk))}).to_netcdf(
            os.path.join(self.outPath, '{}_mask.nc'.format(self.name)))


'''

from blenderize.blenderizator import  BlenderGrid,sortEdges
import numpy as np
grd='/Users/scausio/Desktop/SANIv2_grid/umedbs/atl_Umed_hmin5.grd'
box=[[25.6,30],[39.5,41.6]]
grid=BlenderGrid('transmarmara',grd,box)



f='/Users/scausio/Desktop/SANIv2_grid/umedbs/trans_msk.npy'
#f='/Users/scausio/Desktop/SANIv2_grid/umedbs/transMask_ok.npy'

selection=np.load(f)

bpy.ops.object.select_all(action='DESELECT')
bpy.ops.object.mode_set(mode = 'OBJECT')
mesh = bpy.context.object.data

for e , sel in zip (mesh.edges,selection):
    e.select=sel
bpy.ops.object.mode_set(mode = 'EDIT')

coordsTransect = np.array([i.co.xy for i in mesh.vertices if i.select])

selected_edges = [e for e in mesh.edges if e.select]
selected_nodesId = [tuple((i.vertices[0], i.vertices[1])) for i in selected_edges]
#selected_nodesCo = np.array([tuple((i.vertices[0], i.vertices[1])) for i in selected_edges])
s=[mesh.vertices[i].co.xy for i in selected_nodes]


scurzune= sortEdges(selected_nodes)
s=[mesh.vertices[i].co.xy for i in selected_nodes]

np.save('/Users/scausio/Desktop/SANIv2_grid/umedbs/edges',selected_nodes)




        '''