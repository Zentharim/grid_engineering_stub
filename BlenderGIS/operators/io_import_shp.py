# KEEP

# -*- coding:utf-8 -*-
import os, sys, time
import bpy
from bpy.props import StringProperty, BoolProperty, EnumProperty, IntProperty
from bpy.types import Operator
import bmesh
import math
from mathutils import Vector
from ..core.lib.shapefile import Reader as shpReader

from ..geoscene import GeoScene
from ..core import BBOX
from ..core.proj import Reproj
from ..core.utils import perf_clock

from .utils import adjust3Dview, getBBOX, DropToGround

import logging
log = logging.getLogger(__name__)

PKG, SUBPKG = __package__.split('.', maxsplit=1)

featureType={
0:'Null',
1:'Point',
3:'PolyLine',
5:'Polygon',
8:'MultiPoint',
11:'PointZ',
13:'PolyLineZ',
15:'PolygonZ',
18:'MultiPointZ',
21:'PointM',
23:'PolyLineM',
25:'PolygonM',
28:'MultiPointM',
31:'MultiPatch'
}


"""
dbf fields type:
	C is ASCII characters
	N is a double precision integer limited to around 18 characters in length
	D is for dates in the YYYYMMDD format, with no spaces or hyphens between the sections
	F is for floating point numbers with the same length limits as N
	L is for logical data which is stored in the shapefile's attribute table as a short integer as a 1 (true) or a 0 (false).
	The values it can receive are 1, 0, y, n, Y, N, T, F or the python builtins True and False
"""


# KEEP
class IMPORTGIS_OT_shapefile(Operator):
	"""Import from ESRI shapefile file format (.shp)"""

	bl_idname = "importgis.shapefile" # important since its how bpy.ops.import.shapefile is constructed (allows calling operator from python console or another script)
	#bl_idname rules: must contain one '.' (dot) charactere, no capital letters, no reserved words (like 'import')
	bl_description = 'Import ESRI shapefile (.shp)'
	bl_label = "Import SHP"
	bl_options = {"UNDO"}

	filepath: StringProperty()

	shpCRS: StringProperty(name = "Shapefile CRS", description = "Coordinate Reference System")

	elevSource: StringProperty(name = "Elevation source", description = "Elevation source", default='GEOM') # [NONE, GEOM, OBJ, FIELD]
	objElevName: StringProperty(name = "Elevation object name", description = "")

	fieldElevName: StringProperty(name = "Elevation field", description = "Field name")
	fieldExtrudeName: StringProperty(name = "Extrusion field", description = "Field name")
	fieldObjName: StringProperty(name = "Objects names field", description = "Field name")

	#Extrusion axis
	extrusionAxis: EnumProperty(
			name="Extrude along",
			description="Select extrusion axis",
			items=[ ('Z', 'z axis', "Extrude along Z axis"),
			('NORMAL', 'Normal', "Extrude along normal")]
			)
	#Create separate objects
	separateObjects: BoolProperty(
			name="Separate objects",
			description="Import to separate objects instead one large object",
			default=False
			)

	@classmethod
	def poll(cls, context):
		return context.mode == 'OBJECT'

	def __del__(self):
		bpy.context.window.cursor_set('DEFAULT')

	def execute(self, context):

		# prefs = bpy.context.preferences.addons[PKG].preferences

		#Set cursor representation to 'loading' icon
		w = context.window
		w.cursor_set('WAIT')
		t0 = perf_clock()

		bpy.ops.object.select_all(action='DESELECT')

		#Path
		shpName = os.path.basename(self.filepath)[:-4]

		#Get shp reader
		log.info("Read shapefile...")
		try:
			shp = shpReader(self.filepath)
		except Exception as e:
			log.error("Unable to read shapefile", exc_info=True)
			self.report({'ERROR'}, "Unable to read shapefile, check logs")
			return {'CANCELLED'}

		#Check shape type
		shpType = featureType[shp.shapeType]
		log.info('Feature type : ' + shpType)
		if shpType not in ['Point','PolyLine','Polygon','PointZ','PolyLineZ','PolygonZ']:
			self.report({'ERROR'}, "Cannot process multipoint, multipointZ, pointM, polylineM, polygonM and multipatch feature type")
			return {'CANCELLED'}

		if self.elevSource != 'FIELD':
			self.fieldElevName = ''

		if self.elevSource == 'OBJ':
			scn = bpy.context.scene
			elevObj = scn.objects[self.objElevName]
			rayCaster = DropToGround(scn, elevObj)

		#Get fields
		fields = [field for field in shp.fields if field[0] != 'DeletionFlag'] #ignore default DeletionFlag field
		fieldsNames = [field[0] for field in fields]
		log.debug("DBF fields : "+str(fieldsNames))

		if self.separateObjects or self.fieldElevName or self.fieldObjName or self.fieldExtrudeName:
			self.useDbf = True
		else:
			self.useDbf = False

		if self.fieldObjName and self.separateObjects:
			try:
				nameFieldIdx = fieldsNames.index(self.fieldObjName)
			except Exception as e:
				log.error('Unable to find name field', exc_info=True)
				self.report({'ERROR'}, "Unable to find name field")
				return {'CANCELLED'}

		if self.fieldElevName:
			try:
				zFieldIdx = fieldsNames.index(self.fieldElevName)
			except Exception as e:
				log.error('Unable to find elevation field', exc_info=True)
				self.report({'ERROR'}, "Unable to find elevation field")
				return {'CANCELLED'}

			if fields[zFieldIdx][1] not in ['N', 'F', 'L'] :
				self.report({'ERROR'}, "Elevation field do not contains numeric values")
				return {'CANCELLED'}

		if self.fieldExtrudeName:
			try:
				extrudeFieldIdx = fieldsNames.index(self.fieldExtrudeName)
			except ValueError:
				log.error('Unable to find extrusion field', exc_info=True)
				self.report({'ERROR'}, "Unable to find extrusion field")
				return {'CANCELLED'}

			if fields[extrudeFieldIdx][1] not in ['N', 'F', 'L'] :
				self.report({'ERROR'}, "Extrusion field do not contains numeric values")
				return {'CANCELLED'}

		#Get shp and scene georef infos
		shpCRS = self.shpCRS
		geoscn = GeoScene()
		if geoscn.isBroken:
			self.report({'ERROR'}, "Scene georef is broken, please fix it beforehand")
			return {'CANCELLED'}

		scale = geoscn.scale #TODO

		if not geoscn.hasCRS: #if not geoscn.isGeoref:
			try:
				geoscn.crs = shpCRS
			except Exception as e:
				log.error("Cannot set scene crs", exc_info=True)
				self.report({'ERROR'}, "Cannot set scene crs, check logs for more infos")
				return {'CANCELLED'}

		#Init reprojector class
		if geoscn.crs != shpCRS:
			log.info("Data will be reprojected from {} to {}".format(shpCRS, geoscn.crs))
			try:
				rprj = Reproj(shpCRS, geoscn.crs)
			except Exception as e:
				log.error('Reprojection fails', exc_info=True)
				self.report({'ERROR'}, "Unable to reproject data, check logs for more infos.")
				return {'CANCELLED'}
			if rprj.iproj == 'EPSGIO':
				if shp.numRecords > 100:
					self.report({'ERROR'}, "Reprojection through online epsg.io engine is limited to 100 features. \nPlease install GDAL or pyproj module.")
					return {'CANCELLED'}

		#Get bbox
		bbox = BBOX(shp.bbox)
		if geoscn.crs != shpCRS:
			bbox = rprj.bbox(bbox)

		#Get or set georef dx, dy
		if not geoscn.isGeoref:
			dx, dy = bbox.center
			geoscn.setOriginPrj(dx, dy)
		else:
			dx, dy = geoscn.getOriginPrj()

		#Get reader iterator (using iterator avoids loading all data in memory)
		#warn, shp with zero field will return an empty shapeRecords() iterator
		#to prevent this issue, iter only on shapes if there is no field required
		if self.useDbf:
			#Note: using shapeRecord solve the issue where number of shapes does not match number of table records
			#because it iter only on features with geom and record
			shpIter = shp.iterShapeRecords()
		else:
			shpIter = shp.iterShapes()
		nbFeats = shp.numRecords

		#Create an empty BMesh
		bm = bmesh.new()
		#Extrusion is exponentially slow with large bmesh
		#it's fastest to extrude a small bmesh and then join it to a final large bmesh
		if not self.separateObjects and self.fieldExtrudeName:
			finalBm = bmesh.new()

		progress = -1

		if self.separateObjects:
			layer = bpy.data.collections.new(shpName)
			context.scene.collection.children.link(layer)

		#Main iteration over features
		for i, feat in enumerate(shpIter):

			if self.useDbf:
				shape = feat.shape
				record = feat.record
			else:
				shape = feat

			#Progress infos
			pourcent = round(((i+1)*100)/nbFeats)
			if pourcent in list(range(0, 110, 10)) and pourcent != progress:
				progress = pourcent
				if pourcent == 100:
					print(str(pourcent)+'%')
				else:
					print(str(pourcent), end="%, ")
				sys.stdout.flush() #we need to flush or it won't print anything until after the loop has finished

			#Deal with multipart features
			#If the shape record has multiple parts, the 'parts' attribute will contains the index of
			#the first point of each part. If there is only one part then a list containing 0 is returned
			if (shpType == 'PointZ' or shpType == 'Point'): #point layer has no attribute 'parts'
				partsIdx = [0]
			else:
				try: #prevent "_shape object has no attribute parts" error
					partsIdx = shape.parts
				except Exception as e:
					log.warning('Cannot access "parts" attribute for feature {} : {}'.format(i, e))
					partsIdx = [0]
			nbParts = len(partsIdx)

			#Get list of shape's points
			pts = shape.points
			nbPts = len(pts)

			#Skip null geom
			if nbPts == 0:
				continue #go to next iteration of the loop

			#Reproj geom
			if geoscn.crs != shpCRS:
				pts = rprj.pts(pts)

			#Get extrusion offset
			if self.fieldExtrudeName:
				try:
					offset = float(record[extrudeFieldIdx])
				except Exception as e:
					log.waring('Cannot extract extrusion value for feature {} : {}'.format(i, e))
					offset = 0 #null values will be set to zero

			#Iter over parts
			for j in range(nbParts):

				# EXTRACT 3D GEOM

				geom = [] #will contains a list of 3d points

				#Find first and last part index
				idx1 = partsIdx[j]
				if j+1 == nbParts:
					idx2 = nbPts
				else:
					idx2 = partsIdx[j+1]

				#Build 3d geom
				for k, pt in enumerate(pts[idx1:idx2]):

					if self.elevSource == 'OBJ':
						rcHit = rayCaster.rayCast(x=pt[0]-dx, y=pt[1]-dy)
						z = rcHit.loc.z #will be automatically set to zero if not rcHit.hit

					elif self.elevSource == 'FIELD':
						try:
							z = float(record[zFieldIdx])
						except Exception as e:
							log.warning('Cannot extract elevation value for feature {} : {}'.format(i, e))
							z = 0 #null values will be set to zero

					elif shpType[-1] == 'Z' and self.elevSource == 'GEOM':
						z = shape.z[idx1:idx2][k]

					else:
						z = 0

					geom.append((pt[0], pt[1], z))

				#Shift coords
				geom = [(pt[0]-dx, pt[1]-dy, pt[2]) for pt in geom]


				# BUILD BMESH

				# POINTS
				if (shpType == 'PointZ' or shpType == 'Point'):
					vert = [bm.verts.new(pt) for pt in geom]
					#Extrusion
					if self.fieldExtrudeName and offset > 0:
						vect = (0, 0, offset) #along Z
						result = bmesh.ops.extrude_vert_indiv(bm, verts=vert)
						verts = result['verts']
						bmesh.ops.translate(bm, verts=verts, vec=vect)

				# LINES
				if (shpType == 'PolyLine' or shpType == 'PolyLineZ'):
					verts = [bm.verts.new(pt) for pt in geom]
					edges = []
					for i in range(len(geom)-1):
						edge = bm.edges.new( [verts[i], verts[i+1] ])
						edges.append(edge)
					#Extrusion
					if self.fieldExtrudeName and offset > 0:
						vect = (0, 0, offset) # along Z
						result = bmesh.ops.extrude_edge_only(bm, edges=edges)
						verts = [elem for elem in result['geom'] if isinstance(elem, bmesh.types.BMVert)]
						bmesh.ops.translate(bm, verts=verts, vec=vect)

				# NGONS
				if (shpType == 'Polygon' or shpType == 'PolygonZ'):
					#According to the shapefile spec, polygons points are clockwise and polygon holes are counterclockwise
					#in Blender face is up if points are in anticlockwise order
					geom.reverse() #face up
					geom.pop() #exlude last point because it's the same as first pt
					if len(geom) >= 3: #needs 3 points to get a valid face
						verts = [bm.verts.new(pt) for pt in geom]
						face = bm.faces.new(verts)
						#update normal to avoid null vector
						face.normal_update()
						if face.normal.z < 0: #this is a polygon hole, bmesh cannot handle polygon hole
							pass #TODO
						#Extrusion
						if self.fieldExtrudeName and offset > 0:
							#build translate vector
							if self.extrusionAxis == 'NORMAL':
								normal = face.normal
								vect = normal * offset
							elif self.extrusionAxis == 'Z':
								vect = (0, 0, offset)
							faces = bmesh.ops.extrude_discrete_faces(bm, faces=[face]) #return {'faces': [BMFace]}
							verts = faces['faces'][0].verts
							if self.elevSource == 'OBJ':
								# Making flat roof (TODO add an user input parameter to setup this behaviour)
								z = max([v.co.z for v in verts]) + offset #get max z coord
								for v in verts:
									v.co.z = z
							else:
								##result = bmesh.ops.extrude_face_region(bm, geom=[face]) #return dict {"geom":[BMVert, BMEdge, BMFace]}
								##verts = [elem for elem in result['geom'] if isinstance(elem, bmesh.types.BMVert)] #geom type filter
								bmesh.ops.translate(bm, verts=verts, vec=vect)


			if self.separateObjects:

				if self.fieldObjName:
					try:
						name = record[nameFieldIdx]
					except Exception as e:
						log.warning('Cannot extract name value for feature {} : {}'.format(i, e))
						name = ''
					# null values will return a bytes object containing a blank string of length equal to fields length definition
					if isinstance(name, bytes):
						name = ''
					else:
						name = str(name)
				else:
					name = shpName

				#Calc bmesh bbox
				_bbox = getBBOX.fromBmesh(bm)

				#Calc bmesh geometry origin and translate coords according to it
				#then object location will be set to initial bmesh origin
				#its a work around to bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY')
				ox, oy, oz = _bbox.center
				oz = _bbox.zmin
				bmesh.ops.translate(bm, verts=bm.verts, vec=(-ox, -oy, -oz))

				#Create new mesh from bmesh
				mesh = bpy.data.meshes.new(name)
				bm.to_mesh(mesh)
				bm.clear()

				#Validate new mesh
				mesh.validate(verbose=False)

				#Place obj
				obj = bpy.data.objects.new(name, mesh)
				layer.objects.link(obj)
				context.view_layer.objects.active = obj
				obj.select_set(True)
				obj.location = (ox, oy, oz)

				# bpy operators can be very cumbersome when scene contains lot of objects
				# because it cause implicit scene updates calls
				# so we must avoid using operators when created many objects with the 'separate objects' option)
				##bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY')

				#write attributes data
				for i, field in enumerate(shp.fields):
					fieldName, fieldType, fieldLength, fieldDecLength = field
					if fieldName != 'DeletionFlag':
						if fieldType in ('N', 'F'):
							v = record[i-1]
							if v is not None:
								#cast to float to avoid overflow error when affecting custom property
								obj[fieldName] = float(record[i-1])
						else:
							obj[fieldName] = record[i-1]

			elif self.fieldExtrudeName:
				#Join to final bmesh (use from_mesh method hack)
				buff = bpy.data.meshes.new(".temp")
				bm.to_mesh(buff)
				finalBm.from_mesh(buff)
				bpy.data.meshes.remove(buff)
				bm.clear()

		#Write back the whole mesh
		if not self.separateObjects:

			mesh = bpy.data.meshes.new(shpName)

			if self.fieldExtrudeName:
				bm.free()
				bm = finalBm

			bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.0001)
			bm.to_mesh(mesh)

			#Finish
			#mesh.update(calc_edges=True)
			mesh.validate(verbose=False) #return true if the mesh has been corrected
			obj = bpy.data.objects.new(shpName, mesh)
			context.scene.collection.objects.link(obj)
			context.view_layer.objects.active = obj
			obj.select_set(True)
			bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY')

		#free the bmesh
		bm.free()

		t = perf_clock() - t0
		log.info('Build in %f seconds' % t)

		#Adjust grid size
		bbox.shift(-dx, -dy) #convert shapefile bbox in 3d view space
		adjust3Dview(context, bbox)


		return {'FINISHED'}


classes = [
	IMPORTGIS_OT_shapefile
]


def register():
	for cls in classes:
		try:
			bpy.utils.register_class(cls)
		except ValueError as e:
			log.warning('{} is already registered, now unregister and retry... '.format(cls))
			bpy.utils.unregister_class(cls)
			bpy.utils.register_class(cls)


def unregister():
	for cls in classes:
		bpy.utils.unregister_class(cls)
