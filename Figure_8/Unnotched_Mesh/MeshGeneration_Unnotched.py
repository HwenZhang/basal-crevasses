#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 13:12:15 2021

@author: eart0487
"""

import gmsh
import meshio
import sys
import os
import numpy as np

gmsh.initialize(sys.argv)
gmsh.model.add("mesh")

# half width of patch
w = 1
# half length of strip
L = 10
# Ice thickness
H = 1
# Mesh size around crack tip
Lc1 = 0.002
# Mesh size on the bottom
Lc2 = 0.05
# Mesh Size on sticky patch
Lc3 = 0.2

opengmsh=True
# opengmsh=False

point1 = gmsh.model.geo.addPoint(-L, 0.0, 0, Lc2, 1)
point2 = gmsh.model.geo.addPoint(-w, 0.0, 0, Lc2, 2)
point3 = gmsh.model.geo.addPoint(w, 0.0, 0, Lc2, 3)
point4 = gmsh.model.geo.addPoint(L, 0.0, 0, Lc2, 4)
point5 = gmsh.model.geo.addPoint(L, H, 0, Lc2, 5)
point6 = gmsh.model.geo.addPoint(-L, H, 0, Lc2, 6)

line1 = gmsh.model.geo.addLine(1, 2, 1)
line2 = gmsh.model.geo.addLine(2, 3, 2) # sticky patch
line3 = gmsh.model.geo.addLine(3, 4, 3) # left crack wall
line4 = gmsh.model.geo.addLine(4, 5, 4) # right crack wall
line5 = gmsh.model.geo.addLine(5, 6, 5) # 
line6 = gmsh.model.geo.addLine(6, 1, 6) # left boundary

boundary_loop = gmsh.model.geo.addCurveLoop([1, 2, 3, 4, 5, 6], 7)
plane_surface = gmsh.model.geo.addPlaneSurface([7], 8)

# Neumann boundary condition
boundary0 = gmsh.model.addPhysicalGroup(1, [1], 0)
gmsh.model.setPhysicalName(1, boundary0, "Boundary0")

boundary1 = gmsh.model.addPhysicalGroup(1, [2], 1)
gmsh.model.setPhysicalName(1, boundary1, "Boundary1")

boundary2 = gmsh.model.addPhysicalGroup(1, [3], 2)
gmsh.model.setPhysicalName(1, boundary2, "Boundary2")

boundary3 = gmsh.model.addPhysicalGroup(1, [4], 3)
gmsh.model.setPhysicalName(1, boundary3, "Boundary3")

boundary4 = gmsh.model.addPhysicalGroup(1, [5], 4)
gmsh.model.setPhysicalName(1, boundary4, "Boundary4")

boundary5 = gmsh.model.addPhysicalGroup(1, [6], 5)
gmsh.model.setPhysicalName(1, boundary5, "Boundary5")


gmsh.model.geo.synchronize()


# Mesh Size Initiation
gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(1, "PointsList", [3])

gmsh.model.mesh.field.add("Distance", 2)
gmsh.model.mesh.field.setNumbers(2, "CurvesList", [6])
gmsh.model.mesh.field.setNumber(2, "NumPointsPerCurve", 100)

gmsh.model.mesh.field.add("Box", 3)
gmsh.model.mesh.field.setNumber(3, "VIn", Lc1)
gmsh.model.mesh.field.setNumber(3, "VOut", Lc3)
gmsh.model.mesh.field.setNumber(3, "XMin", w-w/400.0)
gmsh.model.mesh.field.setNumber(3, "XMax", w+w/400.0)
gmsh.model.mesh.field.setNumber(3, "YMin", 0.0)
gmsh.model.mesh.field.setNumber(3, "YMax", 1.0)
gmsh.model.mesh.field.setNumber(3, "Thickness", 0.3)

# gmsh.model.mesh.field.add("MathEval", 4)
# gmsh.model.mesh.field.setString(4, "F",
#                                 "Cos(4*3.14*x) * Sin(4*3.14*y) / 10 + 0.101")


gmsh.model.mesh.field.add("Min", 5)
gmsh.model.mesh.field.setNumbers(5, "FieldsList", [3])

gmsh.model.mesh.field.setAsBackgroundMesh(5)

# quadrilateral elements
# gmsh.model.mesh.setRecombine(2, plane_surface)
# gmsh.option.setNumber("Mesh.Algorithm", 8)
# gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)

gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.SaveAll", 1)
# This will prevent over-refinement due to small mesh sizes on the boundary.

gmsh.model.mesh.generate(2)
gmsh.write("mesh.msh2")

# change file name
os.rename("mesh.msh2", "mesh.msh")

if opengmsh:
    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()
    
gmsh.finalize()

#============================================================
# convert msh file to xdmf file
mesh_from_file = meshio.read("mesh.msh")

prune_z=True
triangle_mesh = meshio.Mesh(points=mesh_from_file.points,
                      cells=[("triangle", mesh_from_file.get_cells_type("triangle"))],
                      cell_data={"Subdomain": [mesh_from_file.cell_data_dict['gmsh:physical']['triangle']]},
                      field_data=mesh_from_file.field_data)
triangle_mesh.prune_z_0()
meshio.write("Unnotched_mesh.xdmf", triangle_mesh)

line_mesh = meshio.Mesh(points=mesh_from_file.points,
                      cells=[("line", mesh_from_file.get_cells_type("line"))],
                      cell_data={"Boundaries": [mesh_from_file.cell_data_dict['gmsh:geometrical']['line']]},
                      field_data=mesh_from_file.field_data)
line_mesh.prune_z_0()
meshio.write("Unnotched_facet_mesh.xdmf", line_mesh)