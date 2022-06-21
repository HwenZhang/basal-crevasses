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

cindex = 0
vlambda1 = np.arange(0.04,0.1,0.02)
vlambda2 = np.arange(0.1,1,0.05)
vlambda = np.append(vlambda1,vlambda2)
Nlambda = np.shape(vlambda)[0]

# folder index
while cindex < Nlambda :
    
    gmsh.initialize(sys.argv)
    gmsh.model.add("mesh")
    
    # half width of patch
    w = 1
    # half length of strip
    L = 10
    # Ice thickness
    H = 1
    # crack length
    C = vlambda[cindex]
    # crack width
    dw = C / 80.0
    # Refined Radius
    Rr = max(C / 10.0, w / 20.0)
    # Mesh size around crack tip
    
    # local refinement around the crack tip
    Lc1 = C / 200.0
    # refienment on the sticky patch
    Lc2 = min([0.05,Lc1*40])
    # far-field mesh size
    Lc3 = 0.2
    
    opengmsh=True
    opengmsh=False
    
    point1 = gmsh.model.geo.addPoint(-L, 0.0, 0, Lc2, 1)
    point2 = gmsh.model.geo.addPoint(-w, 0.0, 0, Lc2, 2)
    point3 = gmsh.model.geo.addPoint(w-dw, 0.0, 0, Lc2, 3)
    point4 = gmsh.model.geo.addPoint(w, C, 0, Lc2, 4) # crack tip
    point5 = gmsh.model.geo.addPoint(w+dw, 0.0, 0, Lc2, 5)
    point6 = gmsh.model.geo.addPoint(L, 0.0, 0, Lc2, 6)
    point7 = gmsh.model.geo.addPoint(L, H, 0, Lc2, 7)
    point8 = gmsh.model.geo.addPoint(-L, H, 0, Lc2, 8)
    
    line1 = gmsh.model.geo.addLine(1, 2, 1)
    line2 = gmsh.model.geo.addLine(2, 3, 2) # sticky patch
    line3 = gmsh.model.geo.addLine(3, 4, 3) # left crack wall
    line4 = gmsh.model.geo.addLine(4, 5, 4) # right crack wall
    line5 = gmsh.model.geo.addLine(5, 6, 5) # 
    line6 = gmsh.model.geo.addLine(6, 7, 6) # right boundary
    line7 = gmsh.model.geo.addLine(7, 8, 7) # top
    line8 = gmsh.model.geo.addLine(8, 1, 8) # left boundary
    
    boundary_loop = gmsh.model.geo.addCurveLoop([1, 2, 3, 4, 5, 6, 7, 8], 9)
    plane_surface = gmsh.model.geo.addPlaneSurface([9], 10)
    
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
    
    boundary6 = gmsh.model.addPhysicalGroup(1, [7], 6)
    gmsh.model.setPhysicalName(1, boundary6, "Boundary6")
    
    boundary7 = gmsh.model.addPhysicalGroup(1, [8], 7)
    gmsh.model.setPhysicalName(1, boundary7, "Boundary7")
    
    
    gmsh.model.geo.synchronize()
    
    
    # Mesh Size Initiation
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "PointsList", [4])

    gmsh.model.mesh.field.add("Distance", 2)
    gmsh.model.mesh.field.setNumbers(2, "CurvesList", [2])
    gmsh.model.mesh.field.setNumber(2, "NumPointsPerCurve", 100)

    gmsh.model.mesh.field.add("Threshold", 3)
    gmsh.model.mesh.field.setNumber(3, "InField", 1)
    gmsh.model.mesh.field.setNumber(3, "SizeMin", Lc1)
    gmsh.model.mesh.field.setNumber(3, "SizeMax", Lc3)
    gmsh.model.mesh.field.setNumber(3, "DistMin", C / 20)
    gmsh.model.mesh.field.setNumber(3, "DistMax", C / 4)

    gmsh.model.mesh.field.add("Threshold", 4)
    gmsh.model.mesh.field.setNumber(4, "InField", 2)
    gmsh.model.mesh.field.setNumber(4, "SizeMin", Lc2)
    gmsh.model.mesh.field.setNumber(4, "SizeMax", Lc3)
    gmsh.model.mesh.field.setNumber(4, "DistMin", w / 20.0)
    gmsh.model.mesh.field.setNumber(4, "DistMax", w / 10.0)

    gmsh.model.mesh.field.add("MathEval", 5)
    gmsh.model.mesh.field.setString(5, "F",
                                    str(Lc1)+"+ F1 / "+str(C)+"*"+str(Lc2-Lc1))
    # gmsh.model.mesh.field.setString(5, "F", str(Lc1)+"*"+"F1"+"/"+str(Rr))

    gmsh.model.mesh.field.add("Min", 6)
    gmsh.model.mesh.field.setNumbers(6, "FieldsList", [3, 4, 5])
    
    gmsh.model.mesh.field.setAsBackgroundMesh(6)
    
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
    meshio.write("Notched_mesh_C_"+str(cindex)+".xdmf", triangle_mesh)
    
    line_mesh = meshio.Mesh(points=mesh_from_file.points,
                          cells=[("line", mesh_from_file.get_cells_type("line"))],
                          cell_data={"Boundaries": [mesh_from_file.cell_data_dict['gmsh:geometrical']['line']]},
                          field_data=mesh_from_file.field_data)
    line_mesh.prune_z_0()
    meshio.write("Notched_facet_mesh_C_"+str(cindex)+".xdmf", line_mesh)
    
    cindex = cindex + 1