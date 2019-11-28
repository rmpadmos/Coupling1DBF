#!/usr/bin/python3
import sys
from os.path import dirname, realpath

sys.path.append(realpath(dirname(__file__)))
import pyacvd
import pyvista
import scipy.spatial
import vtk

import GeneralFunctions
import Patient


def remesh(surfacefile, numbertriangles=40000, output="remeshed.vtp"):
    """
    Remesh a surface mesh using using voronoi clustering. Source and module at https://pypi.org/project/pyacvd/
    :param surfacefile: Surfacefile to be remeshed to a uniform triangulation.
    :param numbertriangles: Number of triangles that the surface will have after the remeshing. Default:40000
    :param output: output file name
    :return: Nothing
    """
    print("Remeshing surface.")
    if surfacefile[-3:] == "vtp":
        reader = vtk.vtkXMLPolyDataReader()
    elif surfacefile[-3:] == "ply":
        reader = vtk.vtkPLYReader()
    else:
        print("Input is not a ply or vtp file.")
        return
    reader.SetFileName(surfacefile)
    reader.Update()

    p = reader.GetOutput()
    surf = pyvista.PolyData(p)
    clus = pyacvd.Clustering(surf)

    clus.subdivide(3)
    clus.cluster(numbertriangles)
    remesh = clus.create_mesh()

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(output)
    writer.SetInputData(remesh)
    writer.Write()


def MapMeshtoMSH(filevtp, filemsh, output="PialSurface.vtp"):
    """
    Apply the mapping on one surface to another surface.
    :param filevtp: The remeshed surface file, see the remesh function.
    :param filemsh: THe file containing the VEASL mapping.
    :param output: Filename of the resulting file.
    :return: Nothing
    """
    print("Mapping msh to vtp.")
    regionsIDs = [4, 21, 22, 23, 24, 25, 26, 30]
    patient = Patient.Patient()
    patient.Perfusion.LoadPrimalGraph(filevtp)
    centroids = patient.Perfusion.PrimalGraph.GetTriangleCentroids()

    msh = GeneralFunctions.MSHfile()
    msh.Loadfile(filemsh)
    positions, elements, indexes = msh.GetSurfaceCentroids(regionsIDs)

    sys.setrecursionlimit(10000)
    KDTree = scipy.spatial.KDTree(positions)
    MinDistance, MinDistanceIndex = KDTree.query(centroids, k=1)

    regiondict = GeneralFunctions.MajorIDdict_inv
    regionsIDs = [regiondict[elements[trianglenumber][3]] for index, trianglenumber in enumerate(MinDistanceIndex)]
    patient.Perfusion.PrimalGraph.PolygonColour = regionsIDs
    patient.Perfusion.PrimalGraph.File = output
    patient.Perfusion.PrimalGraph.GraphToVTP("")


if __name__ == '__main__':
    # remesh the mesh to a uniform triangulation.
    remesh("boundary_4&21&22&23&24&25&26&30.ply", numbertriangles=40000)

    # apply the same mapping to the remeshed file
    MapMeshtoMSH("remeshed.vtp", "labelled_vol_mesh.msh")
