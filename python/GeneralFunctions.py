#!/usr/bin/python3
import contextlib
import os
import sys
from os.path import dirname, realpath

sys.path.append(realpath(dirname(__file__)))
import numpy
import scipy
import scipy.spatial
import vtk
from scipy.interpolate import interpolate
from vtk.util.numpy_support import vtk_to_numpy

import Patient
import Perfusion
import scipy.sparse
import tqdm

MajorIDNames = {0: "CoW and Other",
                2: "R. ACA",
                3: "R. MCA",
                4: "L. MCA",
                5: "L. ACA",
                6: "R. PCA",
                7: "L. PCA",
                8: "Cerebellum",
                9: "Brainstem",
                }

MajorIDdict = {
    5: 21,
    4: 22,
    7: 23,
    2: 24,
    3: 25,
    6: 26,
    8: 4,
    9: 30
}

MajorIDdict_inv = {v: k for k, v in MajorIDdict.items()}

StartClusteringIndex = 20

def is_non_zero_file(fpath):
    """
    Return 1 if file exist and have data.
    """
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def TMatrix(scaling, rotation, translation):
    """
    4x4 rotation and translation matrix.
    Order of transformation is translation, X-Y-Z rotation.
    Scaling is isotropuc
    """
    XCos = numpy.cos(numpy.radians(rotation[0]))
    YCos = numpy.cos(numpy.radians(rotation[1]))
    ZCos = numpy.cos(numpy.radians(rotation[2]))

    XSin = numpy.sin(numpy.radians(rotation[0]))
    YSin = numpy.sin(numpy.radians(rotation[1]))
    ZSin = numpy.sin(numpy.radians(rotation[2]))
    Translate = numpy.array(
        [[scaling, 0, 0, translation[0]], [0, scaling, 0, translation[1]], [0, 0, scaling, translation[2]],
         [0, 0, 0, 1]])
    RotateX = numpy.array([[1, 0, 0, 0], [0, XCos, -XSin, 0], [0, XSin, XCos, 0], [0, 0, 0, 1]])
    RotateY = numpy.array([[YCos, 0, YSin, 0], [0, 1, 0, 0], [-YSin, 0, YCos, 0], [0, 0, 0, 1]])
    RotateZ = numpy.array([[ZCos, -ZSin, 0, 0], [ZSin, ZCos, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    return numpy.dot(RotateZ, numpy.dot(RotateY, numpy.dot(RotateX, Translate)))


def TMatrixNonUniform(scaling, rotation, translation):
    """
    4x4 rotation and translation matrix.
    Order of transformation is translation, X-Y-Z rotation.
    Scaling is isotropuc
    """
    XCos = numpy.cos(numpy.radians(rotation[0]))
    YCos = numpy.cos(numpy.radians(rotation[1]))
    ZCos = numpy.cos(numpy.radians(rotation[2]))

    XSin = numpy.sin(numpy.radians(rotation[0]))
    YSin = numpy.sin(numpy.radians(rotation[1]))
    ZSin = numpy.sin(numpy.radians(rotation[2]))
    Translate = numpy.array(
        [[scaling[0], 0, 0, translation[0]], [0, scaling[1], 0, translation[1]], [0, 0, scaling[2], translation[2]],
         [0, 0, 0, 1]])
    RotateX = numpy.array([[1, 0, 0, 0], [0, XCos, -XSin, 0], [0, XSin, XCos, 0], [0, 0, 0, 1]])
    RotateY = numpy.array([[YCos, 0, YSin, 0], [0, 1, 0, 0], [-YSin, 0, YCos, 0], [0, 0, 0, 1]])
    RotateZ = numpy.array([[ZCos, -ZSin, 0, 0], [ZSin, ZCos, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    return numpy.dot(RotateZ, numpy.dot(RotateY, numpy.dot(RotateX, Translate)))


def TriangleToArea(nodes):
    """
    Take a list of node positions and return the area using Heron's formula.
    """
    a = distancebetweenpoints(nodes[0], nodes[1])
    b = distancebetweenpoints(nodes[1], nodes[2])
    c = distancebetweenpoints(nodes[0], nodes[2])
    s = (a + b + c) / 2
    return numpy.sqrt(s * (s - a) * (s - b) * (s - c))


def meanpos(nodes):
    """
    Take a list of node positions and return the centroid (mean position).
    """
    x = numpy.mean([nodes[0][0], nodes[1][0], nodes[2][0]])
    y = numpy.mean([nodes[0][1], nodes[1][1], nodes[2][1]])
    z = numpy.mean([nodes[0][2], nodes[1][2], nodes[2][2]])
    return [x, y, z]


def distancebetweenpoints(p1, p2):
    """
    Calculate the euclidean distance between two points.
    """
    dsq = sum([numpy.square(p1[i] - p2[i]) for i in range(0, len(p1))])
    return numpy.sqrt(dsq)


def Merge(patient, brava):
    """
    Merge a patient network with a brava network.
    """
    print("Merging patient network with another network.")
    matchpoints = ["R. ACA, A2", "R. MCA", "L. MCA", "L. ACA, A2", "R. PCA, P2", "L. PCA, P2"]
    bravaid = [2, 3, 4, 5, 6, 7]

    vesselnumbers = patient.Topology.Vessels[-1].ID  # last vessel id, renumber appended vessels from here

    # calculate the wk parameters without the merged network
    # make sure that we can find the relevant vessels
    with contextlib.redirect_stdout(None):
        patient.calculate_wk_parameters()
        patient.Topology.RedefineDirection()

    _, _, _, gens = patient.Topology.GetDownstreamVessels("Ascending Aorta")
    for index, gen in enumerate(gens):
        for vessel in gen:
            vessel.GenerationNumber = index

    for index, name in enumerate(matchpoints):
        vessels, _, _, _ = patient.Topology.GetDownstreamVessels(name)
        [vessel.SetMajorVesselID(bravaid[index]) for vessel in vessels]

    # find the relevent vessels and nodes in the patient network!
    # these are by definition the end nodes of the vessel
    nodedict = brava.Topology.MapNodesToVessel()

    # first find all the major branches and their downsteam vessels.
    mainvesselscouplingbif = []
    endvesselsforwkupdate = []
    for loop in range(0, len(bravaid)):
        patientvesselname = bravaid[loop]  # MajorvesselID
        # extract vessels that match name
        vesselsmatch = [brava.Topology.Vessels[i] for i in range(0, len(brava.Topology.Vessels)) if
                        brava.Topology.Vessels[i].MajorVesselID == patientvesselname]

        # The target vessel is connected to other vessels that are not part of the set
        # check for a bifurcation is linked to another vessel type
        comparevesselbrava = []
        vesselnodes = nodedict[patientvesselname]

        for vessel in vesselsmatch:
            # find the bifurcations of the vessels
            bifurcations = [node for node in vessel.Nodes[0].Connections if node in brava.Topology.BifurcationNodes]
            [bifurcations.append(node) for node in vessel.Nodes[-1].Connections if
             node in brava.Topology.BifurcationNodes]

            for bif in bifurcations:
                targetnodes = bif.Connections
                for node in targetnodes:
                    if not (node in vesselnodes):
                        comparevesselbrava = vessel
                        BiflickedtoCow = bif

        firstnode = [node for node in comparevesselbrava.Nodes if node in BiflickedtoCow.Connections][0]
        mainvesselscouplingbif.append((comparevesselbrava, firstnode))
        firstnode.RemoveConnection(BiflickedtoCow)

    # for each major branch, find the downsteam vessels and rename them.
    brava.Topology.InletNodes = [(node[1], "") for node in mainvesselscouplingbif]
    with contextlib.redirect_stdout(None):
        brava.Topology.RedefineDirection()
    for loop in range(0, len(bravaid)):
        vessels, bifurcations, vesselnames, gens = brava.Topology.GetDownstreamVessels(mainvesselscouplingbif[loop][0])
        gens = brava.Topology.ReorderGens(gens)
        for index, gen in enumerate(gens):
            if index == 0:
                name = matchpoints[loop]
            else:
                name = matchpoints[loop] + ", gen:" + str(index)
            for index2, vessel in enumerate(gen):
                if index > 0:
                    vessel.Name = name + "," + str(index2)
                else:
                    vessel.Name = name

        patientvesselname = matchpoints[loop]
        vesselspatient, bifurcationspatient, vesselnamespatient, genspatient = patient.Topology.GetDownstreamVessels(
            patientvesselname)

        # find the equivalent vessels in the bravaset
        # if the patient anatomy contains more vessels beyond the first one.

        matchingvessels = [(genspatient[0][0], gens[0][0])]
        ends = matchingvessels
        while len(vesselspatient) != len(matchingvessels):
            newends = []
            for end in ends:
                patientvessel, bravavessel = end
                _, _, _, newgens = brava.Topology.GetDownstreamVessels(bravavessel)
                _, _, _, newgenspatient = patient.Topology.GetDownstreamVessels(patientvessel)
                newgenspatient = patient.Topology.ReorderGens(newgenspatient)
                newgens = brava.Topology.ReorderGens(newgens)

                if len(newgenspatient) > 1 and len(newgens) > 1:
                    # error if there are anastomoses since these are not present in the dataset
                    if len(newgenspatient[1]) > 2:
                        print("Warning, likely there are anastomoses.")
                        newgenspatient[1] = [vessel for vessel in newgenspatient[1] if
                                             vessel.GenerationNumber > newgenspatient[0][0].GenerationNumber]

                    for index, element in enumerate(newgenspatient[1]):
                        matchingvessels.append((element, newgens[1][index]))
                        newends.append((element, newgens[1][index]))
            ends = newends

        # note that the connectivity in this method might be different from the one above.
        # this seems to lead to some duplications if the ends of a major branch are of different generations
        # patientvessel, bravavessel = (genspatient[0][0], gens[0][0])
        # bravavessels, _, _, newgens = brava.Topology.GetDownstreamVessels(bravavessel)
        # patientvessels, _, _, newgenspatient = patient.Topology.GetDownstreamVessels(patientvessel)
        # matchingvessels = []
        # names = []
        # newgenspatient = patient.Topology.ReorderGens(newgenspatient)
        # newgens = brava.Topology.ReorderGens(newgens)
        #
        # for index, gen in enumerate(newgenspatient):
        #     for indexves, vessel in enumerate(gen):
        #         bravavessel = newgens[index][indexves]
        #         matchingvessels.append((vessel, bravavessel))
        #         names.append((vessel.Name, bravavessel.Name))

        finalends = [i for i in matchingvessels if len(i[0].Nodes[-1].Connections) == 1]
        # for each matching vessel, find the downstream vessel in the brava network
        donorbranches = [brava.Topology.GetDownstreamVessels(vessel[1]) for index, vessel in enumerate(finalends)]
        if len(finalends) > 1:
            for i in finalends:
                print(i[1].Name)
        # Interpolation of the vessels for visualisation
        for index, vessel in enumerate(matchingvessels):
            # Interpolate the connection between the brava and the patient network
            # We only want to change the position in 3d coordinates for visualisation.
            # The actual length used in the simulation should not change!
            if vessel[0].Name in matchpoints:
                # if the vessel is one of the major vessels, we just want to extend it to the next bifurcation.
                nodebrava = vessel[1].Nodes[-1]
                # nodepatient = vessel[0].Nodes[-1]
                x = [node.Position[0] for node in vessel[0].Nodes] + [
                    nodebrava.Position[0]]
                y = [node.Position[1] for node in vessel[0].Nodes] + [
                    nodebrava.Position[1]]
                z = [node.Position[2] for node in vessel[0].Nodes] + [
                    nodebrava.Position[2]]
                s = [node.LengthAlongVessel for node in vessel[0].Nodes] + [
                    vessel[0].Length + vessel[1].Length]
            else:
                # if the vessel is not one of the major vessel, it is an added patient vessel and we want to overwrite all the positions.
                # The original positions are garbage.
                x = [node.Position[0] for node in vessel[1].Nodes]
                y = [node.Position[1] for node in vessel[1].Nodes]
                z = [node.Position[2] for node in vessel[1].Nodes]
                s = [node.LengthAlongVessel for node in vessel[1].Nodes]

            interpolation = 'linear'
            f1 = interpolate.interp1d(s, x, kind=interpolation)
            f2 = interpolate.interp1d(s, y, kind=interpolation)
            f3 = interpolate.interp1d(s, z, kind=interpolation)

            for node in vessel[0].Nodes:
                fractionalongvessellength = node.LengthAlongVessel / vessel[0].Length
                fractionalongvessellength = min(1.0, fractionalongvessellength) * s[-1]
                newx = f1(fractionalongvessellength)
                newy = f2(fractionalongvessellength)
                newz = f3(fractionalongvessellength)
                node.SetPosition([newx, newy, newz])
            # Update the vessel interpolation functions
            vessel[0].UpdateInterpolationFunctions()

        # Extend the vessels of the patient with the vessels from the brava set
        for index, vessel in enumerate(finalends):  # only end vessels!
            # calculate the scaling factor
            # apply scaling to downsteam vessels
            # add vessels to patient
            # add bifurcations to patient
            # add the connections to the bifurcations
            # NOTE: the brava set does not have the ACoA so these scaling factors are off.
            # For this reason, the length ratio is not usable
            scalingsfactor = vessel[0].MeanRadius / vessel[1].MeanRadius
            vessels, bifurcations, names, gens = donorbranches[index]

            connectingvessel = vessel[0]
            endvesselsforwkupdate.append(connectingvessel)
            if len(gens) <= 1:
                # check if patient network continues
                upstream = patient.Topology.GetDownstreamVessels(connectingvessel)
                if len(upstream[0]) > 1:
                    print("Error, patient network continues after the bravaset.")
                # print("No vessels to attach, updating outlet position.")
                # vessel[0].Nodes[-1].SetPosition(vessel[1].Nodes[-1].Position)
                continue

            otherconnetingvessels = gens[1]
            donorvessels = [ves for ves in list(vessels) if ves is not vessel[1]]

            [vessel.ScaleRadius(scalingsfactor) for vessel in donorvessels]

            bif = [bif for bif in bifurcations if bif in otherconnetingvessels[1].Nodes[0].Connections][0]
            bif.ResetConnections()
            bif.AddConnection(connectingvessel.Nodes[-1])
            bif.AddConnection(otherconnetingvessels[1].Nodes[0])
            bif.AddConnection(otherconnetingvessels[0].Nodes[0])
            connectingvessel.Nodes[-1].AddConnection(bif)
            otherconnetingvessels[1].Nodes[0].AddConnection(bif)
            otherconnetingvessels[0].Nodes[0].AddConnection(bif)
            # transfer to patient
            for vessel in donorvessels:
                vesselnumbers += 1
                vessel.SetID(vesselnumbers)
            [patient.Topology.Vessels.append(element) for element in donorvessels]
            [patient.Topology.BifurcationNodes.append(element) for element in bifurcations]
            [patient.Topology.Nodes.append(element) for sublist in donorvessels for element in sublist.Nodes]
            [patient.Topology.Nodes.append(element) for element in bifurcations]

    # add all the other things
    patient.Topology.FindOutletNodes()
    patient.Trees = brava.Trees
    patient.Perfusion = brava.Perfusion
    # The donor networks have outlets that do not get transferred. These need to be removed.
    # If the coupling point node does not have a number, remove it from the list.
    pointstoremove = [point for point in patient.Perfusion.CouplingPoints if not (point.Node in patient.Topology.Nodes)]
    [patient.Perfusion.RemoveCouplingPoint(point) for point in pointstoremove]
    # patient.Topology.RedefineDirection()
    patient.Topology.UpdateTopology()

    # [node.ResetWK() for node in patient.Topology.OutletNodes]
    # patient.calculate_wk_parameters()

    # update the wk_parameters at the ends of the added vessels
    # taking into account that the attached trees have resistance that we have to take into account
    visc = float(patient.ModelParameters["BLOOD_VISC"])
    density = float(patient.ModelParameters["Density"])
    print("Removing resistance added by the merging of networks.")
    for name in endvesselsforwkupdate:
        patient.Topology.DownStreamResistance(name, visc, density)


def TransformFile(file, transformmatrix):
    """
    Transform all points in a vtp file.
    This is not the way that paraview does its transformations.
    """

    if file[-3:] == "vtp":
        reader = vtk.vtkXMLPolyDataReader()
    elif file[-3:] == "ply":
        reader = vtk.vtkPLYReader()
    else:
        print("Error: unreadable file.")
        return 1
    reader.SetFileName(file)
    reader.Update()
    data = reader.GetOutput()

    pos_vtk = reader.GetOutput().GetPoints().GetData()
    pos = vtk_to_numpy(pos_vtk)
    nodes = vtk.vtkPoints()
    for point in pos:
        vec = numpy.array([[point[0]], [point[1]], [point[2]], [1]])
        position = numpy.dot(transformmatrix, vec)
        nodes.InsertNextPoint(position[:-1])
    data.SetPoints(nodes)

    # export to new file
    writer = vtk.vtkXMLPolyDataWriter()
    file = "testingtrans.vtp"
    writer.SetFileName(file)
    writer.SetInputData(data)
    writer.Write()


class MSHfile:
    def __init__(self):
        self.MeshFormat = []
        self.PhysicalNames = []
        self.Nodes = []
        self.Elements = []

    def Loadfile(self, file):
        print("Loading MSH: %s" % file)
        mesh_raw = [i.strip('\n') for i in open(file)]
        mesh = [i.split(' ') for i in mesh_raw]

        startelementsFormat = mesh_raw.index("$MeshFormat")
        endelementsFormat = mesh_raw.index("$EndMeshFormat")
        startelementsNames = mesh_raw.index("$PhysicalNames")
        endelementsNames = mesh_raw.index("$EndPhysicalNames")
        startelementsNodes = mesh_raw.index("$Nodes")
        endelementsNodes = mesh_raw.index("$EndNodes")
        startelements = mesh_raw.index("$Elements")
        endelements = mesh_raw.index("$EndElements")

        self.MeshFormat = [mesh_raw[i] for i in range(startelementsFormat + 1, endelementsFormat)]
        self.PhysicalNames = [[int(mesh[i][0]), int(mesh[i][1]), mesh[i][2]] for i in
                              range(startelementsNames + 2, endelementsNames)]
        self.Nodes = [[int(mesh[i][0]), float(mesh[i][1]), float(mesh[i][2]), float(mesh[i][3])] for i
                      in range(startelementsNodes + 2, endelementsNodes)]
        self.Elements = [[int(x) for x in mesh[i] if x] for i in range(startelements + 2, endelements)]

    def Writefile(self, file):
        print("Writing MSH: %s" % file)
        with open(file, 'w') as f:
            f.write("$MeshFormat\n")
            for i in self.MeshFormat:
                f.write(i + "\n")
            f.write("$EndMeshFormat\n")

            f.write("$PhysicalNames\n")
            f.write(str(len(self.PhysicalNames)) + "\n")
            for i in self.PhysicalNames:
                line = ' '.join(str(x) for x in i)
                f.write(line + "\n")
            f.write("$EndPhysicalNames\n")

            f.write("$Nodes\n")
            f.write(str(len(self.Nodes)) + "\n")
            for i in self.Nodes:
                line = ' '.join(str(x) for x in i)
                f.write(line + "\n")
            f.write("$EndNodes\n")

            f.write("$Elements\n")
            f.write(str(len(self.Elements)) + "\n")
            for i in self.Elements:
                line = ' '.join(str(x) for x in i)
                f.write(line + "\n")
            f.write("$EndElements\n")

    def GetElements(self, ids):
        data = [[index, element] for index, element in enumerate(self.Elements) if
                int(element[4]) in ids or int(element[3]) in ids]
        indexes = [i[0] for i in data]
        elements = [i[1] for i in data]
        return elements, indexes

    def GetSurfaceCentroids(self, ids):
        elements, indexes = self.GetElements(ids)

        triangles = [[i[-1], i[-2], i[-3]] for i in elements]

        trianglespos = [[self.Nodes[triangle[0] - 1][1:],
                         self.Nodes[triangle[1] - 1][1:],
                         self.Nodes[triangle[2] - 1][1:]
                         ] for triangle in triangles]

        positions = [meanpos(i) for i in trianglespos]
        return positions, elements, indexes

    def AreaRegion(self, regionid):
        elements, indexes = self.GetElements([regionid])
        triangles = [[i[-1], i[-2], i[-3]] for i in elements]
        trianglespos = [[self.Nodes[triangle[0] - 1][1:],
                         self.Nodes[triangle[1] - 1][1:],
                         self.Nodes[triangle[2] - 1][1:]
                         ] for triangle in triangles]
        areas = [TriangleToArea(triangle) for triangle in trianglespos]
        totalarea = sum(areas)
        return totalarea, len(triangles)


def VesselMapping(mappingfile, patient):
    # For each vessel/node in patient, map them to the same vessel in the mapping file.
    # For now, only 55 vessels are included.
    mapping = Patient.Patient()
    mapping.LoadVTPFile(mappingfile)
    mapping.Topology.UpdateVesselAtlas()
    for vessel in patient.Topology.Vessels:
        if vessel.ID < 56:  # no mapping beyond these vessels
            mappedvessel = mapping.Topology.VesselAtlas[vessel.ID]
            vessel.MajorVesselID = mappedvessel.MajorVesselID
            for node in vessel.Nodes:
                fractionalongvessellength = node.LengthAlongVessel / vessel.Length
                if vessel.ID == 15:  # starts at the vessel end
                    positioninmap = (1.0 - min(1.0, fractionalongvessellength)) * mappedvessel.Length
                else:
                    positioninmap = min(1.0, fractionalongvessellength) * mappedvessel.Length
                newposx = mappedvessel.InterpolationFunctions[0](positioninmap)
                newposy = mappedvessel.InterpolationFunctions[1](positioninmap)
                newposz = mappedvessel.InterpolationFunctions[2](positioninmap)
                node.SetPosition([newposx, newposy, newposz])
            # update the interpolation functions
            vessel.UpdateInterpolationFunctions()

    # update the bifurcation positions as well
    for bif in patient.Topology.BifurcationNodes:
        bif.SetPosition(next(iter(bif.Connections)).Position)


def MapClusteringToMSH(file, file2, datafolder):
    print("Mapping Clustering to a MSH file.")

    outputfilemsh = datafolder + 'clustered_mesh.msh'
    regionsIDs = [4, 21, 22, 23, 24, 25, 26, 30]  # extracted from the msh file

    Surface1 = Perfusion.Perfusion()
    Surface1.PrimalGraph.LoadSurface(file)
    centers = Surface1.PrimalGraph.GetTriangleCentroids()
    sys.setrecursionlimit(10000)
    KDTree = scipy.spatial.KDTree(centers)

    mshmesh = MSHfile()
    mshmesh.Loadfile(file2)
    positions, elements, indexes = mshmesh.GetSurfaceCentroids(regionsIDs)

    # euclidean distance between outlets and surface
    MinDistance, MinDistanceIndex = KDTree.query(positions, k=1)

    # MinDistance, MinDistanceIndex = KDTree.query([i for i in meshdata2.positions], k=1)
    regiondict = MajorIDdict

    ClusterDict = {}
    for index, trianglenumber in enumerate(MinDistanceIndex):
        region = Surface1.PrimalGraph.map[trianglenumber]
        mshmesh.Elements[indexes[index]][3] = regiondict[region]
        cluster = Surface1.PrimalGraph.PolygonColour[trianglenumber]
        # mshmesh.Elements[indexes[index]][3] = cluster + StartClusteringIndex
        mshmesh.Elements[indexes[index]][4] = cluster + StartClusteringIndex
        ClusterDict[cluster + StartClusteringIndex] = regiondict[region]  # store ids in dict

    # regionsName = ["Cerebellum","R. ACA, A2", "R. MCA", "L. MCA", "L. ACA, A2", "R. PCA, P2", "L. PCA, P2","Brainstem"]
    # regionsIDs = [4, 21, 22, 23, 24, 25, 26, 30]  # region 0, 1, 2, 3, 4, 5

    # list of cluster names
    clusterids = [i for i in range(0, len(set(Surface1.PrimalGraph.PolygonColour)))]
    uniqueclusters = len(set(clusterids))
    namesclusters = [[2, StartClusteringIndex + i, "\"Cluster_" + str(i) + '\"'] for i in range(0, uniqueclusters)]
    mshmesh.PhysicalNames += namesclusters

    indexestoremove = []
    [indexestoremove.append(i) for i, line in enumerate(mshmesh.PhysicalNames) if line[1] in regionsIDs]
    newnames = [line for i, line in enumerate(mshmesh.PhysicalNames) if not (i in indexestoremove)]
    mshmesh.PhysicalNames = newnames

    mshmesh.Writefile(outputfilemsh)




def WriteFlowFilePerfusionModel(Clusterflowdatafile, clusteringfile, datafolder):
    clusterdata = [i.strip('\n') for i in open(clusteringfile)]
    clusterdata = [i.split(',') for i in clusterdata][1:]
    MajorvesselIDs = [int(i[5]) for i in clusterdata]

    flowdata = [i.strip('\n') for i in open(Clusterflowdatafile)]
    flowdata = [i.split(',') for i in flowdata][1:]

    with open(datafolder + 'boundary_condition_file.csv', 'w') as f:
        # f.write("Cluster ID,Major Region ID,Volume flow rate (mL/s),Pressure (Pa)\n")
        f.write("# region I,Q [ml/s],p [Pa],feeding artery ID\n")
        for index, i in enumerate(flowdata):
            # f.write("%d,%d,%f,%f\n" % (StartClusteringIndex + index, MajorvesselIDs[index], float(i[1]), float(i[3])))
            f.write("%d,%f,%f,%d\n" % (MajorIDdict[MajorvesselIDs[index]], float(i[1]), float(i[3]), StartClusteringIndex + index,))
