#!/usr/bin/python3
import collections
import contextlib
import csv
from collections import namedtuple
import sys
from os.path import isfile, dirname, realpath, exists

sys.path.append(realpath(dirname(__file__)))
import networkx as nx
import scipy.sparse
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import math
import numpy
import scipy
import Vessel
import Node
import GeneralFunctions
import BloodflowEquations
import scipy.optimize
from multiprocessing.pool import ThreadPool as Pool

class Topology:
    def __init__(self):
        self.Nodes = []
        self.InletNodes = []
        self.OutletNodes = []
        self.BifurcationNodes = []
        self.Vessels = []
        self.NumberOfVessels = 0
        self.NumberOfNodes = 0
        self.VesselAtlas = dict()
        self.Clots = []
        self.Graph = []

    def GetVesselNameFromNode(self, node):
        for vessel in self.Vessels:
            if node in vessel.Nodes:
                return vessel.Name
        print("Node not found.")
        return None

    def UpdateNodeType(self):
        [node.SetType(0) for node in self.Nodes]
        for node in self.BifurcationNodes:
            node.Type = 1
        for node in self.OutletNodes:
            node.Type = 2
        for node in self.InletNodes:
            node[0].Type = 2

    def UpdateVesselAtlas(self):
        for vessel in self.Vessels:
            self.VesselAtlas[vessel.Name] = vessel
            self.VesselAtlas[vessel.ID] = vessel

    def UpdateTopology(self):
        self.Nodes = []
        [self.Nodes.append(node) for vessel in self.Vessels for node in vessel.Nodes]
        [self.Nodes.append(node) for node in self.BifurcationNodes]
        # [vessel.SetID(id) for id, vessel in enumerate(self.Vessels)]
        [vessel.UpdateVessel() for vessel in self.Vessels]
        self.UpdateVesselAtlas()
        self.UpdateNumbers()
        self.UpdateNodeType()
        self.NumberNodes()
        self.CheckConnectivity()
        # self.FindOutletNodes()

    def WriteNodesCSV(self, filename):
        with open(filename, "w") as f:
            f.write("ID,VesselID,Connectivity,Position,Radius,LengthAlongVessel,Elasticity,R1,R2,C,Thickness,Type\n")
            for node in self.Nodes:
                f.write("%d," % node.Number)
                f.write("%d," % node.VesselID)
                othernodes = ",".join([str(i.Number) for i in list(node.Connections)])
                f.write("\"%s\"," % othernodes)
                pos = ",".join([str(i) for i in node.Position])
                f.write("\"%s\"," % pos)
                f.write("%f," % node.Radius)
                f.write("%f," % node.LengthAlongVessel)
                f.write("%f," % node.YoungsModules)
                f.write("%f," % node.R1) if not (node.R1 is None) else f.write(",")
                f.write("%f," % node.R2) if not (node.R2 is None) else f.write(",")
                f.write("%f," % node.C) if not (node.C is None) else f.write(",")
                f.write("%f," % node.Thickness)
                f.write("%f," % node.Type)
                f.write("\n")

    def BifurcationDict(self):
        '''
        Create a dictionary of the bifurcation nodes.
        Keys are the nodes of the system and return value is the linked bifurcation node
        :return: dict{node:bifnode}
        '''
        duplicatenodes = {}
        for bif in self.BifurcationNodes:
            for node in list(bif.Connections):
                duplicatenodes[node] = bif
        return duplicatenodes

    def UpdateNumbers(self):
        self.NumberOfVessels = len(self.Vessels)
        self.NumberOfNodes = len(self.Nodes)
        [vessel.UpdateNodeNumber() for vessel in self.Vessels]

    def SetThickness(self):
        for node in self.Nodes:
            node.Thickness = BloodflowEquations.thickness(node.Radius)

    def CheckConnectivity(self):
        print("Checking for disconnected nodes.")
        for index, i in enumerate(self.Nodes):
            if not i.Connections:
                raise ValueError("Disconnected node detected at node number %d" % index)
        print("None Found.")

    def SaveVesselAtlas(self, file="Mapping.csv"):
        if self.VesselAtlas is None:
            print("No Atlas defined.")

        print("Writing the vessel Atlas to file: %s" % file)
        # check if there are nodes without a number.
        # if so, renumber all nodes.
        for node in self.Nodes:
            if node.Number is None:
                self.NumberNodes()
                break

        with open(file, 'w') as f:
            for index, vessel in enumerate(self.Vessels):
                vesselname = vessel.Name
                nodes = [node.Number for node in vessel.Nodes]
                nodelist = ",".join(str(x) for x in nodes)
                f.write("\""+vesselname + "\"," + nodelist + "\n")

    def LoadSegmentedVessels_old(self, inputfolder):
        # temp solution

        finalfile = inputfolder + "1-D_Anatomy_Patient.txt"
        readfile = inputfolder + "Feature_Vessel.csv"
        readfile2 = inputfolder + "Image_info.txt"
        scaling = float([i.strip('\n').split(' ') for i in open(readfile2)][0][2])
        file = inputfolder + "1-D_Anatomy.txt"

        vesselscsv = []
        with open(readfile) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=',')
            for row in readCSV:
                vesselscsv.append(row)

        vessels = []
        for vessel in vesselscsv[1:]:
            vesseldata = namedtuple('vessel', 'ID Length Radius')
            vesseldata.ID = int(vessel[10])
            vesseldata.Length = float(vessel[12]) * scaling
            vesseldata.Radius = float(vessel[13]) * scaling
            vessels.append(vesseldata)

        anatomy = [i.strip('\n').split('\t') for i in open(file)]

        # mapping between anatomy and patient segmented vessels
        # VesselID segmentation, name from segmentation, ID in 1-D anatomy
        MapSegmentation = [(1, "ICA L", 18),
                           (2, "ICA R", 21),
                           (3, "M1 L", 23),
                           (4, "M1 R", 24),
                           (5, "M2 L", -1),
                           (6, "M2 R", -1),
                           (7, "A1 L", 25),
                           (8, "A1 R", 26),
                           (9, "A2 L", 29),
                           (10, "A2 R", 30),
                           (11, "AComm", 31),
                           (12, "M3 L", -1),
                           (13, "M3 R", -1),
                           (14, "VA L", -1),  # 17
                           (15, "VA R", -1),  # 14
                           (16, "BA", 22),
                           (17, "P1 L", 27),
                           (18, "P1 R", 28),
                           (19, "P2 L", 32),
                           (20, "P2 R", 33),
                           (21, "PComm L", 19),
                           (22, "PComm R", 20),
                           (23, "OA L", -1),
                           (24, "OA R", -1), ]

        for vessel in MapSegmentation:
            if vessel[2] != -1:
                anatomy[vessel[2]][2] = str(0)

        for vessel in vessels:
            mapped = MapSegmentation[vessel.ID - 1]
            if mapped[2] != -1:
                anatomy[mapped[2]][2] = str(vessel.Length)
                anatomy[mapped[2]][3] = str(vessel.Radius)
                anatomy[mapped[2]][4] = str(vessel.Radius)

        with open(finalfile, 'w') as file:
            for line in anatomy:
                file.write("\t".join(line) + "\n")

    def LoadSegmentedVessels(self, inputfolder):
        # temp solution

        finalfile = inputfolder + "1-D_Anatomy_Patient.txt"
        readfile = inputfolder + "Feature_Vessel.csv"
        readfile2 = inputfolder + "Image_info.txt"
        scaling = float([i.strip('\n').split(' ') for i in open(readfile2)][0][2])

        file = inputfolder + "1-D_Anatomy.txt"
        anatomy = [i.strip('\n').split('\t') for i in open(file)]

        vesselscsv = []
        with open(readfile) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=',')
            for row in readCSV:
                vesselscsv.append(row)

        vessels = []
        for vessel in vesselscsv[1:]:
            vesseldata = namedtuple('vessel', 'ID Length Radius Startmark Endmark vesselnumber StartPos EndPos')
            vesseldata.ID = int(vessel[10])
            vesseldata.StartPos = (float(vessel[4]), float(vessel[5]), float(vessel[6]))
            vesseldata.EndPos = (float(vessel[7]), float(vessel[8]), float(vessel[9]))
            vesseldata.Length = float(vessel[12]) * scaling
            vesseldata.Radius = float(vessel[13]) * scaling
            vesseldata.Startmark = int(vessel[2])
            vesseldata.Endmark = int(vessel[3])
            vessels.append(vesseldata)

        # mapping between anatomy and patient segmented vessels
        # VesselID segmentation, name from segmentation, ID in 1-D anatomy
        MapSegmentation = [(1, "ICA L", -1),
                           (2, "ICA R", -1),
                           (3, "M1 L", 21),
                           (4, "M1 R", 22),
                           (5, "M2 L", -1),  # not in default anatomy
                           (6, "M2 R", -1),  # not in default anatomy
                           (7, "A1 L", 23),
                           (8, "A1 R", 24),
                           (9, "A2 L", 27),
                           (10, "A2 R", 28),
                           (11, "AComm", 29),
                           (12, "M3 L", -1),  # not in default anatomy
                           (13, "M3 R", -1),  # not in default anatomy
                           (14, "VA L", -1),  # 17  # not part of the scan region but is not missing either.
                           (15, "VA R", -1),  # 14 # not part of the scan region but is not missing either.
                           (16, "BA", -1),  ## split into 5!
                           (17, "P1 L", 25),
                           (18, "P1 R", 26),
                           (19, "P2 L", 30),
                           (20, "P2 R", 31),
                           (21, "PComm L", 18),
                           (22, "PComm R", 19),
                           (23, "OA L", -1),  # not part of the scan region but is not missing either.
                           (24, "OA R", -1), ]  # not part of the scan region but is not missing either.

        for vessel in vessels:
            vessel.vesselnumber = MapSegmentation[vessel.ID - 1][2]

        # set missing vessels to zero
        # some vessels are excluded as they have to be present but are not in the scans.
        # example, aortic arch (23,24)
        for vessel in MapSegmentation:
            if not (vessel[2] == -1):
                anatomy[vessel[2]][2] = str(0)

        # overwrite default anatomy
        for vessel in vessels:
            mapped = MapSegmentation[vessel.ID - 1]
            if not (mapped[2] == -1):
                anatomy[mapped[2]][2] = str(vessel.Length)
                anatomy[mapped[2]][3] = str(vessel.Radius)
                anatomy[mapped[2]][4] = str(vessel.Radius)

        # remove the inlet vessels as these are cut off during imaging (incorrect length)
        vessels = [vessel for vessel in vessels if not (vessel.ID in [1, 2, 16])]

        # add vessels not present in default anatomy
        # "M2 L"  # not in default anatomy
        # "M2 R"  # not in default anatomy
        # "M3 L"  # not in default anatomy
        # "M3 R"  # not in default anatomy
        # assumption is that these always come in pairs
        extravessels = []
        newvessels = []
        for vessel in vessels:
            mapped = MapSegmentation[vessel.ID - 1]
            if mapped[2] == -1:
                # print(vessel)
                name = mapped[1]
                if name in extravessels:
                    name = name + ",2"
                else:
                    extravessels.append(name)
                    name = name + ",1"
                length = str(vessel.Length)
                radius = str(vessel.Radius)
                youngsmodulus = str(1.6)
                newves = [['id', name, length, radius, radius, youngsmodulus], vessel]
                newvessels.append(newves)

        # add vessels to the anatomy
        numberoldvessels = anatomy.index(['Bifurcations']) - 1

        for index, vessel in enumerate(newvessels):
            vessel[0][0] = str(numberoldvessels + index)
            vessel[1].vesselnumber = numberoldvessels + index
            anatomy.insert(numberoldvessels + index, vessel[0])

        # extract ends
        endspos = [vessel.StartPos for vessel in vessels] + [vessel.EndPos for vessel in vessels]
        ends = [end for end in list(set(endspos))]
        marks = [[list(end), [], []] for end in ends]
        for vessel in vessels:
            for index, end in enumerate(ends):
                if vessel.StartPos == end:
                    marks[index][1].append(vessel)
                if vessel.EndPos == end:
                    marks[index][2].append(vessel)

        # extract bifurcations
        bifurcations = [mark for mark in marks if len(mark[1]) + len(mark[2]) > 1]

        # assign numbers
        for bif in bifurcations:
            numbers = set()
            for vessel in bif[1]:
                numbers.add(vessel.Startmark)
            for vessel in bif[2]:
                numbers.add(vessel.Endmark)
            if len(numbers) > 1:
                print("Error in segmentation. Bifurcation has more than one ID.")
            bif.append(list(numbers)[0])

        # list the new bifurcations in the 1-D anatomy format
        possiblenewvesselsbif = [7, 8, 13, 14, 25, 26]
        candidatebif = [bif for bif in bifurcations if bif[3] in possiblenewvesselsbif]

        # we assume a certain direction, the segmentation do not seem to follow the same.
        # M1 R, 4->8, 4
        # M2 R, 8->14, 6
        # M3 R, 14->26, 13
        # M1 L, 3->7, 3
        # M2 L, 7->13, 5
        # M3 L, 13->25, 12
        bifends = {7: (3, 5),
                   8: (4, 6),
                   13: (5, 12),
                   14: (6, 13),
                   25: (12,),
                   26: (13,)}

        newbifs = []
        inout = ("o", "i")
        for index, bif in enumerate(candidatebif):
            # ends = [str(i.vesselnumber) + "i" for i in bif[0]]
            # ends += [str(i.vesselnumber) + "o" for i in bif[1]]
            ends = []
            endlist = bifends[bif[3]]
            for vessel in bif[1]:
                index = endlist.index(vessel.ID)
                ends.append(str(vessel.vesselnumber) + inout[index])
            for vessel in bif[2]:
                index = endlist.index(vessel.ID)
                ends.append(str(vessel.vesselnumber) + inout[index])
            newbifs.append(" ".join(ends))

        endsbiflist = anatomy.index(['Boundary Nodes']) - 1
        for index, newbif in enumerate(newbifs):
            anatomy.insert(endsbiflist + index, [newbif])

        ends = ["Outlets Brain:"]
        nodesinbif = {8, 10, 13, 15, 16}  # only possible extra ones atm
        startbiflist = anatomy.index(['Bifurcations'])
        endsbiflist = anatomy.index(['Boundary Nodes']) - 1
        for bif in anatomy[startbiflist:endsbiflist]:
            nodes = bif[0].split(" ")
            [nodesinbif.add(int(node[:-1])) for node in nodes if node[-1] == "o"]

        for i in range(1, startbiflist - 1):
            if i not in nodesinbif:
                ends.append(str(i) + "o")
        anatomy[-1] = [" ".join(ends)]

        with open(finalfile, 'w') as file:
            for line in anatomy:
                file.write("\t".join(line) + "\n")

    def LoadBFSimFiles(self, folder):
        self.LoadTopFile(folder + "System.top")
        self.LoadParFile(folder + "System.par")
        self.LoadRunFile(folder + "Run.txt")
        self.LoadVesselAtlas(folder + "Mapping.csv")
        self.LoadClotFile(folder + "Clots.txt")

    def LoadTopFile(self, file):
        print("Loading topology file.")
        self.ReadNodesFromTopFile(file)
        self.CheckConnectivity()

    def LoadParFile(self, file):
        print("Loading parameter file.")
        pardata = [i.strip('\n').split(' ') for i in open(file)]
        for line in pardata:
            nodenumber = int(line[0])
            node = self.Nodes[nodenumber]
            r1 = float(line[1][3:]) * 1e9
            r2 = float(line[2][3:]) * 1e9
            c = float(line[3][3:]) * 1e-12
            node.SetWK(r1, r2, c)

    def LoadRunFile(self, file):
        print("Loading run file.")
        rundata = [i.strip('\n') for i in open(file)]
        inletdata = rundata[2][11:].split(",")
        for inlet in inletdata:
            nodenumber, inletfile = inlet.split(" ")
            node = self.Nodes[int(nodenumber)]
            self.InletNodes.append((node, inletfile))
            self.OutletNodes.remove(node)

    def LoadClotFile(self, file):
        print("Loading clot file.")
        clotdata = [i.strip('\n').split(",") for i in open(file)][1:]
        clotdata = [[int(i[0]), int(i[1]), float(i[2]), float(i[3])] for i in clotdata]
        Clotids = list(set([i[0] for i in clotdata]))
        for clotid in Clotids:
            clotnodes = [i for i in clotdata if i[0] == clotid]
            nodes = [self.Nodes[i[1]] for i in clotnodes]
            par1 = clotnodes[0][2]
            par2 = clotnodes[0][3]
            self.Clots.append((nodes, par1, par2))

    def LoadVTPFile(self, vtpfile):
        print("Loading %s" % vtpfile)
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(vtpfile)
        reader.Update()

        pos_vtk = reader.GetOutput().GetPoints().GetData()
        pos = vtk_to_numpy(pos_vtk)

        VesselAtlas = None
        vesselids = None
        narray = reader.GetOutput().GetPointData().GetNumberOfArrays()
        for i in range(0, narray):
            arrayname = reader.GetOutput().GetPointData().GetArrayName(i)
            if arrayname == "MaximumInscribedSphereRadius" or arrayname == "Radius":
                radius_vtk = reader.GetOutput().GetPointData().GetArray(i)
            if arrayname == "Type":
                #  load in the vessel atlas
                VesselAtlas = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(i))

        narray = reader.GetOutput().GetCellData().GetNumberOfArrays()
        for i in range(0, narray):
            arrayname = reader.GetOutput().GetCellData().GetArrayName(i)
            if arrayname == "Vessel Ids":
                #  load in the vessel atlas
                vesselids = vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(i))
            if arrayname == "Major Vessel Ids":
                #  load in the vessel atlas
                VesselAtlas = vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(i))

        radius = vtk_to_numpy(radius_vtk)

        for i in range(0, len(pos)):
            node = Node.Node()
            node.Number = int(i)
            node.SetRadius(radius[i])
            node.SetPosition(list(pos[i]))
            self.Nodes.append(node)

        NumberOfVessels = reader.GetNumberOfCells()
        vesselnodes = []
        for i in range(0, NumberOfVessels):
            numberofnodes = reader.GetOutput().GetCell(i).GetNumberOfPoints()
            ids = reader.GetOutput().GetCell(i).GetPointIds()
            vesselnodes.append([self.Nodes[ids.GetId(i)] for i in range(0, numberofnodes)])
            [node.SetVesselID(i) for node in vesselnodes[-1]]
            # these include the bifurcations

        # Connections
        for vessel in vesselnodes:
            for i in range(0, len(vessel) - 1):
                currentnode = vessel[i]
                neighbour = vessel[i + 1]
                currentnode.AddConnection(neighbour)
                neighbour.AddConnection(currentnode)

        if not (VesselAtlas is None):
            for index, id in enumerate(VesselAtlas):
                vessel = vesselnodes[index]
                for node in vessel:
                    node.MajorVesselID = id

        if not (vesselids is None):
            for index, id in enumerate(vesselids):
                vessel = vesselnodes[index]
                for node in vessel:
                    node.VesselID = id

        self.FindBifurcationNodes()

        # extract bifurcations from the vessels
        # redefine nodes
        bifurcations = [Node.Node() for i in range(0, len(self.BifurcationNodes))]
        [bifurcations[i].SetPosition(self.BifurcationNodes[i].Position) for i in range(0, len(self.BifurcationNodes))]
        [bif.SetLengthAlongVessel(0) for bif in bifurcations]
        self.Vessels = [Vessel.Vessel() for vessel in vesselnodes]
        for index, vessel in enumerate(vesselnodes):
            self.Vessels[index].Nodes = vessel

        [vessel.InterpolateVessel3Dto1D() for vessel in self.Vessels]

        for number, node in enumerate(self.BifurcationNodes):
            for index, item in enumerate(vesselnodes):
                if node == item[0]:
                    bifurcations[number].AddConnection(self.Vessels[index].Nodes[0])
                    self.Vessels[index].Nodes[0].AddConnection(bifurcations[number])
                if node == item[-1]:
                    bifurcations[number].AddConnection(self.Vessels[index].Nodes[-1])
                    self.Vessels[index].Nodes[-1].AddConnection(bifurcations[number])

        [bifurcation.SetRadius(max([node.Radius for node in bifurcation.Connections])) for bifurcation in bifurcations]
        self.BifurcationNodes = bifurcations
        self.Nodes = []
        [self.Nodes.append(node) for vessel in self.Vessels for node in vessel.Nodes]
        [self.Nodes.append(node) for node in self.BifurcationNodes]
        self.InletNodes = {}
        self.OutletNodes = []
        self.FindlargestInlets()
        self.CheckConnectivity()
        self.UpdateNodeType()

    def Load1DAnatomy(self, file):
        print("Loading anatomy from %s" % file)
        bfdata = [i.strip('\n').split('\t') for i in open(file)]
        BifurcationsIndex = [x[0] for x in bfdata].index("Bifurcations")
        BCIndex = [x[0] for x in bfdata].index("Boundary Nodes")

        vessels = [i for i in bfdata[1:BifurcationsIndex - 1]]
        bifurcations = [i for i in bfdata[BifurcationsIndex + 1:BCIndex - 1]]
        BC = [i for i in bfdata[BCIndex + 1:]]

        # remove zero length vessels from further processing
        NoVessel = [i for i in vessels if float(i[2]) == 0]

        self.Vessels = [Vessel.Vessel() for element in vessels]
        [self.Vessels[index].SetID(int(element[0])) for index, element in enumerate(vessels)]
        [self.Vessels[index].SetName(element[1]) for index, element in enumerate(vessels)]
        [self.Vessels[index].GenerateVessel(float(element[2]), float(element[3]), float(element[4]),
                                            float(element[5]) * 1e6) for index, element in enumerate(vessels)]
        self.UpdateVesselAtlas()
        self.UpdateNumbers()

        def bifurcationnode(info):
            bifurcationnode = Node.Node()
            bifurcationnode.SetLengthAlongVessel(0)

            nodedata = info[0].split(' ')
            bifurcationinfo = [[int(element[:-1]), element[-1]] for element in nodedata]
            if len(NoVessel) >= 1:
                NoIds = [int(i[0]) for i in NoVessel]
                bifurcationinfo = [node for node in bifurcationinfo if node[0] not in NoIds]

            for node in bifurcationinfo:
                vessel = self.Vessels[node[0] - 1].Nodes  # note that vessel ids start at 1
                if node[1] == "o":
                    link = vessel[-1]
                    bifurcationnode.AddConnection(link)
                    link.AddConnection(bifurcationnode)
                elif node[1] == "i":
                    link = vessel[0]
                    bifurcationnode.AddConnection(link)
                    link.AddConnection(bifurcationnode)
                else:
                    print("Error")

            radius = max([node.Radius for node in bifurcationnode.Connections])
            bifurcationnode.SetRadius(radius)
            return bifurcationnode

        self.BifurcationNodes = [bifurcationnode(element) for element in bifurcations]

        def bcnodes(info):
            inlets = []
            outlets = []

            inletdata = info[0][0].split(' ')
            inletinfo = [[int(element[:-1]), element[-1]] for element in inletdata[1:]]
            if len(NoVessel) >= 1:
                NoIds = [int(i[0]) for i in NoVessel]
                inletinfo = [node for node in inletinfo if node[0] not in NoIds]

            for node in inletinfo:
                vessel = self.Vessels[node[0] - 1].Nodes  # note that vessel ids start at 1
                if node[1] == "o":
                    inlets.append(vessel[-1])
                elif node[1] == "i":
                    inlets.append(vessel[-0])
                else:
                    print("Error")

            for i in range(1, 4):
                outletdata = info[i][0].split(' ')
                outletinfo = [[int(element[:-1]), element[-1]] for element in outletdata[2:]]
                if len(NoVessel) >= 1:
                    NoIds = [int(i[0]) for i in NoVessel]
                    outletinfo = [node for node in outletinfo if node[0] not in NoIds]

                for node in outletinfo:
                    vessel = self.Vessels[node[0] - 1].Nodes  # note that vessel ids start at 1
                    if node[1] == "o":
                        outlets.append(vessel[-1])
                    elif node[1] == "i":
                        outlets.append(vessel[-0])
                    else:
                        print("Error")
            return inlets, outlets

        inletnodes, self.OutletNodes = bcnodes(BC)
        self.InletNodes = [(inlet, "Aorta.txt") for inlet in inletnodes]

        self.Vessels = [vessel for vessel in self.Vessels if not (vessel.Nodes == [])]

        [self.Nodes.append(node) for vessel in self.Vessels for node in vessel.Nodes]
        [self.Nodes.append(node) for node in self.BifurcationNodes]
        self.NumberNodes()
        self.SetThickness()
        self.NumberOfNodes = len(self.Nodes)
        self.NumberOfVessels = len(self.Vessels)
        self.CheckConnectivity()
        self.UpdateNodeType()
        [vessel.UpdateNodeVesselID() for vessel in self.Vessels]
        [vessel.CalculateMeanThickness() for vessel in self.Vessels]
        [vessel.CalculateMeanRadius() for vessel in self.Vessels]
        self.UpdateTopology()

    def NumberNodes(self):
        print("Assigning numbers to the nodes.")
        for node in self.Nodes:
            node.Number = None
        number = 0
        for node in self.Nodes:
            if node.Number is None:
                node.Number = number
                number += 1

    def MapNodesToVessel(self):
        nodedict = dict()
        for i in range(0, len(self.Vessels)):
            [nodedict.setdefault(self.Vessels[i].MajorVesselID, []).append(node) for node in self.Vessels[i].Nodes]
        return nodedict

    def GetDownstreamVessels(self, inputvessel):
        """Get upstream vessels and bifurcation nodes.
        Note that this follows the positive directon of the vessels.
        If vessels are excluded that should be included, check the direction of those vessels.

        The code finds the relevent bifucations and then checks which vessels are connected.


        """
        # print("Finding the upstream vessels.")
        # initiation
        bifset = []
        bifurcations = set()
        if isinstance(inputvessel, str):
            vesselindex = self.VesselAtlas[inputvessel]
            vessel = vesselindex
        else:
            vessel = inputvessel
        # the current vessel and the starting bifurcation
        vesselend = vessel.Nodes[-1]
        for node in vesselend.Connections:
            if node in self.BifurcationNodes:
                bifset = [node]
                bifurcations.add(node)

        gens = []
        gens.append([vessel])
        alreadylisted = [vessel]
        # Continue until we have all upstream vessels
        while len(bifset) > 0:
            connectedvessels = []
            for bifurcationnode in bifset:
                for node in bifurcationnode.Connections:
                    connectedvessels += [vessel for vessel in self.Vessels if node in vessel.Nodes]
            nbifset = set()
            for vessel in connectedvessels:
                ves = vessel.Nodes
                for node in ves[-1].Connections:
                    if node in self.BifurcationNodes and node not in bifurcations:
                        nbifset.add(node)
                        bifurcations.add(node)
            newgen = [vessel for vessel in connectedvessels if vessel not in alreadylisted]
            gens.append(newgen)
            alreadylisted += connectedvessels
            bifset = list(nbifset)

        # only vessels that are connected with the first node should be included.
        # this ensures that we only get the upstream vessels
        vessels = set()
        vessels.add(vessel)
        for bif in bifurcations:
            for node in bif.Connections:
                [vessels.add(vessel) for vessel in self.Vessels if node in vessel.Nodes and vessel not in vessels]

        # export the names of the upstream vessels
        vesselnames = []
        for index, vessel in enumerate(self.Vessels):
            for match in vessels:
                if match == vessel:
                    vesselnames.append(vessel.Name)

        # reorder
        sortedvessels = sorted(list(vessels), key=lambda v: self.Vessels.index(v))

        return sortedvessels, list(bifurcations), vesselnames, gens

    def ReorderGens(self, gens):
        ordering = []
        for generation in gens:
            lengthgen = []
            for vessel in generation:
                v, b, vn, ng = self.GetDownstreamVessels(vessel)
                lengthgen.append(len(ng))
            ordering.append(lengthgen)

        newgens = []
        for index, gen in enumerate(gens):
            newgens.append(sorted(gen, key=lambda g: ordering[index][gen.index(g)]))
        return newgens

    def DownStreamResistance(self, inputvessel, visc, density):
        with contextlib.redirect_stdout(None):
            vessels, bifurcations, vesselnames, gens = self.GetDownstreamVessels(inputvessel)

        C = gens[0][0].Nodes[-1].C
        R = gens[0][0].Nodes[-1].R1 + gens[0][0].Nodes[-1].R2
        gens[0][0].Nodes[-1].ResetWK()

        vesselres = [vessel.VesselResistance(visc) for vessel in vessels]
        vesselres[vessels.index(gens[0][0])] = 0  # only get the downstream resistance

        vesselcomp = [vessel.VesselCompliance() for vessel in vessels]
        vesselcomp[vessels.index(gens[0][0])] = 0  # only get the downstream resistance

        for i in range(len(gens) - 1, -1, -1):
            for vessel in gens[i]:
                for index, vesselmatch in enumerate(self.Vessels):
                    if vesselmatch == vessel:
                        vesselnames = vesselmatch.Name
                with contextlib.redirect_stdout(None):
                    vesselstemp, bifurcationstemp, vesselnamestemp, genstemp = self.GetDownstreamVessels(vesselnames)
                if len(genstemp) > 1:
                    downstreamres = [vesselres[vessels.index(vesselt)] for vesselt in genstemp[1]]
                    vesselres[vessels.index(vessel)] = self.sumresvessels(vesselres[vessels.index(vessel)],
                                                                          downstreamres)
                    # vesselc = sum([vesselcomp[vessels.index(vesselt)] for vesselt in genstemp[1]])
                    # if vesselcomp[vessels.index(vessel)] > 0:
                    #     vesselcomp[vessels.index(vessel)] = 1/(1/vesselc + 1/vesselcomp[vessels.index(vessel)])
                    # else:
                    #     vesselcomp[vessels.index(vessel)] = vesselc

        outlets = []
        for vessel in vessels:
            for out in self.OutletNodes:
                if out in vessel.Nodes:
                    outlets.append(out)

        rcubed = 0
        scaling = 1e-3
        RT = R - vesselres[vessels.index(gens[0][0])]
        # Ct = 1/(1/C - 1/vesselcomp[vessels.index(gens[0][0])])
        if RT < 0:
            print("Found RT to be negative!?")
            RT = 0.01 * R
        # RT = R
        print("Resistance scaled: %f" % (RT / R))
        # print("Compliance scaled: %f" % (Ct / C))
        for i in outlets:
            radius = i.Radius * scaling
            rcubed += math.pow(radius, 3)

        C1D = 0
        for vessel in vessels:
            # mean radius and thickness
            lengthvessel = vessel.Length * scaling
            meanradius = vessel.MeanRadius * scaling
            meanh = vessel.MeanThickness * scaling
            C1D += 2 * numpy.power(numpy.pi * meanradius * meanradius, 1.5) * lengthvessel / (
                    (4 / 3) * numpy.sqrt(numpy.pi) * vessel.YoungsModules * meanh)

        for i in outlets:
            radius = i.Radius * scaling
            h = i.CalculateThickness() * scaling
            rt = RT * rcubed / math.pow(radius, 3)
            A = math.pi * math.pow(radius, 2)
            beta = (4 / 3) * math.sqrt(math.pi) * i.YoungsModules * h / A
            c0 = abs(math.sqrt(beta / (2 * density))) * math.pow(A, 0.25)
            r1 = (density / A) * c0
            r2 = rt - r1
            if r2 < 0:
                r2 = 0.1e9
                r1 = rt - r2
            c = max(0, (C - C1D) * RT / rt)
            i.SetWK(r1, r2, c)

    # def UpstreamResistancesyms(self, inputvessel, visc, density):
    #     rt = sympy.symbols("rt")
    #     vessels,  bifurcations, vesselnames, gens = self.GetUpstreamVessels(inputvessel)
    #
    #
    #     C = gens[0][0][-1].C
    #     R = gens[0][0][-1].R1 + gens[0][0][-1].R2
    #     gens[0][0][-1].ResetWK()
    #
    #     vesselres = [VesselResistance(vessel, visc) for vessel in vessels]
    #     # vesselres[vessels.index(gens[0][0])] = 0 #assumption that the major vessels do not have meaningful resistance.
    #
    #     # add resistance of the outlets
    #     outlets = []
    #     for vessel in vessels:
    #         for out in self.OutletNodes:
    #             if out in vessel:
    #                 outlets.append(out)
    #     rcubed = 0
    #     scaling = 1e-3
    #     for i in outlets:
    #         radius = i.Radius * scaling
    #         rcubed += math.pow(radius, 3)
    #
    #     for vessel in vessels:
    #         for out in self.OutletNodes:
    #             if out in vessel:
    #                 vesselres[vessels.index(vessel)] += rt*rcubed / math.pow(out.Radius*scaling, 3)
    #
    #     for i in range(len(gens)-1,-1,-1):
    #         for vessel in gens[i]:
    #             for index, vesselmatch in enumerate(self.Vessels):
    #                 if vesselmatch == vessel:
    #                     vesselnames = self.VesselAtlas[index]
    #             with contextlib.redirect_stdout(None):
    #                 vesselstemp, bifurcationstemp, vesselnamestemp, genstemp = self.GetUpstreamVessels(vesselnames)
    #             if len(genstemp)>1:
    #                 downstreamres = [vesselres[vessels.index(vesselt)] for vesselt in genstemp[1]]
    #                 vesselres[vessels.index(vessel)] = self.sumresvessels(vesselres[vessels.index(vessel)],downstreamres)
    #
    #     expression =sympy.lambdify(rt, R - vesselres[vessels.index(gens[0][0])],'numpy')
    #     RT = float(fsolve(expression,expression(0)))
    #     # print((RT)/R)
    #
    #     C1D = 0
    #     for vessel in vessels:
    #         # mean radius and thickness
    #         lengthvessel = max([vessel[0].Position, vessel[-1].Position]) * scaling
    #         meanradius = numpy.trapz([i.Radius * scaling for i in vessel],
    #                                  [i.Position * scaling for i in vessel]) / lengthvessel
    #         meanh = numpy.trapz([thickness(i.Radius * scaling) for i in vessel],
    #                                  [i.Position * scaling for i in vessel]) / lengthvessel
    #         C1D += 2 * numpy.power(numpy.pi * meanradius * meanradius, 1.5) * lengthvessel / (
    #                 (4 / 3) * numpy.sqrt(numpy.pi) * vessel[0].E * meanh)
    #
    #     for i in outlets:
    #         radius = i.Radius * scaling
    #         h = thickness(radius)
    #         rt = RT * rcubed / math.pow(radius, 3)
    #         A = math.pi * math.pow(radius, 2)
    #         beta = (4 / 3) * math.sqrt(math.pi) * i.E * h / A
    #         c0 = abs(math.sqrt(beta / (2 * density))) * math.pow(A, 0.25)
    #         r1 = (density / A) * c0
    #         r2 = rt - r1
    #         if r2 < 0:
    #             r2 = 0.1e9
    #             r1 = rt - r2
    #         c = (C-C1D) * RT / rt
    #         i.set_WK(r1, r2, c)

    def sumresvessels(self, upstream, downstream):
        parallelvesselsres = 1 / sum([1 / vesselres for vesselres in downstream])
        serialvesselres = upstream
        totalres = parallelvesselsres + serialvesselres
        return totalres

    def RedefineDirection(self):
        print("Redefining the direction of the vessels.")
        direction = [0 for vessel in self.Vessels]

        for inlet in self.InletNodes:
            for index, vessel in enumerate(self.Vessels):
                if inlet[0] is vessel.Nodes[0]:
                    direction[index] = 1
                    vessel.GenerationNumber = 0
                elif inlet[0] is vessel.Nodes[-1]:
                    print("Vessel reversed.")
                    self.Vessels[index].Nodes = self.Vessels[index].Nodes[::-1]
                    for nodenumber in range(0, len(self.Vessels[index].Nodes)):
                        if nodenumber < len(self.Vessels[index].Nodes) - nodenumber:
                            # self.Vessels[index].Nodes[nodenumber].Position, self.Vessels[index].Nodes[
                            #     -1 - nodenumber].Position = \
                            #     self.Vessels[index].Nodes[-1 - nodenumber].Position, self.Vessels[index].Nodes[
                            #         nodenumber].Position
                            self.Vessels[index].Nodes[nodenumber].LengthAlongVessel, self.Vessels[index].Nodes[
                                -1 - nodenumber].LengthAlongVessel = \
                                self.Vessels[index].Nodes[-1 - nodenumber].LengthAlongVessel, self.Vessels[index].Nodes[
                                    nodenumber].LengthAlongVessel
                    direction[index] = -1

        inletvessels = []
        for inlet in self.InletNodes:
            for index, vessel in enumerate(self.Vessels):
                if inlet[0] is vessel.Nodes[0]:
                    inletvessels.append(vessel.Nodes)

        bifs = set()
        for inves in inletvessels:
            for bif in inves[-1].Connections:
                if bif in self.BifurcationNodes:
                    bifs.add(bif)

        while len(bifs) > 0:
            newbifs = set()
            for bif in bifs:
                vessels = [self.Vessels.index(vessel) for vessel in self.Vessels if
                           vessel.Nodes[0] in bif.Connections or vessel.Nodes[-1] in bif.Connections]
                vesselundef = [vessel for vessel in vessels if direction[vessel] == 0]

                for vessel in vesselundef:
                    # print(self.VesselAtlas[vessel])
                    if self.Vessels[vessel].Nodes[0] in bif.Connections:
                        direction[vessel] = 1
                    else:
                        print("Vessel reversed.")
                        self.Vessels[vessel].Nodes = self.Vessels[vessel].Nodes[::-1]

                        for nodenumber in range(0, len(self.Vessels[vessel].Nodes)):
                            if nodenumber < len(self.Vessels[vessel].Nodes) - nodenumber:
                                # self.Vessels[vessel].Nodes[nodenumber].Position, self.Vessels[vessel].Nodes[
                                #     -1 - nodenumber].Position = \
                                #     self.Vessels[vessel].Nodes[-1 - nodenumber].Position, self.Vessels[vessel].Nodes[
                                #         nodenumber].Position
                                self.Vessels[vessel].Nodes[nodenumber].LengthAlongVessel, self.Vessels[vessel].Nodes[
                                    -1 - nodenumber].LengthAlongVessel = \
                                    self.Vessels[vessel].Nodes[-1 - nodenumber].LengthAlongVessel, \
                                    self.Vessels[vessel].Nodes[
                                        nodenumber].LengthAlongVessel
                        direction[vessel] = -1

                    newcandidates = [bift for bift in self.Vessels[vessel].Nodes[-1].Connections if
                                     bift in self.BifurcationNodes]
                    if newcandidates:
                        newbifs.add(newcandidates[0])
            bifs = newbifs

    def GetVesselsFromType(self, inputtype):
        print("Extracting all vessels of a type.")
        # return all vessels that are labelled with the input type
        vessels = [self.Vessels[i] for i in range(0, len(self.Vessels)) if self.VesselAtlas[i] == inputtype]
        return vessels

    def ReadNodesFromTopFile(self, file):
        print("Load nodes from topology file.")
        bfdata = [i.strip('\n').split(' ') for i in open(file)]
        bonds_index = [x[0] for x in bfdata].index("Bonds:")
        self.Nodes = [Node.Node() for i in bfdata[2:bonds_index - 1]]

        [self.Nodes[index].SetNodeFromTopLine(i) for index, i in enumerate(bfdata[2:bonds_index - 1])]

        links = bfdata[bonds_index + 1:len(bfdata)]
        for link in links:
            [self.Nodes[int(link[0])].AddConnection(self.Nodes[int(link[i])]) for i in
             range(1, len(link))]

        self.BifurcationNodes = [node for node in self.Nodes if node.Type == 1]
        self.OutletNodes = [node for node in self.Nodes if node.Type == 2]

    def Read3DNodesFromTopFile(self, file):
        print("Load nodes from topology file.")
        bfdata = [i.strip('\n').split(' ') for i in open(file)]
        bonds_index = [x[0] for x in bfdata].index("Bonds:")
        self.Nodes = [Node.Node() for i in bfdata[2:bonds_index - 1]]

        for index, i in enumerate(bfdata[2:bonds_index - 1]):
            self.Nodes[index].SetNumber(int(i[0]))
            self.Nodes[index].SetPosition([1e3 * float(i[1][2:]), 1e3 * float(i[2][2:]), 1e3 * float(i[3][2:])])
            self.Nodes[index].SetRadius(float(i[4][2:]))

        links = bfdata[bonds_index + 1:len(bfdata)]
        for link in links:
            [self.Nodes[int(link[0])].AddConnection(self.Nodes[int(link[i])]) for i in
             range(1, len(link) - 1)]

        self.BifurcationNodes = [node for node in self.Nodes if node.Type == 1]
        self.OutletNodes = [node for node in self.Nodes if node.Type == 2]

    def LoadVesselAtlas(self, file):
        print("Loading vessel atlas from file.")

        if file[-3:] == "csv":
            with open(file, 'r') as csvfile:
                reader = csv.reader(csvfile, delimiter=',')
                vesselatlas = [i for i in reader]
        else:
            with open(file, 'r') as csvfile:
                reader = csv.reader(csvfile, delimiter='\t')
                vessel = [i for i in reader]
                vesselatlas = [i[1].split(" ") for i in vessel]
                vesselnames = [i[0] for i in vessel]
                [vesselatlas[index].insert(0, name) for index, name in enumerate(vesselnames)]

        # vesselatlas = [i.strip('\n').split('\t') for i in open(file)]
        self.Vessels = [Vessel.Vessel() for i in vesselatlas]

        for index, line in enumerate(vesselatlas):
            self.Vessels[index].SetName(line[0])
            self.Vessels[index].SetID(index)
            nodenumbers = [int(i) for i in line[1:]]
            nodes = [self.Nodes[i] for i in nodenumbers]
            self.Vessels[index].SetLength(nodes[-1].LengthAlongVessel)
            self.Vessels[index].SetGridSize(nodes[-1].LengthAlongVessel - nodes[-2].LengthAlongVessel)
            self.Vessels[index].SetNodes(nodes)
            self.Vessels[index].UpdateVessel()
            self.Vessels[index].CalculateMeanRadius()
            self.Vessels[index].CalculateMeanThickness()

        self.UpdateVesselAtlas()

    def SetInletNodes(self, inlets):
        print("Setting inlet nodes.")
        self.InletNodes = inlets

    def SetOutletNodes(self, outlets):
        print("Setting outlet nodes.")
        self.OutletNodes = outlets

    def SetBifurcationNodes(self, bifurcationnodes=None):
        print("Setting bifurcation nodes.")
        if bifurcationnodes is None:
            self.FindBifurcationNodes()
        else:
            self.BifurcationNodes = bifurcationnodes

    def FindBifurcationNodes(self):
        print("Finding the bifurcation nodes.")
        bifurcationnodes = []
        for nodes in self.Nodes:
            if len(nodes.Connections) > 2:
                bifurcationnodes.append(nodes)
        self.BifurcationNodes = bifurcationnodes

    def FindOutletNodes(self):
        print("Finding the outlet nodes.")
        InletNodes = [i[0] for i in self.InletNodes]
        OutletsNodes = [i for i in self.Nodes if len(i.Connections) == 1 and i not in InletNodes]
        self.SetOutletNodes(OutletsNodes)

    def FindlargestInlets(self, InletNumber=3):
        print("Finding the largest %d inlets." % InletNumber)
        # InletNodes = [i.strip('\n').split(',') for i in open(inletfilename)]
        Ends = [i for i in self.Nodes if len(i.Connections) == 1]

        InletNodes = sorted(range(len(Ends)), key=lambda k: Ends[k].Radius, reverse=True)
        inletnames = ["left_carotid_inlet.txt", "right_carotid_inlet.txt", "basilar.txt"]
        Inlets = [(Ends[line], inletnames.pop(0)) for line in InletNodes[0:InletNumber]]  # Dictionary
        self.SetInletNodes(Inlets)
        self.FindOutletNodes()

    def ApplyTransformation(self, matrix):
        for node in self.Nodes:
            pos = node.Position
            vec = numpy.array([[pos[0]], [pos[1]], [pos[2]], [1]])
            position = numpy.dot(matrix, vec)
            posnew = [position[0][0], position[1][0], position[2][0]]
            node.SetPosition(posnew)

    def AnatomyToVessels(self):
        """
        Get vessels from topology data
        start at an end point
        get the neighbour, if connection of neighbour is 2 then get next neighbour and repeat.
        If connection is 1, stop. If connection is more than 2, stop and start new vessel.
        repeat for every end point, keep track of nodes that are already assigned to a vesse (except bif nodes)
        bifurcation nodes are nodes with more than 2 connections
        """

        print("Extracting vessels from topology data.")

        if not self.BifurcationNodes:
            self.FindBifurcationNodes()

        if not self.OutletNodes:
            self.FindOutletNodes()

        ProcessedNodes = []
        vessels = []
        # Note: gives an error is there are no vessel nodes
        for BifurcationNode in self.BifurcationNodes:
            # check connected vessels
            newvessels = list(set(BifurcationNode.Connections) - set(ProcessedNodes))
            for node in newvessels:
                vessel = []
                vessel.append(BifurcationNode)
                ProcessedNodes.append(BifurcationNode)
                currentnode = node
                while 1:
                    numberconnection = len(currentnode.Connections)
                    if numberconnection == 2:
                        vessel.append(currentnode)
                        ProcessedNodes.append(currentnode)
                        currentnode = list(set(currentnode.Connections) - set(vessel))[0]
                    elif numberconnection == 1:
                        vessel.append(currentnode)
                        ProcessedNodes.append(currentnode)
                        vessels.append(vessel)
                        break
                    else:
                        vessel.append(currentnode)
                        ProcessedNodes.append(currentnode)
                        vessels.append(vessel)
                        break
        self.Vessels = [Vessel.Vessel() for i in range(0, len(vessels))]
        [self.Vessels[index].SetNodes(elements) for index, elements in enumerate(vessels)]
        [self.Vessels[index].SetMajorVesselID(vessel.Nodes[1].MajorVesselID) for index, vessel in
         enumerate(self.Vessels)]
        self.UpdateNodeType()

    def removebifurcation(self):
        for bif in self.BifurcationNodes:
            nodes = list(bif.Connections)
            linked = []
            for node in nodes:
                for vessel in self.Vessels:
                    if node in vessel.Nodes:
                        vessel.Nodes[vessel.Nodes.index(node)] = bif
                othernodes = list(node.Connections)
                [linked.append(n) for n in othernodes if not (n is bif)]
                for lastmostvesselend in othernodes:
                    lastmostvesselend.RemoveConnection(node)
                self.Nodes.remove(node)
            bif.Connections = set()
            for link in linked:
                bif.AddConnection(link)
                link.AddConnection(bif)
        self.BifurcationNodes = []
        self.UpdateTopology()

    def TopologyToVTP(self, filename):
        """
        Write the network to a .vtp file.
        This only works with 3-D coordinates.
        """
        if not self.Nodes[0].Number:
            self.NumberNodes()

        print("Writing topology to %s." % filename)
        nodes = vtk.vtkPoints()
        vessels = vtk.vtkCellArray()
        radius = vtk.vtkFloatArray()

        radius.SetNumberOfComponents(1)
        radius.SetName("Radius")

        nodetype = vtk.vtkIntArray()
        nodetype.SetName("Node Type ")

        # Add radius and position to data array
        for i in self.Nodes:
            nodes.InsertNextPoint(i.Position)
            radius.InsertNextValue(float(i.Radius))
            nodetype.InsertNextValue(i.Type)

        # Add vessels to cell array
        for vessel in self.Vessels:
            line = vtk.vtkLine()
            line.GetPointIds().SetNumberOfIds(len(vessel.Nodes))
            for i in range(0, len(vessel.Nodes)):
                line.GetPointIds().SetId(i, vessel.Nodes[i].Number)
            vessels.InsertNextCell(line)

        # Create a polydata to store everything in
        VesselsPolyData = vtk.vtkPolyData()

        # Add the nodes to the polydata
        VesselsPolyData.SetPoints(nodes)

        # Add the vessels to the polydata
        VesselsPolyData.SetLines(vessels)

        # Assign radii to the nodes
        VesselsPolyData.GetPointData().SetScalars(radius)
        VesselsPolyData.GetPointData().AddArray(nodetype)

        atlas = vtk.vtkIntArray()
        [atlas.InsertNextValue(i.ID) for index, i in enumerate(self.Vessels)]
        atlas.SetName("Vessel Ids")
        VesselsPolyData.GetCellData().AddArray(atlas)

        majoratlas = vtk.vtkIntArray()
        [majoratlas.InsertNextValue(i.MajorVesselID) for index, i in enumerate(self.Vessels)]
        majoratlas.SetName("Major Vessel Ids")
        VesselsPolyData.GetCellData().AddArray(majoratlas)

        # Save everyting in a vtp file
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(filename)
        writer.SetInputData(VesselsPolyData)
        writer.Write()

    def GetMapping(self):
        Vessel = collections.namedtuple('Vessel', ['Vesselname', 'Numberofnodes', 'Startnode', 'Endnode'])
        vessellist = []
        if self.VesselAtlas is None:
            vessellist.append(Vessel(Vesselname='System',
                                     Numberofnodes=len(self.Nodes),
                                     Startnode=0,
                                     Endnode=len(self.Nodes) - 1))
        else:
            for index, vessel in enumerate(self.VesselAtlas):
                vessellist.append(Vessel(Vesselname=vessel,
                                         Numberofnodes=len(self.Vessels[index]),
                                         Startnode=self.Vessels[index][0].Number,
                                         Endnode=self.Vessels[index][-1].Number))
        return vessellist

    def CalculateMaximumTimestep(self, density):
        wavespeed = [2 * vessel.CalculateMaxWaveSpeed(density) for vessel in self.Vessels]
        gridsize = [vessel.GridSize * 1e-3 for vessel in self.Vessels]
        CourantNumber = 1
        timestep = [CourantNumber * gridsize[i] / speed for i, speed in enumerate(wavespeed)]
        dt = min(timestep)

        print("Maximum timestep: %1.1e s" % dt)
        # newgridsize = [1e3 * speed * dt / CourantNumber for speed in wavespeed]  # minimum gridsize, there is no max
        # max gridsize is set by the vessel length and by the need for detail
        return dt

    def GetAdjacencyMatrix(self, bloodvisc):
        """
        Get the adjacency matrix of the network

        :return: scipy sparse matrix, nodes in order
        """
        # for each vessel, get the bifurcations at the ends
        # then check at what parts of other vessels these bifurcations connect.
        # matrix = numpy.zeros((len(self.Vessels),len(self.Vessels)))
        #
        # for index, vessel in enumerate(self.Vessels):
        #     bif = vessel.GetDistalBifurcation()
        #     if not (bif is None):
        #         vesselsnumbers = bif.GetConnectedVesselIDs()
        #         vesselsnumbers.remove(vessel.ID)
        #         for number in vesselsnumbers:
        #             matrix[vessel.ID][number] = 1
        #
        #     bif = vessel.GetProximalBifurcation()
        #     if not (bif is None):
        #         vesselsnumbers = bif.GetConnectedVesselIDs()
        #         vesselsnumbers.remove(vessel.ID)
        #         for number in vesselsnumbers:
        #             matrix[vessel.ID][number] = -1

        network = nx.Graph()
        inlets = [node[0] for node in self.InletNodes]
        bifs = self.BifurcationNodes
        outlets = self.OutletNodes  # vessel ends
        newoutletnodes = [Node.Node() for node in self.OutletNodes]
        allnodes = inlets + bifs + outlets + newoutletnodes

        network.add_nodes_from(allnodes)

        self.GetVesselResistance(bloodvisc)
        edges = [[vessel.GetEndNodes(), vessel.Length, vessel.Resistance] for vessel in self.Vessels]

        outletedges = [[[node, newoutletnodes[i]], 0, node.R1 + node.R2] for i, node in enumerate(outlets)]
        edges += outletedges

        for edge in edges:
            network.add_edge(edge[0][0], edge[0][1], length=edge[1], weight=edge[2])

        labels = {}
        for idx, node in enumerate(network.nodes()):
            labels[node] = idx

        # import matplotlib.pylab as plt
        # pos = nx.spring_layout(network, weight='length')
        # nx.draw_networkx_nodes(network, pos)
        # nx.draw_networkx_edges(network, pos)
        # nx.draw_networkx_labels(network, pos, labels, font_size=16)
        # plt.show()
        matrix = nx.adjacency_matrix(network)
        return matrix, allnodes

    def GetAdjacencyMatrix2(self, bloodvisc):
        """
        Get the adjacency matrix of the network

        :return: scipy sparse matrix, nodes in order
        """
        network = nx.Graph()
        inlets = [node[0] for node in self.InletNodes]
        bifs = self.BifurcationNodes
        outlets = self.OutletNodes  # vessel ends

        allnodes = inlets + bifs + outlets

        network.add_nodes_from(allnodes)

        self.GetVesselResistance(bloodvisc)
        edges = [[vessel.GetEndNodes(), vessel.Length, vessel.Resistance] for vessel in self.Vessels]

        for edge in edges:
            network.add_edge(edge[0][0], edge[0][1], length=edge[1], weight=edge[2])

        labels = {}
        for idx, node in enumerate(network.nodes()):
            labels[node] = idx

        # import matplotlib.pylab as plt
        # pos = nx.spring_layout(network, weight='length')
        # nx.draw_networkx_nodes(network, pos)
        # nx.draw_networkx_edges(network, pos)
        # nx.draw_networkx_labels(network, pos, labels, font_size=16)
        # plt.show()
        matrix = nx.adjacency_matrix(network)
        return matrix, allnodes

    def GetVesselResistance(self, bloodvisc):
        resistance = [vessel.VesselResistance(bloodvisc) for vessel in self.Vessels]
        return resistance

    def SegmentConductance(self, node1, node2, bloodvisc):
        length = abs(node1.LengthAlongVessel - node2.LengthAlongVessel)
        if length == 0:
            length = GeneralFunctions.distancebetweenpoints(node1.Position, node2.Position)

        # a = node1.Radius
        # b = (node2.Radius - node1.Radius) / length
        # radiusint = scipy.integrate.quad(lambda l: 1 / scipy.power(a + b * l, 4), 0, length)
        radius = (node1.Radius + node2.Radius) / 2
        G = 1e-9 * (numpy.pi * numpy.power(radius, 4)) / (22 * bloodvisc * length)
        # G = 1e-9*numpy.pi / (8 * bloodvisc*radiusint[0])
        return G

    def TopologyToGraph(self):
        network = nx.Graph()
        for vessel in self.Vessels:
            for i in range(1, len(vessel.Nodes)):
                network.add_edge(vessel.Nodes[i - 1], vessel.Nodes[i],
                                 weight=vessel.Nodes[i].LengthAlongVessel - vessel.Nodes[i - 1].LengthAlongVessel)

        for bif in self.BifurcationNodes:
            for connection in bif.Connections:
                network.add_edge(bif, connection, weight=0)
        self.Graph = network

        # pos = nx.spring_layout(network)
        # nx.draw(network, pos)
        # plt.show()

    def _Get1DsteadyNetwork(self, bloodvisc):
        print("Calculating steady state model matrix.")
        # calculate the conductance per segment
        edges = []
        for vessel in self.Vessels:
            for i in range(1, len(vessel.Nodes)):
                edges.append([vessel.Nodes[i - 1], vessel.Nodes[i],
                              self.SegmentConductance(vessel.Nodes[i - 1], vessel.Nodes[i], bloodvisc)])

        # outlet nodes, fixed pressure no windkessels
        outletnodes = [[outlet, next(iter(outlet.Connections))] for outlet in self.OutletNodes if outlet.R1 is None]

        # overwrite with wk nodes
        # outletnodes = []
        for outlet in self.OutletNodes:
            if not (outlet.R1 is None):
                g = 1 / (outlet.R1 + outlet.R2)
                wknode = Node.Node()
                wknode.OutPressure = outlet.OutPressure
                outletnodes.append([wknode, outlet])
                edges.append([outlet, wknode, g])

        # map the nodes at the bifurcations correctly
        for bif in self.BifurcationNodes:
            for node in bif.Connections:
                for index, edge in enumerate(edges):
                    if edge[0] == node:
                        edges[index][0] = bif
                    if edge[1] == node:
                        edges[index][1] = bif

        # list of all nodes in the system
        nodes = [i[0] for i in edges]
        [nodes.append(i[1]) for i in edges]
        nodes = list(set(nodes))

        # create a sparse matrix and fill it with the conductances
        mtx = scipy.sparse.lil_matrix((len(nodes), len(nodes)))
        for edge in edges:
            index1 = nodes.index(edge[0])
            index2 = nodes.index(edge[1])
            mtx[index1, index2] = -1 * edge[2]
            mtx[index2, index1] = -1 * edge[2]

        rowsum = mtx.sum(axis=1)
        for index in range(0, len(nodes)):
            mtx[index, index] = -1 * rowsum[index]

        # create the input vector
        inputvector = numpy.zeros((len(nodes),))

        for inletnode in self.InletNodes:
            inputvector[nodes.index(inletnode[0])] = inletnode[0].InletFlowRate  # standard units m^3/s

        for index, outlet in enumerate(outletnodes):
            index1 = nodes.index(outlet[0])
            rownumber = index1
            mtx[rownumber, :] = 0
            mtx[rownumber, rownumber] = 1
            inputvector[rownumber] = outlet[0].OutPressure  # pressure at that wk node, outpressure

        return mtx, inputvector, nodes

    def Get1DsteadyNetwork(self, bloodvisc, clotactive=False):
        edges = []
        # bloodvisc = patient.ModelParameters["BLOOD_VISC"]
        for vessel in self.Vessels:
            for i in range(1, len(vessel.Nodes)):
                edges.append([vessel.Nodes[i - 1],
                              vessel.Nodes[i],
                              self.SegmentConductance(vessel.Nodes[i - 1], vessel.Nodes[i], bloodvisc)])

        # generate more nodes for the wk elements
        newoutlets = []
        for outlet in self.OutletNodes:
            if not (outlet.R1 is None):
                g = 1 / (outlet.R1 + outlet.R2)
                wknode = Node.Node()
                wknode.OutPressure = outlet.OutPressure
                wknode.InputPressure = wknode.OutPressure
                edges.append([outlet, wknode, g])
                    # outlet.AddConnection(wknode)
                wknode.AddConnection(outlet)
                newoutlets.append(wknode)
            else:
                newoutlets.append(outlet)
        OutletNodes = newoutlets

        # modelling clots
        # if we want to use this code as a replacement for the current implementation
        # we can simply update the edges
        # set global parameter ClotActive to true
        # set conductance to zero.
        # find edge in list of edges, set conductance to zero
        if clotactive:
            for clot in self.Clots:
                clotnodes = clot[0]
                for index, node in enumerate(clotnodes[:-1]):
                    for edge in edges:
                        if (edge[0] == node or edge[0] == clotnodes[index + 1]) and (
                                edge[1] == node or edge[1] == clotnodes[index + 1]):
                            edge[2] = 0
                            print("Set conductance to zero.")

        duplicatenodes = self.BifurcationDict()

        nodestoreplace1 = []
        nodestoreplace2 = []
        oldnodes1 = []
        oldnodes2 = []
        for index, edge in enumerate(edges):
            if edge[0] in duplicatenodes:
                nodestoreplace1.append((index, duplicatenodes[edge[0]]))
                oldnodes1.append((index, edge[0]))
            if edge[1] in duplicatenodes:
                nodestoreplace2.append((index, duplicatenodes[edge[1]]))
                oldnodes2.append((index, edge[1]))
        for replace in nodestoreplace1:
            edges[replace[0]][0] = replace[1]
        for replace in nodestoreplace2:
            edges[replace[0]][1] = replace[1]

        nodeset = set()
        for edge in edges:
            nodeset.add(edge[0])
            nodeset.add(edge[1])

        nodesetlist = list(nodeset)
        for index, node in enumerate(nodesetlist):
            node.ssindex = index

        mtx = scipy.sparse.lil_matrix((len(nodesetlist), len(nodesetlist)))

        def updatemat(edge):
            index1 = edge[0].ssindex
            index2 = edge[1].ssindex
            mtx[index1, index2] = -1 * edge[2]
            mtx[index2, index1] = -1 * edge[2]

        with Pool() as pool:
            pool.map(updatemat, edges)

        rowsum = mtx.sum(axis=1)
        for index in range(0, len(nodesetlist)):
            mtx[index, index] = -1 * rowsum[index]

        inputvector = scipy.sparse.lil_matrix((len(nodesetlist), 1))
        for inletnode in self.InletNodes:
            inputvector[inletnode[0].ssindex] = inletnode[0].InletFlowRate  # standard units m^3/s
        for outlet in OutletNodes:
            inputvector[outlet.ssindex] = outlet.OutPressure  # pressure at that wk node, outpressure

        for outlet in OutletNodes:
            index1 = outlet.ssindex
            mtx[index1, next(iter(outlet.Connections)).ssindex] = 0
            mtx[index1, index1] = 1


        return mtx.tocsr(), inputvector, nodesetlist, edges ,(oldnodes1, oldnodes2)

    def nonlinearmodel(self, bloodvisc, density):
        # for each vessel node, create a list of conservation of flow
        # use the matrix from the linear model but with the bifurcations removed

        # we can use the same code for most
        # be sure to subtract the input vector and only return the pressure of the correct nodes

        # the matrix is fine for the linear parts
        # use matrix notation for the nonlinear parts as well
        # code needs to be fast in the future.

        # matrix for the linear system with removed bifurcations replacement
        # return residuals for the vessel nodes inc outlets and inlets

        # calculate the velocity at the bifurcation node
        # calcualte the residuals for the bifurcations

        # probably better to use a loop
        # make a list of edges to store the connectivity information
        # loop over bifurcations, inlet and outlets.

        # store a list of information about which node the pressure belongs to
        nodeslength = len(self.Nodes) - len(self.BifurcationNodes)
        nodes = self.Nodes[:nodeslength]
        # assumption is that the pressure input is in the same order as the nodes

        conductancematrix = numpy.zeros((nodeslength, nodeslength))
        for vessel in self.Vessels:
            # first and last nodes
            conductancematrix[vessel.Nodes[0].Number, vessel.Nodes[1].Number] = -1 * self.SegmentConductance(
                vessel.Nodes[0],
                vessel.Nodes[1],
                bloodvisc)

            conductancematrix[vessel.Nodes[-1].Number, vessel.Nodes[-2].Number] = -1 * self.SegmentConductance(
                vessel.Nodes[-1],
                vessel.Nodes[-2],
                bloodvisc)

            for index, node in enumerate(vessel.Nodes[1:-1]):
                conductancematrix[node.Number, vessel.Nodes[index].Number] = -1 * self.SegmentConductance(
                    node,
                    vessel.Nodes[index],
                    bloodvisc)
                conductancematrix[node.Number, vessel.Nodes[index + 2].Number] = -1 * self.SegmentConductance(
                    node,
                    vessel.Nodes[index + 2],
                    bloodvisc)

        rowsum = conductancematrix.sum(axis=1)
        for index in range(0, len(nodes)):
            conductancematrix[index, index] = -1 * rowsum[index]

        # add the conductance from the WK outlets
        for index, outlet in enumerate(self.OutletNodes):
            conductancematrix[outlet.Number, outlet.Number] += 1 / (outlet.R1 + outlet.R2)

        for vessel in self.Vessels:
            vessel.Nodes[0].NearestVesselNode = vessel.Nodes[1]
            vessel.Nodes[-1].NearestVesselNode = vessel.Nodes[-2]

        # add the vesselnodes that connect to vesselbifurcationnodes for easy access.
        vesselbifnodes = [list(bif.Connections) for bif in self.BifurcationNodes]

        sources = numpy.zeros((nodeslength,))

        # inlets
        for index, inlet in enumerate(self.InletNodes):
            sources[inlet[0].Number] = inlet[0].InletFlowRate

        # outlets
        for index, outlet in enumerate(self.OutletNodes):
            sources[outlet.Number] = outlet.OutPressure / (outlet.R1 + outlet.R2)

        def func(pressures):
            residuals = conductancematrix.dot(pressures) - sources
            for bif in vesselbifnodes:
                # flow = sum([(pressures[node.Number]-pressures[node.NearestVesselNode.Number])*conductancematrix[node.Number,node.NearestVesselNode.Number] for node in bif])
                residuals[bif[0].Number] = (pressures[bif[0].Number] - pressures[bif[0].NearestVesselNode.Number]) * \
                                           conductancematrix[bif[0].Number, bif[0].NearestVesselNode.Number]
                for index, node in enumerate(bif[1:]):
                    residuals[bif[0].Number] += (pressures[node.Number] - pressures[node.NearestVesselNode.Number]) * \
                                                conductancematrix[node.Number, node.NearestVesselNode.Number]
                    velocity = (pressures[node.Number] - pressures[node.NearestVesselNode.Number]) * conductancematrix[
                        node.Number, node.NearestVesselNode.Number] / (numpy.pi * node.Radius * node.Radius)
                    bernoulli = pressures[node.Number] + (density / 2) * velocity * velocity
                    velocity2 = (pressures[bif[index - 1].Number] - pressures[
                        bif[index - 1].NearestVesselNode.Number]) * conductancematrix[
                                    bif[index - 1].Number, bif[index - 1].NearestVesselNode.Number] / (
                                            numpy.pi * bif[index - 1].Radius * bif[index - 1].Radius)
                    bernoulli2 = pressures[bif[index - 1].Number] + (density / 2) * velocity2 * velocity2
                    residuals[node.Number] = bernoulli - bernoulli2

            return residuals

        # guess = func(numpy.zeros((nodeslength,)))
        # guess = [node.Pressure for vessel in self.Vessels for node in vessel.Nodes]
        guess = []
        for vessel in self.Vessels:
            for node in vessel.Nodes:
                guess.append(node.Pressure)
        #

        solution = scipy.optimize.fsolve(func, guess)  # ,full_output=True)
        # print(solution[1])
        # solution = scipy.linalg.solve(conductancematrix, sources)
        return solution, nodes


class Tree:
    def __init__(self, node):
        self.NumberOfGenerations = 0
        node.SetDirectionVector()
        self.Direction = [-1 * i for i in node.DirectionVector]
        self.InitialNode = node
        self.StartingTreeNodes = None
        self.EndNodes = [node]
        self.Nodes = []
        self.Vessels = []
        self.BifurcationNodes = []

    def GenerateTree(self, cutoff):
        while 1:
            Newendnodes = []
            for endnode in self.EndNodes:
                if endnode.Radius > cutoff:
                    # calculate node properties
                    r1 = BloodflowEquations.murraylaw(endnode.Radius)
                    r2 = BloodflowEquations.murraylaw(endnode.Radius)
                    l1 = BloodflowEquations.RadiusToLength(r1)
                    l2 = BloodflowEquations.RadiusToLength(r2)

                    direction = [-1 * i for i in endnode.DirectionVector]

                    angle = math.atan2(direction[1], direction[0])
                    theta = math.atan(1) + angle
                    theta2 = math.atan(-1) + angle

                    posvec = [l1 * math.cos(theta),
                              l1 * math.sin(theta),
                              0]
                    posvec2 = [l2 * math.cos(theta2),
                               l2 * math.sin(theta2),
                               0]

                    pos1end = [endnode.Position[0] + posvec[0],
                               endnode.Position[1] + posvec[1],
                               endnode.Position[2] + posvec[2]]
                    pos2end = [endnode.Position[0] + posvec2[0],
                               endnode.Position[1] + posvec2[1],
                               endnode.Position[2] + posvec2[2]]

                    pos1mid = [endnode.Position[0] + 0.5 * posvec[0],
                               endnode.Position[1] + 0.5 * posvec[1],
                               endnode.Position[2] + 0.5 * posvec[2]]
                    pos2mid = [endnode.Position[0] + 0.5 * posvec2[0],
                               endnode.Position[1] + 0.5 * posvec2[1],
                               endnode.Position[2] + 0.5 * posvec2[2]]

                    # create new nodes
                    bifurcationnode = Node.Node()
                    bifurcationnode.SetPosition(endnode.Position)
                    bifurcationnode.SetMajorVesselID(endnode.MajorVesselID)

                    vessel1node1 = Node.Node()
                    vessel1node1.SetPosition(endnode.Position)
                    vessel1node1.SetLengthAlongVessel(0.0)
                    vessel1node1.SetRadius(r1)
                    vessel1node1.SetMajorVesselID(endnode.MajorVesselID)
                    vessel1node2 = Node.Node()
                    vessel1node2.SetPosition(pos1mid)
                    vessel1node2.SetLengthAlongVessel(0.5 * l1)
                    vessel1node2.SetRadius(r1)
                    vessel1node2.SetMajorVesselID(endnode.MajorVesselID)
                    vessel1node3 = Node.Node()
                    vessel1node3.SetPosition(pos1end)
                    vessel1node3.SetLengthAlongVessel(l1)
                    vessel1node3.SetRadius(r1)
                    vessel1node3.SetMajorVesselID(endnode.MajorVesselID)

                    vessel2node1 = Node.Node()
                    vessel2node1.SetPosition(endnode.Position)
                    vessel2node1.SetLengthAlongVessel(0.0)
                    vessel2node1.SetRadius(r2)
                    vessel2node1.SetMajorVesselID(endnode.MajorVesselID)
                    vessel2node2 = Node.Node()
                    vessel2node2.SetPosition(pos2mid)
                    vessel2node2.SetLengthAlongVessel(0.5 * l2)
                    vessel2node2.SetRadius(r2)
                    vessel2node2.SetMajorVesselID(endnode.MajorVesselID)
                    vessel2node3 = Node.Node()
                    vessel2node3.SetPosition(pos2end)
                    vessel2node3.SetLengthAlongVessel(l2)
                    vessel2node3.SetRadius(r2)
                    vessel2node3.SetMajorVesselID(endnode.MajorVesselID)

                    vessel1 = Vessel.Vessel()
                    vessel2 = Vessel.Vessel()
                    vessel1.SetNodes([vessel1node1, vessel1node2, vessel1node3])
                    vessel2.SetNodes([vessel2node1, vessel2node2, vessel2node3])

                    vessel1.SetLength(l1)
                    vessel2.SetLength(l2)
                    vessel1.SetMeanRadius(r1)
                    vessel2.SetMeanRadius(r2)
                    vessel1.SetGridSize(0.5 * l1)
                    vessel2.SetGridSize(0.5 * l2)
                    vessel1.SetMajorVesselID(endnode.MajorVesselID)
                    vessel2.SetMajorVesselID(endnode.MajorVesselID)
                    # Set connections
                    # do not add them to the main topology for now
                    bifurcationnode.AddConnection(vessel1node1)
                    bifurcationnode.AddConnection(vessel2node1)

                    bifurcationnode.SetRadius((r1+endnode.Radius+r2)/3)

                    vessel1node1.AddConnection(bifurcationnode)
                    vessel1node1.AddConnection(vessel1node2)
                    vessel1node2.AddConnection(vessel1node1)
                    vessel1node2.AddConnection(vessel1node3)
                    vessel1node3.AddConnection(vessel1node2)

                    vessel2node1.AddConnection(bifurcationnode)
                    vessel2node1.AddConnection(vessel2node2)
                    vessel2node2.AddConnection(vessel2node1)
                    vessel2node2.AddConnection(vessel2node3)
                    vessel2node3.AddConnection(vessel2node2)

                    if self.NumberOfGenerations != 0:
                        endnode.AddConnection(bifurcationnode)
                        bifurcationnode.AddConnection(endnode)
                    else:
                        self.StartingTreeNodes = bifurcationnode

                    self.BifurcationNodes.append(bifurcationnode)
                    Newendnodes.append(vessel1node3)
                    Newendnodes.append(vessel2node3)

                    self.Vessels.append(vessel1)
                    self.Vessels.append(vessel2)
                    self.Nodes.append(bifurcationnode)
                    self.Nodes.append(vessel1node1)
                    self.Nodes.append(vessel1node2)
                    self.Nodes.append(vessel1node3)
                    self.Nodes.append(vessel2node1)
                    self.Nodes.append(vessel2node2)
                    self.Nodes.append(vessel2node3)
                    vessel1node3.SetDirectionVector()
                    vessel2node3.SetDirectionVector()

            if len(Newendnodes) == 0:
                break
            self.EndNodes = Newendnodes
            self.NumberOfGenerations += 1
        # rotation to match the direction of the original vessel

        Vo = [-1 * self.Direction[0], -1 * self.Direction[1], -1 * self.Direction[2]]
        Vn = [-1 * self.Direction[0], -1 * self.Direction[1], 0]

        Vcross = numpy.cross(Vn, Vo)
        C = numpy.dot(Vn, Vo)
        Vx = [[0, -1 * Vcross[2], Vcross[1]],
              [Vcross[2], 0, -1 * Vcross[0]],
              [-1 * Vcross[1], Vcross[0], 0]]

        R = numpy.eye(3) + Vx + numpy.matmul(Vx, Vx) * (1 / (1 + C))
        Rnew = numpy.eye(4)
        Rnew[0:3, 0:3] = R

        Translate1 = numpy.array(
            [[1, 0, 0, -1 * self.InitialNode.Position[0]],
             [0, 1, 0, -1 * self.InitialNode.Position[1]],
             [0, 0, 1, -1 * self.InitialNode.Position[2]],
             [0, 0, 0, 1]])
        Translate2 = numpy.array(
            [[1, 0, 0, 1 * self.InitialNode.Position[0]],
             [0, 1, 0, 1 * self.InitialNode.Position[1]],
             [0, 0, 1, 1 * self.InitialNode.Position[2]],
             [0, 0, 0, 1]])

        TMatrix = numpy.dot(Translate2, numpy.dot(Rnew, Translate1))

        for node in self.Nodes:
            pos = node.Position
            vec = numpy.array([[pos[0]], [pos[1]], [pos[2]], [1]])
            position = numpy.dot(TMatrix, vec)
            posnew = [position[0][0], position[1][0], position[2][0]]
            node.SetPosition(posnew)

    def ConnectTree(self):
        self.InitialNode.AddConnection(self.StartingTreeNodes)
        self.StartingTreeNodes.AddConnection(self.InitialNode)
