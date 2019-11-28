#!/usr/bin/python3
import sys
import pathos.multiprocessing

import os
import random
import shutil
import csv
import time

from os.path import dirname, realpath

sys.path.append(realpath(dirname(__file__)))

import math
import numpy
import scipy
import scipy.spatial
import vtk
from scipy.optimize import newton
import scipy.integrate
from vtk.util.numpy_support import vtk_to_numpy
from scipy.sparse.linalg import spsolve as spsolve
import networkx as nx
import BloodflowEquations
import GeneralFunctions
import Perfusion
import Results
import Topology
import Metadata
import Vessel
import Node
import Remesh

import copy

fractions = [0.65, 0.05, 0.05, 0.25]
bodyparts = ['Thoracic aorta', 'R. brachial', 'L. brachial']

class Patient:
    def __init__(self, folder=os.path.dirname(os.path.realpath(sys.argv[0]))):
        folder = folder + "/"
        self.Folders = Metadata.Folders(folder)
        self.ModelParameters = Metadata.ModelParameter()
        self.PatientData = Metadata.PatientData()
        self.Topology = Topology.Topology()
        self.Results = Results.Results()
        self.Perfusion = Perfusion.Perfusion()

        self.Trees = []
        print("Patient folder: " + self.Folders.PatientFolder)

    def RemoveOldSimFiles(self):
        folder = self.Folders.ModellingFolder
        deletelist = ["Results.dyn", "ResultsPrev.dyn", "ResultsTotal.dyn", "Convergence.csv", "Conv.txt"]
        for name in deletelist:
            if os.path.isfile(folder + name):
                try:
                    os.remove(folder + name)
                except:
                    pass
        deletefolder = ["TimeSeries"]
        for name in deletefolder:
            try:
                shutil.rmtree(folder + name)
            except:
                pass

    def ResetModellingFolder(self):
        folder = self.Folders.ModellingFolder
        filenames = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]
        blacklist = ["boundary_4&21&22&23&24&25&26&30.h5", "boundary_4&21&22&23&24&25&26&30.ply",
                     "Boundary.ply", "labelled_vol_mesh_high.msh", "labelled_vol_mesh.msh",
                     "permeability_tensors.h5", "Model_parameters.txt", "PialSurface.vtp"]
        for name in filenames:
            if name not in blacklist:
                try:
                    os.remove(folder + name)
                    print("Deleting file: ", folder + name)
                except:
                    print("Error while deleting file: ", folder + name)

        folders = [f for f in os.listdir(folder) if not os.path.isfile(os.path.join(folder, f))]
        for name in folders:
            try:
                shutil.rmtree(folder + name)
            except:
                print("Error while deleting file ", folder + name)

    def LoadBFSimFiles(self):
        self.Topology.LoadBFSimFiles(self.Folders.ModellingFolder)
        self.LoadPatientData()
        self.LoadModelParameters()
        self.LoadRegionMapping()

    def LoadResults(self, file="Results.dyn", correct=True):
        self.Results.LoadResults(self.Folders.ModellingFolder + file)
        if correct == False:
            self.CorrectforDirectionVectorOutlets()
        self.Results.CalculateVelocity()
        self.Results.CalculateMeanResultsPerNode()

    def GetMeanResults(self):
        self.Results.GetMeanResults(self.Topology.Vessels)

    def ExportMeanResults(self,file="ResultsPerVessel.csv"):
        self.Results.ExportMeanResults(self.Folders.ModellingFolder,file)

    def LoadPatientData(self, file="config.xml"):
        if file[-3:] == "txt":
            self.PatientData.LoadPatientData(self.Folders.InputFolder + file)
        elif file[-3:] == "xml":
            self.PatientData.LoadPatientDataXML(self.Folders.InputFolder + file)
        else:
            print("File open error")
            exit(1)

    def LoadClusteringMapping(self, file):
        print("Load Clustering map to file: %s" % file)
        data = [line.strip('\n').split(',') for line in open(file)][1:]
        for line in data:
            cp = Perfusion.CouplingPoint(self.Topology.Nodes[int(line[0])])
            cp.Node.Position = [float(i) for i in line[1].split(" ")]
            cp.NumberOfTriangles = int(line[2])
            cp.Area = float(line[4])
            cp.Node.MajorVesselID = int(line[5])
            self.Perfusion.CouplingPoints.append(cp)

    def LoadModelParameters(self, file="Model_parameters.txt"):
        self.ModelParameters.LoadModelParameters(self.Folders.ModellingFolder + file)

    def LoadSegmentedVessels(self):
        self.Topology.LoadSegmentedVessels(self.Folders.InputFolder)

    def Load1DAnatomy(self, file="1-D_Anatomy_Patient.txt"):
        self.Topology.Load1DAnatomy(self.Folders.InputFolder + file)

    def Load3DTopFile(self, file="0_2_cutoff.top"):
        self.Topology.Read3DNodesFromTopFile(file)
        self.Topology.CheckConnectivity()
        self.Topology.AnatomyToVessels()

    def ClusteringByRegion(self, dualgraph, method=""):
        self.Perfusion.ClusteringByRegion(dualgraph, method)

    def SetFolder(self, folder):
        self.Folders.SetPatientFolder(folder)

    def TopologyToVTP(self, filename="Topology.vtp"):
        self.Topology.TopologyToVTP(self.Folders.ModellingFolder + filename)

    def SaveVesselAtlas(self, filename="Mapping.csv"):
        self.Topology.SaveVesselAtlas(self.Folders.ModellingFolder + filename)

    def LoadVTPFile(self, vtpfile):
        self.Topology.LoadVTPFile(vtpfile)

    def WriteClusteringMapping(self, filename="Clusters.csv"):
        self.Perfusion.WriteClusteringMapping(self.Folders.ModellingFolder + filename)

    def WriteModelParameters(self):
        file = self.Folders.ModellingFolder + "Model_parameters.txt"
        print("Writing the Model parameters to %s" % file)
        self.ModelParameters.WriteModelParameters(file)

    def CalculateMaximumTimestep(self):
        timestep = self.Topology.CalculateMaximumTimestep(self.ModelParameters["Density"])
        t = float(format(timestep, '1.1e'))
        # iterations = numpy.ceil(self.ModelParameters["Beat_Duration"] / t)
        # t= self.ModelParameters["Beat_Duration"] / iterations
        # if timestep < self.ModelParameters["TIMESTEP"]:
        #     print("Warning: Simulation can be unstable!")
        #     print("Timestep exceeds the maximum timestep.")
        print("Timestep set to maximum timestep: %1.1e s" % t)
        self.ModelParameters["TIMESTEP"] = t

    def CalculateDistanceFromTheHeart(self):
        """
        Returns the distance from the heart in the shortest sense bansed on Dijkstra's distance.
        Direction is determined by the positive pressure direction.
        :return: Nothing
        """
        if not self.Topology.Graph:
            self.Topology.TopologyToGraph()

        for edge in self.Topology.Graph.edges:
            node1 = edge[0]
            node2 = edge[1]
            p1 = self.Results.MeanPressurePerNode[-1][node1.Number]
            p2 = self.Results.MeanPressurePerNode[-1][node2.Number]

            dp = p1 - p2
            dx = abs(self.Topology.Graph[node1][node2]['weight'])
            if dp == 0:
                continue

            if dp < 0:
                self.Topology.Graph[node2][node1]['direction'] = dx
            elif dp > 0:
                self.Topology.Graph[node1][node2]['direction'] = dx

        length, _ = nx.single_source_dijkstra(self.Topology.Graph, self.Topology.Nodes[0], weight='direction')
        for node in self.Topology.Nodes:
            node.DistancefromHeart = length[node]
        self.Results.DistancefromHeart = [node.DistancefromHeart for node in self.Topology.Nodes]

    def CalculateTimedelayFromTheHeart(self):
        """
        Returns the timedelay from the heart in the shortest sense bansed on Dijkstra's distance.
        Direction is determined by the positive pressure direction.
        :return: Nothing
        """
        if not self.Topology.Graph:
            self.Topology.TopologyToGraph()

        for edge in self.Topology.Graph.edges:
            node1 = edge[0]
            node2 = edge[1]
            v1 = self.Results.MeanVelocityPerNode[-1][node1.Number]
            v2 = self.Results.MeanVelocityPerNode[-1][node2.Number]

            meanv = (v1+v2)/2
            dx = self.Topology.Graph[node1][node2]['weight']*1e-3
            t = dx/meanv
            if meanv == 0:
                continue

            if t < 0:
                self.Topology.Graph[node2][node1]['timedelay'] = abs(t)
            elif t > 0:
                self.Topology.Graph[node1][node2]['timedelay'] = t
            else:
                self.Topology.Graph[node2][node1]['timedelay'] = 0
                self.Topology.Graph[node1][node2]['timedelay'] = 0

        length, _ = nx.single_source_dijkstra(self.Topology.Graph, self.Topology.Nodes[0], weight='timedelay')
        for node in self.Topology.Nodes:
            node.TimeDelay = length[node]
        self.Results.DistancefromHeart = [node.DistancefromHeart for node in self.Topology.Nodes]

    def UpdateNodeType(self):
        self.Topology.UpdateNodeType()

    def CorrectforDirectionVectorOutlets(self):
        numbers = [vessel.Nodes[-1].Number for vessel in self.Topology.Vessels]
        self.Results.CorrectForDirectionEnds(numbers)

    def CerebellumBrainstemMapping(self):
        # add labels and positions to the vessels.
        cerebellum = 8
        brainstem = 9

        cerebellumnames = ["R. SCA",
                           "L. SCA",
                           "R. AICA",
                           "L. AICA",
                           "R. PICA",
                           "L. PICA", ]

        brainstemnames = ["Pontine I",
                          "Pontine II",
                          "Pontine III",
                          "Pontine IV",
                          "Pontine V",
                          "Pontine VI",
                          "Pontine VII",
                          "Pontine VIII",
                          "Pontine IX",
                          "Pontine X",
                          "Pontine XI",
                          "Pontine XII",
                          "Pontine XIII"]

        cerebellumpos = [[25, -74, -36],
                         [-25, -74, -34],
                         [28, -55, -48],
                         [-28, -55, -48],
                         [32, -40, -60],
                         [-30, -48, -60]]

        # starting at the bottom of the stem
        pontinepos = [[-3.1, -16.3, -58.5],
                      [1.6, -15.7, -57.2],
                      [-6.1, -14.0, -50.4],
                      [2.2, -12.8, -51.7],
                      [-8.2, -10.9, -48.2],
                      [2.39, -9.2, -48.9],
                      [-8.8, -5.81, -46.1],
                      [2.7, -3.2, -46.6],
                      [-8.4, -0.6, -40.7],
                      [0.3, 1.3, -41.3],
                      [-6.8, 3.0, -34.6],
                      [-3.0, 1.4, -34.0],
                      [0, 0, 0]]

        for vessel in self.Topology.Vessels:
            if vessel.Name in cerebellumnames:
                vessel.SetMajorVesselID(cerebellum)
                vessel.Nodes[-1].SetPosition(cerebellumpos[cerebellumnames.index(vessel.Name)])
                vessel.Nodes[-1].SetMajorVesselID(cerebellum)
            if vessel.Name in brainstemnames:
                vessel.SetMajorVesselID(brainstem)
                vessel.Nodes[-1].SetPosition(pontinepos[brainstemnames.index(vessel.Name)])
                vessel.Nodes[-1].SetMajorVesselID(brainstem)

    def calculate_wk_parameters(self):
        """
        Calculate the windkessel parameters for the outlets of the network.
        """
        print("Calculating WindKessel Parameters.")
        scaling = 1e-3  # default for mm



        C1D = 0
        for vessel in self.Topology.Vessels:
            # mean radius and thickness
            lengthvessel = vessel.Length * scaling
            meanradius = vessel.MeanRadius * scaling
            meanh = vessel.MeanThickness * scaling
            C1D += 2 * numpy.power(numpy.pi * meanradius * meanradius, 1.5) * lengthvessel / (
                    (4 / 3) * numpy.sqrt(numpy.pi) * vessel.YoungsModules * meanh)

        Cperipheral = self.ModelParameters["CTotal"] - C1D
        self.ModelParameters["Cperipheral"] = Cperipheral
        if Cperipheral < 0:
            raise Exception('Cperipheral should not be negative. The value of Cperipheral was: {}'.format(Cperipheral))

        # identify outlets
        outletmap = []
        for i in self.Topology.OutletNodes:
            for index, vessel in enumerate(self.Topology.Vessels):
                if i in vessel.Nodes:
                    outletmap.append(vessel.Name)

        bodymap = []
        for out in outletmap:
            part = len(fractions) - 1
            for index, body in enumerate(bodyparts):
                if out == body:
                    part = index
            bodymap.append(part)

        for bodynumber, fraction in enumerate(fractions):
            outlets = [out for index, out in enumerate(self.Topology.OutletNodes) if bodymap[index] == bodynumber]
            rcubed = 0
            for i in outlets:
                radius = i.Radius * scaling
                rcubed += math.pow(radius, 3)

            for i in outlets:
                radius = i.Radius * scaling
                h = i.Thickness * scaling
                rt = (self.ModelParameters["RTotal"] / fraction) * rcubed / math.pow(radius, 3)
                A = math.pi * math.pow(radius, 2)
                beta = (4 / 3) * math.sqrt(math.pi) * i.YoungsModules * h / A
                c0 = abs(math.sqrt(beta / (2 * self.ModelParameters["Density"]))) * math.pow(A, 0.25)
                r1 = (self.ModelParameters["Density"] / A) * c0
                r2 = rt - r1
                if r2 < 0:
                    r2 = 0.1e9
                    r1 = rt - r2
                c = self.ModelParameters["Cperipheral"] * self.ModelParameters["RTotal"] / rt
                i.SetWK(r1, r2, c)

    def calculate_wk_parameters_evenly(self):
        """
        Calculate the windkessel parameters for the outlets of the network.
        """
        print("Calculating WindKessel Parameters.")
        scaling = 1e-3  # default for mm

        C1D = 0
        for vessel in self.Topology.Vessels:
            # mean radius and thickness
            lengthvessel = vessel.Length * scaling
            meanradius = vessel.MeanRadius * scaling
            meanh = vessel.MeanThickness * scaling
            C1D += 2 * numpy.power(numpy.pi * meanradius * meanradius, 1.5) * lengthvessel / (
                    (4 / 3) * numpy.sqrt(numpy.pi) * vessel.YoungsModules * meanh)

        Cperipheral = self.ModelParameters["CTotal"] - C1D
        self.ModelParameters["Cperipheral"] = Cperipheral
        if Cperipheral < 0:
            raise Exception('Cperipheral should not be negative. The value of Cperipheral was: {}'.format(Cperipheral))

        # identify outlets
        outletmap = []
        for i in self.Topology.OutletNodes:
            for index, vessel in enumerate(self.Topology.Vessels):
                if i in vessel.Nodes:
                    outletmap.append(vessel.Name)

        outlets = [out for index, out in enumerate(self.Topology.OutletNodes)]
        rcubed = 0
        for i in outlets:
            radius = i.Radius * scaling
            rcubed += math.pow(radius, 3)

        for i in outlets:
            radius = i.Radius * scaling
            h = i.Thickness * scaling
            rt = (self.ModelParameters["RTotal"]) * rcubed / math.pow(radius, 3)
            A = math.pi * math.pow(radius, 2)
            beta = (4 / 3) * math.sqrt(math.pi) * i.YoungsModules * h / A
            c0 = abs(math.sqrt(beta / (2 * self.ModelParameters["Density"]))) * math.pow(A, 0.25)
            r1 = (self.ModelParameters["Density"] / A) * c0
            r2 = rt - r1
            if r2 < 0:
                r2 = 0.1e9
                r1 = rt - r2
            c = self.ModelParameters["Cperipheral"] * self.ModelParameters["RTotal"] / rt
            i.SetWK(r1, r2, c)

    def WriteSimFiles(self):
        print("Writing simulation files.")
        self.Topology.UpdateTopology()
        self.SaveVesselAtlas()
        self.WriteTopologyFile()
        self.WriteParFile()
        self.WriteFlowProfiles()
        self.WriteRunFile()
        self.WriteModelParameters()
        self.WriteClotFile()
        self.WriteRegionMapping()

    def ExportSurface(self, file, graph):
        self.Perfusion.ExportSurface(file, graph)

    def SelectPatient(self):
        # Load input files
        # Randomly select a patient
        segmentationfolder = self.Folders.ScriptFolder + "/Segmentations/"
        segfolders = [f for f in os.listdir(segmentationfolder)]
        selectedfolder = random.sample(segfolders, k=1)[0]

        segmentationfile = segmentationfolder + selectedfolder + "/Feature_Vessel.csv"
        infofile = segmentationfolder + selectedfolder + "/Image_info.txt"

        shutil.copy(segmentationfile, self.Folders.InputFolder)
        shutil.copy(infofile, self.Folders.InputFolder)
        self.ModelParameters["Patient_Segmentation"] = selectedfolder
        print("Selected segmentation: %s" % selectedfolder)

    def SelectDonorNetwork(self):
        """
        Select a network from the Brava set to extend the patient network.
        For now, the donor is randomly selected from the folder.
        """
        scriptfolder = self.Folders.ScriptFolder
        bravafolder = scriptfolder + "/Brava/"
        bravafiles = [f for f in os.listdir(bravafolder) if os.path.isfile(os.path.join(bravafolder, f)) if
                      f[-3:] == "vtp"]
        selectedfile = random.sample(bravafiles, k=1)[0]
        # selectedfile = "BH0023_ColorCoded.CNG.vtp"
        self.ModelParameters["Donor_Network"] = selectedfile

    def UpdateModelParameters(self):
        self.ModelParameters["Beat_Duration"] = 60 / self.PatientData["HeartRate"]
        self.ModelParameters["SISTOLIC_PRESSURE"] = self.PatientData["SystolePressure"]
        self.ModelParameters["DIASTOLIC_PRESSURE"] = self.PatientData["DiastolePressure"]
        self.ModelParameters["RTotal"] = ((self.PatientData["SystolePressure"] * (1 / 3) + self.PatientData[
            "DiastolePressure"] * (2 / 3)) - self.PatientData["MeanRightAtrialPressure"]) / (
                                                 self.PatientData["StrokeVolume"] * (
                                                 self.PatientData["HeartRate"] / 60) * 1e-6)
        self.ModelParameters["TimeConstantDiastolicDecay"] = 1.34
        self.ModelParameters["CTotal"] = self.ModelParameters["TimeConstantDiastolicDecay"] / self.ModelParameters[
            "RTotal"]
        self.ModelParameters["ScriptLocation"] = os.path.dirname(os.path.realpath(__file__)) + "/Check_Convergence.py"

    def WriteFlowProfilesAlastruey2007(self):
        """
        Take patient parameters and generate flow profiles.
        Assumption is that the two inlets with the largest radius are the carotid arteries.
        Flowrate given here are ml/s
        """
        print("Writing Flow Profiles.")
        Aortafile = self.Folders.ModellingFolder + "Aorta.txt"
        period = 60 / self.PatientData["HeartRate"]
        time = numpy.linspace(0, period, 1000)

        VolumeFlowRate = self.PatientData["StrokeVolume"] * (self.PatientData["HeartRate"] / 60)
        # AortaFlow = BloodflowEquations.FlowRateAorta(time, period, VolumeFlowRate)
        scaling = scipy.integrate.simps(BloodflowEquations.FlowRateAlastruey2007(time,1), time)
        AortaFlow = BloodflowEquations.FlowRateAlastruey2007(time, VolumeFlowRate/scaling)
        # AortaFlow = [self.PatientData["StrokeVolume"] for i in AortaFlow]
        # totalflow = numpy.trapz(AortaFlow,time)
        # totalflows = simps(AortaFlow, time)
        with open(Aortafile, 'w') as f:
            f.write("1.0e-6\n")
            for i in range(0, len(time)):
                f.write('{}\t{}\n'.format(time[i], AortaFlow[i]))

    def WriteFlowProfiles(self):
        """
        Take patient parameters and generate flow profiles.
        Assumption is that the two inlets with the largest radius are the carotid arteries.
        Flowrate given here are ml/s
        """
        print("Writing Flow Profiles.")
        Aortafile = self.Folders.ModellingFolder + "Aorta.txt"
        period = 60 / self.PatientData["HeartRate"]
        time = numpy.linspace(0, period, 1000)

        VolumeFlowRate = self.PatientData["StrokeVolume"] * (self.PatientData["HeartRate"] / 60)
        AortaFlow = BloodflowEquations.FlowRateAorta(time, period, VolumeFlowRate)
        # scaling = scipy.integrate.simps(BloodflowEquations.FlowRateAlastruey2007(time,1), time)
        # AortaFlow = BloodflowEquations.FlowRateAlastruey2007(time, VolumeFlowRate/scaling)
        # AortaFlow = [self.PatientData["StrokeVolume"] for i in AortaFlow]
        # totalflow = numpy.trapz(AortaFlow,time)
        # totalflows = simps(AortaFlow, time)
        with open(Aortafile, 'w') as f:
            f.write("1.0e-6\n")
            for i in range(0, len(time)):
                f.write('{}\t{}\n'.format(time[i], AortaFlow[i]))

    # def WriteResultsToVTP(self):
    #     print("Writing Results to vtp files")
    #     resultsfolder = self.Folders.ModellingFolder
    #     if not os.path.exists(resultsfolder):
    #         os.mkdir(resultsfolder)
    #
    #     TimeFile = open(resultsfolder + "Time.txt", 'w', encoding='utf-8')
    #     for i in range(0, len(self.Results.TimePoints)):
    #         TimeFile.write(str(i) + " " + str(self.Results.TimePoints[i].WT) + '\n')
    #     TimeFile.close()
    #
    #     for i in range(0, len(self.Results.TimePoints)):
    #         self.FrameToVTP(resultsfolder + "TimePoint" + str(i) + ".vtp", self.Results.TimePoints[i])
    #
    # def FrameToVTP(self, filename, timepoint):
    #     print("Writing frame at time:%f" % timepoint.WT)
    #     nodes = vtk.vtkPoints()
    #     vessels = vtk.vtkCellArray()
    #
    #     radius = vtk.vtkFloatArray()
    #     radius.SetNumberOfComponents(1)
    #     radius.SetName("Radius")
    #
    #     pressure = vtk.vtkFloatArray()
    #     pressure.SetNumberOfComponents(1)
    #     pressure.SetName("Pressure")
    #
    #     Flow = vtk.vtkFloatArray()
    #     Flow.SetNumberOfComponents(1)
    #     Flow.SetName("Flow Rate")
    #
    #     Radius = vtk.vtkFloatArray()
    #     Radius.SetNumberOfComponents(1)
    #     Radius.SetName("Radius")
    #
    #     for i in range(0, len(timepoint.Pressure)):
    #         pressure.InsertNextValue(timepoint.Pressure[i])
    #         Flow.InsertNextValue(timepoint.Flow[i])
    #         Radius.InsertNextValue(timepoint.Radius[i])
    #
    #     # Add radius and position to data array
    #     for node in self.Topology.Nodes:
    #         nodes.InsertNextPoint(node.Position)
    #         radius.InsertNextValue(node.Radius)
    #
    #     # Create a polydata to store everything in
    #     VesselsPolyData = vtk.vtkPolyData()
    #
    #     # Add the points to the dataset
    #     VesselsPolyData.SetPoints(nodes)
    #
    #     # Add vessels to cell array
    #     for vessel in self.Topology.Vessels:
    #         line = vtk.vtkLine()
    #         line.GetPointIds().SetNumberOfIds(len(vessel))
    #         for i in range(0, len(vessel)):
    #             line.GetPointIds().SetId(i, vessel[i])
    #         vessels.InsertNextCell(line)
    #     # Add the lines to the dataset
    #     VesselsPolyData.SetLines(vessels)
    #
    #     # Assign radii to the nodes
    #     VesselsPolyData.GetPointData().SetScalars(radius)
    #     VesselsPolyData.GetPointData().AddArray(pressure)
    #     VesselsPolyData.GetPointData().AddArray(Flow)
    #     VesselsPolyData.GetPointData().AddArray(Radius)
    #
    #     # Save everyting in a vtk file
    #     writer = vtk.vtkXMLPolyDataWriter()
    #     writer.SetFileName(filename)
    #     writer.SetInputData(VesselsPolyData)
    #     writer.Write()

    def VesselToMeshAllignmentSides(self, graph):
        print("Optimizing alignment between the mesh and the outlets.")

        pos = []
        majornodesID = []
        for index, node in enumerate(self.Topology.OutletNodes):
            if (not (node.Position is None)) and node.MajorVesselID >= 0:
                pos.append(node.Position)
                majornodesID.append(node.MajorVesselID)

        # left and right sides of the pial surface
        rightpialsurface = [index for index, node in enumerate(graph.PialSurface) if graph.SideInfo[index] > 0]
        leftpialsurface = [index for index, node in enumerate(graph.PialSurface) if graph.SideInfo[index] < 0]

        # [cow, - , "R. ACA, A2", "R. MCA", "L. MCA", "L. ACA, A2", "R. PCA, P2", "L. PCA, P2"]
        # leftvessels = [4, 5, 7]
        # rightvessels = [2, 3, 6]
        leftvessels = self.Perfusion.leftregionids
        rightvessels = self.Perfusion.rightregionids

        leftnodes = [index for index, node in enumerate(pos) if majornodesID[index] in leftvessels]
        rightnodes = [index for index, node in enumerate(pos) if majornodesID[index] in rightvessels]

        sys.setrecursionlimit(10000)
        KDTreeleft = scipy.spatial.KDTree([graph.PialSurface[i] for i in leftpialsurface])
        KDTreeright = scipy.spatial.KDTree([graph.PialSurface[i] for i in rightpialsurface])

        def optimfun(args):
            scaling = args[0]
            translation = [args[1], args[2], args[3]]
            rotation = [args[4], args[5], args[6]]
            # Calculate the transform matrix
            tmatrix = GeneralFunctions.TMatrix(scaling, rotation, translation)
            # Apply matrix to the nodes
            tpos = []
            for node in pos:
                vec = numpy.array([[node[0]], [node[1]], [node[2]], [1]])
                position = numpy.dot(tmatrix, vec)
                nodenew = [position[0][0], position[1][0], position[2][0]]
                tpos.append(nodenew)

            meshdis = [-1 for node in tpos]

            MinDistanceR, MinDistanceIndexR = KDTreeright.query([tpos[i] for i in rightnodes], k=1)
            MinDistanceL, MinDistanceIndexL = KDTreeleft.query([tpos[i] for i in leftnodes], k=1)

            for index, nodeindex in enumerate(rightnodes):
                meshdis[nodeindex] = MinDistanceR[index]

            for index, nodeindex in enumerate(leftnodes):
                meshdis[nodeindex] = MinDistanceL[index]

            return sum(meshdis)

        # initial guess for the allignment
        result = scipy.optimize.minimize(optimfun, numpy.array([1.5, 0, 0, 0, 0, 0, 0]),
                                         method='Nelder-Mead', options={'disp': True})
        values = result.x
        self.Perfusion.AllignmentResult = result
        print(values)
        scaling = values[0]
        translation = [values[1], values[2], values[3]]
        rotation = [values[4], values[5], values[6]]
        self.Perfusion.TMatrix = GeneralFunctions.TMatrix(scaling, rotation, translation)
        self.ApplyTransformation(self.Perfusion.TMatrix)

    def ApplyTransformation(self, matrix):
        for vesselnode in self.Topology.Nodes:
            if not (vesselnode.Position is None):
                node = vesselnode.Position
                vec = numpy.array([[node[0]], [node[1]], [node[2]], [1]])
                position = numpy.dot(matrix, vec)
                nodenew = [position[0][0], position[1][0], position[2][0]]
                vesselnode.SetPosition(nodenew)

    def FindCouplingPoints(self):
        """
        Map the outlets to the nearest point on the surface.
        """
        surfaceoutlets = [i for i in self.Topology.OutletNodes if i.MajorVesselID >= 2]
        print("Finding outlets close to the surface.")
        for index, node in enumerate(surfaceoutlets):
            couplingpoint = Perfusion.CouplingPoint(node)
            self.Perfusion.AddCouplingPoint(couplingpoint)

    def RemoveCoWOutletsFromPerfusion(self):
        """
        The donor networks have outlets that do not get transferred.
        These need to be removed.
        If the coupling point node does not have a number, remove it from the list.
        """
        print("Removing coupling points outside topology.")
        pointstoremove = [point for point in self.Perfusion.CouplingPoints if not (point.Node in self.Topology.Nodes)]
        [self.Perfusion.RemoveCouplingPoint(point) for point in pointstoremove]

    def WriteCouplingPointsFile(self):
        print("Writing coupling points to file.")
        if not self.Perfusion.CouplingPoints:
            self.FindCouplingPoints()

        with open(self.Folders.ModellingFolder + 'Coupling.csv', 'w') as f:
            f.write("ID,Position,Direction,Radius,MeshPoint\n")
            for couplingpoint in self.Perfusion.CouplingPoints:
                f.write(str(couplingpoint.NodeID)
                        + "," + str(couplingpoint.Node.Position)
                        + "," + str(couplingpoint.Node.DirectionVector)
                        + "," + str(couplingpoint.Node.Position.Radius)
                        + "," + str(couplingpoint.PialSurfacePoint) + "\n")

    def WriteClotFile(self):
        print("Writing Clots to file.")
        with open(self.Folders.ModellingFolder + 'Clots.txt', 'w') as f:
            f.write("ClotID,NodeID,Permeability,Porosity\n")
            for index, clot in enumerate(self.Topology.Clots):
                for node in clot[0]:
                    f.write(str(index) + "," + str(node.Number) + "," + str(clot[1]) + "," + str(clot[2]) + "\n")

    def WriteTopologyFile(self):
        file = self.Folders.ModellingFolder + "System.top"
        print("Writing Topology file to %s" % file)

        OutputFile = open(file, 'w', encoding='utf-8')
        OutputFile.write('Name: System_0' + '\nCoordinates:\n')

        # Output Vessel position and radius
        for node in self.Topology.Nodes:
            OutputLine = str(node.Number) \
                         + " L:" + "{0:0.5f}".format(node.LengthAlongVessel) \
                         + " R:" + "{0:0.5f}".format(node.Radius) \
                         + " E:" + "{0:0.1f}".format(node.YoungsModules * 1e-6) \
                         + " T:" + "{0:d}".format(node.Type) + '\n'

            OutputFile.write(OutputLine)
        # Output Connections
        OutputFile.write('\nBonds:\n')
        for node in self.Topology.Nodes:
            BondsLine = ''
            for bond in node.Connections:
                BondsLine += ' ' + str(bond.Number)
            OutputFile.write(str(node.Number) + BondsLine + '\n')
        OutputFile.close()

    def WriteParFile(self):
        file = self.Folders.ModellingFolder + "System.par"
        print("Writing parameter file to %s" % file)
        # Output outlets
        OutputFilePar = open(file, 'w', encoding='utf-8')
        if self.Topology.OutletNodes[0].R1 is None:
            self.calculate_wk_parameters()

        for node in self.Topology.OutletNodes:
            OutputFilePar.write(str(node.Number)
                                + ' R1:' + "{0:0.5f}".format(node.R1 * 1e-9)
                                + ' R2:' + "{0:0.5f}".format(node.R2 * 1e-9)
                                + ' C:' + "{0:0.5f}".format(node.C * 1e12) + '\n')
        OutputFilePar.close()

    def WriteRunFile(self):
        file = self.Folders.ModellingFolder + "Run.txt"
        print("Writing run file to %s" % file)
        # Output Run file
        OutputRunFile = open(file, 'w', encoding='utf-8')
        OutputRunFile.write('<System>\n')
        OutputRunFile.write('Topology: ' + "System.top" + '\n')
        InletLine = ""
        for node in self.Topology.InletNodes:
            InletLine += str(node[0].Number) + ' ' + node[1] + ', '
        OutputRunFile.write('InletFlux: ' + InletLine[0:-2] + '\n')
        OutputRunFile.write('OutletParams: ' + "System.par" + '\n')
        OutputRunFile.write('Task: ')
        OutputRunFile.write('\n' + 'NumberCoupledPoints: 0' + '\n')
        modelparameterfile = "Model_parameters.txt"
        OutputRunFile.write('Modelparameters: ' + modelparameterfile + '\n')
        OutputRunFile.close()

    def GenerateTreesAtCps(self, cutoff):
        print("Generating bifurcating arterial trees at the outlet nodes.")
        for id in self.Perfusion.CouplingPoints:
            tree = self.GenerateTree(id.Node, cutoff)
            id.Tree = tree

        print("Generated %d bifurcating trees." % len(self.Trees))

    def GenerateTrees(self):
        print("Generating bifurcating arterial trees at the outlet nodes.")
        for id in self.Topology.OutletNodes:
            self.GenerateTree(id)
        print("Generated %d bifurcating trees." % len(self.Trees))

    def getTotalEndsofTrees(self):
        numberofterminals = sum([len(tree.EndNodes) for tree in self.Trees])
        print("Total number of terminal points of the arterial trees: %d" % numberofterminals)

    def PerfusionEndsNumber(self):
        print("Calculating the number of terminal points per outlet.")
        for i in range(0, len(self.Perfusion.CouplingPoints)):
            self.Perfusion.CouplingPoints[i].NumberOfPialPoints = max(1, len(self.Trees[i].EndNodes))

    def GenerateTree(self, node, cutoff=0.125):
        tree = Topology.Tree(node)
        tree.GenerateTree(cutoff)
        self.Trees.append(tree)
        return tree

    def TreesToTopology(self):
        print("Creating a tree topology.")
        treetop = Topology.Topology()
        for tree in self.Trees:
            treetop.Nodes += tree.Nodes
        return treetop

    def AddTreesToTop(self):
        print("Adding trees to the vessel topology.")
        for tree in self.Trees:
            if len(tree.Nodes) > 0:
                tree.ConnectTree()
                self.Topology.Nodes.extend(tree.Nodes)
                # self.Topology.OutletNodes += tree.EndNodes
                # self.Topology.OutletNodes.remove(tree.InitialNode)
                self.Topology.BifurcationNodes.extend(tree.BifurcationNodes)
                name = self.Topology.GetVesselNameFromNode(tree.InitialNode)
                currentvesselnumber = self.Topology.Vessels[-1].ID
                [t.SetID(currentvesselnumber + i + 1) for i, t in enumerate(tree.Vessels)]
                self.Topology.Vessels.extend(tree.Vessels)
                [vessel.SetName(name + "Tree_" + str(i)) for i, vessel in enumerate(tree.Vessels)]

        self.Topology.UpdateTopology()

    def TreeEndsToEndNodes(self):
        for tree in self.Trees:
            self.Topology.OutletNodes.remove(tree.InitialNode)
            self.Topology.OutletNodes.extend(tree.EndNodes)

        print("Removing resistance added by the merging of networks.")
        visc = float(self.ModelParameters["BLOOD_VISC"])
        density = float(self.ModelParameters["Density"])
        for tree in self.Trees:
            name = self.Topology.GetVesselNameFromNode(tree.InitialNode)
            self.Topology.DownStreamResistance(name, visc, density)

    def _Add1TreeToTop(self):  # todo fix
        print("Adding trees to the vessel topology.")
        tree = self.Trees[0]
        if len(tree.Nodes) > 0:
            tree.ConnectTree()
            self.Topology.Nodes += tree.Nodes
            self.Topology.OutletNodes += tree.EndNodes
            self.Topology.OutletNodes.remove(tree.InitialNode)
            self.Topology.BifurcationNodes += tree.BifurcationNodes
            name = self.Topology.GetVesselNameFromNode(tree.InitialNode)
            self.Topology.Vessels += tree.Vessels
            [vessel.SetName(name + "Tree_" + str(i)) for i, vessel in enumerate(tree.Vessels)]
        self.Topology.UpdateTopology()

    # def WriteMeanVesselResultsToVTP(self):
    #     print("Writing mean results to Mean.vtp.")
    #     self.FrameToVTP(self.Folders.ModellingFolder + 'Mean.vtp', self.Results.CalcMean())

    def LoadClusters(self):
        file = self.Folders.ModellingFolder + "Clusters.csv"
        clusters = [i.strip('\n').split(',') for i in open(file)]

        NodeID = [int(i[0]) for i in clusters[1:]]
        Position = [i[1].split(' ') for i in clusters[1:]]
        Position = [[float(i[0]), float(i[1]), float(i[2])] for i in Position]
        CouplingPointsN = [int(i[2]) for i in clusters[1:]]
        ClusterID = [int(i[3]) for i in clusters[1:]]
        areas = [float(i[4]) for i in clusters[1:]]
        majorvesselid = [int(i[5]) for i in clusters[1:]]
        # nodes = [i[5] for i in clusters[1:]]

        return NodeID, Position, CouplingPointsN, ClusterID, areas, majorvesselid

    def LoadSurfaceMapping(self,filename="SurfaceNodesMapping.csv"):
        print("Loading surface mapping.")
        file = self.Folders.ModellingFolder + filename
        if file[-3:] == "csv":
            with open(file, 'r') as csvfile:
                spamreader = csv.reader(csvfile, delimiter=',')
                clusters = [i for i in spamreader]
        else:
            with open(file, 'r') as csvfile:
                spamreader = csv.reader(csvfile, delimiter='\t')
                clusters = [i for i in spamreader]


        # clusters = [i.strip('\n').split(',') for i in open(file)]

        NodeID = [int(i[0]) for i in clusters[1:]]
        ClusterID = [[int(j) for j in i[1].split(',')] for i in clusters[1:]]
        Areas = [float(i[2]) for i in clusters[1:]]
        # print(len(self.Perfusion.CouplingPoints) == len(set([y for x in ClusterID for y in x])))
        if len(self.Perfusion.CouplingPoints) == len(set([y for x in ClusterID for y in x])):
            for node, cluster in zip(NodeID, ClusterID):
                for id in cluster:
                    self.Perfusion.CouplingPoints[id].SurfaceNodes.append(node)

        return NodeID, ClusterID, Areas

    def DistributeFlowVertex(self):
        print("Distributing flow to the coupling points")
        # use this when the clustering is done with the primal graph, i.e. the vertices.
        # the only difference with the dual graph clustering is the type of data, cell vs point and the clustersize might be different
        nodeids, _, _, clusterids, Areas, majorvesselid = self.LoadClusters()
        surfacenodes, clustermap, triangleArea = self.LoadSurfaceMapping()

        pointspercluster = [sum([1 for i in clustermap if i == id]) for id in clusterids]

        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(self.Folders.ModellingFolder + "Clustering_mapped.vtp")
        reader.Update()
        data = reader.GetOutput()
        resultsfolder = self.Folders.ModellingFolder
        if not os.path.exists(resultsfolder):
            os.mkdir(resultsfolder)

        for timeindex, timepoint in enumerate(self.Results.TimePoints):
            flow = vtk.vtkFloatArray()
            flow.SetNumberOfComponents(1)
            flow.SetName("Flow")

            for surfacenode in surfacenodes:
                cluster = clustermap[surfacenode]
                fractionflow = 0
                for i in cluster:
                    outlet = nodeids[i]
                    outletflow = timepoint.Flow[outlet]
                    clustersize = pointspercluster[i]
                    fractionflow += outletflow / clustersize

                flow.InsertNextValue(fractionflow)

            writer = vtk.vtkXMLPolyDataWriter()
            filename = resultsfolder + "FlowDistributed" + str(timeindex) + ".vtp"
            writer.SetFileName(filename)
            data.GetPointData().AddArray(flow)
            writer.SetInputData(data)
            writer.Write()

        meanflow = [numpy.mean([tp.Flow[node] for tp in self.Results.TimePoints]) for node in nodeids]
        flow = vtk.vtkFloatArray()
        flow.SetNumberOfComponents(1)
        flow.SetName("Flow")
        for surfacenode in surfacenodes:
            cluster = clustermap[surfacenode]
            fractionflow = 0
            for i in cluster:
                outletflow = meanflow[i]
                clustersize = pointspercluster[i]
                fractionflow += outletflow / clustersize
            flow.InsertNextValue(fractionflow)

        writer = vtk.vtkXMLPolyDataWriter()
        filename = resultsfolder + "FlowDistributedMean.vtp"
        writer.SetFileName(filename)
        data.GetPointData().AddArray(flow)
        writer.SetInputData(data)
        writer.Write()

    def GetMeanFlowRates(self):
        nodeids, _, _, clusterids, areas, majorvesselid = self.LoadClusters()
        time = [tp.WT for tp in self.Results.TimePoints]
        totalflow = 0
        regionsflowvolume = [0 for i in range(0, 8)]
        clusterflowrate = []

        for index, node in enumerate(nodeids):
            flowdatanode = [tp.Flow[node] for tp in self.Results.TimePoints]
            flowvolume = numpy.trapz(flowdatanode, time)
            flowrate = flowvolume / time[-1]
            totalflow += flowvolume
            major = majorvesselid[index] - 2
            regionsflowvolume[major] += flowvolume
            clusterflowrate.append(flowrate)

        regionsflowrateregion = [element / time[-1] for element in regionsflowvolume]
        sumflow = sum([flow for index, flow in enumerate(clusterflowrate) if majorvesselid[index] == 3])

        return regionsflowrateregion, clusterflowrate

    def DistributeFlowTriangles(self):
        print("Distributing flow to the coupling points")
        # use this when the clustering is done with the primal graph, i.e. the vertices.
        # the only difference with the dual graph clustering is the type of data, cell vs point and the clustersize might be different
        nodeids, positions, couplingpointsN, clusterids, areas, majorvesselid = self.LoadClusters()
        surfacenodes, clustermap, triangleAreas = self.LoadSurfaceMapping()

        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(self.Folders.ModellingFolder + "Clustering.vtp")
        reader.Update()
        data = reader.GetOutput()
        outputfolder = self.Folders.ModellingFolder + "TimeSeries/"
        if not os.path.exists(outputfolder):
            os.mkdir(outputfolder)

        for timeindex, timepoint in enumerate(self.Results.TimePoints):
            flow = vtk.vtkFloatArray()
            flow.SetNumberOfComponents(1)
            flow.SetName("Volume Flow Rate (mL/s)")

            pressure = vtk.vtkFloatArray()
            pressure.SetNumberOfComponents(1)
            pressure.SetName("Pressure (Pa)")

            flowrate = vtk.vtkFloatArray()
            flowrate.SetNumberOfComponents(1)
            flowrate.SetName("Flow Velocity (m/s)")

            for clusters, surfacenode,triangleArea in zip(clustermap, surfacenodes, triangleAreas):
                FlowRatePerOutlet = [timepoint.Flow[nodeids[i]] / areas[i] for i in clusters]
                PressurePerCluster = [timepoint.Pressure[nodeids[i]]for i in clusters]
                FlowPerTriangle = sum(FlowRatePerOutlet)*triangleArea

                flow.InsertNextValue(FlowPerTriangle)
                flowrate.InsertNextValue(sum(FlowRatePerOutlet))
                pressure.InsertNextValue(numpy.mean(PressurePerCluster))

            writer = vtk.vtkXMLPolyDataWriter()
            filename = outputfolder + "FlowDistributed" + str(timeindex) + ".vtp"
            writer.SetFileName(filename)
            data.GetCellData().AddArray(flow)
            data.GetCellData().AddArray(flowrate)
            data.GetCellData().AddArray(pressure)
            writer.SetInputData(data)
            writer.Write()

        time = [tp.WT for tp in self.Results.TimePoints]
        _, meanflowrate = self.GetMeanFlowRates()
        flow = vtk.vtkFloatArray()
        flow.SetNumberOfComponents(1)
        flow.SetName("Volume Flow Rate (mL/s)")
        flowrate = vtk.vtkFloatArray()
        flowrate.SetNumberOfComponents(1)
        flowrate.SetName("Flow Velocity (m/s)")
        pressure = vtk.vtkFloatArray()
        pressure.SetNumberOfComponents(1)
        pressure.SetName("Pressure (Pa)")

        with open(self.Folders.ModellingFolder + "Results.pvd", "w") as f:
            f.write("<?xml version=\"1.0\"?>\n")
            f.write(
                "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n")
            f.write("<Collection>\n")
            for timeindex, timepoint in enumerate(self.Results.TimePoints):
                name = "TimeSeries/FlowDistributed" + str(timeindex) + ".vtp"
                time = str(timepoint.WT)
                f.write("<DataSet timestep=\"" + time + "\" group=\"\" part=\"0\" file=\"" + name + "\"/>\n")
            f.write("</Collection>\n")
            f.write("</VTKFile>")

        for clusters, surfacenode, triangleArea in zip(clustermap, surfacenodes, triangleAreas):
            FlowRatePerOutlet = [self.Results.MeanVolumeFlowRatePerNode[-1][nodeids[i]] / areas[i] for i in clusters]
            PressurePerCluster = [self.Results.MeanPressurePerNode[-1][nodeids[i]] for i in clusters]
            FlowPerTriangle = sum(FlowRatePerOutlet) * triangleArea

            flow.InsertNextValue(FlowPerTriangle)
            flowrate.InsertNextValue(sum(FlowRatePerOutlet))
            pressure.InsertNextValue(numpy.mean(PressurePerCluster))

        writer = vtk.vtkXMLPolyDataWriter()
        filename = self.Folders.ModellingFolder + "FlowDistributedMean.vtp"
        writer.SetFileName(filename)
        data.GetCellData().AddArray(flow)
        data.GetCellData().AddArray(flowrate)
        data.GetCellData().AddArray(pressure)
        writer.SetInputData(data)
        writer.Write()

    def WriteTimeseriesVessels(self):
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(self.Folders.ModellingFolder + "Topology.vtp")
        reader.Update()
        data = reader.GetOutput()
        outputfolder = self.Folders.ModellingFolder + "TimeSeriesBF/"
        if not os.path.exists(outputfolder):
            os.mkdir(outputfolder)

        for timeindex, timepoint in enumerate(self.Results.Time[-1]):
            flow = vtk.vtkFloatArray()
            flow.SetNumberOfComponents(1)
            flow.SetName("Volume Flowrate (mL/s)")

            velocity = vtk.vtkFloatArray()
            velocity.SetNumberOfComponents(1)
            velocity.SetName("Velocity (m/s)")

            pressure = vtk.vtkFloatArray()
            pressure.SetNumberOfComponents(1)
            pressure.SetName("Pressure (pa)")

            radius = vtk.vtkFloatArray()
            radius.SetNumberOfComponents(1)
            radius.SetName("Radius (mm)")

            for nodeindex in range(0, len(self.Results.Pressure[-1])):
                flow.InsertNextValue(self.Results.VolumeFlowRate[-1][nodeindex][timeindex])
                velocity.InsertNextValue(self.Results.Velocity[-1][nodeindex][timeindex])
                pressure.InsertNextValue(self.Results.Pressure[-1][nodeindex][timeindex])
                radius.InsertNextValue(self.Results.Radius[-1][nodeindex][timeindex])

            writer = vtk.vtkXMLPolyDataWriter()
            filename = outputfolder + "Bloodflow" + str(timeindex) + ".vtp"
            writer.SetFileName(filename)
            data.GetPointData().AddArray(flow)
            data.GetPointData().AddArray(velocity)
            data.GetPointData().AddArray(pressure)
            data.GetPointData().AddArray(radius)
            writer.SetInputData(data)
            writer.Write()

        with open(self.Folders.ModellingFolder + "ResultsBF.pvd", "w") as f:
            f.write("<?xml version=\"1.0\"?>\n")
            f.write(
                "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n")
            f.write("<Collection>\n")
            for timeindex, timepoint in enumerate(self.Results.TimePoints):
                name = "TimeSeriesBF/Bloodflow" + str(timeindex) + ".vtp"
                time = str(timepoint.WT)
                f.write("<DataSet timestep=\"" + time + "\" group=\"\" part=\"0\" file=\"" + name + "\"/>\n")
            f.write("</Collection>\n")
            f.write("</VTKFile>")

    def ImportClotNodeFile(self):
        file = self.Folders.ModellingFolder + "Clots.txt"
        clotdata = [i.strip('\n').split(',') for i in open(file)][1:]
        ClotID = [int(i[0]) for i in clotdata]
        NodeID = [int(i[1]) for i in clotdata]
        Permeability = [float(i[2]) for i in clotdata]
        Porosity = [float(i[3]) for i in clotdata]
        uniqueclots = set(ClotID)
        clots = [[] for _ in uniqueclots]
        for id in uniqueclots:
            for index, clotid in enumerate(ClotID):
                if clotid == id:
                    clots[id].append((NodeID[index], Permeability[index], Porosity[index]))
        print(1)

    def ImportClots(self, file=""):
        if not file:
            file = self.Folders.InputFolder + "Clots.txt"
        # first add new nodes with at the location of the clot
        # then map the clot to the node
        if not GeneralFunctions.is_non_zero_file(file):
            print("No clotfile found.")
            return
        print("Importing Clots.")
        Clots = [i.strip('\n').split('\t') for i in open(file)]

        for clot in Clots[1:]:
            vesselname, location, length, permeability, porosity = clot
            vessel = self.Topology.VesselAtlas[str(vesselname)]
            numbergridnodes = int(numpy.ceil(vessel.Length / 2.5))  # update resolution to about 2.5mm
            vessel.UpdateResolution(numbergridnodes)
            self.Topology.UpdateTopology()
            location, length, permeability, porosity = float(location), float(length), float(permeability), float(
                porosity)

            # nearest node in the vessel
            clotnodes = set()
            # boundary
            boundary1 = numpy.argmin([abs(node.LengthAlongVessel - location) for node in vessel.Nodes])
            clotnodes.add(vessel.Nodes[boundary1])
            # internal nodes
            [clotnodes.add(node) for node in vessel.Nodes if
             node.LengthAlongVessel >= location and node.LengthAlongVessel <= (location + length)]
            # other boundary
            boundary2 = numpy.argmin(
                [abs(node.LengthAlongVessel - (location + length)) for node in vessel.Nodes])
            clotnodes.add(vessel.Nodes[boundary2])

            clotnodes = list(clotnodes)
            # if clot extends past a bifurcation, map the clot nodes to the same clot
            self.Topology.Clots.append((clotnodes, permeability, porosity, length))

    def ShowResults(self):
        # wait for results to be saved
        # while not is_non_zero_file(self.Results.File):
        #     print("Waiting for data.")
        #     time.sleep(1)
        # self.Results.LoadResults(self.Results.File)
        # open the results window
        app = Results.QtWidgets.QApplication(sys.argv)
        win = Results.ResultsWindow(self)
        app.exec()
        # sys.exit(app.exec_())

    def MappingSixMajorArteries(self):
        """
        Follow the major cerebreal arteries upstream and return the outlets.
        """
        print("Mapping the outlets of the six major cerebral arteries.")
        self.Topology.RedefineDirection()
        arteries = ["R. ACA, A2",
                    "R. MCA",
                    "L. MCA",
                    "L. ACA, A2",
                    "R. PCA, P2",
                    "L. PCA, P2", ]

        upstream = []
        for vessel in arteries:
            upstream.append(self.Topology.GetDownstreamVessels(vessel))

        for index, region in enumerate(upstream):
            vessels = region[0]
            for vessel in vessels:
                for node in vessel.Nodes:
                    node.MajorVesselID = index + 2

        artertydict = dict()
        for index, node in enumerate(self.Topology.Nodes):
            artertydict[node] = node.MajorVesselID
        self.Perfusion.MajorArteries = artertydict

        outletarteries = [artertydict[node] for node in self.Topology.OutletNodes]
        for indexregion, region in enumerate(arteries):
            numberofoutlets = len([out for out in outletarteries if out == indexregion + 2])
            print(region + ": " + str(numberofoutlets) + " outlets")

    def ExportDualClustering(self, primal="GMsurface.vtp", dual="Clustering_mapped.vtp"):
        # import the primal graph
        if primal[-3:] == "vtp":
            reader = vtk.vtkXMLPolyDataReader()
        elif primal[-3:] == "ply":
            reader = vtk.vtkPLYReader()
        else:
            print("Error: unreadable file.")
            return 1

        reader.SetFileName(primal)
        reader.Update()
        primaldata = reader.GetOutput()

        # import the dual graph with the triangle colour
        if dual[-3:] == "vtp":
            reader = vtk.vtkXMLPolyDataReader()
        elif dual[-3:] == "ply":
            reader = vtk.vtkPLYReader()
        reader.SetFileName(dual)
        reader.Update()

        narray = reader.GetOutput().GetPointData().GetNumberOfArrays()
        # vertex colour of the dual graph is the triangle colour of the primal graph
        for i in range(0, narray):
            arrayname = reader.GetOutput().GetPointData().GetArrayName(i)
            if arrayname == "Colour":
                clusteringdata = reader.GetOutput().GetPointData().GetArray(i)
                primaldata.GetCellData().AddArray(clusteringdata)
            if arrayname == "Major Cerebral Artery":
                clusteringdata = reader.GetOutput().GetPointData().GetArray(i)
                primaldata.GetCellData().AddArray(clusteringdata)

        # export to new file
        writer = vtk.vtkXMLPolyDataWriter()
        filename = self.Folders.ModellingFolder + "ClusteringByTriangle.vtp"
        writer.SetFileName(filename)
        writer.SetInputData(primaldata)
        writer.Write()

    def ExportTriangleFlowData(self):
        print("Exporting Flow Data.")
        file = self.Folders.ModellingFolder + "FlowDistributedMean.vtp"

        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(file)
        reader.Update()

        narray = reader.GetOutput().GetCellData().GetNumberOfArrays()
        for i in range(0, narray):
            arrayname = reader.GetOutput().GetCellData().GetArrayName(i)
            if arrayname == "Major Cerebral Artery":
                mapCA = vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(i))
            if arrayname == "Volume Flow Rate (mL/s)":
                flowdata = vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(i))
            if arrayname == "Pressure (Pa)":
                pressuredata = vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(i))
            if arrayname == "Colour":
                flowpercluster = vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(i))
                # todo change for multiple elements per cluster
                flowpercluster = [int(i) for i in flowpercluster]

        flowrateperregion, clusterflowrate = self.GetMeanFlowRates()

        trianglesperregion = [sum([1 for i in mapCA if i == cluster + 2]) for cluster in range(0, 8)]
        surfacenodes, clustermap, triangleAreas = self.LoadSurfaceMapping()

        numberofclusters = len(list(set(flowpercluster)))
        trianglespercluster = [sum([1 for i in flowpercluster if i == cluster]) for cluster in
                               range(0, numberofclusters)]

        totalAreaPerCluster = [0 for _ in range(numberofclusters)]
        for clusterid in range(0,numberofclusters):
            for triangle, area in zip(clustermap,triangleAreas):
                if clusterid in triangle:
                    totalAreaPerCluster[clusterid] +=area

        flowPerAreaPerCluster = [clusterflowrate[i]/totalAreaPerCluster[i] for i in range(0,numberofclusters)]

        flowPerTriangle = [0 for _ in surfacenodes]
        for triandleID, trianglemapping, area in zip(surfacenodes, clustermap,triangleAreas):
            flowPerTriangle[triandleID] = sum([flowPerAreaPerCluster[id]*area for id in trianglemapping])

        with open(self.Folders.ModellingFolder + "TriangleToRegion.csv", "w") as f:
            f.write("Triangle ID,Cluster ID,Region ID,Flowrate over the triangle(mL/s),Pressure(Pa)\n")
            for index, flow in zip(surfacenodes,flowPerTriangle):
                clustersPerTriangle = ','.join(map(str, clustermap[index]))
                f.write("%d,\"%s\",%d,%.10f,%.10f\n" % (index, clustersPerTriangle, mapCA[index], flow, pressuredata[index]))

        with open(self.Folders.ModellingFolder + "RegionFlowData.csv", "w") as f:
            f.write("Region ID,Flow rate over the region(mL/s),number of triangles\n")
            for index in range(0, 8):
                f.write("%d,%.10f,%d\n" % (index, flowrateperregion[index], trianglesperregion[index]))

        with open(self.Folders.ModellingFolder + "ClusterFlowData.csv", "w") as f:
            f.write("Cluster ID,Flow rate over the region(mL/s),number of triangles,Pressure(Pa)\n")
            for index in range(0, numberofclusters):
                f.write("%d,%.10f,%d,%.10f\n" % (
                index, clusterflowrate[index], trianglespercluster[index], pressuredata[index]))

    def MajorArteriesToPialSurfaceNN(self, graph):
        print("Mapping the surface to the nearest major artery.")
        # map to nearest node
        # colour based on nearest outlet node.
        majornodes = []
        majornodesID = []
        for index, node in enumerate(self.Topology.Nodes):
            if (not (node.Position is None)) and node.MajorVesselID >= 0:
                majornodes.append(node.Position)
                majornodesID.append(node.MajorVesselID)

        # left and right sides of the pial surface
        rightpialsurface = [index for index, node in enumerate(graph.PialSurface) if
                            graph.SideInfo[index] > 0]
        leftpialsurface = [index for index, node in enumerate(graph.PialSurface) if
                           graph.SideInfo[index] < 0]

        # [cow, - , "R. ACA, A2", "R. MCA", "L. MCA", "L. ACA, A2", "R. PCA, P2", "L. PCA, P2"]
        leftvessels = [4, 5, 7, 8, 9]
        rightvessels = [2, 3, 6, 8, 9]

        leftnodes = [index for index, node in enumerate(majornodes) if majornodesID[index] in leftvessels]
        rightnodes = [index for index, node in enumerate(majornodes) if majornodesID[index] in rightvessels]

        sys.setrecursionlimit(10000)
        KDTreeleft = scipy.spatial.KDTree([majornodes[i] for i in leftnodes])
        KDTreeright = scipy.spatial.KDTree([majornodes[i] for i in rightnodes])

        meshcolour = [-1 for node in graph.PialSurface]

        MinDistanceR, MinDistanceIndexR = KDTreeright.query([graph.PialSurface[i] for i in rightpialsurface], k=1)
        MinDistanceL, MinDistanceIndexL = KDTreeleft.query([graph.PialSurface[i] for i in leftpialsurface], k=1)
        for index, nodeindex in enumerate(rightpialsurface):
            meshcolour[nodeindex] = majornodesID[rightnodes[MinDistanceIndexR[index]]]
        for index, nodeindex in enumerate(leftpialsurface):
            meshcolour[nodeindex] = majornodesID[leftnodes[MinDistanceIndexL[index]]]

        graph.MajorVesselID = meshcolour

    def MajorArteriesToPialSurfaceOutlets(self, graph):
        print("Mapping the surface to the nearest major artery.")
        # map to nearest node
        # colour based on nearest outlet node.
        majornodes = []
        majornodesID = []
        for index, node in enumerate(self.Topology.OutletNodes):
            if (not (node.Position is None)) and node.MajorVesselID >= 0:
                majornodes.append(node.Position)
                majornodesID.append(node.MajorVesselID)

        # left and right sides of the pial surface
        rightpialsurface = [index for index, node in enumerate(graph.PialSurface) if
                            graph.SideInfo[index] > 0]
        leftpialsurface = [index for index, node in enumerate(graph.PialSurface) if
                           graph.SideInfo[index] < 0]

        # [cow, - , "R. ACA, A2", "R. MCA", "L. MCA", "L. ACA, A2", "R. PCA, P2", "L. PCA, P2"]
        leftvessels = [4, 5, 7, 8, 9]
        rightvessels = [2, 3, 6, 8, 9]

        leftnodes = [index for index, node in enumerate(majornodes) if majornodesID[index] in leftvessels]
        rightnodes = [index for index, node in enumerate(majornodes) if majornodesID[index] in rightvessels]

        sys.setrecursionlimit(10000)
        KDTreeleft = scipy.spatial.KDTree([majornodes[i] for i in leftnodes])
        KDTreeright = scipy.spatial.KDTree([majornodes[i] for i in rightnodes])

        meshcolour = [-1 for node in graph.PialSurface]

        MinDistanceR, MinDistanceIndexR = KDTreeright.query([graph.PialSurface[i] for i in rightpialsurface], k=1)
        MinDistanceL, MinDistanceIndexL = KDTreeleft.query([graph.PialSurface[i] for i in leftpialsurface], k=1)
        for index, nodeindex in enumerate(rightpialsurface):
            meshcolour[nodeindex] = majornodesID[rightnodes[MinDistanceIndexR[index]]]
        for index, nodeindex in enumerate(leftpialsurface):
            meshcolour[nodeindex] = majornodesID[leftnodes[MinDistanceIndexL[index]]]

        graph.MajorVesselID = meshcolour

    def WriteRegionMapping(self,filename="RegionMap.csv"):
        file = self.Folders.ModellingFolder + filename
        print("Writing the mapping per region to file: %s" % file)

        with open(file, 'w') as f:
            f.write("NodeID,MajorVesselID\n")
            for index, node in enumerate(self.Topology.Nodes):
                regionid = node.MajorVesselID
                f.write(str(node.Number) + "," + str(regionid) + "\n")

    def LoadRegionMapping(self,filename="RegionMap.csv"):
        file = self.Folders.ModellingFolder + filename
        regiondata = [i.strip('\n').split(',') for i in open(file)]
        regiondict = dict()
        for nodedata in regiondata[1:]:
            regiondict[self.Topology.Nodes[int(nodedata[0])]] = int(nodedata[1])
        self.Perfusion.MajorArteries = regiondict

    def Initiate1DSteadyStateModel(self):
        print("Initiating steady state model.")
        referencepressure = self.PatientData["DiastolePressure"]
        outpressure = self.ModelParameters["OUT_PRESSURE"]

        for inlet in self.Topology.InletNodes:
            inlet[0].InletFlowRate = self.PatientData["StrokeVolume"] * (self.PatientData["HeartRate"] / 60) * 1e-6

        for node in self.Topology.Nodes:
            node.RefRadius = copy.deepcopy(node.Radius)
            node.RefPressure = copy.deepcopy(referencepressure)
            node.SetPressureAreaEquation()

        # Set pressure at the outlets (without windkessel)
        for outlet in self.Topology.OutletNodes:
            if outlet.R1 is None:
                outlet.OutPressure = outlet.InputPressure  # set pressure here per outlet
            else:
                outlet.OutPressure = outpressure

    def Solve1DSteadyState(self, clotactive=False):
        visc = self.ModelParameters["BLOOD_VISC"]

        matrix, inputvector, nodes, edges, mapping = self.Topology.Get1DsteadyNetwork(visc,clotactive)
        oldnodes1, oldnodes2 = mapping

        # solve system
        # solution = scipy.linalg.solve(matrix.todense(), inputvector)
        print(f"Solving system of size: {matrix.shape}")
        solution = spsolve(matrix, inputvector)
        # solved = numpy.allclose(matrix.dot(solution), inputvector)

        # set pressure to all nodes
        for i in range(0, len(nodes)):
            nodes[i].Pressure = solution[i]

        for bif in self.Topology.BifurcationNodes:
            for node in bif.Connections:
                node.Pressure = bif.Pressure

        # update edges to include the dublicate ends again.
        for replace in oldnodes1:
            edges[replace[0]][0] = replace[1]
        for replace in oldnodes2:
            edges[replace[0]][1] = replace[1]

        for edge in edges:
            dp = edge[0].Pressure - edge[1].Pressure
            flowrate = dp * edge[2] * 1e6
            edge[0].FlowRate = flowrate
            edge[1].FlowRate = flowrate
            radius = edge[0].Radius * 0.5 + 0.5 * edge[1].Radius
            velocity = flowrate / (scipy.pi * radius * radius)
            edge[0].Velocity = velocity
            edge[1].Velocity = velocity

        # calculate new radius
        # for vessel in self.Topology.Vessels:
        for node in self.Topology.Nodes:
            node.UpdateRadius()
            # vessel.CalculateMeanRadius()

        solution = [node.Pressure for node in self.Topology.Nodes]
        return solution

    def Solve1DSteadyStateNonLinear(self):
        visc = self.ModelParameters["BLOOD_VISC"]

        solution, nodes = self.Topology.nonlinearmodel(visc, self.ModelParameters["Density"])

        # solve system
        # solution = scipy.linalg.solve(matrix.todense(), inputvector)
        # solved = numpy.allclose(matrix.dot(solution), inputvector)

        # set pressure to all nodes
        for i in range(0, len(nodes)):
            nodes[i].Pressure = solution[i]

        # for bif in self.Topology.BifurcationNodes:
        #     for node in bif.Connections:
        #         node.Pressure = bif.Pressure

        # calculate new radius
        for vessel in self.Topology.Vessels:
            for node in vessel.Nodes:
                node.UpdateRadius()
            vessel.CalculateMeanRadius()

        solution = [node.Pressure for node in self.Topology.Nodes]
        # self.Topology.nonlinearmodel(visc, self.ModelParameters["Density"])
        return solution

    def UpdateOutletResistanceToPialSurface(self, PialSurfacePressure=8000):
        for index, cp in enumerate(self.Perfusion.CouplingPoints):
            newresistancetotal = (cp.Node.Pressure - PialSurfacePressure) / (cp.Node.FlowRate * 1e-6)
            cp.Node.R2 = newresistancetotal - cp.Node.R1
            if cp.Node.R2 < 0:
                cp.Node.R2 = 0.1e9
                cp.Node.R1 = newresistancetotal - cp.Node.R2
                print("Error: R2<0. Setting R2 to 0.1e9 and updating R1.")
            cp.Node.OutPressure = PialSurfacePressure  # pressure now fixed to this value at this node

    def UpdatePressureCouplingPoints(self, pressures):
        for index, cp in enumerate(self.Perfusion.CouplingPoints):
            cp.Node.R1 = None  # remove wk elements
            cp.Node.OutPressure = pressures[index]  # pressure now fixed to this value at this node

    def Run1DSteadyStateModel(self, model="Linear", clotactive=False):
        print("Solving 1D steady state system.")
        # initialise all values
        start = time.time()

        solutionold = self.Solve1DSteadyState(clotactive)
        if model == "NonLinear":
            solutionold = self.Solve1DSteadyStateNonLinear()
        # iterate system again
        iter = 0
        while True:
            if model == "Linear":
                solutionnew = self.Solve1DSteadyState(clotactive)
            elif model == "NonLinear":
                solutionnew = self.Solve1DSteadyStateNonLinear()
            else:
                print("Please specify if model is Linear or NonLinear.")
                return
            # difference = numpy.sqrt(sum([(solutionold[i] - solution) * (solutionold[i] - solution) for i, solution in
            #                   enumerate(solutionnew)]))

            difference = numpy.linalg.norm(
                [solutionold[i] - solution for i, solution in enumerate(solutionnew)]) / numpy.linalg.norm(solutionnew)
            print("Convergence criteria %f" % difference)
            solutionold = solutionnew
            iter += 1
            if difference < 1e-6:
                print("1D steady state convergence in %d iterations." % iter)
                end = time.time()
                print(end - start)
                break

        # calculate flowrate
        # visc = self.ModelParameters["BLOOD_VISC"]
        # for vessel in self.Topology.Vessels:
        #     for i in range(1, len(vessel.Nodes)):
        #         dp = vessel.Nodes[i - 1].Pressure - vessel.Nodes[i].Pressure
        #         flowrate = dp * self.Topology.SegmentConductance(vessel.Nodes[i - 1], vessel.Nodes[i],
        #                                                          visc) * 1e6
        #         vessel.Nodes[i - 1].FlowRate = flowrate
        #         vessel.Nodes[i - 1].Velocity = flowrate / (
        #                 scipy.pi * vessel.Nodes[i - 1].Radius * vessel.Nodes[i - 1].Radius)

        # # vessel ends
        # for vessel in self.Topology.Vessels:
        #     dp = vessel.Nodes[-2].Pressure - vessel.Nodes[-1].Pressure
        #     flowrate = dp * self.Topology.SegmentConductance(vessel.Nodes[-1], vessel.Nodes[-2],
        #                                                      visc) * 1e6
        #     vessel.Nodes[- 1].FlowRate = flowrate
        #     vessel.Nodes[- 1].Velocity = flowrate / (
        #             scipy.pi * vessel.Nodes[- 1].Radius * vessel.Nodes[- 1].Radius)

        for bif in self.Topology.BifurcationNodes:
            bif.Pressure = numpy.mean([node.Pressure for node in list(bif.Connections)])
            bif.FlowRate = numpy.mean([node.FlowRate for node in list(bif.Connections)])
            bif.Velocity = numpy.mean([node.Velocity for node in list(bif.Connections)])

    def Results1DSteadyStateModel(self):
        # export to results
        # we get mean results per node by definition.
        self.Results.MeanPressurePerNode = [[node.Pressure for node in self.Topology.Nodes]]
        self.Results.MeanVelocityPerNode = [[node.Velocity for node in self.Topology.Nodes]]
        self.Results.MeanVolumeFlowRatePerNode = [[node.FlowRate for node in self.Topology.Nodes]]
        self.Results.MeanRadiusPerNode = [[node.Radius for node in self.Topology.Nodes]]

        self.Results.PulsatilityIndexPressure = [[0 for _ in self.Topology.Nodes]]
        self.Results.PulsatilityIndexVelocity = [[0 for _ in self.Topology.Nodes]]

        vesselnames = [vessel.Name for vessel in self.Topology.Vessels]
        vesselnodes = [[node.Number for node in vessel.Nodes] for vessel in self.Topology.Vessels]

        self.Results.MeanResults = []
        for index, vessel in enumerate(vesselnodes):
            meanvolumeflowrate = numpy.mean(
                [self.Results.MeanVolumeFlowRatePerNode[-1][nodenumber] for nodenumber in vessel])
            meanpressure = numpy.mean([self.Results.MeanPressurePerNode[-1][nodenumber] for nodenumber in vessel])
            meanradius = numpy.mean([self.Results.MeanRadiusPerNode[-1][nodenumber] for nodenumber in vessel])
            meanvelocity = numpy.mean([self.Results.MeanVelocityPerNode[-1][nodenumber] for nodenumber in vessel])
            meanvalues = (meanvolumeflowrate, meanpressure, meanradius, meanvelocity, 0, 0)
            self.Results.MeanResults.append((vesselnames[index], meanvalues))

    def Export1DSteadyClusterFlowRate(self,file="ClusterFlowData.csv"):
        with open(self.Folders.ModellingFolder + file, "w") as f:
            f.write("Cluster ID,Flow rate over the region(mL/s),number of triangles,Pressure(Pa)\n")
            for index, cp in enumerate(self.Perfusion.CouplingPoints):
                f.write("%d,%.10f,%d,%.10f\n" % (index, cp.Node.FlowRate, cp.NumberOfTriangles, cp.Node.Pressure))
