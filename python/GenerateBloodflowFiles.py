#!/usr/bin/python3
# -*- coding:utf-8 -*-
"""1D Blood Flow Simulator

Usage:
  tavernaBloodFlow.py <patient_folder>
  tavernaBloodFlow.py (-h | --help)

Options:
  -h --help     Show this screen.
"""
from datetime import datetime

import sys
from os.path import isfile, dirname, realpath
sys.path.append(realpath(dirname(__file__)))

import GeneralFunctions
import Patient as PatientClass
import transcript
import Remesh
from docopt import docopt

"""
Script to run a blood flow simulation.
Input argument is the patient folder with a folder for input files and a folder for modelling files
The input folder contains patient data and the modelling file folder contains files for the models such as the parameters and the surface mesh.
"""

def generatebloodflowfiles(patient_folder):
    start_time = datetime.now()
    transcript.start(patient_folder + 'logfile.log')

    # Load input files
    Patient = PatientClass.Patient(patient_folder)
    Patient.ResetModellingFolder()
    Patient.SelectPatient()
    Patient.LoadSegmentedVessels()  # This creates a new file, this is temporary until there is a module 2.
    Patient.Load1DAnatomy()

    # Patient.Load1DAnatomy("1-D_Anatomy.txt")
    Patient.LoadPatientData()
    Patient.LoadModelParameters()
    Patient.UpdateModelParameters()
    Patient.SelectDonorNetwork()
    # Mapping to system vessels
    # The other vessels are mapped in the merge function
    mappingfile = Patient.Folders.ScriptFolder+"/DefaultFiles/FinalSystemVessels.vtp"
    GeneralFunctions.VesselMapping(mappingfile, Patient)

    # Load selected donor network
    Brava = PatientClass.Patient(patient_folder)
    DonorFile = Brava.Folders.ScriptFolder + "/Brava/" + Patient.ModelParameters["Donor_Network"]
    Brava.LoadVTPFile(DonorFile)

    # Load the pial surface
    surfacefile = Patient.Folders.ModellingFolder + "PialSurface.vtp"
    if not GeneralFunctions.is_non_zero_file(surfacefile):
        # if file does not exist, create a new one from the .ply file.
        # remesh the mesh to a uniform triangulation.
        print("Pial Surface file not found.")
        Remesh.remesh(Patient.Folders.ModellingFolder+"boundary_4&21&22&23&24&25&26&30.ply", numbertriangles=40000, output=Patient.Folders.ModellingFolder+"remeshed.vtp")
        # apply the same mapping to the remeshed file
        Remesh.MapMeshtoMSH(Patient.Folders.ModellingFolder+"remeshed.vtp", Patient.Folders.ModellingFolder+"labelled_vol_mesh.msh",output=Patient.Folders.ModellingFolder + "PialSurface.vtp")

    Brava.Perfusion.PrimalGraph.LoadSurface(surfacefile)
    Brava.Perfusion.regionids = [2, 3, 4, 5, 6, 7, 8, 9]
    Brava.VesselToMeshAllignmentSides(Brava.Perfusion.PrimalGraph)
    Brava.Perfusion.SetDualGraph(method="vertices")
    Brava.Perfusion.RemapMajorRegions()

    GeneralFunctions.Merge(Patient, Brava)
    del Brava
    Patient.TopologyToVTP()
    Patient.CerebellumBrainstemMapping()
    # Calculate mapping
    # Patient.MajorArteriesToPialSurfaceNN(Patient.Perfusion.DualGraph)
    # mapping from VEASL data
    Patient.Perfusion.DualGraph.MajorVesselID = Patient.Perfusion.PrimalGraph.map

    Patient.FindCouplingPoints()
    # Patient.GenerateTreesAtCps(0.20)  # cutoff strongly scales the number of coupling points.
    # Patient.PerfusionEndsNumber()
    # Patient.MappingSixMajorArteries()

    Patient.Perfusion.ClusteringByRegionScalingLaw(Patient.Perfusion.DualGraph, method="", fractiontriangles=(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.5), debug=False)  # default is the Dijkstra distance
    Patient.Perfusion.MapDualGraphToPrimalGraph()
    Patient.ImportClots()
    Patient.CalculateMaximumTimestep()

    # write all files needed for the BF simulation
    Patient.WriteModelParameters()
    Patient.WriteSimFiles()
    Patient.Perfusion.PrimalGraph.ExportTriangleColour(Patient.Folders.ModellingFolder)
    Patient.Perfusion.WriteClusteringMapping(Patient.Folders.ModellingFolder)
    Patient.ExportSurface(Patient.Folders.ModellingFolder + "Clustering.vtp", Patient.Perfusion.PrimalGraph)

    # optional export functions
    Patient.Perfusion.PrimalGraph.GraphToVTP(Patient.Folders.ModellingFolder)
    Patient.Perfusion.DualGraph.GraphToVTP(Patient.Folders.ModellingFolder)

    # map the clustering back to the msh file.
    file1 = Patient.Folders.ModellingFolder + "Clustering.vtp"
    filehighres = Patient.Folders.ModellingFolder + "labelled_vol_mesh.msh"
    GeneralFunctions.MapClusteringToMSH(file1, filehighres, Patient.Folders.ModellingFolder)

    time_elapsed = datetime.now() - start_time
    print('Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))
    transcript.stop()


if __name__ == '__main__':
    arguments = docopt(__doc__, version='0.1')
    patient_folder = arguments["<patient_folder>"]
    # patient_folder = "/drive/1d-blood-flow/Generated_Patients/patient_test/"
    generatebloodflowfiles(patient_folder)
