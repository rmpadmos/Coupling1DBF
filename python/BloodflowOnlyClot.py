#!/usr/bin/python3
# -*- coding:utf-8 -*-
"""1D Blood Flow Simulator

Usage:
  tavernaBloodFlow.py <executable_location> <patient_folder> <clot_present>
  tavernaBloodFlow.py (-h | --help)

Options:
  -h --help     Show this screen.
"""
import os
import subprocess
import sys
from os.path import isfile, dirname, realpath, exists
sys.path.append(realpath(dirname(__file__) +"/scripts") )
from datetime import datetime

import GeneralFunctions
import Patient as PatientModule
import transcript
from docopt import docopt
import Results
import BloodflowEquations
"""
Script to run a blood flow simulation.
Input argument is the patient folder with a folder for input files and a folder for modelling files
The input folder contains patient data and the modelling file folder contains files for the models such as the parameters and the surface mesh.
"""


def blood_flow_script(executable, patient_folder, clot_present):
    Patient = PatientModule.Patient(patient_folder)
    # remove result files if these were created during a previous run
    Patient.RemoveOldSimFiles()
    Patient.LoadBFSimFiles()

    # if not os.path.exists(resultsfolder):
    #     os.mkdir(resultsfolder)
    model = "Pulsatile"
    # model = "Steady"
    resultsfolder = Patient.Folders.ModellingFolder
    if model == "Pulsatile":
        if clot_present == "true" or clot_present == "True":
            subprocess.call(["mono", executable, resultsfolder + "Run.txt",
                             resultsfolder + "Results.dyn", resultsfolder + "Clots.txt"])
        else:
            subprocess.call(["mono", executable, resultsfolder + "Run.txt",
                             resultsfolder + "Results.dyn"])
        Patient.LoadResults("Results.dyn", correct=False)
        # BloodflowEquations.Bloodflow1D(Patient)
    else:
        Patient.Initiate1DSteadyStateModel()
        if clot_present == "True" or clot_present == "true":
            Patient.Run1DSteadyStateModel(model="Linear", clotactive=True)
        else:
            Patient.Run1DSteadyStateModel(model="Linear", clotactive=False)

        # export data
        TimePoint = Results.TimePoint(0)
        TimePoint.Flow = [node.FlowRate for node in Patient.Topology.Nodes]
        TimePoint.Pressure = [node.Pressure for node in Patient.Topology.Nodes]
        TimePoint.Radius = [node.Radius for node in Patient.Topology.Nodes]

        TimePoint2 = Results.TimePoint(Patient.ModelParameters['Beat_Duration'])
        TimePoint2.Flow = TimePoint.Flow
        TimePoint2.Pressure = TimePoint.Pressure
        TimePoint2.Radius = TimePoint.Radius

        Patient.Results.TimePoints = [TimePoint, TimePoint2]
        Patient.Results.ExportResults(resultsfolder+"Results.dyn")

        Patient.LoadResults("Results.dyn")

    Patient.GetMeanResults()
    Patient.ExportMeanResults()

    # depending on triangles or vertexes
    # Patient.DistributeFlowVertex()
    Patient.DistributeFlowTriangles()
    Patient.ExportTriangleFlowData()

    Patient.WriteTimeseriesVessels()
    Patient.Results.AddResultsPerNodeToFile(Patient.Folders.ModellingFolder + "Topology.vtp")
    Patient.Results.AddResultsPerVesselToFile(Patient.Folders.ModellingFolder + "Topology.vtp")

    # add results to the mesh used in the perfusion model
    # file1 = Patient.Folders.ModellingFolder + "FlowDistributedMean.vtp"
    # clusteringflowdata = Patient.Folders.ModellingFolder + "ClusterFlowData.csv"
    # filehighres = Patient.Folders.ModellingFolder + "labelled_vol_mesh.msh"
    # GeneralFunctions.MapClusteringToMSH(file1, filehighres, clusteringflowdata, Patient.Folders.ModellingFolder)

    # update the flow to the perfusion model
    clusteringflowdata = Patient.Folders.ModellingFolder + "ClusterFlowData.csv"
    clusteringfile = Patient.Folders.ModellingFolder + "Clusters.csv"
    datafolder = Patient.Folders.ModellingFolder
    GeneralFunctions.WriteFlowFilePerfusionModel(clusteringflowdata, clusteringfile, datafolder)

if __name__ == '__main__':
    arguments = docopt(__doc__, version='0.1')
    executable = arguments["<executable_location>"]
    patient_folder = arguments["<patient_folder>"]
    clot_present = arguments["<clot_present>"]

    start_time = datetime.now()
    transcript.start(patient_folder + 'logfile.log')
    blood_flow_script(executable, patient_folder, clot_present)
    time_elapsed = datetime.now() - start_time
    print('Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))
    transcript.stop()
