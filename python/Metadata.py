#!/usr/bin/python3
import os
import sys
from os.path import isfile, dirname, realpath, exists
from lxml import etree
# sys.path.insert(0, "/drive/in-silico-trial/software/1d-blood-flow/scripts/")

class Folders:
    def __init__(self, patientfolder):
        self.PatientFolder = patientfolder
        self.InputFolder = patientfolder
        self.ModellingFolder = patientfolder + "bf_sim/"
        self.ScriptFolder = os.path.dirname(realpath(dirname(__file__))) + "/scripts/"

    def SetPatientFolder(self, patientfolder):
        self.PatientFolder = patientfolder
        self.InputFolder = patientfolder
        self.ModellingFolder = patientfolder + "bf_sim/"


class PatientData(dict):
    def __init__(self):
        dict.__init__(self)
        self.File = ""
        self["Age"] = 50  # years
        self["Sex"] = "Male"
        # self["BMI"] = 20  # kg/m^2
        self["HeartRate"] = 60  # per minute
        self["SystolePressure"] = 17300 # Pa
        self["DiastolePressure"] = 10100 # Pa
        self["MeanRightAtrialPressure"] = 0 # Pa
        self["StrokeVolume"] = 95 #mL/s

    def LoadPatientData(self, file):
        print("Loading %s" % file)
        self.File = file
        data = [line.strip('\n').split('=') for line in open(file)]
        for line in data:
            try:
                datavalue = int(line[1])
            except ValueError:
                try:
                    datavalue = float(line[1])
                except ValueError:
                    datavalue = line[1]
            self[line[0]] = datavalue

    def LoadPatientDataXML(self, xml_file):
        print("Loading %s" % xml_file)
        try:
            f = open(xml_file, "rb+")
        except IOError:
            print("File open error")
            exit(1)

        root_xml = etree.parse(f)
        patient_xml = root_xml.find("Patient")
        self.File = xml_file
        self["Age"] = patient_xml.find("Age").text  # years
        self["Sex"] = patient_xml.find("Sex").text
        self["HeartRate"] = float(patient_xml.find("HeartRate").text)  # per minute
        self["SystolePressure"] = float(patient_xml.find("SystolePressure").text)  # Pa
        self["DiastolePressure"] = float(patient_xml.find("DiastolePressure").text)  # Pa
        self["MeanRightAtrialPressure"] = float(patient_xml.find("MeanRightAtrialPressure").text)  # Pa
        self["StrokeVolume"] = float(patient_xml.find("StrokeVolume").text)  # mL/s

    def WritePatientData(self, file="Patient_parameters.txt"):
        print("Writing: %s" % file)
        with open(file, 'w') as f:
            for key, val in self.items():
                f.write(key + "=" + str(val) + "\n")


class ModelParameter(dict):
    def __init__(self):
        dict.__init__(self)
        self["RTotal"] = 1.19e8  # N S^-1 m^-5
        self["CTotal"] = 2.38e-9  # m^5 N^-1
        self["Density"] = 1040  # kg M^-3
        # self["YoungsModules"] = 225e3  # Pa
        # self["Alpha"] = 1.1

    def LoadModelParameters(self, file):
        print("Loading %s" % file)
        data = [line.strip('\n').split('=') for line in open(file)]
        for line in data:
            try:
                datavalue = int(line[1])
            except ValueError:
                try:
                    datavalue = float(line[1])
                except ValueError:
                    datavalue = line[1]
            self[line[0]] = datavalue

    def WriteModelParameters(self, file="Model_parameters.txt"):
        print("Writing: %s" % file)
        with open(file, 'w') as f:
            for key, val in self.items():
                f.write(key + "=" + str(val) + "\n")
