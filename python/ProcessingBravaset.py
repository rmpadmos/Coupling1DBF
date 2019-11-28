#!/usr/bin/python3
import contextlib
import os
#import sys
# sys.path.insert(0, "/drive/in-silico-trial/software/1d-blood-flow/scripts/")
import numpy

import GeneralFunctions
import Patient as PatientModule
from Node import *


def SWC_Processing1d(filename):
    """
    Convert a .swc file to a .vtp file
    Note that this assumes a particular order and meaning of the columns
    """
    print("Converting swc file.")
    newname = filename[:-3] + "vtp"
    dataset = [text.split() for text in open(filename)]

    number = [int(line[0]) for line in dataset]
    colour = [int(line[1]) for line in dataset]
    positions = [[float(line[2]), float(line[3]), float(line[4])] for line in dataset]
    radius = [float(line[5]) for line in dataset]
    connection = [int(line[6]) for line in dataset]

    links = [[] for node in range(0, len(number))]
    for i in range(0, len(number)):
        if connection[i] > 0:
            nodenumber = number[i] - 1
            linkedto = connection[i] - 1
            links[nodenumber].append(linkedto)
            links[linkedto].append(nodenumber)

    nodes = [Node() for node in range(0, len(number))]
    [nodes[index].SetRadius(rad) for index, rad in enumerate(radius)]
    [nodes[index].SetPosition(pos) for index, pos in enumerate(positions)]
    [nodes[index].SetMajorVesselID(c) for index, c in enumerate(colour)]

    for index, link in enumerate(links):
        for con in link:
            nodes[index].AddConnection(nodes[con])

    patient = PatientModule.Patient()
    patient.Topology.Nodes = nodes
    patient.Topology.AnatomyToVessels()
    vesseltype = [colour[patient.Topology.Nodes.index(vessel.Nodes[1])] for vessel in patient.Topology.Vessels]
    patient.Topology.VesselAtlas = vesseltype

    pos = [node.Position for node in patient.Topology.Nodes]
    meanx = numpy.mean([p[0] for p in pos])
    meany = numpy.mean([p[1] for p in pos])
    meanz = numpy.mean([p[2] for p in pos])
    centermatrix = GeneralFunctions.TMatrix(1, [90, 0, 0], [-meanx, -meany, -meanz])
    patient.Topology.ApplyTransformation(centermatrix)

    patient.Topology.TopologyToVTP(newname)
    return patient


def Dataset_SWC_To_VTP(path):
    """
    Convert all files in a folder to .vtp files.
    This assumes all files are .swc
    """
    print("Converting all files from swc to vtp.")
    filenames = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

    for name in filenames:
        if name[-3:] == "swc":
            print("Converting: %s" % name)
            with contextlib.redirect_stdout(None):
                SWC_Processing1d(path + name)


def ConvertBravaSet(path):
    Dataset_SWC_To_VTP(path)

    # one file seems to be rotated with respect to the others
    with contextlib.redirect_stdout(None):
        file = path + "Set8_ColorCoded.CNG.vtp"
        patient = PatientModule.Patient()
        patient.LoadVTPFile(file)
        totalmatrix = GeneralFunctions.TMatrix(1, [0, 0, 180], [0, 0, 0])
        GeneralFunctions.TransformFile(file, totalmatrix)


if __name__ == '__main__':
    path = "/home/raymond/temp/"
    ConvertBravaSet(path)
