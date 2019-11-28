#!/usr/bin/python3
import os
import sys

# sys.path.insert(0, "/drive/in-silico-trial/software/1d-blood-flow/scripts/")
import numpy

import Results as ResultClass

if len(sys.argv) > 1:
    Folder = str(sys.argv[1])
else:
    Folder = "/drive/in-silico-trial/software/1d-blood-flow/patient1/Modelling_files/"

ResultsFile = "Results.dyn"
ResultsFile2 = "ResultsPrev.dyn"
ResultsFile3 = "ResultsTotal.dyn"
Conv = 0


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


if is_non_zero_file(Folder + ResultsFile) and is_non_zero_file(Folder + ResultsFile2):
    Results = ResultClass.Results()
    Results.LoadResults(Folder + ResultsFile)

    Results2 = ResultClass.Results()
    Results2.LoadResults(Folder + ResultsFile2)

    numberofnodes = len(Results.TimePoints[0].Flow)
    # pressureres = 0
    # flowres = 0
    # radiusres = 0
    #
    # pressuremaxres = 0.0
    # flowresmaxres = 0.0
    # radiusresmaxres = 0.0
    #
    # maxrtol = 0.0
    #
    # for index, timepoint in enumerate(Results.TimePoints):
    #     comparetimepoint = Results2.TimePoints[index]
    #
    #     for i in range(0, numberofnodes):
    #         pressureres = abs((timepoint.Pressure[i] - comparetimepoint.Pressure[i])) / (
    #                 1 + min(abs(timepoint.Pressure[i]), abs(comparetimepoint.Pressure[i])))
    #         flowres = abs((timepoint.Flow[i] - comparetimepoint.Flow[i])) / (
    #                 1 + min(abs(timepoint.Flow[i]), abs(comparetimepoint.Flow[i])))
    #         radiusres = abs((timepoint.Radius[i] - comparetimepoint.Radius[i])) / (
    #                 1 + min(abs(timepoint.Radius[i]), abs(comparetimepoint.Radius[i])))
    #         pressuremaxres = max(pressuremaxres, pressureres)
    #         flowresmaxres = max(flowresmaxres, flowres)
    #         radiusresmaxres = max(radiusresmaxres, radiusres)

    p1 = numpy.array([item for sublist in Results.Pressure[-1] for item in sublist])
    p2 = numpy.array([item for sublist in Results2.Pressure[-1] for item in sublist])
    p1_norm = numpy.linalg.norm(p1)
    pressurenorm = numpy.linalg.norm(p1 - p2) / p1_norm

    v1 = numpy.array([item for sublist in Results.VolumeFlowRate[-1] for item in sublist])
    v2 = numpy.array([item for sublist in Results2.VolumeFlowRate[-1] for item in sublist])
    v1_norm = numpy.linalg.norm(v1)
    flownorm = numpy.linalg.norm(v1 - v2) / v1_norm

    r1 = numpy.array([item for sublist in Results.Radius[-1] for item in sublist])
    r2 = numpy.array([item for sublist in Results2.Radius[-1] for item in sublist])
    r1_norm = numpy.linalg.norm(r1)
    radiusnorm = numpy.linalg.norm(r1 - r2) / r1_norm

    print("Max Residual: Pressure: %f  Flow rate: %f Radius: %f" % (pressurenorm, flownorm, radiusnorm))

    # print("Max Residual: Pressure: %f  Flow rate: %f Radius: %f" % (pressuremaxres, flowresmaxres, radiusresmaxres))
    # if pressureres < 1e-3 and flowres < 1e-3 and velores < 1e-3:
    with open(Folder + "Convergence.csv", "a") as f:
        f.write("%f,%f,%f\n" % (pressurenorm, flownorm, radiusnorm))

    if pressurenorm < 1e-3:
        Conv = 1

with open(Folder + ResultsFile) as f:
    with open(Folder + ResultsFile3, "a") as f1:
        for line in f:
            f1.write(line)

if Conv == 0:
    with open(Folder + ResultsFile) as f:
        with open(Folder + ResultsFile2, "w") as f1:
            for line in f:
                f1.write(line)

with open(Folder + "Conv.txt", "w") as f:
    f.write(str(Conv))
exit()
