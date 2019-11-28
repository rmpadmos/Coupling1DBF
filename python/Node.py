#!/usr/bin/python3
import math
import sys
from os.path import isfile, dirname, realpath, exists
sys.path.append(realpath(dirname(__file__)) )
import BloodflowEquations
import numpy

class Node:
    def __init__(self):
        self.Number = None
        self.Position = [0, 0, 0]
        self.Radius = 0
        self.Connections = set()
        self.R1 = None
        self.R2 = None
        self.C = None
        self.YoungsModules = 1.6e6
        self.Type = 0
        self.VesselID = -1
        self.MajorVesselID = -1
        self.LengthAlongVessel = 0
        self.Thickness = None
        self.DirectionVector = [0, 0, 0]

        self.Pressure = 0
        self.Velocity = 0
        self.FlowRate = 0
        self.PressureRadiusEquation = None
        self.RefRadius = 0
        self.RefPressure = 0

    def flux_function(self, Area):
        # chrt_back = -4 * numpy.power(self.Lumen0, 0.25) * numpy.sqrt(self.Beta / 2.0 / self.density)
        return Area * self.charback + 4 * numpy.power(Area, 1.25) * numpy.sqrt(
            self.Beta / 2.0 / self.Density) - self.CurrentFlowrate

    def chrt_function(self, Area):
        chrt_frw_right = self.CurrentVolumeFlowRate / Area + 4 * numpy.power(Area, 0.25) * numpy.sqrt(
            self.Beta / 2.0 / self.Density)
        return self.chrt_frw_left - chrt_frw_right

    def SetPressureAreaEquation(self):
        function = lambda pressure: (pressure - self.RefPressure) * self.RefRadius * self.RefRadius * 3 / (
                    4 * self.YoungsModules * self.Thickness) + self.RefRadius
        self.PressureRadiusEquation = function

    def UpdateRadius(self):
        self.Radius = self.PressureRadiusEquation(self.Pressure)

    def SetNumber(self, number):
        self.Number = number

    def SetYoungsModules(self, e):
        self.YoungsModules = e

    def SetLengthAlongVessel(self, l):
        self.LengthAlongVessel = l

    def SetRadius(self, r):
        self.Radius = r
        self.CalculateThickness()

    def SetPosition(self, pos):
        self.Position = []
        for p in pos:
            self.Position.append(float(p))

    def SetVesselID(self, _id):
        self.VesselID = _id

    def SetType(self, type):
        self.Type = type

    def CalculateThickness(self):
        self.Thickness = BloodflowEquations.thickness(self.Radius)
        return self.Thickness

    def SetMajorVesselID(self, _id):
        self.MajorVesselID = _id

    def GetConnectedVesselIDs(self):
        listIDS = []
        for node in self.Connections:
            listIDS.append(node.VesselID)
        return listIDS

    def GetConnectedMajorVesselIDs(self):
        listIDS = []
        for node in self.Connections:
            listIDS.append(node.MajorVesselID)
        return listIDS

    def AddConnection(self, node):
        self.Connections.add(node)

    def RemoveConnection(self, node):
        self.Connections.remove(node)

    def ResetConnections(self):
        self.Connections.clear()

    def SetWK(self, r1, r2, c):
        self.R1 = r1
        self.R2 = r2
        self.C = c

    def ResetWK(self):
        self.R1 = None
        self.R2 = None
        self.C = None

    def SetNodeFromTopLine(self, line):
        self.SetNumber(int(line[0]))
        self.SetLengthAlongVessel(float(line[1][2:]))
        self.SetRadius(float(line[2][2:]))
        self.SetYoungsModules(float(line[3][2:]) * 1e6)
        self.SetType(int(line[4][2:]))

    def SetDirectionVector(self):
        """
        Calculate the direction vector of a node
        Only valid for vessel nodes
        """
        if len(self.Connections) > 2:
            print("Error: Direction Vector not defined for bifurcations.")
            return 1
        if len(self.Connections) == 0:
            print("Error: Node has no connections.")
            return 1

        linkednode = max(list(self.Connections), key=lambda x: x.LengthAlongVessel)
        nbpos = linkednode.Position
        pos = self.Position
        direction = [nbpos[0] - pos[0], nbpos[1] - pos[1], nbpos[2] - pos[2]]
        length = math.sqrt((direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]))
        directionvector = [direction[0] / length, direction[1] / length, direction[2] / length]
        self.DirectionVector = directionvector
