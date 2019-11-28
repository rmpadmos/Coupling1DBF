#!/usr/bin/python3
import math
import numpy
from scipy.interpolate import interpolate
import sys
from os.path import isfile, dirname, realpath, exists
sys.path.append(realpath(dirname(__file__)) )
import Node
import BloodflowEquations

class Vessel:
    def __init__(self):
        self.Nodes = []
        self.NumberOfNodes = 0
        self.Length = None
        self.MeanRadius = None
        self.MeanThickness = None
        self.GridSize = 1
        self.Name = ""
        self.ID = None
        self.MajorVesselID = -1
        self.YoungsModules = 1.6e6
        self.Resistance = None
        self.Compliance = None
        self.InterpolationFunctions = []
        self.GenerationNumber = None

    def GetDistalBifurcation(self):
        """
        :return: Bifurcation node at the distal end of the vessel.
        """
        nodes = self.Nodes[-1].Connections
        for node in nodes:
            if node.Type == 1:
                return node
        return None

    def GetProximalBifurcation(self):
        """
        :return: Bifurcation node at the proximal end of the vessel.
        """
        nodes = self.Nodes[0].Connections
        for node in nodes:
            if node.Type == 1:
                return node
        return None

    def GetEndNodes(self):
        proximal = self.GetProximalBifurcation()
        if proximal is None:
            proximal = self.Nodes[0]

        distal = self.GetDistalBifurcation()
        if distal is None:
            distal = self.Nodes[-1]

        return [proximal, distal]

    def UpdateNodeVesselID(self):
        [node.SetVesselID(self.ID) for node in self.Nodes]
        [node.SetMajorVesselID(self.MajorVesselID) for node in self.Nodes]

    def UpdateNodeNumber(self):
        self.NumberOfNodes = len(self.Nodes)

    def SetID(self, idnumber):
        self.ID = idnumber

    def SetName(self, name):
        self.Name = name

    def UpdateNodeDirectionVectors(self):
        [node.SetDirectionVector() for node in self.Nodes]

    def SetNodes(self, nodes):
        self.Nodes = nodes
        self.NumberOfNodes = len(nodes)

    def SetLength(self, length):
        self.Length = length

    def SetMeanRadius(self, rad):
        self.MeanRadius = rad
        self.MeanThickness = BloodflowEquations.thickness(rad)

    def SetMeanThickness(self, h):
        self.MeanThickness = h

    def SetGridSize(self, gridsize):
        self.GridSize = gridsize

    def SetMajorVesselID(self, majorid):
        self.MajorVesselID = majorid
        for node in self.Nodes:
            node.SetMajorVesselID(majorid)

    def ScaleRadius(self, r):
        [node.SetRadius(node.Radius * r) for node in self.Nodes]
        self.MeanRadius *= r
        self.InterpolationFunctions[3].y *= r
        self.MeanThickness = BloodflowEquations.thickness(self.MeanRadius)

    def UpdateVessel(self):
        self.UpdateNodeVesselID()
        self.UpdateNodeNumber()

    def CalculateMeanRadius(self):
        self.MeanRadius = numpy.trapz([i.Radius for i in self.Nodes],
                                      [i.LengthAlongVessel for i in self.Nodes]) / self.Length
        return self.MeanRadius

    def CalculateMeanThickness(self):
        self.MeanThickness = numpy.trapz([i.Thickness for i in self.Nodes],
                                         [i.LengthAlongVessel for i in self.Nodes]) / self.Length

    def VesselResistance(self, bloodvisc):
        length = self.Length * 1e-3
        radius = self.MeanRadius * 1e-3
        # R = 8 * bloodvisc * length / (numpy.pi * numpy.power(radius, 4))
        R = 22 * bloodvisc * length / (numpy.pi * numpy.power(radius, 4))
        # func = [22 * bloodvisc * self.GridSize * 1e-3 / (numpy.pi * numpy.power(n.Radius* 1e-3, 4)) for n in self.Nodes]
        # R2 = sum(func)

        self.Resistance = R
        return R

    def VesselCompliance(self):
        C1D = 2 * numpy.power(numpy.pi * self.MeanRadius * 1e-3 * self.MeanRadius * 1e-3, 1.5) * self.Length * 1e-3 / (
                (4 / 3) * numpy.sqrt(numpy.pi) * self.YoungsModules * self.MeanThickness * 1e-3)
        self.Compliance = C1D
        return C1D

    def GenerateVessel(self, length, inletradius, outletradius, Elastic):
        numbernodes = int(max(3, math.ceil(length / self.GridSize) + 1))
        nodelist = [Node.Node() for i in range(0, numbernodes)]

        dx = length / (numbernodes - 1)
        if length == 0:
            self.Length = 0
            return
        drdl = (outletradius - inletradius) / length

        for i in range(0, numbernodes):
            nodelist[i].SetLengthAlongVessel(dx * i)
            nodelist[i].SetPosition([0, 0, dx * i])
            nodelist[i].SetRadius(inletradius + drdl * dx * i)

        for i in range(0, numbernodes - 1):
            nodelist[i].AddConnection(nodelist[i + 1])
            nodelist[i + 1].AddConnection(nodelist[i])

        [node.SetYoungsModules(Elastic) for node in nodelist]
        self.Nodes = nodelist
        self.NumberOfNodes = numbernodes
        self.YoungsModules = Elastic
        self.GridSize = length / (numbernodes - 1)
        self.Length = length
        self.MeanRadius = (inletradius + outletradius) / 2
        self.MeanThickness = BloodflowEquations.thickness(self.MeanRadius)

    def InterpolateVessel3Dto1D(self):
        # print("Interpolating Vessel.")
        x = [i.Position[0] for i in self.Nodes]
        y = [i.Position[1] for i in self.Nodes]
        z = [i.Position[2] for i in self.Nodes]
        r = [i.Radius for i in self.Nodes]
        id = self.Nodes[1].MajorVesselID
        vesselid = self.Nodes[1].VesselID
        npts = len(x)
        s = numpy.zeros(npts, dtype=float)
        for j in range(1, npts):
            dx = x[j] - x[j - 1]
            dy = y[j] - y[j - 1]
            dz = z[j] - z[j - 1]
            vec = numpy.array([dx, dy, dz])
            s[j] = s[j - 1] + numpy.linalg.norm(vec)

        interpolation = 'linear'
        # Create new interpolation function for each dimension against the norm
        f1 = interpolate.interp1d(s, x, kind=interpolation)
        f2 = interpolate.interp1d(s, y, kind=interpolation)
        f3 = interpolate.interp1d(s, z, kind=interpolation)
        f4 = interpolate.interp1d(s, r, kind=interpolation)
        meanradius = numpy.trapz(r, s) / s[-1]

        gridnodes = max(3, math.ceil(s[-1] / self.GridSize) + 1)
        xvec = numpy.linspace(s[0], s[-1], gridnodes)
        radius = f4(xvec)

        position = [xvec[i] for i in range(0, gridnodes)]
        for r in radius:
            if r < 0:
                print("Radius below zero, interpolation is off.")
        position3d = [[f1(pos), f2(pos), f3(pos)] for pos in xvec]

        nodes = [Node.Node() for i in range(0, gridnodes)]
        [nodes[i].SetLengthAlongVessel(position[i]) for i in range(0, gridnodes)]
        [nodes[i].SetRadius(meanradius) for i in range(0, gridnodes)]
        [nodes[i].SetYoungsModules(self.YoungsModules) for i in range(0, gridnodes)]
        [nodes[i].SetPosition(position3d[i]) for i in range(0, gridnodes)]
        [nodes[i].SetMajorVesselID(id) for i in range(0, gridnodes)]
        [nodes[i].SetVesselID(vesselid) for i in range(0, gridnodes)]
        self.InterpolationFunctions = [f1, f2, f3, f4]
        self.Nodes = nodes
        self.Length = s[-1]
        self.MeanRadius = meanradius
        self.NumberOfNodes = len(nodes)
        self.GridSize = self.Length / (self.NumberOfNodes - 1)
        self.MajorVesselID = id
        self.ID = vesselid
        for i in range(0, len(self.Nodes) - 1):
            self.Nodes[i].AddConnection(self.Nodes[i + 1])
            self.Nodes[i + 1].AddConnection(self.Nodes[i])
        self.MeanThickness = BloodflowEquations.thickness(self.MeanRadius)

    def CreateInterpolationFunctions(self):
        x = [i.Position[0] for i in self.Nodes]
        y = [i.Position[1] for i in self.Nodes]
        z = [i.Position[2] for i in self.Nodes]
        r = [i.Radius for i in self.Nodes]

        npts = len(x)
        s = numpy.zeros(npts, dtype=float)
        for j in range(1, npts):
            dx = x[j] - x[j - 1]
            dy = y[j] - y[j - 1]
            dz = z[j] - z[j - 1]
            vec = numpy.array([dx, dy, dz])
            s[j] = s[j - 1] + numpy.linalg.norm(vec)

        self.Length = s[-1]

        interpolation = 'linear'
        # Create new interpolation function for each dimension against the norm
        f1 = interpolate.interp1d(s, x, kind=interpolation)
        f2 = interpolate.interp1d(s, y, kind=interpolation)
        f3 = interpolate.interp1d(s, z, kind=interpolation)
        f4 = interpolate.interp1d(s, r, kind=interpolation)
        self.InterpolationFunctions = [f1, f2, f3, f4]

    def UpdateInterpolationFunctions(self):
        x = [i.Position[0] for i in self.Nodes]
        y = [i.Position[1] for i in self.Nodes]
        z = [i.Position[2] for i in self.Nodes]
        r = [i.Radius for i in self.Nodes]
        s = [i.LengthAlongVessel for i in self.Nodes]

        interpolation = 'linear'
        # Create new interpolation function for each dimension against the norm
        f1 = interpolate.interp1d(s, x, kind=interpolation)
        f2 = interpolate.interp1d(s, y, kind=interpolation)
        f3 = interpolate.interp1d(s, z, kind=interpolation)
        f4 = interpolate.interp1d(s, r, kind=interpolation)
        self.InterpolationFunctions = [f1, f2, f3, f4]

    def Interpolate3D(self, lengthalongvessel):
        if len(self.InterpolationFunctions) == 0:
            self.CreateInterpolationFunctions()

        [f1, f2, f3, f4] = self.InterpolationFunctions
        return (
            numpy.vstack((f1(lengthalongvessel), f2(lengthalongvessel), f3(lengthalongvessel))), f4(lengthalongvessel))

    def CalculateMaxWaveSpeed(self, density):
        # find the largest wavespeed in the vessel
        wavespeeds = []
        for node in self.Nodes:
            A = numpy.pi * node.Radius * 1e-3 * node.Radius * 1e-3
            beta = (4 / 3) * math.sqrt(math.pi) * node.YoungsModules * node.Thickness * 1e-3 / A
            c0 = abs(math.sqrt(beta / (2 * density))) * math.pow(A, 0.25)
            wavespeeds.append(c0)
        maxwavespeed = max(wavespeeds)
        return maxwavespeed  # m/s

    def UpdateResolution(self, numbernodes):
        # note that 2 is the minimum
        if numbernodes < 2:
            print("Error: Minimum number of nodes allowed per vessel is 2.")
            return 0

        self.GridSize = self.Length / (numbernodes - 1)
        xvec = numpy.linspace(0, self.Length, numbernodes)
        xvec = [min(i, self.Length) for i in xvec]  # 12.00000001 is above the interpolation limit if length is 12.
        Interpolationvalues = self.Interpolate3D(xvec)

        # remove old connections
        for i in range(0, len(self.Nodes) - 1):
            self.Nodes[i].RemoveConnection(self.Nodes[i + 1])
            self.Nodes[i + 1].RemoveConnection(self.Nodes[i])

        # create new nodes
        newnodes = [self.Nodes[0]]
        for i in range(1, numbernodes - 1):
            newnode = Node.Node()
            newnode.SetRadius(Interpolationvalues[1][i])
            newnode.SetLengthAlongVessel(xvec[i])
            newnode.SetPosition(Interpolationvalues[0][:, i])
            newnode.SetYoungsModules(self.YoungsModules)
            newnode.SetVesselID(self.ID)
            newnode.SetMajorVesselID(self.MajorVesselID)
            newnodes.append(newnode)
        newnodes.append(self.Nodes[-1])
        self.SetNodes(newnodes)

        # add connections
        for i in range(0, len(self.Nodes) - 1):
            self.Nodes[i].AddConnection(self.Nodes[i + 1])
            self.Nodes[i + 1].AddConnection(self.Nodes[i])
