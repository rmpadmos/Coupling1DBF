#!/usr/bin/python3
import copy
import itertools
import sys
from os.path import dirname, realpath

sys.path.append(realpath(dirname(__file__)))
import networkx as nx
import numpy
import scipy
import scipy.sparse
import scipy.spatial
import vtk
from vtk.util.numpy_support import vtk_to_numpy

import GeneralFunctions
import Node



class Surface:
    def __init__(self):
        self.PialSurface = []
        self.Graph = nx.Graph()
        self.Links = []
        self.Polygons = []
        self.Triangles = []
        self.SideInfo = []
        self.Areas = []  # per polygon
        self.File = None
        self.Weights = []

        self.Labels = []  # of points
        self.Points = []
        self.Roots = []
        self.NodeColour = []
        self.PolygonColour = []
        self.MajorVesselID = []

    def CalculateScipyGraph(self):
        self.ScipyGraph = nx.to_scipy_sparse_matrix(self.Graph)
        # d=nx.adjacency_matrix(self.Graph)
        # a = scipy.sparse.csgraph.dijkstra(d, directed=False, limit=300)
        # m=self.ScipyGraph.todense()
        # print(self.ScipyGraph.todense())
        # self.CalculateDijkstraMatrix()

    def GetSurfaceRegion(self, regionids):
        trianglesids = [i for i, element in enumerate(self.MajorVesselID) if element in regionids]
        triangles = [self.Triangles[i] for i in trianglesids]
        pointids = [y for x in triangles for y in x]
        points = list(set(pointids))
        pointspos = [self.PialSurface[i] for i in points]
        newtriangles = [[points.index(p) for p in triangle] for triangle in triangles]

        return newtriangles, pointspos

    def SeparatePialNetwork(self):
        file = "Test.vtp"
        vesselspolydata = vtk.vtkPolyData()
        nodes = vtk.vtkPoints()
        for i in self.PialSurface:
            nodes.InsertNextPoint(i)
        vesselspolydata.SetPoints(nodes)

        NodeColour = vtk.vtkIntArray()
        NodeColour.SetNumberOfComponents(1)
        NodeColour.SetName("Node Colour")
        for colour in self.NodeColour:
            NodeColour.InsertNextValue(colour)
        vesselspolydata.GetPointData().AddArray(NodeColour)

        # removing edges
        edges = []
        for index, link in enumerate(self.Links):
            for linkednode in link:
                if self.NodeColour[index] == self.NodeColour[linkednode]:
                    edges.append((index, linkednode))

        lines = vtk.vtkCellArray()
        for index, element in enumerate(edges):
            line0 = vtk.vtkLine()
            line0.GetPointIds().SetNumberOfIds(2)
            line0.GetPointIds().SetId(0, element[0])
            line0.GetPointIds().SetId(1, element[1])
            lines.InsertNextCell(line0)
        vesselspolydata.SetLines(lines)

        # removing polygons
        # polygons = vtk.vtkCellArray()
        # for poly in self.PolygonsList:
        #     nodeides = [self.NodeColour[p] for p in poly]
        #     if len(set(nodeides)) == 1:
        #         polygon = vtk.vtkPolygon()
        #         polygon.GetPointIds().SetNumberOfIds(len(poly))
        #         for i in range(0, len(poly)):
        #             polygon.GetPointIds().SetId(i, poly[i])
        #         polygons.InsertNextCell(polygon)
        # vesselspolydata.SetPolys(polygons)

        print("Saving to %s" % file)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(file)
        writer.SetInputData(vesselspolydata)
        writer.Write()

    def CalculateDijkstraMatrix(self):
        # Warning, this can require TBs of space.
        # a = numpy.memmap('test.mymemmap', dtype='float32', mode='w+', shape=(len(self.PialSurface), len(self.PialSurface)))
        # self.CalculateScipyGraph()
        # for i in range(0, len(self.PialSurface)):
        #     a[i, :] = scipy.sparse.csgraph.dijkstra(self.ScipyGraph, i, limit=300)
        a = scipy.sparse.csgraph.dijkstra(self.ScipyGraph, directed=False, limit=300)
        print(a)
        print("Matrix Done!")

    def FindNearestSurfaceNode(self, couplingpoints, threshold=1e10):
        """
        Map the outlets to the nearest point on the surface.
        A threshold can be used to remove outlets far away.
        """
        KDTree = scipy.spatial.KDTree(self.PialSurface)
        pos = [i.Node.Position for i in couplingpoints]
        # euclidean distance between outlets and surface
        MinDistance, MinDistanceIndex = KDTree.query(pos, k=1)
        print("Finding outlets close to the surface.")
        for index, node in enumerate(couplingpoints):
            if MinDistance[index] < threshold:
                node.PialSurfacePointID = MinDistanceIndex[index]
                self.Roots.append(MinDistanceIndex[index])
            else:
                print("Outlet not close to the surface.")

    def ExportTriangleColour(self, folder):
        file = folder + "SurfaceNodesMapping.csv"
        print("Writing surface node mapping to file: %s" % file)
        with open(file, 'w') as f:
            f.write("TriangleNumber,ClusterID,Area\n")
            for index, item in enumerate(self.Triangles):
                f.write(str(index) + "," + str(self.PolygonColour[index]) + "," +str(self.Areas[index]) + "\n")

    def CalculateAreas(self):
        print("Calculating area per triangle.")
        trianglenodes = []
        for triangle in self.Triangles:
            nodes = [self.PialSurface[i].tolist() for i in triangle]
            trianglenodes.append(nodes)
        self.Areas = [GeneralFunctions.TriangleToArea(triangle) for triangle in trianglenodes]

    def DefineSides(self):
        for index, node in enumerate(self.PialSurface):
            if node[0] > 0:
                self.SideInfo.append(1)
            else:
                self.SideInfo.append(-1)

    def WeightedSampling(self, region, n):
        weights = numpy.array([self.Weights[i] for i in region])
        prob = weights / weights.sum()
        sample = numpy.random.choice(region, n, replace=False, p=prob)
        return list(sample)

    def CalculateDistanceMatRegion(self, method, roots, surfacepoints):
        print("Calculating Distance between roots and the pial surface.")
        if method == "euclidean":
            rootspos = scipy.array([self.PialSurface[point] for point in roots])
            pointpos = scipy.array([self.PialSurface[point] for point in surfacepoints])
            return scipy.spatial.distance.cdist(rootspos, pointpos, 'euclidean')
        else:
            results = scipy.sparse.csgraph.dijkstra(self.ScipyGraph, directed=False, indices=roots)
            results = results[:, surfacepoints]
            return results

    def ToGraph(self):
        print("Converting pialSurface to a weighted Graph.")
        for number in range(0, len(self.Links)):
            node1 = number
            for othernodes in self.Links[node1]:
                if othernodes > node1:
                    self.Graph.add_edge(node1, othernodes,
                                        weight=GeneralFunctions.distancebetweenpoints(self.PialSurface[node1],
                                                                                      self.PialSurface[othernodes]))

    def LoadSurface(self, file="Surface.vtp"):
        print("Loading surface: %s" % file)
        if file[-3:] == "vtp":
            reader = vtk.vtkXMLPolyDataReader()
        elif file[-3:] == "ply":
            reader = vtk.vtkPLYReader()
        elif file[-3:] == "msh":
            self.LoadSurfaceMSH(file)
            return
        else:
            return
        reader.SetFileName(file)
        reader.Update()

        pos_vtk = reader.GetOutput().GetPoints().GetData()
        pos = vtk_to_numpy(pos_vtk)

        NumberOfPoints = reader.GetOutput().GetNumberOfPoints()

        pialsurface = []
        for i in range(0, NumberOfPoints):
            pialsurface.append(pos[i])

        Connections_vtk = reader.GetOutput().GetPolys().GetData()
        Connections = vtk_to_numpy(Connections_vtk)
        NumberOfPolys = reader.GetOutput().GetNumberOfPolys()

        sideinfo = []
        narray = reader.GetOutput().GetPointData().GetNumberOfArrays()
        for i in range(0, narray):
            arrayname = reader.GetOutput().GetPointData().GetArrayName(i)
            if arrayname == "Result":
                sideinfo = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(i))
            if arrayname == "Major Cerebral Artery":
                self.map = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(i))
            if arrayname == "Volume Flow Rate (mL/s)":
                self.flowdata = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(i))
            if arrayname == "Colour":
                flowpercluster = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(i))
                self.PolygonColour = [int(i) for i in flowpercluster]
            if arrayname == "Node Colour":
                flowpercluster = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(i))
                self.NodeColour = [int(i) for i in flowpercluster]

        narray = reader.GetOutput().GetCellData().GetNumberOfArrays()
        for i in range(0, narray):
            arrayname = reader.GetOutput().GetCellData().GetArrayName(i)
            if arrayname == "Major Cerebral Artery":
                self.map = vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(i))
            if arrayname == "Volume Flow Rate (mL/s)":
                self.flowdata = vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(i))
            if arrayname == "Colour":
                flowpercluster = vtk_to_numpy(reader.GetOutput().GetCellData().GetArray(i))
                self.PolygonColour = [int(i) for i in flowpercluster]
        #
        # Triangles = []  # triangles
        # for i in range(0, NumberOfPolys):
        #     # 3 points per polygon
        #     # every 4 int from connections
        #     Triangles.append(Connections[i * 4 + 1:i * 4 + 4])

        currentpos = 0
        Polygons = []
        while currentpos < len(Connections):
            currentpoly = Connections[currentpos]
            poly = Connections[currentpos + 1:currentpos + currentpoly + 1]
            currentpos += currentpoly + 1
            Polygons.append(poly)

        links = [set() for i in range(0, NumberOfPoints)]
        for t in Polygons:
            for i in range(0, len(t)):
                node1 = t[i]
                node2 = t[i - 1]
                links[node1].add(node2)
                links[node2].add(node1)

        # for t in Triangles:
        #     t1 = t[0]
        #     t2 = t[1]
        #     t3 = t[2]
        #     links[t1].add(t2)
        #     links[t1].add(t3)
        #     links[t2].add(t1)
        #     links[t2].add(t3)
        #     links[t3].add(t1)
        #     links[t3].add(t2)

        self.PialSurface = pialsurface
        if len(sideinfo) == 0:
            self.DefineSides()
        else:
            self.SideInfo = list(sideinfo)

        self.Triangles = [i for i in Polygons if len(i) == 3]
        self.CalculateAreas()
        self.Polygons = reader.GetOutput().GetPolys()  # polygons
        self.PolygonsList = Polygons
        self.Links = links
        self.ToGraph()

    def LoadSurfaceMSH(self, file="labeled_vol_mesh.msh"):
        self.MSH = GeneralFunctions.MSHfile()
        self.MSH.Loadfile(file)
        regionsIDs = [4, 21, 22, 23, 24, 25, 26, 30]
        elements, indexes = self.MSH.GetElements(regionsIDs)

        Trianglesmsh = [[i[-1], i[-2], i[-3]] for i in elements]

        pointids = set()
        for i in Trianglesmsh:
            pointids.add(i[0])
            pointids.add(i[1])
            pointids.add(i[2])

        nodeids = sorted(list(pointids))
        NumberOfPoints = len(nodeids)

        Triangles = [[nodeids.index(i) for i in triangle] for triangle in Trianglesmsh]

        pos = [self.MSH.Nodes[i - 1][1:] for i in nodeids]

        links = [set() for i in range(0, NumberOfPoints)]
        for t in Triangles:
            t1 = t[0]
            t2 = t[1]
            t3 = t[2]
            links[t1].add(t2)
            links[t1].add(t3)
            links[t2].add(t1)
            links[t2].add(t3)
            links[t3].add(t1)
            links[t3].add(t2)

        self.PialSurface = pos
        sideinfo = []
        if len(sideinfo) == 0:
            self.DefineSides()
        else:
            self.SideInfo = list(sideinfo)
        self.Triangles = Triangles
        self.CalculateAreas()
        polygons = vtk.vtkCellArray()
        for poly in Triangles:
            polygon = vtk.vtkPolygon()
            polygon.GetPointIds().SetNumberOfIds(len(poly))
            for i in range(0, len(poly)):
                polygon.GetPointIds().SetId(i, poly[i])
            polygons.InsertNextCell(polygon)
        self.Polygons = polygons
        self.Links = links
        self.ToGraph()

        regiondict = {
            21: 5,
            22: 4,
            23: 7,
            24: 2,
            25: 3,
            26: 6,
            4: 8,
            30: 9
        }

        self.MajorVesselID = [regiondict[i[3]] for i in elements]

    def GetTriangleCentroids(self):
        centers = [[] for i in range(0, len(self.Triangles))]
        for index, triangle in enumerate(self.Triangles):
            x = numpy.mean([self.PialSurface[i][0] for i in triangle])
            y = numpy.mean([self.PialSurface[i][1] for i in triangle])
            z = numpy.mean([self.PialSurface[i][2] for i in triangle])
            centers[index] = [x, y, z]
        return centers

    def GetCentersGraph(self, **args):
        print("Calculating the triangle centers.")
        centers = self.GetTriangleCentroids()

        centergraph = nx.Graph()
        [centergraph.add_node(i) for i in range(0, len(self.Triangles))]

        nodestotriagles = [[] for i in self.PialSurface]
        for index, triangle in enumerate(self.Triangles):
            for node in triangle:
                nodestotriagles[node].append(index)
        print("Calculating the connections.")
        # for each edge find the triangles
        if args["method"] == "vertices":
            # method 2: connect center by shared vertex
            links = []
            for index, node in enumerate(nodestotriagles):
                link = list(itertools.permutations(node, 2))
                [links.append(i) for i in link]
            [centergraph.add_edge(i[0], i[1],
                                  weight=GeneralFunctions.distancebetweenpoints(centers[i[0]], centers[i[1]])) for i in
             links]
        else:
            # method 1: connect center by shared edges
            edges = [tuple(sorted(s)) for s in list(self.Graph.edges())]
            edgesdict = dict()
            for index, edge in enumerate(edges):
                edgesdict[edge] = index

            triangletoedge = [[] for i in self.Triangles]
            for index, triangle in enumerate(self.Triangles):
                sortedtriangle = sorted(triangle)
                triangletoedge[index] = [edgesdict[i] for i in
                                         [(sortedtriangle[0], sortedtriangle[1]),
                                          (sortedtriangle[0], sortedtriangle[2]),
                                          (sortedtriangle[1], sortedtriangle[2])]]

            edgestotriangle = [[] for i in edges]
            for index, triangle in enumerate(triangletoedge):
                for edge in triangle:
                    edgestotriangle[edge].append(index)
            [centergraph.add_edge(i[0], i[1],
                                  weight=GeneralFunctions.distancebetweenpoints(centers[i[0]], centers[i[1]])) for i in
             edgestotriangle if len(i) > 1]

        # get the links for each point
        # print(list(centergraph.edges()))
        lines = [tuple(sorted(s)) for s in list(centergraph.edges())]
        linkspercenter = [set() for i in range(0, len(centers))]
        for t in lines:
            t1 = t[0]
            t2 = t[1]
            linkspercenter[t1].add(t2)
            linkspercenter[t2].add(t1)

        sideinfotriangles = []
        if len(self.SideInfo) > 0:
            # for each triangle determine the side, since every node of the triangle is on the same side, this is easy
            for index, triangle in enumerate(self.Triangles):
                side = self.SideInfo[triangle[0]]
                sideinfotriangles.append(side)
        else:
            self.DefineSides()

        print("Determining the new polygons.")
        # Calculate the polygons
        polygons = vtk.vtkCellArray()
        # the order of the polygon is imported
        for numberpoly, element in enumerate(nodestotriagles):
            if element:
                # first point is starting point
                np = len(element)
                newelement = [element[0]]
                element.remove(element[0])
                for index in range(1, np):
                    distanceto = [GeneralFunctions.distancebetweenpoints(centers[i], centers[newelement[-1]]) for i in
                                  element]
                    closest = min(range(len(distanceto)), key=lambda k: distanceto[k])
                    newelement.append(element[closest])
                    element.remove(element[closest])
                nodestotriagles[numberpoly] = newelement

        for poly in nodestotriagles:
            polygon = vtk.vtkPolygon()
            polygon.GetPointIds().SetNumberOfIds(len(poly))
            for i in range(0, len(poly)):
                polygon.GetPointIds().SetId(i, poly[i])
            polygons.InsertNextCell(polygon)

        return centergraph, centers, linkspercenter, polygons, sideinfotriangles

    def GraphToVTP(self, folder):
        file = folder + self.File
        vesselspolydata = vtk.vtkPolyData()
        # Add all the points
        nodes = vtk.vtkPoints()
        for i in self.PialSurface:
            nodes.InsertNextPoint(i)
        vesselspolydata.SetPoints(nodes)

        vesselspolydata.SetPolys(self.Polygons)
        if len(self.SideInfo) > 0:
            nodetype = vtk.vtkIntArray()
            nodetype.SetName("Side")
            for index, triangle in enumerate(self.PialSurface):
                nodetype.InsertNextValue(self.SideInfo[index])
            vesselspolydata.GetPointData().SetScalars(nodetype)

        if len(self.PolygonColour) > 0:
            Colour = vtk.vtkFloatArray()
            Colour.SetNumberOfComponents(1)
            Colour.SetName("Major Cerebral Artery")
            for colour in self.PolygonColour:
                Colour.InsertNextValue(colour)
            vesselspolydata.GetCellData().AddArray(Colour)

        if len(self.NodeColour) > 0:
            NodeColour = vtk.vtkIntArray()
            NodeColour.SetNumberOfComponents(1)
            NodeColour.SetName("Node Colour")
            for colour in self.NodeColour:
                NodeColour.InsertNextValue(colour)
            vesselspolydata.GetPointData().AddArray(NodeColour)

        print("Saving to %s" % file)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(file)
        writer.SetInputData(vesselspolydata)
        writer.Write()


class Perfusion:
    def __init__(self):
        self.CouplingPoints = []  # these are the nodes that connect to the pial surface with their 3d position
        self.PrimalGraph = Surface()
        self.PrimalGraph.File = "PrimalGraph.vtp"
        self.DualGraph = Surface()
        self.DualGraph.File = "DualGraph.vtp"

        self.regionids = [2, 3, 4, 5, 6, 7]
        self.leftregionids = [4, 5, 7]
        self.rightregionids = [2, 3, 6]

        self.TMatrix = None
        self.AllignmentResult = None

    def AddCouplingPoint(self, couplingpoint):
        self.CouplingPoints.append(couplingpoint)

    def RemoveCouplingPoint(self, couplingpoint):
        self.CouplingPoints.remove(couplingpoint)

    def LoadPrimalGraph(self, file):
        self.PrimalGraph.LoadSurface(file)

    def ClusterArea(self, graph):
        for index, cp in enumerate(self.CouplingPoints):
            cp.SurfaceNodes = [i for i, label in enumerate(self.DualGraph.NodeColour) if label == index]
            cp.Area = sum([graph.Areas[i] for i in cp.SurfaceNodes])
            cp.AreaSampling = sum([graph.Areas[i] for i in cp.Pialpoints])

    def WriteClusteringMapping(self, folder):
        file = folder + "Clusters.csv"
        print("Writing Clustering map to file: %s" % file)
        with open(file, 'w') as f:
            f.write(
                "NodeID,Position,Number of CouplingPoints,ClusterID,Area,MajorVesselID\n")
            for index, couplingpoint in enumerate(self.CouplingPoints):
                number = couplingpoint.Node.Number
                position = couplingpoint.Node.Position
                numbernodes = couplingpoint.NumberOfPialPoints
                line = ",".join((str(number),
                                 (str(position[0]) + " " + str(position[1]) + " " + str(position[2])),
                                 str(numbernodes),
                                 str(index),
                                 str(couplingpoint.Area),
                                 str(couplingpoint.Node.MajorVesselID),
                                 # ";".join([str(i) for i in couplingpoint.Pialpoints])))
                                 ))
                f.write(line)
                f.write("\n")

    def UpdateMappedRegionsFlowdata(self, file):
        # data = [line.strip('\n').split(' ') for line in open(file)][1:]

        # with open(file, 'w+') as f:
        #     f.write("Cluster ID,Major Region ID,Volume flow rate (mL/s),Pressure (Pa)\n")
        #     for index, i in enumerate(self.CouplingPoints):
        #         # f.write("%d %.10f %d %f\n" % (GeneralFunctions.StartClusteringIndex + index, i.Node.FlowRate, i.NumberOfTriangles, i.Area))
        #         f.write("%d,%df,%f,%f\n" % (GeneralFunctions.StartClusteringIndex + index, i.Node.MajorVesselID, i.Node.FlowRate, i.Node.Pressure))

        with open(file, 'w') as f:
            # f.write("Cluster ID,Major Region ID,Volume flow rate (mL/s),Pressure (Pa)\n")
            f.write("# region I,Q [ml/s],p [Pa],feeding artery ID\n")
            for index, i in enumerate(self.CouplingPoints):
                # f.write("%d,%d,%f,%f\n" % (GeneralFunctions.StartClusteringIndex + index, MajorvesselIDs[index], float(i[1]), float(i[3])))
                f.write("%d,%f,%f,%d\n" % (GeneralFunctions.StartClusteringIndex + index , i.Node.FlowRate, i.Node.Pressure, GeneralFunctions.MajorIDdict[i.Node.MajorVesselID]))

    def SetDualGraph(self, method="vertices"):
        print("Calculating the dual graph.")
        graph, centers, links, poly, sideinfotriangles = self.PrimalGraph.GetCentersGraph(method=method)
        self.DualGraph.Graph = graph
        self.DualGraph.PialSurface = centers
        self.DualGraph.Links = links
        self.DualGraph.Polygons = poly
        self.DualGraph.SideInfo = sideinfotriangles
        self.DualGraph.Weights = self.PrimalGraph.Areas

    def MapDualGraphToPrimalGraph(self):
        self.PrimalGraph.PolygonColour = self.DualGraph.NodeColour
        self.ClusterArea(self.PrimalGraph)

    def SquareLatticeCP(self, clusters):
        self.PrimalGraph.MajorVesselID = [2 for node in self.PrimalGraph.PialSurface]
        self.CouplingPoints = [CouplingPoint(Node.Node()) for node in range(0, clusters)]
        nodes = numpy.random.choice(len(self.PrimalGraph.PialSurface), size=clusters, replace=False)
        for index, node in enumerate(self.CouplingPoints):
            node.Node.SetPosition(self.PrimalGraph.PialSurface[nodes[index]])
            node.NumberOfPialPoints = int(1 * len(self.PrimalGraph.PialSurface) / clusters)
            node.Node.MajorVesselID = 2
        self.regionids = [2]

    def SquareLattice(self, graph, nn):
        print("generating square lattice for testing")
        lattice = nx.generators.lattice.grid_2d_graph(nn, nn)
        # graph.PialSurface = [[0, numpy.random.normal()*0.1+i // nn, numpy.random.normal()*0.1+i % nn] for i in range(0, nn * nn)]
        graph.PialSurface = [[0, i // nn, i % nn] for i in
                             range(0, nn * nn)]
        # graph.Graph = nx.Graph()
        links = [set() for i in range(0, nn * nn)]
        for edge in lattice.edges():
            u, v = edge
            node1 = u[0] * nn + u[1]
            node2 = v[0] * nn + v[1]
            links[node1].add(node2)
            links[node2].add(node1)

        d = numpy.sqrt(2)
        graph.Graph = lattice
        for n, p in enumerate(graph.PialSurface):
            graph.Graph.node[tuple(p[1:])]['pos'] = p[1:]

        for index, node in enumerate(graph.PialSurface):
            if not (index <= nn or index > nn * nn - nn or index % nn == 0 or index % nn == nn - 1):
                # graph.Graph.add_edge(index, index + nn + 1, weight=d)
                # graph.Graph.add_edge(index, index + nn - 1, weight=d)
                # graph.Graph.add_edge(index, index - nn + 1, weight=d)
                # graph.Graph.add_edge(index, index - nn - 1, weight=d)
                links[index].add(index + nn + 1)
                links[index].add(index + nn - 1)
                links[index].add(index - nn + 1)
                links[index].add(index - nn - 1)
                links[index + nn + 1].add(index)
                links[index + nn - 1].add(index)
                links[index - nn + 1].add(index)
                links[index - nn - 1].add(index)
                # graph.Graph.add_edge(index + nn + 1, index, weight=d)
                # graph.Graph.add_edge(index + nn - 1, index, weight=d)
                # graph.Graph.add_edge(index - nn + 1, index , weight=d)
                # graph.Graph.add_edge(index - nn - 1, index , weight=d)

        # graph.Graph.add_edge(1, nn, weight=d)
        # graph.Graph.add_edge(nn - 2, 2 * nn - 1, weight=d)
        # graph.Graph.add_edge(nn * nn - 2, nn * nn - 1 - nn, weight=d)
        # graph.Graph.add_edge(nn * nn - 2 * nn, nn * nn - 1 * nn + 1, weight=d)

        # graph.Graph.add_edge(nn, 1, weight=d)
        # graph.Graph.add_edge(2 * nn - 1, nn - 2, weight=d)
        # graph.Graph.add_edge(nn * nn - 1 - nn,nn * nn - 2, weight=d)
        # graph.Graph.add_edge(nn * nn - 1 * nn + 1, nn * nn - 2 * nn, weight=d)

        links[1].add(nn)
        links[nn].add(1)
        links[nn - 2].add(2 * nn - 1)
        links[2 * nn - 1].add(nn - 2)
        links[nn * nn - 1 - nn].add(nn * nn - 2)
        links[nn * nn - 2].add(nn * nn - 1 - nn)
        links[nn * nn - 1 * nn + 1].add(nn * nn - 2 * nn)
        links[nn * nn - 2 * nn].add(nn * nn - 1 * nn + 1)

        graph.Links = [list(i) for i in links]
        graph.Areas = [1 for i in graph.PialSurface]  # equal area for all
        graph.Weights = [1 for i in graph.PialSurface]  # equal area for all

        for index, element in enumerate(links):
            for index2, element2 in enumerate(element):
                node1 = tuple(graph.PialSurface[index][1:])
                node2 = tuple(graph.PialSurface[element2][1:])
                d = GeneralFunctions.distancebetweenpoints(node1, node2)
                graph.Graph.add_edge(node1, node2, weight=d)

        # for index, element in enumerate(graph.PialSurface):
        #     for index2, element2 in enumerate(graph.PialSurface):
        #         node1 = tuple(graph.PialSurface[index][1:])
        #         node2 = tuple(graph.PialSurface[index2][1:])
        #         d = distancebetweenpoints(node1,node2)
        #         graph.Graph.add_edge(node1,node2, weight=d)

        # import matplotlib.pyplot as plt
        # labels = nx.get_edge_attributes(graph.Graph, 'weight')
        # pos = nx.get_node_attributes(graph.Graph, 'pos')
        # nx.draw(graph.Graph,pos)
        # # nx.draw_networkx_edge_labels(graph.Graph, edge_labels=labels)
        # plt.show()

    def ExportSurface(self, file, graph):
        nodes = vtk.vtkPoints()
        for index, item in enumerate(graph.PialSurface):
            nodes.InsertNextPoint(item)

        Colour = vtk.vtkFloatArray()
        Colour.SetNumberOfComponents(1)
        Colour.SetName("Colour")
        for colour in graph.PolygonColour:
            Colour.InsertNextValue(colour)

        MajorArteries = vtk.vtkIntArray()
        MajorArteries.SetNumberOfComponents(1)
        MajorArteries.SetName("Major Cerebral Artery")
        for i in range(0, len(graph.PolygonColour)):
            regionid = graph.PolygonColour[i]
            couplingpoint = self.CouplingPoints[regionid]
            MajorArteries.InsertNextValue(couplingpoint.Node.MajorVesselID)

        VesselsPolyData = vtk.vtkPolyData()
        VesselsPolyData.SetPoints(nodes)
        VesselsPolyData.SetPolys(graph.Polygons)
        VesselsPolyData.GetCellData().AddArray(Colour)
        VesselsPolyData.GetCellData().AddArray(MajorArteries)

        writer = vtk.vtkXMLPolyDataWriter()
        print("Writing Clustering to file: %s" % file)
        writer.SetFileName(file)
        writer.SetInputData(VesselsPolyData)
        writer.Write()

    def ExportSurfacePoints(self, file, graph):
        nodes = vtk.vtkPoints()
        for index, item in enumerate(graph.PialSurface):
            nodes.InsertNextPoint(item)

        Colour = vtk.vtkFloatArray()
        Colour.SetNumberOfComponents(1)
        Colour.SetName("Colour")
        for index, colour in enumerate(graph.PolygonColour):
            # if index in graph.Roots:
            #     Colour.InsertNextValue(-1)
            # else:
            Colour.InsertNextValue(colour)

        MajorArteries = vtk.vtkIntArray()
        MajorArteries.SetNumberOfComponents(1)
        MajorArteries.SetName("Major Cerebral Artery")
        for i in range(0, len(graph.PolygonColour)):
            regionid = graph.PolygonColour[i]
            couplingpoint = self.CouplingPoints[regionid]
            MajorArteries.InsertNextValue(couplingpoint.Node.MajorVesselID)

        lines = vtk.vtkCellArray()
        for index, element in enumerate(graph.Links):
            for otherelement in element:
                line0 = vtk.vtkLine()
                line0.GetPointIds().SetId(0, index)
                line0.GetPointIds().SetId(1, otherelement)
                lines.InsertNextCell(line0)

        VesselsPolyData = vtk.vtkPolyData()
        VesselsPolyData.SetPoints(nodes)
        # VesselsPolyData.SetPolys(graph.Polygons)
        VesselsPolyData.GetPointData().AddArray(Colour)
        VesselsPolyData.GetPointData().AddArray(MajorArteries)
        VesselsPolyData.SetLines(lines)

        writer = vtk.vtkXMLPolyDataWriter()
        print("Writing Clustering to file: %s" % file)
        writer.SetFileName(file)
        writer.SetInputData(VesselsPolyData)
        writer.Write()

    def ClusteringByRegion(self, dualgraph, method):
        print("Start clustering.")
        resultcluster = [[] for i in self.regionids]

        dualgraph.FindNearestSurfaceNode(self.CouplingPoints)
        dualgraph.CalculateScipyGraph()

        for regionindex, currentregion in enumerate(self.regionids):
            regionsum = len([i for i, x in enumerate(dualgraph.MajorVesselID) if x == currentregion])
            regionroots = [node for node in self.CouplingPoints if node.Node.MajorVesselID == currentregion]
            rootsum = sum([node.NumberOfPialPoints for node in regionroots])
            if rootsum == len(regionroots):
                print("Setting number of surface elements to 50 per root.")
                for node in regionroots:
                    node.NumberOfPialPoints = 50
                rootsum = len(regionroots) * 50
            print("RegionID: " + str(currentregion) + ", Number of coupling points: " + str(
                rootsum) + ", number of surface elements: " + str(regionsum))
            if rootsum > regionsum:
                errorline = "Error: Number of coupling points: " + str(
                    rootsum) + " exceeds the number of surface elements: " + str(regionsum)
                raise Exception(errorline)

        for regionindex, currentregion in enumerate(self.regionids):
            region = [i for i, x in enumerate(dualgraph.MajorVesselID) if x == currentregion]
            regionroots = [node for node in self.CouplingPoints if node.Node.MajorVesselID == currentregion]
            NRootsLabels = [node.NumberOfPialPoints for node in regionroots]
            PointsN = sum(NRootsLabels)  # total points
            RootsN = len(regionroots)
            regionpoints = dualgraph.WeightedSampling(region, PointsN)

            # calculate distance matrix for the pial surface.
            roots = [point.PialSurfacePointID for point in regionroots]
            PialDMat = scipy.sparse.csgraph.dijkstra(dualgraph.ScipyGraph, directed=False, indices=regionpoints)
            PialDMat = PialDMat[:, regionpoints]

            print("Optimizing clustering.")
            DistanceMat = dualgraph.CalculateDistanceMatRegion(method, roots, regionpoints)
            iter = 0
            rootpositions = []
            while 1:
                # save the roots during this iteration to plot their trajectory.
                rootpositions.append(copy.deepcopy(roots))

                flat_list = [item for sublist in DistanceMat for item in sublist]
                AssignedLabelsN = [0 for i in range(0, RootsN)]
                Labels = [-1 for i in range(0, PointsN)]

                sortedindex = sorted(range(len(flat_list)), key=lambda k: flat_list[k])

                for index in sortedindex:
                    ClusterIndex = index // PointsN
                    PointIndex = index % PointsN

                    if AssignedLabelsN[ClusterIndex] < NRootsLabels[ClusterIndex] and Labels[PointIndex] == -1:
                        # if Labels[PointIndex] == -1:
                        Labels[PointIndex] = ClusterIndex
                        AssignedLabelsN[ClusterIndex] += 1

                swapiter = 0
                while 1:
                    swaplist = []
                    for point in range(0, PointsN):
                        possiblelabels = DistanceMat[:, point]
                        sortedlabels = sorted(range(len(possiblelabels)), key=lambda k: possiblelabels[k])
                        currentlabel = Labels[point]
                        if currentlabel != sortedlabels[0]:
                            swaplist.append(point)
                    changecounter = 0
                    for swap in swaplist:
                        differenceafterswap = [(DistanceMat[Labels[swap], swap] + DistanceMat[Labels[i], i])
                                               - (DistanceMat[Labels[i], swap] + DistanceMat[Labels[swap], i])
                                               for i in range(0, PointsN)]

                        differencebyindex = max(range(len(differenceafterswap)), key=lambda k: differenceafterswap[k])
                        if differenceafterswap[differencebyindex] > 0:
                            Labels[differencebyindex], Labels[swap] = Labels[swap], Labels[differencebyindex]
                            changecounter += 1
                    print("Total changes: %d" % changecounter)
                    swapiter += 1
                    if changecounter == 0 or swapiter == 20:
                        break

                stopcrit = 0

                print("Recalculating roots.")
                clusterelements = [[i for i in range(0, len(Labels)) if Labels[i] == element] for element
                                   in range(0, RootsN)]

                for index, cluster in enumerate(clusterelements):
                    # if method == "euclidean":
                    # newpos = numpy.mean([dualgraph.PialSurface[regionpoints[i]] for i in cluster], axis=0)
                    # diff = distancebetweenpoints(roots[index], newpos)
                    # stopcrit += diff
                    # roots[index] = newpos

                    distancematforcluster = numpy.empty([len(cluster), len(cluster)])
                    for num, point in enumerate(cluster):
                        for num2, point2 in enumerate(cluster):
                            distancematforcluster[num, num2] = PialDMat[point, point2]
                    distances = numpy.amax(distancematforcluster, axis=0)
                    newguess = min(range(0, len(distances)), key=lambda k: distances[k])
                    newguessPoint = cluster[newguess]
                    # new guess for root node is the Jardan Center of the cluster
                    # calculate the change in position of the root
                    newpos = dualgraph.PialSurface[regionpoints[newguessPoint]]
                    oldpos = dualgraph.PialSurface[roots[index]]
                    diff = GeneralFunctions.distancebetweenpoints(oldpos, newpos)
                    stopcrit += diff
                    roots[index] = regionpoints[newguessPoint]
                    # print(roots[index])
                print("Root displacement: %f" % stopcrit)
                iter += 1
                # print(rootpositions)
                if stopcrit == 0 or iter == 20:
                    resultcluster[regionindex] = (regionpoints, Labels, rootpositions)
                    break

                # if method == "euclidean":
                #     DistanceMat = dualgraph.CalculateDistanceMatRegion(method, roots, regionpoints)
                # else:
                DistanceMat = dualgraph.CalculateDistanceMatRegion(method, roots, regionpoints)
                # DistanceMat = PialDMat[roots, :]

        dualgraph.Points = []
        dualgraph.Labels = []
        CurrentID = 0
        reorderedcouplingpoints = []
        for index, region in enumerate(resultcluster):
            dualgraph.Points += list(region[0])
            labels = [i + CurrentID for i in region[1]]
            dualgraph.Labels += labels
            regionroots = [node for node in self.CouplingPoints if node.Node.MajorVesselID == self.regionids[index]]
            # rootposition = [[iteration[i] for iteration in region[2]] for i in range(0, len(regionroots))]  # per region
            for index2, root in enumerate(regionroots):
                rootpos = [iteration[index2] for iteration in region[2]]
                root.RootPos = rootpos
            reorderedcouplingpoints += regionroots
            CurrentID += len(regionroots)

        self.CouplingPoints = reorderedcouplingpoints

        # For each outlet, save the connected pial points
        for iter in range(0, len(self.CouplingPoints)):
            points = [dualgraph.Points[i] for i in range(0, len(dualgraph.Labels)) if dualgraph.Labels[i] == iter]
            self.CouplingPoints[iter].Pialpoints = points

        # Colour the mesh with the nearest coupling point
        # Colouring the mesh
        meshcolour = [-1 for node in dualgraph.PialSurface]
        for i in range(0, len(dualgraph.Points)):
            point = dualgraph.Points[i]
            meshcolour[point] = dualgraph.Labels[i]

        iter = 0
        while 1:
            nocolour = 0
            colourupdate = []
            for index, item in enumerate(meshcolour):
                if item == -1:
                    colouredlinkednodes = [i for i in dualgraph.Links[index] if meshcolour[i] >= 0]
                    if len(colouredlinkednodes) > 0:
                        # calculate distance to each linked node and take the nearest node's colour
                        distancetocolourednodes = [
                            GeneralFunctions.distancebetweenpoints(dualgraph.PialSurface[index],
                                                                   dualgraph.PialSurface[p])
                            for p in colouredlinkednodes]

                        mindifferencebyindex = min(range(len(distancetocolourednodes)),
                                                   key=lambda k: distancetocolourednodes[k])
                        # update after the for loop
                        colourupdate.append((index, colouredlinkednodes[mindifferencebyindex]))
                    else:
                        nocolour += 1
            for element in colourupdate:
                meshcolour[element[0]] = meshcolour[element[1]]
            iter += 1
            if nocolour == 0 or iter == 1000:
                break
        dualgraph.NodeColour = meshcolour

    def ClusteringByRegionNoConstraints(self, dualgraph, method):
        print("Start clustering.")
        resultcluster = [[] for i in self.regionids]

        dualgraph.FindNearestSurfaceNode(self.CouplingPoints)
        dualgraph.CalculateScipyGraph()

        for regionindex, currentregion in enumerate(self.regionids):
            regionsum = len([i for i, x in enumerate(dualgraph.MajorVesselID) if x == currentregion])
            regionroots = [node for node in self.CouplingPoints if node.Node.MajorVesselID == currentregion]
            rootsum = sum([node.NumberOfPialPoints for node in regionroots])
            # if rootsum == len(regionroots):
            #     print("Setting number of surface elements to 50 per root.")
            #     for node in regionroots:
            #         node.NumberOfPialPoints = 50
            #     rootsum = len(regionroots) * 50
            # print("RegionID: " + str(currentregion) + ", Number of coupling points: " + str(
            #     rootsum) + ", number of surface elements: " + str(regionsum))
            # if rootsum > regionsum:
            #     errorline = "Error: Number of coupling points: " + str(
            #         rootsum) + " exceeds the number of surface elements: " + str(regionsum)
            #     raise Exception(errorline)

        for regionindex, currentregion in enumerate(self.regionids):
            region = [i for i, x in enumerate(dualgraph.MajorVesselID) if x == currentregion]
            regionroots = [node for node in self.CouplingPoints if node.Node.MajorVesselID == currentregion]
            NRootsLabels = [node.NumberOfPialPoints for node in regionroots]
            PointsN = sum(NRootsLabels)  # total points
            RootsN = len(regionroots)
            regionpoints = dualgraph.WeightedSampling(region, PointsN)

            # calculate distance matrix for the pial surface.
            roots = [point.PialSurfacePointID for point in regionroots]
            PialDMat = scipy.sparse.csgraph.dijkstra(dualgraph.ScipyGraph, directed=False, indices=regionpoints)
            PialDMat = PialDMat[:, regionpoints]

            print("Optimizing clustering.")
            DistanceMat = dualgraph.CalculateDistanceMatRegion(method, roots, regionpoints)
            iter = 0
            rootpositions = []
            while 1:
                # save the roots during this iteration to plot their trajectory.
                rootpositions.append(copy.deepcopy(roots))

                flat_list = [item for sublist in DistanceMat for item in sublist]
                AssignedLabelsN = [0 for i in range(0, RootsN)]
                Labels = [-1 for i in range(0, PointsN)]

                sortedindex = sorted(range(len(flat_list)), key=lambda k: flat_list[k])

                for index in sortedindex:
                    ClusterIndex = index // PointsN
                    PointIndex = index % PointsN

                    # if AssignedLabelsN[ClusterIndex] < NRootsLabels[ClusterIndex] and Labels[PointIndex] == -1:
                    if Labels[PointIndex] == -1:
                        Labels[PointIndex] = ClusterIndex
                        AssignedLabelsN[ClusterIndex] += 1

                swapiter = 0
                while 1:
                    swaplist = []
                    for point in range(0, PointsN):
                        possiblelabels = DistanceMat[:, point]
                        sortedlabels = sorted(range(len(possiblelabels)), key=lambda k: possiblelabels[k])
                        currentlabel = Labels[point]
                        if currentlabel != sortedlabels[0]:
                            swaplist.append(point)
                    changecounter = 0
                    for swap in swaplist:
                        differenceafterswap = [(DistanceMat[Labels[swap], swap] + DistanceMat[Labels[i], i])
                                               - (DistanceMat[Labels[i], swap] + DistanceMat[Labels[swap], i])
                                               for i in range(0, PointsN)]

                        differencebyindex = max(range(len(differenceafterswap)), key=lambda k: differenceafterswap[k])
                        if differenceafterswap[differencebyindex] > 0:
                            Labels[differencebyindex], Labels[swap] = Labels[swap], Labels[differencebyindex]
                            changecounter += 1
                    print("Total changes: %d" % changecounter)
                    swapiter += 1
                    if changecounter == 0 or swapiter == 20:
                        break

                stopcrit = 0

                print("Recalculating roots.")
                clusterelements = [[i for i in range(0, len(Labels)) if Labels[i] == element] for element
                                   in range(0, RootsN)]

                for index, cluster in enumerate(clusterelements):
                    # if method == "euclidean":
                    # newpos = numpy.mean([dualgraph.PialSurface[regionpoints[i]] for i in cluster], axis=0)
                    # diff = distancebetweenpoints(roots[index], newpos)
                    # stopcrit += diff
                    # roots[index] = newpos

                    distancematforcluster = numpy.empty([len(cluster), len(cluster)])
                    for num, point in enumerate(cluster):
                        for num2, point2 in enumerate(cluster):
                            distancematforcluster[num, num2] = PialDMat[point, point2]
                    distances = numpy.amax(distancematforcluster, axis=0)
                    newguess = min(range(0, len(distances)), key=lambda k: distances[k])
                    newguessPoint = cluster[newguess]
                    # new guess for root node is the Jardan Center of the cluster
                    # calculate the change in position of the root
                    newpos = dualgraph.PialSurface[regionpoints[newguessPoint]]
                    oldpos = dualgraph.PialSurface[roots[index]]
                    diff = GeneralFunctions.distancebetweenpoints(oldpos, newpos)
                    stopcrit += diff
                    roots[index] = regionpoints[newguessPoint]
                    # print(roots[index])
                print("Root displacement: %f" % stopcrit)
                iter += 1
                # print(rootpositions)
                if stopcrit == 0 or iter == 20:
                    resultcluster[regionindex] = (regionpoints, Labels, rootpositions)
                    break

                # if method == "euclidean":
                #     DistanceMat = dualgraph.CalculateDistanceMatRegion(method, roots, regionpoints)
                # else:
                DistanceMat = dualgraph.CalculateDistanceMatRegion(method, roots, regionpoints)
                # DistanceMat = PialDMat[roots, :]

        dualgraph.Points = []
        dualgraph.Labels = []
        CurrentID = 0
        reorderedcouplingpoints = []
        for index, region in enumerate(resultcluster):
            dualgraph.Points += list(region[0])
            labels = [i + CurrentID for i in region[1]]
            dualgraph.Labels += labels
            regionroots = [node for node in self.CouplingPoints if node.Node.MajorVesselID == self.regionids[index]]
            # rootposition = [[iteration[i] for iteration in region[2]] for i in range(0, len(regionroots))]  # per region
            for index2, root in enumerate(regionroots):
                rootpos = [iteration[index2] for iteration in region[2]]
                root.RootPos = rootpos
            reorderedcouplingpoints += regionroots
            CurrentID += len(regionroots)

        self.CouplingPoints = reorderedcouplingpoints

        # For each outlet, save the connected pial points
        for iter in range(0, len(self.CouplingPoints)):
            points = [dualgraph.Points[i] for i in range(0, len(dualgraph.Labels)) if dualgraph.Labels[i] == iter]
            self.CouplingPoints[iter].Pialpoints = points

        # Colour the mesh with the nearest coupling point
        # Colouring the mesh
        meshcolour = [-1 for node in dualgraph.PialSurface]
        for i in range(0, len(dualgraph.Points)):
            point = dualgraph.Points[i]
            meshcolour[point] = dualgraph.Labels[i]

        iter = 0
        while 1:
            nocolour = 0
            colourupdate = []
            for index, item in enumerate(meshcolour):
                if item == -1:
                    colouredlinkednodes = [i for i in dualgraph.Links[index] if meshcolour[i] >= 0]
                    if len(colouredlinkednodes) > 0:
                        # calculate distance to each linked node and take the nearest node's colour
                        distancetocolourednodes = [
                            GeneralFunctions.distancebetweenpoints(dualgraph.PialSurface[index],
                                                                   dualgraph.PialSurface[p])
                            for p in colouredlinkednodes]

                        mindifferencebyindex = min(range(len(distancetocolourednodes)),
                                                   key=lambda k: distancetocolourednodes[k])
                        # update after the for loop
                        colourupdate.append((index, colouredlinkednodes[mindifferencebyindex]))
                    else:
                        nocolour += 1
            for element in colourupdate:
                meshcolour[element[0]] = meshcolour[element[1]]
            iter += 1
            if nocolour == 0 or iter == 1000:
                break
        dualgraph.NodeColour = meshcolour

    def ClusteringByRegionScalingLaw(self, dualgraph, method="",
                                     fractiontriangles=(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.5), debug=False):
        """
        Algorithm to estimate perfusion territories.
        :param dualgraph: Graph to use for the clustering. Clustering happens with the vertices of the graph, clustering of triangles requires the dual graph.
        :param method: Default option: Dijkstra, Other option: 'euclidean'
        :param fractiontriangles: Sets the fraction of triangles to use for the initial clustering. Mainly for performance reasons.
        :param debug: Default:False. True sets the number of points per cluster to one. Use this to get a quickyclustering to the nearst outlet.
        :return: Clustering result is stored under dualgraph.NodeColour. See MapDualGraphToPrimalGraph() for mapping back to the primal graph.
        """
        print("Start clustering.")
        resultcluster = [[] for i in self.regionids]

        dualgraph.FindNearestSurfaceNode(self.CouplingPoints)
        dualgraph.CalculateScipyGraph()

        for regionindex, currentregion in enumerate(self.regionids):
            regionsum = len([i for i, x in enumerate(dualgraph.MajorVesselID) if x == currentregion])
            regionroots = [node for node in self.CouplingPoints if node.Node.MajorVesselID == currentregion]
            # rootsum = sum([node.NumberOfPialPoints for node in regionroots])

            TotalTriangles = int(regionsum * fractiontriangles[regionindex])
            # print("Using %d of the %d triangles for region %d" % (TotalTriangles, regionsum, currentregion))
            inputRadii = [i.Node.Radius for i in regionroots]
            RtCubed = sum([numpy.power(i, 3) for i in inputRadii])

            fractionArea = [numpy.power(i, 3) / RtCubed for i in inputRadii]

            numberPialSurfaceElements = [i * TotalTriangles for i in fractionArea]

            roundedNumberPialSurfaceElements = [int(round(i)) for i in numberPialSurfaceElements]
            for index, node in enumerate(regionroots):
                node.NumberOfPialPoints = roundedNumberPialSurfaceElements[index]

            rootsum = sum(roundedNumberPialSurfaceElements)
            # if sum(roundedNumberPialSurfaceElements) > regionsum:
            #     print("Error in the number of PialSurfaceElements: too many assigned")

            # if rootsum == len(regionroots):
            #     print("Setting number of surface elements to 50 per root.")
            #     for node in regionroots:
            #         node.NumberOfPialPoints = 50
            #     rootsum = len(regionroots) * 50
            if debug:
                print("Debug Mode.")
                for index, node in enumerate(regionroots):
                    node.NumberOfPialPoints = 1

            print("RegionID: " + str(currentregion) + ", Number of coupling points: " + str(
                rootsum) + ", number of surface elements: " + str(regionsum))
            if rootsum > regionsum:
                errorline = "Error: Number of coupling points: " + str(
                    rootsum) + " exceeds the number of surface elements: " + str(regionsum)
                raise Exception(errorline)

        for regionindex, currentregion in enumerate(self.regionids):
            region = [i for i, x in enumerate(dualgraph.MajorVesselID) if x == currentregion]
            regionroots = [node for node in self.CouplingPoints if node.Node.MajorVesselID == currentregion]
            NRootsLabels = [node.NumberOfPialPoints for node in regionroots]
            PointsN = sum(NRootsLabels)  # total points
            RootsN = len(regionroots)
            regionpoints = dualgraph.WeightedSampling(region, PointsN)

            # calculate distance matrix for the pial surface.
            roots = [point.PialSurfacePointID for point in regionroots]
            PialDMat = scipy.sparse.csgraph.dijkstra(dualgraph.ScipyGraph, directed=False, indices=regionpoints)
            PialDMat = PialDMat[:, regionpoints]

            print("Optimizing clustering.")
            DistanceMat = dualgraph.CalculateDistanceMatRegion(method, roots, regionpoints)
            iter = 0
            rootpositions = []
            while 1:
                # save the roots during this iteration to plot their trajectory.
                rootpositions.append(copy.deepcopy(roots))

                flat_list = [item for sublist in DistanceMat for item in sublist]
                AssignedLabelsN = [0 for i in range(0, RootsN)]
                Labels = [-1 for i in range(0, PointsN)]

                sortedindex = sorted(range(len(flat_list)), key=lambda k: flat_list[k])

                for index in sortedindex:
                    ClusterIndex = index // PointsN
                    PointIndex = index % PointsN

                    if AssignedLabelsN[ClusterIndex] < NRootsLabels[ClusterIndex] and Labels[PointIndex] == -1:
                        # if Labels[PointIndex] == -1:
                        Labels[PointIndex] = ClusterIndex
                        AssignedLabelsN[ClusterIndex] += 1

                swapiter = 0
                while 1:
                    swaplist = []
                    for point in range(0, PointsN):
                        possiblelabels = DistanceMat[:, point]
                        sortedlabels = sorted(range(len(possiblelabels)), key=lambda k: possiblelabels[k])
                        currentlabel = Labels[point]
                        if currentlabel != sortedlabels[0]:
                            swaplist.append(point)
                    changecounter = 0
                    for swap in swaplist:
                        differenceafterswap = [(DistanceMat[Labels[swap], swap] + DistanceMat[Labels[i], i])
                                               - (DistanceMat[Labels[i], swap] + DistanceMat[Labels[swap], i])
                                               for i in range(0, PointsN)]

                        differencebyindex = max(range(len(differenceafterswap)), key=lambda k: differenceafterswap[k])
                        if differenceafterswap[differencebyindex] > 0:
                            Labels[differencebyindex], Labels[swap] = Labels[swap], Labels[differencebyindex]
                            changecounter += 1
                    print("Total changes: %d" % changecounter)
                    swapiter += 1
                    if changecounter == 0 or swapiter == 20:
                        break

                stopcrit = 0

                print("Recalculating roots.")
                clusterelements = [[i for i in range(0, len(Labels)) if Labels[i] == element] for element
                                   in range(0, RootsN)]

                for index, cluster in enumerate(clusterelements):
                    # if method == "euclidean":
                    # newpos = numpy.mean([dualgraph.PialSurface[regionpoints[i]] for i in cluster], axis=0)
                    # diff = distancebetweenpoints(roots[index], newpos)
                    # stopcrit += diff
                    # roots[index] = newpos

                    distancematforcluster = numpy.empty([len(cluster), len(cluster)])
                    for num, point in enumerate(cluster):
                        for num2, point2 in enumerate(cluster):
                            distancematforcluster[num, num2] = PialDMat[point, point2]
                    distances = numpy.amax(distancematforcluster, axis=0)
                    newguess = min(range(0, len(distances)), key=lambda k: distances[k])
                    newguessPoint = cluster[newguess]
                    # new guess for root node is the Jardan Center of the cluster
                    # calculate the change in position of the root
                    newpos = dualgraph.PialSurface[regionpoints[newguessPoint]]
                    oldpos = dualgraph.PialSurface[roots[index]]
                    diff = GeneralFunctions.distancebetweenpoints(oldpos, newpos)
                    stopcrit += diff
                    roots[index] = regionpoints[newguessPoint]
                    # print(roots[index])
                print("Root displacement: %f" % stopcrit)
                iter += 1
                # print(rootpositions)
                if stopcrit == 0 or iter == 20:
                    resultcluster[regionindex] = (regionpoints, Labels, rootpositions)
                    break

                # if method == "euclidean":
                #     DistanceMat = dualgraph.CalculateDistanceMatRegion(method, roots, regionpoints)
                # else:
                DistanceMat = dualgraph.CalculateDistanceMatRegion(method, roots, regionpoints)
                # DistanceMat = PialDMat[roots, :]

        dualgraph.Points = []
        dualgraph.Labels = []
        CurrentID = 0
        reorderedcouplingpoints = []
        for index, region in enumerate(resultcluster):
            dualgraph.Points += list(region[0])
            labels = [i + CurrentID for i in region[1]]
            dualgraph.Labels += labels
            regionroots = [node for node in self.CouplingPoints if node.Node.MajorVesselID == self.regionids[index]]
            # rootposition = [[iteration[i] for iteration in region[2]] for i in range(0, len(regionroots))]  # per region
            for index2, root in enumerate(regionroots):
                rootpos = [iteration[index2] for iteration in region[2]]
                root.RootPos = rootpos
            reorderedcouplingpoints += regionroots
            CurrentID += len(regionroots)

        self.CouplingPoints = reorderedcouplingpoints

        # For each outlet, save the connected pial points
        for iter in range(0, len(self.CouplingPoints)):
            points = [dualgraph.Points[i] for i in range(0, len(dualgraph.Labels)) if dualgraph.Labels[i] == iter]
            self.CouplingPoints[iter].Pialpoints = points

        # Colour the mesh with the nearest coupling point
        # Colouring the mesh
        meshcolour = [-1 for node in dualgraph.PialSurface]
        for i in range(0, len(dualgraph.Points)):
            point = dualgraph.Points[i]
            meshcolour[point] = dualgraph.Labels[i]

        iter = 0
        while 1:
            nocolour = 0
            colourupdate = []
            for index, item in enumerate(meshcolour):
                if item == -1:
                    colouredlinkednodes = [i for i in dualgraph.Links[index] if meshcolour[i] >= 0]
                    if len(colouredlinkednodes) > 0:
                        # calculate distance to each linked node and take the nearest node's colour
                        distancetocolourednodes = [
                            GeneralFunctions.distancebetweenpoints(dualgraph.PialSurface[index],
                                                                   dualgraph.PialSurface[p])
                            for p in colouredlinkednodes]

                        mindifferencebyindex = min(range(len(distancetocolourednodes)),
                                                   key=lambda k: distancetocolourednodes[k])
                        # update after the for loop
                        colourupdate.append((index, colouredlinkednodes[mindifferencebyindex]))
                    else:
                        nocolour += 1
            for element in colourupdate:
                meshcolour[element[0]] = meshcolour[element[1]]
            iter += 1
            if nocolour == 0 or iter == 1000:
                break
        dualgraph.NodeColour = meshcolour

    # def GeneratePAs(self, graph, n=250000):
    #
    #     weights = numpy.array(self.PrimalGraph.Areas)
    #     prob = weights / weights.sum()
    #     # n = 250000
    #     sample = numpy.random.choice(range(0, len(self.PrimalGraph.Triangles)), n, replace=True, p=prob)
    #     # sample = numpy.random.choice(range(0,10000,100), n, replace=True)
    #
    #     uniformdistribution = numpy.random.uniform(0, 1, (len(sample), 2))
    #
    #     penetratingartery = []  # list of positions within a triangle
    #     for index, sampledtriagle in enumerate(sample):
    #         triangle = self.PrimalGraph.Triangles[sampledtriagle]
    #         v1 = self.PrimalGraph.PialSurface[triangle[1]] - self.PrimalGraph.PialSurface[
    #             triangle[0]]
    #         v2 = self.PrimalGraph.PialSurface[triangle[2]] - self.PrimalGraph.PialSurface[
    #             triangle[0]]
    #
    #         if uniformdistribution[index][0] > 0 and uniformdistribution[index][1] > 0 and uniformdistribution[index][
    #             0] + uniformdistribution[index][1] < 1:
    #             x = self.PrimalGraph.PialSurface[triangle[0]] + uniformdistribution[index][0] * v1 + \
    #                 uniformdistribution[index][1] * v2
    #         else:
    #             x = self.PrimalGraph.PialSurface[triangle[0]] + (1 - uniformdistribution[index][0]) * v1 + \
    #                 (1 - uniformdistribution[index][1]) * v2
    #
    #         penetratingartery.append(x)
    #
    #     # save the connecting vessels
    #     connectionsPA = []
    #     for index, sampledtriagle in enumerate(sample):
    #         numb = index + len(self.PrimalGraph.Triangles)  # number of the generated point
    #         # sampled triangle is the node that connects to the generated point
    #         connectionsPA.append((sampledtriagle, numb))
    #
    #     return penetratingartery, connectionsPA

    # def SplitRegionInTwo(self, dualgraph, points, PialDMat):
    #     print("Splitting region in two.")
    #
    #     RootsN = 2
    #     PointsN = len(points)
    #     NRootsLabels = [numpy.floor(PointsN / 2), numpy.ceil(PointsN / 2)]
    #     regionpoints = points
    #     roots = numpy.random.choice(points, 2, replace=False)
    #
    #     sampledict = {sample: index for index, sample in enumerate(points)}
    #
    #     # DistanceMat = dualgraph.CalculateDistanceMatRegion("", roots, regionpoints)
    #     # DistanceMat2 = PialDMat[[sampledict[root] for root in roots],:]
    #     iter = 0
    #     result = []
    #     while 1:
    #         DistanceMat = PialDMat[[sampledict[root] for root in roots], :]
    #
    #         flat_list = [item for sublist in DistanceMat for item in sublist]
    #         AssignedLabelsN = [0 for i in range(0, RootsN)]
    #         Labels = [-1 for i in range(0, PointsN)]
    #
    #         sortedindex = sorted(range(len(flat_list)), key=lambda k: flat_list[k])
    #
    #         for index in sortedindex:
    #             ClusterIndex = index // PointsN
    #             PointIndex = index % PointsN
    #
    #             if AssignedLabelsN[ClusterIndex] < NRootsLabels[ClusterIndex] and Labels[PointIndex] == -1:
    #                 # if Labels[PointIndex] == -1:
    #                 Labels[PointIndex] = ClusterIndex
    #                 AssignedLabelsN[ClusterIndex] += 1
    #
    #         swapiter = 0
    #         while 1:
    #             swaplist = []
    #             for point in range(0, PointsN):
    #                 possiblelabels = DistanceMat[:, point]
    #                 sortedlabels = sorted(range(len(possiblelabels)), key=lambda k: possiblelabels[k])
    #                 currentlabel = Labels[point]
    #                 if currentlabel != sortedlabels[0]:
    #                     swaplist.append(point)
    #             changecounter = 0
    #             for swap in swaplist:
    #                 differenceafterswap = [(DistanceMat[Labels[swap], swap] + DistanceMat[Labels[i], i])
    #                                        - (DistanceMat[Labels[i], swap] + DistanceMat[Labels[swap], i])
    #                                        for i in range(0, PointsN)]
    #
    #                 differencebyindex = max(range(len(differenceafterswap)), key=lambda k: differenceafterswap[k])
    #                 if differenceafterswap[differencebyindex] > 0:
    #                     Labels[differencebyindex], Labels[swap] = Labels[swap], Labels[differencebyindex]
    #                     changecounter += 1
    #             print("Total changes: %d" % changecounter)
    #             swapiter += 1
    #             if changecounter == 0 or swapiter == 20:
    #                 break
    #
    #         stopcrit = 0
    #
    #         print("Recalculating roots.")
    #         clusterelements = [[i for i in range(0, len(Labels)) if Labels[i] == element] for element
    #                            in range(0, RootsN)]
    #
    #         for index, cluster in enumerate(clusterelements):
    #             # if method == "euclidean":
    #             # newpos = numpy.mean([dualgraph.PialSurface[regionpoints[i]] for i in cluster], axis=0)
    #             # diff = distancebetweenpoints(roots[index], newpos)
    #             # stopcrit += diff
    #             # roots[index] = newpos
    #
    #             distancematforcluster = numpy.empty([len(cluster), len(cluster)])
    #             for num, point in enumerate(cluster):
    #                 for num2, point2 in enumerate(cluster):
    #                     distancematforcluster[num, num2] = PialDMat[point, point2]
    #             print(distancematforcluster)
    #             distances = numpy.amax(distancematforcluster, axis=0)
    #             newguess = min(range(0, len(distances)), key=lambda k: distances[k])
    #             newguessPoint = cluster[newguess]
    #             # new guess for root node is the Jardan Center of the cluster
    #             # calculate the change in position of the root
    #             newpos = dualgraph.PialSurface[regionpoints[newguessPoint]]
    #             oldpos = dualgraph.PialSurface[roots[index]]
    #             diff = GeneralFunctions.distancebetweenpoints(oldpos, newpos)
    #             stopcrit += diff
    #             roots[index] = regionpoints[newguessPoint]
    #             # print(roots[index])
    #         print("Root displacement: %f" % stopcrit)
    #         iter += 1
    #         # print(rootpositions)
    #         if stopcrit == 0 or iter == 20:
    #             result = Labels
    #             break
    #
    #         # DistanceMat = dualgraph.CalculateDistanceMatRegion("", roots, regionpoints)
    #
    #     region1 = [points[i] for i in range(0, len(result)) if result[i] == 0]
    #     region2 = [points[i] for i in range(0, len(result)) if result[i] == 1]
    #     return region1, region2

    def RemapMajorRegions(self):
        """
        Remap wrongly mapped regions to the nearst major region.
        :return: Updated mappings
        """
        print("Remapping the Major Regions.")
        regionsMajorID = []

        uniqueids = set(self.PrimalGraph.map)
        for id in uniqueids:
            points = [i for i, idcomp in enumerate(self.PrimalGraph.map) if idcomp == id]
            regions = []
            while len(points) > 0:
                region = [points[0]]
                candicates = [points[0]]
                points.remove(points[0])
                while len(candicates) > 0:
                    newcandicates = []
                    for index, item in enumerate(candicates):
                        for index2, item2 in enumerate(self.DualGraph.Links[item]):
                            if item2 in points:
                                region.append(item2)
                                newcandicates.append(item2)
                                points.remove(item2)
                    candicates = newcandicates
                regions.append(region)
            regionsMajorID.append(regions)

        LargestRegions = [sorted(regions, key=lambda x: len(x), reverse=True) for regions in regionsMajorID]
        # print(LargestRegions)
        # keptregions = [i[0] for i in LargestRegions]

        nodes = []
        wronglabel = [i[1:] for i in LargestRegions]
        for i, item in enumerate(wronglabel):
            for element in item:
                nodes.extend(element)

        for node in nodes:
            self.PrimalGraph.map[node] = -1

        iter = 0
        while 1:
            nocolour = 0
            colourupdate = []
            for _, item in enumerate(nodes):
                if self.PrimalGraph.map[item] == -1:
                    colouredlinkednodes = [i for i in self.DualGraph.Links[item] if self.PrimalGraph.map[i] >= 0]
                    if len(colouredlinkednodes) > 0:
                        # calculate distance to each linked node and take the nearest node's colour
                        distancetocolourednodes = [
                            GeneralFunctions.distancebetweenpoints(self.DualGraph.PialSurface[item],
                                                                   self.DualGraph.PialSurface[p])
                            for p in colouredlinkednodes]

                        mindifferencebyindex = min(range(len(distancetocolourednodes)),
                                                   key=lambda k: distancetocolourednodes[k])
                        # update after the for loop
                        colourupdate.append((item, colouredlinkednodes[mindifferencebyindex]))
                    else:
                        nocolour += 1
            for element in colourupdate:
                self.PrimalGraph.map[element[0]] = self.PrimalGraph.map[element[1]]
            iter += 1
            if nocolour == 0 or iter == 1000:
                break


class CouplingPoint:
    def __init__(self, node):
        self.Node = node
        self.PialSurfacePointID = None
        self.Pialpoints = []  # sampling used in the clustering
        self.NumberOfPialPoints = 0
        self.SurfaceNodes = []  # mapping after the clusting (triangles in primal graph)
        self.Area = 0
        self.Tree = None


def SplitRegionInTwo(points, PialDMat, debug=False):
    print("Splitting region in two.")

    RootsN = 2
    PointsN = len(points)
    NRootsLabels = [numpy.floor(PointsN / 2), numpy.ceil(PointsN / 2)]
    regionpoints = points
    roots = numpy.random.choice(points, 2, replace=False)
    sampledict = {sample: index for index, sample in enumerate(points)}

    if debug:
        region1 = numpy.random.choice(points, int(NRootsLabels[0]), replace=False)
        region2 = [point for point in points if point not in region1]
        return region1, region2

    # DistanceMat = dualgraph.CalculateDistanceMatRegion("", roots, regionpoints)
    # DistanceMat2 = PialDMat[[sampledict[root] for root in roots],:]
    iter = 0
    result = []
    while 1:
        DistanceMat = PialDMat[[sampledict[root] for root in roots], :]

        flat_list = [item for sublist in DistanceMat for item in sublist]
        AssignedLabelsN = [0 for _ in range(0, RootsN)]
        Labels = [-1 for _ in range(0, PointsN)]

        sortedindex = sorted(range(len(flat_list)), key=lambda k: flat_list[k])

        for index in sortedindex:
            ClusterIndex = index // PointsN
            PointIndex = index % PointsN

            if AssignedLabelsN[ClusterIndex] < NRootsLabels[ClusterIndex] and Labels[PointIndex] == -1:
                # if Labels[PointIndex] == -1:
                Labels[PointIndex] = ClusterIndex
                AssignedLabelsN[ClusterIndex] += 1

        swapiter = 0
        while 1:
            swaplist = []
            for point in range(0, PointsN):
                possiblelabels = DistanceMat[:, point]
                sortedlabels = sorted(range(len(possiblelabels)), key=lambda k: possiblelabels[k])
                currentlabel = Labels[point]
                if currentlabel != sortedlabels[0]:
                    swaplist.append(point)
            changecounter = 0
            for swap in swaplist:
                differenceafterswap = [(DistanceMat[Labels[swap], swap] + DistanceMat[Labels[i], i])
                                       - (DistanceMat[Labels[i], swap] + DistanceMat[Labels[swap], i])
                                       for i in range(0, PointsN)]
                differencebyindex = max(range(len(differenceafterswap)), key=lambda k: differenceafterswap[k])
                if differenceafterswap[differencebyindex] > 0:
                    Labels[differencebyindex], Labels[swap] = Labels[swap], Labels[differencebyindex]
                    changecounter += 1
            print("Total changes: %d" % changecounter)
            swapiter += 1
            if changecounter == 0 or swapiter == 20:
                break

        stopcrit = 0

        print("Recalculating roots.")
        clusterelements = [[i for i in range(0, len(Labels)) if Labels[i] == element] for element
                           in range(0, RootsN)]

        for index, cluster in enumerate(clusterelements):
            # if method == "euclidean":
            # newpos = numpy.mean([dualgraph.PialSurface[regionpoints[i]] for i in cluster], axis=0)
            # diff = distancebetweenpoints(roots[index], newpos)
            # stopcrit += diff
            # roots[index] = newpos

            distancematforcluster = numpy.empty([len(cluster), len(cluster)])
            for num, point in enumerate(cluster):
                for num2, point2 in enumerate(cluster):
                    distancematforcluster[num, num2] = PialDMat[point, point2]
            print(distancematforcluster)
            distances = numpy.amax(distancematforcluster, axis=0)
            newguess = min(range(0, len(distances)), key=lambda k: distances[k])
            newguessPoint = cluster[newguess]
            # new guess for root node is the Jardan Center of the cluster
            # calculate the change in position of the root
            # newpos = PialSurface[regionpoints[newguessPoint]]
            # oldpos = PialSurface[roots[index]]
            # diff = GeneralFunctions.distancebetweenpoints(oldpos, newpos)
            diff = abs(regionpoints[newguessPoint] - roots[index])
            stopcrit += diff
            roots[index] = regionpoints[newguessPoint]
            # print(roots[index])
        print("Root displacement: %f" % stopcrit)
        iter += 1
        # print(rootpositions)
        if stopcrit == 0 or iter == 20:
            result = Labels
            break

        # DistanceMat = dualgraph.CalculateDistanceMatRegion("", roots, regionpoints)

    region1 = [points[i] for i in range(0, len(result)) if result[i] == 0]
    region2 = [points[i] for i in range(0, len(result)) if result[i] == 1]
    return region1, region2
