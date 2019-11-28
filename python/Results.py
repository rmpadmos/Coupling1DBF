#!/usr/bin/python3
import time

# import sys
# sys.path.insert(0, "/drive/in-silico-trial/software/1d-blood-flow/scripts/")
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy
import vtk
from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from scipy import integrate


class Results:
    def __init__(self):
        self.TimePoints = []
        self.File = None

        self.Time = []
        self.VolumeFlowRate = []
        self.Velocity = []
        self.Pressure = []
        self.Radius = []

        self.MeanResults = []  # per vessel

        self.MeanPressurePerNode = []
        self.MeanVelocityPerNode = []
        self.MeanVolumeFlowRatePerNode = []
        self.MeanRadiusPerNode = []
        self.PulsatilityIndexPressure = []
        self.PulsatilityIndexVelocity = []

        self.ConvergenceFile = None
        self.Convergence = []

    def ClearResults(self):
        self.TimePoints = []
        self.File = None

        self.Time = []
        self.VolumeFlowRate = []
        self.Velocity = []
        self.Pressure = []
        self.Radius = []

        self.MeanResults = []  # per vessel

        self.MeanPressurePerNode = []
        self.MeanVelocityPerNode = []
        self.MeanVolumeFlowRatePerNode = []
        self.MeanRadiusPerNode = []
        self.PulsatilityIndexPressure = []
        self.PulsatilityIndexVelocity = []

        self.ConvergenceFile = None
        self.Convergence = []

    def addtimepoint(self, time):
        self.TimePoints.append(TimePoint(time))

    def LoadResults(self, datafile):
        print("Loading results.")
        self.File = datafile
        bfdatatemp = [i.strip('\n').split('\t') for i in open(datafile)]
        self.TimePoints = []

        for line in bfdatatemp:
            if "WT" in line[0]:
                self.addtimepoint(float(line[0][3:]))
            else:
                result = [float(line[1]), float(line[2]), float(line[3])]
                self.TimePoints[-1].addresults(result)

        self.SplitHeartbeats()

    def ExportResults(self, datafile="ExportedResults.dyn"):
        lengthnodes = len(self.TimePoints[0].Flow)
        with open(datafile, "w") as f:
            for index, time in enumerate(self.TimePoints):
                f.write("WT: %f\n" % time.WT)
                for i in range(0, lengthnodes):
                    f.write("%d\t%f\t%f\t%f\n" % (i, time.Flow[i], time.Pressure[i], time.Radius[i]))

    def get_results(self):
        time = [timestep.WT for timestep in self.TimePoints]

        pressurepernode = []
        radiuspernode = []
        flowpernode = []
        for index in range(0, len(self.TimePoints[0].Pressure)):
            pressurepernode.append([tp.Pressure[index] for tp in self.TimePoints])
            radiuspernode.append([tp.Radius[index] for tp in self.TimePoints])
            flowpernode.append([tp.Flow[index] for tp in self.TimePoints])

        return time, flowpernode, pressurepernode, radiuspernode

    def SplitHeartbeats(self):
        # returns results split into separate heartbeats
        time, flowpernode, pressurepernode, radiuspernode = self.get_results()

        hbends = [i for i, e in enumerate(time) if e == 0.0]
        numberhbs = len(hbends)
        hbends.append(len(time))

        hbstime = [[] for i in range(numberhbs)]
        hbsflowpernode = [[] for i in range(numberhbs)]
        hbspressurepernode = [[] for i in range(numberhbs)]
        hbsradiuspernode = [[] for i in range(numberhbs)]

        for hb in range(numberhbs):
            hbstime[hb] = time[hbends[hb]:hbends[hb + 1]]
            hbsflowpernode[hb] = [node[hbends[hb]:hbends[hb + 1]] for node in flowpernode]
            hbspressurepernode[hb] = [node[hbends[hb]:hbends[hb + 1]] for node in pressurepernode]
            hbsradiuspernode[hb] = [node[hbends[hb]:hbends[hb + 1]] for node in radiuspernode]

        self.Time = hbstime
        self.VolumeFlowRate = hbsflowpernode
        self.Pressure = hbspressurepernode
        self.Radius = hbsradiuspernode
        return hbstime, hbsflowpernode, hbspressurepernode, hbsradiuspernode

    def CorrectForDirectionEnds(self, ids):
        for index in ids:
            for heartbeat in self.VolumeFlowRate:
                heartbeat[index] = [volumeflowrate * -1 for volumeflowrate in heartbeat[index]]
            for tp in self.TimePoints:
                tp.Flow[index] *= -1

    def CalculateMeanVessel(self, ids):
        # take the last heartbeat
        positions = [node.LengthAlongVessel for node in ids]
        lengthvessel = max(positions)
        meanvolumeflowrate = integrate.simps([self.MeanVolumeFlowRatePerNode[-1][i.Number] for i in ids],
                                             positions) / lengthvessel
        meanpressure = integrate.simps([self.MeanPressurePerNode[-1][i.Number] for i in ids], positions) / lengthvessel
        meanradius = integrate.simps([self.MeanRadiusPerNode[-1][i.Number] for i in ids], positions) / lengthvessel
        # meanvelocity =integrate.simps([self.MeanVelocityPerNode[-1][i.Number] for i in ids],positions)/lengthvessel
        meanvelocity = meanvolumeflowrate / (numpy.pi * meanradius * meanradius)
        meanpipressure = integrate.simps([self.PulsatilityIndexPressure[-1][i.Number] for i in ids],
                                         positions) / lengthvessel
        meanpivelocity = integrate.simps([self.PulsatilityIndexVelocity[-1][i.Number] for i in ids],
                                         positions) / lengthvessel
        return (meanvolumeflowrate, meanpressure, meanradius, meanvelocity, meanpipressure, meanpivelocity)

    def GetMeanResults(self, vessels):
        vesselnames = [vessel.Name for vessel in vessels]
        vesselnodes = [[node for node in vessel.Nodes] for vessel in vessels]

        for index, vessel in enumerate(vesselnodes):
            meanvalues = self.CalculateMeanVessel(vessel)
            self.MeanResults.append((vesselnames[index], meanvalues, vessel))

    def ExportMeanResults(self, folder, file="ResultsPerVessel.csv"):
        with open(folder + file, 'w') as f:
            f.write(
                "VesselName,VolumeFlowrate(mL/s),Pressure(Pa),Radius(mm),Velocity(m/s),PulsatilityIndexPressure,PulsatilityIndexVelocity,VolumeFlowRate(mL/min),CardiacOutputFraction\n")
            for result in self.MeanResults:
                f.write("\"%s\"," % result[0])
                f.write("%s," % result[1][0])
                f.write("%s," % result[1][1])
                f.write("%s," % result[1][2])
                f.write("%s," % result[1][3])
                f.write("%s," % result[1][4])
                f.write("%s," % result[1][5])
                f.write("%s," % (result[1][0] * 60))
                # assume node 0 of vessel 0 is the inlet
                f.write("%s\n" % (result[1][0] / (self.MeanResults[0][1][0])))

    def SimTime(self):
        return [timepoint.WT for timepoint in self.TimePoints]

    def CalculateVelocity(self):
        for index, hb in enumerate(self.VolumeFlowRate):
            velocity = []
            for index2, flowrate in enumerate(hb):
                velres = [(flowrate[i]) / (numpy.pi * self.Radius[index][index2][i] * self.Radius[index][index2][i])
                          for i in range(0, len(flowrate))]
                velocity.append(velres)
            self.Velocity.append(velocity)

    def CalculateMeanResultsPerNode(self):
        for index, heartbeat in enumerate(self.Pressure):
            pulseindex = [(numpy.max(measurementpernode) - numpy.min(measurementpernode)) /
                          (integrate.simps(measurementpernode, self.Time[index]) / self.Time[index][-1])
                          for i, measurementpernode in enumerate(heartbeat)]
            pressure = [integrate.simps(measurementpernode, self.Time[index]) / self.Time[index][-1]
                        for i, measurementpernode in enumerate(heartbeat)]
            # velocity = [integrate.simps(measurementpernode, self.Time[index]) / self.Time[index][-1]
            #             for i, measurementpernode in enumerate(self.Velocity[index])]
            volumeflowrate = [integrate.simps(measurementpernode, self.Time[index]) / self.Time[index][-1]
                              for i, measurementpernode in enumerate(self.VolumeFlowRate[index])]
            radius = [integrate.simps(measurementpernode, self.Time[index]) / self.Time[index][-1]
                      for i, measurementpernode in enumerate(self.Radius[index])]

            velocity = [volumeflowrate[i] / (numpy.pi * radius[i] * radius[i]) for i in range(0, len(radius))]

            timeindexDiastolicPressure = [p.index(min(p)) for p in heartbeat]
            meanvelocity = [integrate.simps(measurementpernode, self.Time[index]) / self.Time[index][-1] for
                            i, measurementpernode in enumerate(self.Velocity[index])]
            pulseindexvelocity = [
                (numpy.max(measurementpernode) - measurementpernode[timeindexDiastolicPressure[i]]) / meanvelocity[i] if
                meanvelocity[i] > 0 else 0
                for i, measurementpernode in enumerate(self.Velocity[index])]

            self.MeanPressurePerNode.append(pressure)
            self.MeanVelocityPerNode.append(velocity)
            self.MeanVolumeFlowRatePerNode.append(volumeflowrate)
            self.MeanRadiusPerNode.append(radius)
            self.PulsatilityIndexPressure.append(pulseindex)
            self.PulsatilityIndexVelocity.append(pulseindexvelocity)

    def PlotPIRadius(self, figname="",figuredpi=72):
        bravasetvessels = [i for i in self.MeanResults if "gen" in i[0] and not ("Tree" in i[0])]
        Treevessels = [i for i in self.MeanResults if "Tree" in i[0]]
        Largevessels = [i for i in self.MeanResults if i not in bravasetvessels and i not in Treevessels]

        coloursvessels = numpy.zeros(len(self.MeanResults))
        nodessum = sum([len(i[2]) for i in self.MeanResults])
        coloursnodes = numpy.zeros(nodessum)
        pi = numpy.zeros(nodessum)
        pi2 = numpy.zeros(nodessum)
        r = numpy.zeros(nodessum)

        for vessel in Largevessels:
            # coloursvessels[self.MeanResults.index(vessel)] = 0
            for node in vessel[2]:
                # coloursnodes[node] = 0
                r[node.Number] = self.MeanRadiusPerNode[-1][node.Number]
                pi[node.Number] = self.PulsatilityIndexPressure[-1][node.Number]
                pi2[node.Number] = self.PulsatilityIndexVelocity[-1][node.Number]

        for vessel in Treevessels:
            coloursvessels[self.MeanResults.index(vessel)] = .5
            for node in vessel[2]:
                coloursnodes[node.Number] = 0.5
                r[node.Number] = self.MeanRadiusPerNode[-1][node.Number]
                pi[node.Number] = self.PulsatilityIndexPressure[-1][node.Number]
                pi2[node.Number] = self.PulsatilityIndexVelocity[-1][node.Number]

        for vessel in bravasetvessels:
            coloursvessels[self.MeanResults.index(vessel)] = 1
            for node in vessel[2]:
                coloursnodes[node.Number] = 1
                r[node.Number] = self.MeanRadiusPerNode[-1][node.Number]
                pi[node.Number] = self.PulsatilityIndexPressure[-1][node.Number]
                pi2[node.Number] = self.PulsatilityIndexVelocity[-1][node.Number]

        fig = plt.figure()
        DPI = fig.get_dpi()
        fig.set_size_inches(1500.0 / float(DPI), 1500.0 / float(DPI))

        axis = fig.add_subplot(221)
        axis.set_title("Pressure Pulsatility Index over Radius Per Node", fontsize=30)
        axis.set_xlabel("Radius (mm)", fontsize=30)
        axis.set_ylabel("Pulsatility Index", fontsize=30)
        axis.xaxis.set_tick_params(labelsize=25)
        axis.yaxis.set_tick_params(labelsize=25)
        axis.scatter(r, pi, c=cm.rainbow(coloursnodes))
        axis.grid()
        axis.set_xscale('log')

        axis2 = fig.add_subplot(222)
        axis2.set_title("Pressure Pulsatility Index over Radius Per Vessel", fontsize=30)
        axis2.set_xlabel("Radius (mm)", fontsize=30)
        axis2.set_ylabel("Pulsatility Index", fontsize=30)
        axis2.xaxis.set_tick_params(labelsize=25)
        axis2.yaxis.set_tick_params(labelsize=25)
        axis2.grid()
        axis2.set_xscale('log')

        axis3 = fig.add_subplot(223)
        axis3.set_title("Velocity Pulsatility Index over Radius Per Vessel", fontsize=30)
        axis3.set_xlabel("Radius (mm)", fontsize=30)
        axis3.set_ylabel("Pulsatility Index", fontsize=30)
        axis3.xaxis.set_tick_params(labelsize=25)
        axis3.yaxis.set_tick_params(labelsize=25)
        axis3.grid()
        axis3.set_xscale('log')

        axis4 = fig.add_subplot(224)
        axis4.set_title("Velocity Pulsatility Index over Radius Per Vessel", fontsize=30)
        axis4.set_xlabel("Radius (mm)", fontsize=30)
        axis4.set_ylabel("Pulsatility Index", fontsize=30)
        axis4.xaxis.set_tick_params(labelsize=25)
        axis4.yaxis.set_tick_params(labelsize=25)
        axis4.grid()
        axis4.set_xscale('log')

        meanr = [i[1][2] for i in self.MeanResults]
        meanpi = [i[1][4] for i in self.MeanResults]
        meanpi2 = [i[1][5] for i in self.MeanResults]

        axis.scatter(r, pi, c=cm.rainbow(coloursnodes))
        axis2.scatter(meanr, meanpi, c=cm.rainbow(coloursvessels))
        axis3.scatter(r, pi2, c=cm.rainbow(coloursnodes))
        axis4.scatter(meanr, meanpi2, c=cm.rainbow(coloursvessels))

        fig.canvas.draw_idle()
        fig.tight_layout()
        if figname:
            print("Saving figure.")
            fig.savefig(figname,dpi=figuredpi)
        plt.show()


    def PlotPIRadiusVessel(self, figname="",figuredpi=72):
        bravasetvessels = [i for i in self.MeanResults if "gen" in i[0] and not ("Tree" in i[0])]
        Treevessels = [i for i in self.MeanResults if "Tree" in i[0]]
        Largevessels = [i for i in self.MeanResults if i not in bravasetvessels and i not in Treevessels]

        coloursvessels = numpy.zeros(len(self.MeanResults))
        nodessum = sum([len(i[2]) for i in self.MeanResults])
        coloursnodes = numpy.zeros(nodessum)
        pi = numpy.zeros(nodessum)
        pi2 = numpy.zeros(nodessum)
        r = numpy.zeros(nodessum)

        for vessel in Largevessels:
            # coloursvessels[self.MeanResults.index(vessel)] = 0
            for node in vessel[2]:
                # coloursnodes[node] = 0
                r[node.Number] = self.MeanRadiusPerNode[-1][node.Number]
                pi[node.Number] = self.PulsatilityIndexPressure[-1][node.Number]
                pi2[node.Number] = self.PulsatilityIndexVelocity[-1][node.Number]

        for vessel in Treevessels:
            coloursvessels[self.MeanResults.index(vessel)] = .5
            for node in vessel[2]:
                coloursnodes[node.Number] = 0.5
                r[node.Number] = self.MeanRadiusPerNode[-1][node.Number]
                pi[node.Number] = self.PulsatilityIndexPressure[-1][node.Number]
                pi2[node.Number] = self.PulsatilityIndexVelocity[-1][node.Number]

        for vessel in bravasetvessels:
            coloursvessels[self.MeanResults.index(vessel)] = 1
            for node in vessel[2]:
                coloursnodes[node.Number] = 1
                r[node.Number] = self.MeanRadiusPerNode[-1][node.Number]
                pi[node.Number] = self.PulsatilityIndexPressure[-1][node.Number]
                pi2[node.Number] = self.PulsatilityIndexVelocity[-1][node.Number]

        fig = plt.figure()
        DPI = fig.get_dpi()
        fig.set_size_inches(2000.0 / float(DPI), 1000.0 / float(DPI))

        axis1 = fig.add_subplot(121)
        axis1.set_title("Pressure Pulsatility Index over Radius Per Tree Node", fontsize=30)
        axis1.set_xlabel("Radius (mm)", fontsize=30)
        axis1.set_ylabel("Pulsatility Index", fontsize=30)
        axis1.xaxis.set_tick_params(labelsize=25)
        axis1.yaxis.set_tick_params(labelsize=25)
        axis1.set_xscale('log')
        axis1.grid()

        axis2 = fig.add_subplot(122)
        axis2.set_title("Velocity Pulsatility Index over Radius Per Tree Node", fontsize=30)
        axis2.set_xlabel("Radius (mm)", fontsize=30)
        axis2.set_ylabel("Pulsatility Index", fontsize=30)
        axis2.xaxis.set_tick_params(labelsize=25)
        axis2.yaxis.set_tick_params(labelsize=25)
        axis2.set_xscale('log')
        axis2.grid()

        treenames = set([i[0].split("Tree_")[0] for i in Treevessels])
        colours = numpy.linspace(0, 1, len(treenames))
        radius = []
        pulseindex = []
        pulseindex2 = []
        pointcolour = []
        for index, name in enumerate(treenames):
            colour = colours[index]
            vessels = [i for i in Treevessels if name in i[0]]
            for vessel in vessels:
                for node in vessel[2]:
                    radius.append(self.MeanRadiusPerNode[-1][node.Number])
                    pulseindex.append(self.PulsatilityIndexPressure[-1][node.Number])
                    pulseindex2.append(self.PulsatilityIndexVelocity[-1][node.Number])
                    pointcolour.append(colour)

        axis1.scatter(radius, pulseindex, c=cm.rainbow(pointcolour))
        axis2.scatter(radius, pulseindex2, c=cm.rainbow(pointcolour))

        fig.canvas.draw_idle()
        fig.tight_layout()
        if figname:
            print("Saving figure.")
            fig.savefig(figname,dpi=figuredpi)
        plt.show()

    def plotallprofiles(self, ids, name="",figuredpi=72):
        self.fig = plt.figure()
        DPI = self.fig.get_dpi()
        self.fig.set_size_inches(2000.0 / float(DPI), 1000.0 / float(DPI))
        self.ax1 = self.fig.add_subplot(131)
        self.ax2 = self.fig.add_subplot(132)
        self.ax2_c = self.ax2.twinx()
        self.ax3 = self.fig.add_subplot(133)
        # self.ax4 = self.fig.add_subplot(144)

        self.ax1.set_aspect('auto')
        self.ax2.set_aspect('auto')
        self.ax3.set_aspect('auto')
        # self.ax4.set_aspect('auto')
        self.ax1.set_title("Flow rate over time", fontsize=30)
        self.ax2.set_title("Pressure over time", fontsize=30)
        self.ax3.set_title("Radius over time", fontsize=30)
        # self.ax4.set_title("Velocity over time", fontsize=30)
        self.ax1.set_xlabel("Time (s)", fontsize=30)
        self.ax2.set_xlabel("Time (s)", fontsize=30)
        self.ax3.set_xlabel("Time (s)", fontsize=30)
        # self.ax4.set_xlabel("Time (s)", fontsize=30)
        self.ax1.set_ylabel("Flow rate (mL/s)", fontsize=30)
        self.ax2.set_ylabel("Pressure (Pa)", fontsize=30)
        self.ax2_c.set_ylabel("Pressure (mmHg)", fontsize=30)
        # self.ax3.set_ylabel("Radius Change (%)", fontsize=30)
        self.ax3.set_ylabel("Radius (mm)", fontsize=30)
        # self.ax4.set_ylabel("Velocity (m/s)", fontsize=30)

        self.ax1.xaxis.set_tick_params(labelsize=25)
        self.ax2.xaxis.set_tick_params(labelsize=25)
        self.ax3.xaxis.set_tick_params(labelsize=25)
        # self.ax4.xaxis.set_tick_params(labelsize=25)
        self.ax1.yaxis.set_tick_params(labelsize=25)
        self.ax2.yaxis.set_tick_params(labelsize=25)
        self.ax2_c.yaxis.set_tick_params(labelsize=25)
        self.ax3.yaxis.set_tick_params(labelsize=25)
        # self.ax4.yaxis.set_tick_params(labelsize=25)

        linestyles = ['-', '--', '-.', ':']
        for linestyle, hb in zip(linestyles,ids):
            self.ax1.plot(self.Time[0], self.VolumeFlowRate[-1][hb], linewidth=4.0,linestyle=linestyle)
            self.ax2.plot(self.Time[0], self.Pressure[-1][hb], linewidth=4.0,linestyle=linestyle)
            meanr = integrate.simps(self.Radius[-1][hb], self.Time[0]) / self.Time[0][-1]
            radius = [100 * (i - meanr) / meanr for i in self.Radius[-1][hb]]
            # self.ax3.plot(self.Time[0], radius)

            self.ax3.plot(self.Time[0], self.Radius[-1][hb], linewidth=4.0,linestyle=linestyle)
            # self.ax4.plot(self.Time[0], self.Velocity[-1][hb])

        if len(ids) > 1:
            self.ax1.set_yscale('log')
            # self.ax4.set_yscale('log')
        y1, y2 = self.ax2.get_ylim()
        self.ax2_c.set_ylim(y1 * 0.007500617, y2 * 0.007500617)
        self.ax1.grid()
        self.ax2.grid()
        self.ax3.grid()
        # self.ax4.grid()
        self.fig.canvas.draw_idle()
        self.fig.tight_layout()
        if name:
            self.fig.savefig(name,dpi=figuredpi)
        plt.show()

    def PlotDistancePI(self, figname="",figuredpi=72):
        bravasetvessels = [i for i in self.MeanResults if "gen" in i[0] and not ("Tree" in i[0])]
        Treevessels = [i for i in self.MeanResults if "Tree" in i[0]]
        Largevessels = [i for i in self.MeanResults if i not in bravasetvessels and i not in Treevessels]

        coloursvessels = numpy.zeros(len(self.MeanResults))
        nodessum = sum([len(i[2]) for i in self.MeanResults])
        coloursnodes = numpy.zeros(nodessum)
        pi = numpy.zeros(nodessum)
        pi2 = numpy.zeros(nodessum)
        DistancefromHeart = numpy.zeros(nodessum)

        for vessel in Largevessels:
            for node in vessel[2]:
                # coloursnodes[node] = 0
                DistancefromHeart[node.Number] = self.DistancefromHeart[node.Number]
                pi[node.Number] = self.PulsatilityIndexPressure[-1][node.Number]
                pi2[node.Number] = self.PulsatilityIndexVelocity[-1][node.Number]

        for vessel in Treevessels:
            coloursvessels[self.MeanResults.index(vessel)] = .5
            for node in vessel[2]:
                coloursnodes[node.Number] = 0.5
                DistancefromHeart[node.Number] = self.DistancefromHeart[node.Number]
                pi[node.Number] = self.PulsatilityIndexPressure[-1][node.Number]
                pi2[node.Number] = self.PulsatilityIndexVelocity[-1][node.Number]

        for vessel in bravasetvessels:
            coloursvessels[self.MeanResults.index(vessel)] = 1
            for node in vessel[2]:
                coloursnodes[node.Number] = 1
                DistancefromHeart[node.Number] = self.DistancefromHeart[node.Number]
                pi[node.Number] = self.PulsatilityIndexPressure[-1][node.Number]
                pi2[node.Number] = self.PulsatilityIndexVelocity[-1][node.Number]

                # treenames = set([i[0].split("Tree_")[0] for i in Treevessels])
        # colours = numpy.linspace(0, 1, len(treenames))
        # DistancefromHeart = []
        # pulseindex = []
        # pulseindex2 = []
        # pointcolour = []
        # for index, name in enumerate(treenames):
        #     colour = colours[index]
        #     vessels = [i for i in Treevessels if name in i[0]]
        #     for vessel in vessels:
        #         for node in vessel[2]:
        #             DistancefromHeart.append(self.MeanRadiusPerNode[-1][node.Number])
        #             pulseindex.append(self.PulsatilityIndexPressure[-1][node.Number])
        #             pulseindex2.append(self.PulsatilityIndexVelocity[-1][node.Number])
        #             pointcolour.append(colour)

        fig = plt.figure()
        DPI = fig.get_dpi()
        fig.set_size_inches(2000.0 / float(DPI), 1000.0 / float(DPI))

        axis1 = fig.add_subplot(121)
        axis1.set_title("Pressure Pulsatility Index as a function of distance", fontsize=30)
        axis1.set_xlabel("Distance from the heart (mm)", fontsize=30)
        axis1.set_ylabel("Pulsatility Index", fontsize=30)
        axis1.xaxis.set_tick_params(labelsize=25)
        axis1.yaxis.set_tick_params(labelsize=25)
        # axis1.set_xscale('log')
        axis1.grid()

        axis1.scatter(DistancefromHeart, pi, c=cm.rainbow(coloursnodes))

        axis2 = fig.add_subplot(122)
        axis2.set_title("Velocity Pulsatility Index as a function of distance", fontsize=30)
        axis2.set_xlabel("Distance from the heart (mm)", fontsize=30)
        axis2.set_ylabel("Pulsatility Index", fontsize=30)
        axis2.xaxis.set_tick_params(labelsize=25)
        axis2.yaxis.set_tick_params(labelsize=25)
        # axis1.set_xscale('log')
        axis2.grid()
        axis2.scatter(DistancefromHeart, pi2, c=cm.rainbow(coloursnodes))

        fig.canvas.draw_idle()
        fig.tight_layout()
        if figname:
            fig.savefig(figname,dpi=figuredpi)
        plt.show()
        return fig

    def AddResultsPerNodeToFile(self, file):
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(file)
        reader.Update()
        data = reader.GetOutput()

        names = ["Mean VolumeFlowrate (mL/s)",
                 "Mean Pressure (Pa)",
                 "Mean Radius (mm)",
                 "Mean Velocity (m/s)",
                 "Mean Pressure Pulsatility Index",
                 "Mean Velocity Pulsatility Index"]

        results = [self.MeanVolumeFlowRatePerNode[-1],
                   self.MeanPressurePerNode[-1],
                   self.MeanRadiusPerNode[-1],
                   self.MeanVelocityPerNode[-1],
                   self.PulsatilityIndexPressure[-1],
                   self.PulsatilityIndexVelocity[-1]]

        for resultnodes, name in zip(results, names):
            resultdata = vtk.vtkFloatArray()
            resultdata.SetNumberOfComponents(1)
            resultdata.SetName(name)
            [resultdata.InsertNextValue(result) for result in resultnodes]
            data.GetPointData().AddArray(resultdata)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(file)
        writer.SetInputData(data)
        writer.Write()

    def AddResultsPerVesselToFile(self, file):
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(file)
        reader.Update()
        data = reader.GetOutput()

        names = ["Mean VolumeFlowrate (mL/s)",
                 "Mean Pressure (Pa)",
                 "Mean Radius (mm)",
                 "Mean Velocity (m/s)",
                 "Mean Pressure Pulsatility Index",
                 "Mean Velocity Pulsatility Index"]

        for index, name in enumerate(names):
            resultdata = vtk.vtkFloatArray()
            resultdata.SetNumberOfComponents(1)
            resultdata.SetName(name)
            [resultdata.InsertNextValue(result[1][index]) for result in self.MeanResults]
            data.GetCellData().AddArray(resultdata)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(file)
        writer.SetInputData(data)
        writer.Write()


class TimePoint:
    def __init__(self, time):
        self.WT = time
        self.Flow = []
        self.Pressure = []
        self.Radius = []

    def addresults(self, results):
        self.Flow.append(results[0])
        self.Pressure.append(results[1])
        self.Radius.append(results[2])


class ResultsWindow(QtWidgets.QMainWindow):
    def __init__(self, patient):
        super(ResultsWindow, self).__init__()

        self.Patient = patient
        self.Results = patient.Results
        self.Vessels = patient.Topology.Vessels

        self.main_widget = QtWidgets.QWidget(self)
        self.fs_watcher = QtCore.QFileSystemWatcher([self.Results.File])
        self.fs_watcher.fileChanged.connect(self.file_changed)

        self.fig = Figure()
        self.ax1 = self.fig.add_subplot(141)
        self.ax2 = self.fig.add_subplot(142)
        self.ax2_c = self.ax2.twinx()
        self.ax3 = self.fig.add_subplot(143)
        self.ax4 = self.fig.add_subplot(144)
        self.canvas = FigureCanvasQTAgg(self.fig)

        self.ax1.set_aspect('auto')
        self.ax2.set_aspect('auto')
        self.ax3.set_aspect('auto')

        self.canvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                  QtWidgets.QSizePolicy.Expanding)
        self.canvas.updateGeometry()

        # Vessel selector
        self.VesselSelector = QtWidgets.QComboBox()
        self.VesselSelector.setStyleSheet("QComboBox { combobox-popup: 0; }")
        self.VesselSelector.setMaxVisibleItems(10)
        self.VesselSelector.addItems([vessel.Name for vessel in self.Vessels])

        # Node selector
        self.NodeSelector = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.NodeSelector.setTickPosition(QtWidgets.QSlider.TicksBelow)

        # Save figure button
        self.SaveButton = QtWidgets.QPushButton("Save")
        # Exit Button
        self.ExitButton = QtWidgets.QPushButton("Exit")

        # Linking stuff to actions
        self.VesselSelector.currentIndexChanged.connect(self.updateVessel)
        self.NodeSelector.valueChanged.connect(self.updatefig)
        self.SaveButton.clicked.connect(self.savefig)
        self.ExitButton.clicked.connect(self.exit)

        # Add widgets to the layout
        self.layout = QtWidgets.QGridLayout(self.main_widget)
        self.layout.addWidget(QtWidgets.QLabel("Select Vessel"))
        self.layout.addWidget(self.VesselSelector)
        self.layout.addWidget(self.NodeSelector)
        self.layout.addWidget(self.canvas, 0, 0)

        self.groupBox = QtWidgets.QGroupBox()
        self.Buttons = QtWidgets.QGridLayout()

        # self.Buttons.addWidget(self.NodeSelector)
        self.Buttons.addWidget(self.SaveButton)
        self.Buttons.addWidget(self.ExitButton)
        self.groupBox.setLayout(self.Buttons)

        self.layout.addWidget(self.groupBox, 0, 1)
        self.setCentralWidget(self.main_widget)
        self.resize(2000, 800)
        self.updateVessel()
        self.updatefig()
        self.show()

    def updateVessel(self):
        # print("Changed Vessel")
        CurrentVessel = self.VesselSelector.currentText()
        StartNode, EndNode = [[vessel.Nodes[0], vessel.Nodes[-1]]
                              for vessel in self.Vessels if vessel.Name == CurrentVessel][0]
        self.NodeSelector.setMinimum(StartNode.Number)
        self.NodeSelector.setMaximum(EndNode.Number)

    @QtCore.pyqtSlot(str)
    def file_changed(self):
        print('File Changed.')
        time.sleep(.100)  # wait for writing to finish
        self.Patient.Results.LoadResults(self.Results.File)
        self.Patient.CorrectforDirectionVectorOutlets()
        self.Patient.Results.CalculateVelocity()
        self.Patient.Results.CalculateMeanResultsPerNode()
        self.updatefig()

    def updatefig(self):
        self.ax1.clear()
        self.ax2.clear()
        self.ax2_c.clear()
        self.ax3.clear()
        self.ax4.clear()

        CurrentVessel = self.VesselSelector.currentText()
        CurrentNode = self.NodeSelector.value()
        # print('Vessel: ' + CurrentVessel + ', Node: ' + str(CurrentNode))

        self.ax1.set_title("Flow rate over time", fontsize=30)
        self.ax2.set_title("Pressure over time", fontsize=30)
        self.ax3.set_title("Radius over time", fontsize=30)
        self.ax4.set_title("Velocity over time", fontsize=30)

        self.ax1.set_xlabel("Time (s)", fontsize=30)
        self.ax2.set_xlabel("Time (s)", fontsize=30)
        self.ax3.set_xlabel("Time (s)", fontsize=30)
        self.ax4.set_xlabel("Time (s)", fontsize=30)

        self.ax1.set_ylabel("Flow rate (mL/s)", fontsize=30)
        self.ax2.set_ylabel("Pressure (Pa)", fontsize=30)
        self.ax2_c.set_ylabel("Pressure (mmHg)", fontsize=30)
        self.ax3.set_ylabel("Radius (mm)", fontsize=30)
        self.ax4.set_ylabel("Velocity (m/s)", fontsize=30)

        self.ax2.callbacks.connect("ylim_changed", self.PressureTommhg)

        self.ax1.xaxis.set_tick_params(labelsize=25)
        self.ax2.xaxis.set_tick_params(labelsize=25)
        self.ax3.xaxis.set_tick_params(labelsize=25)
        self.ax4.xaxis.set_tick_params(labelsize=25)

        self.ax1.yaxis.set_tick_params(labelsize=25)
        self.ax2.yaxis.set_tick_params(labelsize=25)
        self.ax2_c.yaxis.set_tick_params(labelsize=25)
        self.ax3.yaxis.set_tick_params(labelsize=25)
        self.ax4.yaxis.set_tick_params(labelsize=25)

        # print(CurrentNode)
        for hb in range(0, len(self.Results.VolumeFlowRate)):
            self.ax1.plot(self.Results.Time[0], self.Results.VolumeFlowRate[hb][CurrentNode])
            self.ax2.plot(self.Results.Time[0], self.Results.Pressure[hb][CurrentNode])
            self.ax3.plot(self.Results.Time[0], self.Results.Radius[hb][CurrentNode])
            self.ax4.plot(self.Results.Time[0], self.Results.Velocity[hb][CurrentNode])

        self.ax1.text(0.0, 0.0, str(self.Results.MeanVolumeFlowRatePerNode[-1][CurrentNode]),
                      fontsize=25, transform=self.ax1.transAxes)
        self.ax2.text(0.0, 0.0, str(self.Results.MeanPressurePerNode[-1][CurrentNode]),
                      fontsize=25, transform=self.ax2.transAxes)
        self.ax3.text(0.0, 0.0, str(self.Results.MeanRadiusPerNode[-1][CurrentNode]),
                      fontsize=25, transform=self.ax3.transAxes)
        self.ax4.text(0.0, 0.0, str(self.Results.MeanVelocityPerNode[-1][CurrentNode]),
                      fontsize=25, transform=self.ax4.transAxes)

        self.ax1.grid()
        self.ax2.grid()
        self.ax3.grid()
        self.ax4.grid()
        self.fig.canvas.draw_idle()
        self.fig.tight_layout()

    def savefig(self,figuredpi=72):
        # print("Save Button Pressed")
        node = self.NodeSelector.value()
        path = str(QtWidgets.QFileDialog.getExistingDirectory(directory="Select Directory"))

        self.fig.savefig(path + '/Figure Node:' + str(node) + '.png',dpi=figuredpi)

        with open(path + '/Node: ' + str(node) + '.csv', 'w') as f:
            f.write("time(s),Flowrate(mL/s),Pressure(Pa),Radius(m/s),Velocity(m/s)\n")
            for i in range(len(self.Results.Time[-1])):
                f.write("%s," % self.Results.Time[-1][i])
                f.write("%s," % self.Results.VolumeFlowRate[-1][node][i])
                f.write("%s," % self.Results.Pressure[-1][node][i])
                f.write("%s," % self.Results.Radius[-1][node][i])
                f.write("%s\n" % self.Results.Velocity[-1][node][i])

    def PressureTommhg(self, ax2):
        y1, y2 = ax2.get_ylim()
        self.ax2_c.set_ylim(y1 * 0.007500617, y2 * 0.007500617)

    def exit(self):
        # print("Exit Button Pressed")
        self.close()
