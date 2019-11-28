#!/usr/bin/python3
import copy
from multiprocessing.dummy import Pool as ThreadPool
import sys
from os.path import isfile, dirname, realpath, exists
sys.path.append(realpath(dirname(__file__)) )
import math
import numpy
import scipy
from scipy.optimize import newton, fsolve

import Results


def flux_function(self, Area):
    return Area * self.charback + 4 * numpy.power(Area, 1.25) * numpy.sqrt(
        self.Beta / 2.0 / self.Density) - self.CurrentFlowrate


def chrt_function(self, Area):
    chrt_frw_right = self.CurrentVolumeFlowRate / Area + 4 * numpy.power(Area, 0.25) * numpy.sqrt(
        self.Beta / 2.0 / self.Density)
    return self.chrt_frw_left - chrt_frw_right


def thickness(r):  # radius defined in mm
    """
    Emperical function to calculate vessel thickness given the radius.
    """
    radius = r * 1e-3
    h1 = 0.2802
    h2 = -5.053 * 1e2
    h3 = 0.1324
    h4 = -0.1114 * 1e2
    h = radius * (h1 * math.exp(h2 * radius) + h3 * math.exp(h4 * radius))
    return h * 1e3  # mm


def SolveOutletRCR(args):
    Outletnode, PreviousNode = args

    Outletnode.CurrentVolumeFlowRate = Outletnode.CurrentLumen * Outletnode.CurrentVelocity
    dQdt = (Outletnode.CurrentVolumeFlowRate - Outletnode.PreviousVolumeFlowRate) / Outletnode.TimeStep

    Outletnode.CurrentPressure = (Outletnode.PreviousVolumeFlowRate * (1 + Outletnode.R1 / Outletnode.R2) +
                                  Outletnode.C * Outletnode.R1 * dQdt + Outletnode.OutPressure / Outletnode.R2 +
                                  Outletnode.PreviousPressure * Outletnode.C / Outletnode.TimeStep) / \
                                 (Outletnode.C / Outletnode.TimeStep + 1 / Outletnode.R2)

    Outletnode.chrt_frw_left = Outletnode.PreviousVelocity + 4 * numpy.power(PreviousNode.CurrentLumen, 0.25) * \
                               numpy.sqrt(Outletnode.Beta / 2.0 / Outletnode.Density)
    dfdx = lambda x: (Outletnode.chrt_function(x + 1e-10) - Outletnode.chrt_function(x)) / 1e-10
    Outletnode.CurrentLumen = newton(Outletnode.chrt_function, PreviousNode.CurrentLumen, fprime=dfdx, tol=1e-12)
    Outletnode.CurrentVelocity = Outletnode.CurrentVolumeFlowRate / Outletnode.CurrentLumen
    Outletnode.PreviousVolumeFlowRate = Outletnode.CurrentVolumeFlowRate


def SolveInlet(args):
    Inlet, NextNode = args
    Inlet.CurrentFlowrate = Inlet.FlowRateFunction(Inlet.TimeStep * Inlet.Iter)
    dfdx = lambda x: (Inlet.flux_function(x + Inlet.Lumen0 * 1e-6) - Inlet.flux_function(x)) / (Inlet.Lumen0 * 1e-6)

    Inlet.CurrentLumen = newton(Inlet.flux_function, Inlet.CurrentLumen, fprime=dfdx, tol=1e-12)
    U = Inlet.CurrentFlowrate / Inlet.CurrentLumen

    dx = abs(NextNode.LengthAlongVessel - Inlet.LengthAlongVessel)
    dUdt = (U - Inlet.CurrentVelocity) / Inlet.TimeStep
    dUdZ = (NextNode.CurrentVelocity - U) / dx
    visc_term = (Inlet.Friction * Inlet.Viscosity * numpy.pi * U / Inlet.Density / Inlet.CurrentLumen)
    Inlet.CurrentPressure = NextNode.CurrentPressure + (dUdt + U * dUdZ + visc_term) * dx * Inlet.Density

    Inlet.CurrentVelocity = U
    Inlet.charback = Inlet.CurrentVelocity - 4 * numpy.power(NextNode.CurrentLumen, 0.25) * numpy.sqrt(
        Inlet.Beta / 2.0 / Inlet.Density)


class Bifurcation:
    def __init__(self, nodes, signs):
        self.Nodes = nodes
        self.Signs = signs

    def funcs(self, values):
        functionvalues = []
        v = values[:len(self.Nodes)]
        l = values[len(self.Nodes):]
        # print(l)
        sum_flux = 0
        for index, node in enumerate(self.Nodes):
            sum_flux += v[index] * l[index] * self.Signs[index]
            wave_speed = 4 * (
                    numpy.power(node.CurrentLumen, 0.25) * numpy.sqrt(node.Beta / 2.0 / node.Density) - node.c_0)
            if v[index] * self.Signs[index] > 0:
                chrt_f = abs(node.CurrentVelocity) + wave_speed
                function = abs(v[index]) + 4 * (
                        numpy.power(l[index], 0.25) * numpy.sqrt(node.Beta / 2.0 / node.Density) - node.c_0) - chrt_f
            else:
                chrt_b = abs(node.CurrentVelocity) - wave_speed
                function = abs(v[index]) - 4 * (
                        numpy.power(l[index], 0.25) * numpy.sqrt(node.Beta / 2.0 / node.Density) - node.c_0) - chrt_b
            functionvalues.append(function)

        functionvalues.append(sum_flux)

        v0 = v[0]
        p0 = self.Nodes[0].Beta * (numpy.sqrt(l[0]) - numpy.sqrt(self.Nodes[0].Lumen0)) + self.Nodes[
            0].DiastolicPressure
        for index in range(1, len(self.Nodes)):
            p = self.Nodes[index].Beta * (
                    numpy.sqrt(l[index]) - numpy.sqrt(self.Nodes[index].Lumen0)) + self.Nodes[index].DiastolicPressure
            function = self.Nodes[index].Density * (v0 * v0 - v[index] * v[index]) / 2 + p0 - p
            # function = p0 - p
            functionvalues.append(function)
        return numpy.array(functionvalues).reshape(len(values))

    def Jac(self, values):
        current = self.funcs(values)
        dx = 1e-12
        N = len(values)
        jac = scipy.zeros((N, N))

        for index in range(0, N):
            currentval = current[index]
            for index2 in range(0, N):
                values[index2] = values[index2] + dx
                currentvalue2 = self.funcs(values)[index]
                jac[index, index2] = (currentvalue2 - currentval) / dx
                values[index2] = values[index2] - dx
        return jac

    def SolveBif(self):
        init_X = [node.CurrentVelocity for node in self.Nodes]
        init_X += [node.CurrentLumen for node in self.Nodes]
        init_X = numpy.array(init_X).reshape(2 * len(self.Nodes), 1)
        solution = fsolve(self.funcs, init_X, xtol=1e-12, epsfcn=1e-12)

        for index, node in enumerate(self.Nodes):
            node.CurrentVelocity = solution[index]
            node.CurrentLumen = solution[index + len(self.Nodes)]
            node.CurrentPressure = node.Beta * (
                    numpy.sqrt(node.CurrentLumen) - numpy.sqrt(node.Lumen0)) + node.RefPressure


def SolvingBif(args):
    Bif = args
    Bif.SolveBif()


def SolveVessel(vessel):
    # vessel, dt, friction, density, viscosity = args
    length = len(vessel.Nodes)

    velocity_pred = copy.deepcopy(vessel.CurrentVelocity)
    lumen_sq_pred = copy.deepcopy(vessel.CurrentLumen)
    pressure_pred = copy.deepcopy(vessel.CurrentPressure)

    for number in range(1, length - 1):
        dx = abs(vessel.Nodes[number].LengthAlongVessel - vessel.Nodes[number - 1].LengthAlongVessel)

        velocity_pred[number] = vessel.CurrentVelocity[number] - vessel.Nodes[number].TimeStep / dx * \
                                ((vessel.CurrentVelocity[number] * vessel.CurrentVelocity[number] / 2 +
                                  vessel.CurrentPressure[number] / vessel.Nodes[number].Density) -
                                 (vessel.CurrentVelocity[number - 1] * vessel.CurrentVelocity[number - 1] / 2 +
                                  vessel.CurrentPressure[number - 1] / vessel.Nodes[number].Density))

        velocity_pred[number] = velocity_pred[number] - 1.0 / vessel.Nodes[number].Density / vessel.CurrentLumen[
            number] * \
                                (vessel.Nodes[number].Friction * vessel.Nodes[number].Viscosity * numpy.pi *
                                 velocity_pred[
                                     number]) * vessel.Nodes[number].TimeStep

        lumen_sq_pred[number] = vessel.CurrentLumen[number] - vessel.Nodes[number].TimeStep / dx * (
                vessel.CurrentLumen[number] *
                vessel.CurrentVelocity[number] -
                vessel.CurrentLumen[number - 1] *
                vessel.CurrentVelocity[number - 1])

        pressure_pred[number] = vessel.Nodes[number].RefPressure + (
                vessel.Nodes[number].Beta * (
                numpy.sqrt(lumen_sq_pred[number]) - numpy.sqrt(vessel.Nodes[number].Lumen0)))
        # print(1)

    for number in range(1, length - 1):
        dx = abs(vessel.Nodes[number].LengthAlongVessel - vessel.Nodes[number + 1].LengthAlongVessel)
        vessel.CurrentVelocity[number] = (vessel.CurrentVelocity[number] + velocity_pred[number]) / 2 - \
                                         vessel.Nodes[
                                             number].TimeStep / dx / 2 * \
                                         ((velocity_pred[number + 1] * velocity_pred[number + 1] / 2 +
                                           pressure_pred[number + 1] / vessel.Nodes[number].Density) -
                                          (velocity_pred[number] * velocity_pred[number] / 2 +
                                           pressure_pred[number] / vessel.Nodes[number].Density))

        vessel.CurrentVelocity[number] = vessel.CurrentVelocity[number] - 1.0 / 2.0 / vessel.Nodes[
            number].Density / \
                                         vessel.CurrentLumen[number] * (
                                                 vessel.Nodes[number].Friction * vessel.Nodes[
                                             number].Viscosity * numpy.pi *
                                                 vessel.CurrentVelocity[number]) * vessel.Nodes[
                                             number].TimeStep
        vessel.CurrentLumen[number] = (vessel.CurrentLumen[number] + lumen_sq_pred[number]) / 2 - \
                                      vessel.Nodes[
                                          number].TimeStep / dx / 2 * (
                                              lumen_sq_pred[number + 1] * velocity_pred[number + 1] -
                                              lumen_sq_pred[number] * velocity_pred[number])

        vessel.CurrentPressure[number] = vessel.Nodes[number].RefPressure + (
                vessel.Nodes[number].Beta * (
                numpy.sqrt(vessel.CurrentLumen[number]) - numpy.sqrt(vessel.Nodes[number].Lumen0)))

    vessel.CurrentVelocity[0] = vessel.CurrentVelocity[1]
    vessel.CurrentLumen[0] = vessel.CurrentLumen[1]
    vessel.CurrentPressure[0] = vessel.CurrentPressure[1]

    vessel.CurrentVelocity[-1] = vessel.CurrentVelocity[-2]
    vessel.CurrentLumen[-1] = vessel.CurrentLumen[-2]
    vessel.CurrentPressure[-1] = vessel.CurrentPressure[-2]

    for number, node in enumerate(vessel.Nodes):
        node.CurrentVelocity = vessel.CurrentVelocity[number]
        node.CurrentLumen = vessel.CurrentLumen[number]
        node.CurrentPressure = vessel.CurrentPressure[number]


def flowrate(t, Period, scaling):
    """
    Function from Boileau2015 (scaled inflow from common carotid)
    scaling 1 found give 1ml/s
    """
    T = Period
    flow = (6.5 + 3.294 * numpy.sin(2 * numpy.pi * t / T - 0.023974) + 1.9262 * numpy.sin(
        4 * numpy.pi * t / T - 1.1801) - 1.4219 * numpy.sin(6 * numpy.pi * t / T + 0.92701) - 0.66627 * numpy.sin(
        8 * numpy.pi * t / T - 0.24118) - 0.33933 * numpy.sin(10 * numpy.pi * t / T - 0.27471) - 0.37914 * numpy.sin(
        12 * numpy.pi * t / T - 1.0557) + 0.22396 * numpy.sin(14 * numpy.pi * t / T + 1.22) + 0.1507 * numpy.sin(
        16 * numpy.pi * t / T + 1.0984) + 0.18735 * numpy.sin(18 * numpy.pi * t / T + 0.067483) + 0.038625 * numpy.sin(
        20 * numpy.pi * t / T + 0.22262) + 0.012643 * numpy.sin(
        22 * numpy.pi * t / T - 0.10093) - 0.0042453 * numpy.sin(24 * numpy.pi * t / T - 1.1044) - 0.012781 * numpy.sin(
        26 * numpy.pi * t / T - 1.3739) + 0.014805 * numpy.sin(28 * numpy.pi * t / T + 1.2797) + 0.012249 * numpy.sin(
        30 * numpy.pi * t / T + 0.80827) + 0.0076502 * numpy.sin(
        32 * numpy.pi * t / T + 0.40757) + 0.0030692 * numpy.sin(34 * numpy.pi * t / T + 0.195) - 0.0012271 * numpy.sin(
        36 * numpy.pi * t / T - 1.1371) - 0.0042581 * numpy.sin(
        38 * numpy.pi * t / T - 0.92102) - 0.0069785 * numpy.sin(
        40 * numpy.pi * t / T - 1.2364) + 0.0085652 * numpy.sin(42 * numpy.pi * t / T + 1.4539) + 0.0081881 * numpy.sin(
        44 * numpy.pi * t / T + 0.89599) + 0.0056549 * numpy.sin(
        46 * numpy.pi * t / T + 0.17623) + 0.0026358 * numpy.sin(
        48 * numpy.pi * t / T - 1.3003) - 0.0050868 * numpy.sin(
        50 * numpy.pi * t / T - 0.011056) - 0.0085829 * numpy.sin(52 * numpy.pi * t / T - 0.86463)) * (scaling / 6.5)
    return flow


def FlowRateAorta(t, period, scaling):
    """
    Function from Boileau2015 (scaled inflow )
    Benchmark6_ADAN56
    scaling 1 found give 1ml/s
    """
    time = numpy.array(
        [0, 0.02193, 0.041305, 0.06068, 0.080055, 0.09943, 0.118805, 0.128492, 0.13818, 0.157555, 0.17693, 0.196305,
         0.244742, 0.322242, 0.351305, 0.390055, 0.415117, 0.437555, 0.48693, 0.515992, 0.583805, 1])
    flow = numpy.array(
        [0, 0, 105.089, 273.23, 404.592, 535.952, 567.479, 572.734, 567.479, 535.952, 488.663, 451.882, 315.266,
         178.651, 136.615, 0, -162.888, -21.0178, 0, 5.25444, 0, 0])
    meanflow = numpy.trapz(flow, time)
    currenttime = t / period - (t // period)
    flowrate = numpy.interp(currenttime, time, flow) * scaling / meanflow

    return flowrate


def FlowRateAlastruey2007(timevalues, scaling):
    returnvalues = [0 for i in timevalues]
    for index, time in enumerate(timevalues):
        if time <= 0.3:
            # returnvalues[index] = 485*numpy.sin(numpy.pi*time/0.3)
            returnvalues[index] = scaling * numpy.sin(numpy.pi * time / 0.3)
    return returnvalues

def FlowRateAlastruey2007function(t, period, scaling):
    time = numpy.linspace(0,1,1000)
    flow = [0 for i in time]
    for index, tt in enumerate(time):
        if tt <= 0.3:
            # returnvalues[index] = 485*numpy.sin(numpy.pi*time/0.3)
            flow[index] = numpy.sin(numpy.pi * tt / 0.3)

    meanflow = numpy.trapz(flow, time)
    currenttime = t / period - (t // period)
    flowrate = numpy.interp(currenttime, time, flow) * scaling / meanflow
    return flowrate




def FlowRateAorta2(t, period, scaling):
    """
    Function from Boileau2015 (scaled inflow )
    Benchmark4_AorticBifurcation
    scaling 1 found give 1ml/s
    """
    inflow = (7.9853e-06 + 2.6617e-05 * numpy.sin(
        2 * numpy.pi * t / period + 0.29498) + 2.3616e-05 * numpy.sin(
        4 * numpy.pi * t / period - 1.1403) - 1.9016e-05 * numpy.sin(
        6 * numpy.pi * t / period + 0.40435) - 8.5899e-06 * numpy.sin(
        8 * numpy.pi * t / period - 1.1892) - 2.436e-06 * numpy.sin(
        10 * numpy.pi * t / period - 1.4918) + 1.4905e-06 * numpy.sin(
        12 * numpy.pi * t / period + 1.0536) + 1.3581e-06 * numpy.sin(
        14 * numpy.pi * t / period - 0.47666) - 6.3031e-07 * numpy.sin(
        16 * numpy.pi * t / period + 0.93768) - 4.5335e-07 * numpy.sin(
        18 * numpy.pi * t / period - 0.79472) - 4.5184e-07 * numpy.sin(
        20 * numpy.pi * t / period - 1.4095) - 5.6583e-07 * numpy.sin(
        22 * numpy.pi * t / period - 1.3629) + 4.9522e-07 * numpy.sin(
        24 * numpy.pi * t / period + 0.52495) + 1.3049e-07 * numpy.sin(
        26 * numpy.pi * t / period - 0.97261) - 4.1072e-08 * numpy.sin(
        28 * numpy.pi * t / period - 0.15685) - 2.4182e-07 * numpy.sin(
        30 * numpy.pi * t / period - 1.4052) - 6.6217e-08 * numpy.sin(
        32 * numpy.pi * t / period - 1.3785) - 1.5511e-07 * numpy.sin(
        34 * numpy.pi * t / period - 1.2927) + 2.2149e-07 * numpy.sin(
        36 * numpy.pi * t / period + 0.68178) + 6.7621e-08 * numpy.sin(
        38 * numpy.pi * t / period - 0.98825) + 1.0973e-07 * numpy.sin(
        40 * numpy.pi * t / period + 1.4327) - 2.5559e-08 * numpy.sin(
        42 * numpy.pi * t / period - 1.2372) - 3.5079e-08 * numpy.sin(
        44 * numpy.pi * t / period + 0.2328)) * (scaling / 7.9853e-06)
    return inflow


def Bloodflow1D(self):
    print("Start Blood Flow Simulations.")
    # timestep = self.ModelParameters["TIMESTEP"]
    totaliter = numpy.ceil(self.ModelParameters["Beat_Duration"] / self.ModelParameters["TIMESTEP"])
    timestep = self.ModelParameters["Beat_Duration"] / totaliter
    pool = ThreadPool(8)
    OutputPeriod = int(self.ModelParameters["OUTPUT_PERIOD"])

    for node in self.Topology.Nodes:
        node.Radius *= 1e-3
        node.LengthAlongVessel *= 1e-3
        node.Thickness *= 1e-3
        node.RefRadius = copy.deepcopy(node.Radius)
        node.RefPressure = copy.deepcopy(self.PatientData["DiastolePressure"])
        node.SetPressureAreaEquation()

    for vessel in self.Topology.Vessels:
        for node in vessel.Nodes:
            node.PreviousPressure = self.PatientData["DiastolePressure"]
            node.PreviousVelocity = 0
            node.PreviousLumen = node.Radius * node.Radius * numpy.pi
            node.CurrentPressure = self.PatientData["DiastolePressure"]
            node.CurrentVelocity = 0
            node.CurrentLumen = node.Radius * node.Radius * numpy.pi
            node.Connections = [node.Number for node in node.Connections]
            node.Lumen0 = copy.deepcopy(node.Radius * node.Radius * numpy.pi)
            node.Density = self.ModelParameters["Density"]
            node.Iter = 0
            node.Friction = 22
            node.Viscosity = self.ModelParameters["BLOOD_VISC"]
            node.TimeStep = timestep
            node.DiastolicPressure = self.ModelParameters["DIASTOLIC_PRESSURE"]
            node.Beta = 4.0 / 3.0 * numpy.sqrt(numpy.pi) * node.YoungsModules * node.Thickness / node.Lumen0

        vessel.CurrentVelocity = [0 for node in vessel.Nodes]
        vessel.CurrentLumen = [node.Radius * node.Radius * numpy.pi for node in vessel.Nodes]
        vessel.CurrentPressure = [self.PatientData["DiastolePressure"] for node in vessel.Nodes]

    Bifurcations = []
    for bifnode in self.Topology.BifurcationNodes:
        bifnode.directions = []
        conlist = []
        for node in bifnode.Connections:
            conlist.append(node)
            direction = 1
            if node.LengthAlongVessel <= 1e-9:
                direction = -1
            bifnode.directions.append(direction)
        bf = Bifurcation(conlist, bifnode.directions)
        Bifurcations.append(bf)

    for node in self.Topology.OutletNodes:
        node.CurrentVolumeFlowRate = 0
        node.PreviousVolumeFlowRate = 0
        node.OutPressure = self.ModelParameters["OUT_PRESSURE"]
        node.chrt_frw_left = node.PreviousVelocity + 4 * numpy.power(node.PreviousLumen, 0.25) * numpy.sqrt(
            node.Beta / 2.0 / node.Density)

    def flowrate(t):
        #FlowRateAorta
        return FlowRateAorta(t, 60 / self.PatientData["HeartRate"],
                             self.PatientData["StrokeVolume"] * (self.PatientData["HeartRate"] / 60)) * 1e-6
        # return self.PatientData["StrokeVolume"] * (self.PatientData["HeartRate"] / 60) * 1e-6
    #
    # time = numpy.linspace(0,1,totaliter)
    # volume = numpy.trapz(flowrate(time),time)*1e6
    # print(volume)

    for element in self.Topology.InletNodes:
        node = element[0]
        node.FlowRateFunction = flowrate
        node.CurrentFlowrate = node.FlowRateFunction(0)
        node.charback = -4 * numpy.power(node.Lumen0, 0.25) * numpy.sqrt(node.Beta / 2.0 / node.Density)

    for vessel in self.Topology.Vessels:
        vessel.Nodes[0].c_0 = numpy.power(vessel.Nodes[0].Lumen0, 0.25) * numpy.sqrt(
            vessel.Nodes[0].Beta / 2.0 / vessel.Nodes[0].Density)
        vessel.Nodes[-1].c_0 = numpy.power(vessel.Nodes[-1].Lumen0, 0.25) * numpy.sqrt(
            vessel.Nodes[-1].Beta / 2.0 / vessel.Nodes[-1].Density)

    simvessels = self.Topology.Vessels
    simoutlets = [(node, self.Topology.Nodes[node.Connections[0]]) for node in self.Topology.OutletNodes]
    siminlets = [(node[0], self.Topology.Nodes[node[0].Connections[0]]) for node in self.Topology.InletNodes]

    def Update(vessel):
        for node in vessel.Nodes:
            node.PreviousPressure = node.CurrentPressure
            node.PreviousVelocity = node.CurrentVelocity
            node.PreviousLumen = node.CurrentLumen
            node.Iter += 1

    loops = 0
    maxloops = 20
    while loops < maxloops:
        loops += 1
        currenttime = 0.0
        currentiter = 0
        for node in self.Topology.Nodes:
            node.Iter = 0
        while 1:  # todo add convergence check
            # print("Time: %f" % currenttime)
            if (currentiter % (OutputPeriod) == 0 or currentiter == totaliter):# and loops == maxloops:
                tp = Results.TimePoint(currenttime)
                tp.velocity = [node.CurrentVelocity for vessel in self.Topology.Vessels for node in vessel.Nodes]
                tp.Pressure = [node.CurrentPressure for vessel in self.Topology.Vessels for node in vessel.Nodes]
                tp.Flow = [node.CurrentLumen * node.CurrentVelocity * 1e6 for vessel in self.Topology.Vessels for node
                           in vessel.Nodes]
                tp.Radius = [node.PressureRadiusEquation(node.CurrentPressure) * 1e3 for vessel in self.Topology.Vessels
                             for node in vessel.Nodes]
                print("Writing at Loop %d, Time: %f, Iter %d" % (loops, currenttime, currentiter))
                self.Results.TimePoints.append(tp)

            if currentiter == totaliter:
                break

            currentiter += 1
            currenttime = currentiter * timestep
            # print("Current Iter: %d" % currentiter)

            pool.map(SolveVessel, simvessels)
            pool.map(SolveOutletRCR, simoutlets)
            [SolveInlet(inlet) for inlet in siminlets]
            pool.map(SolvingBif, Bifurcations)

            # with mp.Pool(numberofcores) as pool:
            #     simvessels = list(pool.map(SolveVessel, simvessels))

            # with mp.Pool() as pool:
            #     pool.map_async(SolveVessel, simvessels)
            # for vessel in simvessels:
            #     t = Process(target=SolveVessel, args=(vessel,))
            #     t.start()
            # t.join()
            # Parallel(n_jobs=-1, verbose=0, backend="threading")(
            #     map(delayed(SolveVessel), simvessels))

            # for vessel in simvessels:
            #     SolveVessel(vessel)
            #
            # for vessel in simoutlets:
            #     SolveOutletRCR(vessel)
            #
            # for vessel in Bifurcations:
            #     SolvingBif(vessel)

            pool.map(Update, simvessels)
            for vessel in self.Topology.Vessels:
                vessel.CurrentPressure = [node.CurrentPressure for node in vessel.Nodes]
                vessel.CurrentLumen = [node.CurrentLumen for node in vessel.Nodes]
                vessel.CurrentVelocity = [node.CurrentVelocity for node in vessel.Nodes]

    pool.close()
    pool.join()
    file = self.Folders.ModellingFolder + "Results.dyn"
    OutputFile = open(file, 'w', encoding='utf-8')
    for res in self.Results.TimePoints:
        OutputFile.write("WT: %f\n" % res.WT)
        for index, data in enumerate(res.Pressure):
            OutputFile.write("%d\t%f\t%f\t%f\n" % (index, res.Flow[index], data, res.Radius[index]))

    print("Sim Done.")


def murraylaw(radius):
    return ((radius ** 3) / 2) ** (1 / 3)


def RadiusToLength(radius):
    return radius * 10


def BloodFlowSimR(patient):
    adjacencymatrix, nodes = patient.Topology.GetAdjacencyMatrix(patient.ModelParameters['BLOOD_VISC'])

    # build the matrix
    length = len(nodes)
    resistancematrix = numpy.zeros((length, length))
    for index, node in enumerate(nodes):
        data = adjacencymatrix[index, :]
        for indice in data.indices:
            resistancematrix[index][index] += 1e6 / (adjacencymatrix[index, indice])
            resistancematrix[index][indice] = -1e6 / (adjacencymatrix[index, indice])
            # update to call resistance of vessel itself

    startlengthoutlets = len(patient.Topology.InletNodes) + \
                         len(patient.Topology.BifurcationNodes) + \
                         len(patient.Topology.OutletNodes)
    # remove the wk nodes themselves, we already know the pressure there.
    resm = resistancematrix[:startlengthoutlets, :startlengthoutlets]

    inputvector = numpy.zeros((startlengthoutlets, 1))
    inputflowrate = patient.PatientData["StrokeVolume"] * (patient.PatientData["HeartRate"] / 60)
    inputvector[0, 0] = inputflowrate  # for all inlets

    for index, node in enumerate(patient.Topology.OutletNodes):
        inputvector[len(patient.Topology.InletNodes) + len(patient.Topology.BifurcationNodes) + index] = \
            patient.ModelParameters["OUT_PRESSURE"] / (
                    node.R1 + node.R2)

    solution = scipy.linalg.solve(resm, inputvector)

    # set pressure for each node
    for index in range(0, startlengthoutlets):
        nodes[index].Pressure = solution[index, 0]

    for bif in patient.Topology.BifurcationNodes:
        for node in bif.Connections:
            node.Pressure = bif.Pressure

    # determine flowrate for each vessel
    # flow in a vessel, Q=G*dP
    # we have R and P
    for vessel in patient.Topology.Vessels:
        vessel.FlowRate = 1e6 * (vessel.Nodes[0].Pressure - vessel.Nodes[-1].Pressure) / vessel.Resistance
        print(vessel.Name)
        print(vessel.FlowRate)
    # totalvolume = 0
    # for node in patient.Topology.OutletNodes:
    #     for vessel in patient.Topology.Vessels:
    #         if node in vessel.Nodes:
    #             totalvolume += vessel.FlowRate

    # plt.imshow(resistancematrix, cmap="gray")
    # plt.show()
    # print(1)


def BloodFlowSimDirect(patient):
    adjacencymatrix, nodes = patient.Topology.GetAdjacencyMatrix2(patient.ModelParameters['BLOOD_VISC'])
    outletpressure = [i.Pressure for i in patient.Topology.OutletNodes]
    # build the matrix
    length = len(nodes)
    resistancematrix = numpy.zeros((length, length))
    for index, node in enumerate(nodes):
        data = adjacencymatrix[index, :]
        for indice in data.indices:
            resistancematrix[index][index] += 1e6 / (adjacencymatrix[index, indice])
            resistancematrix[index][indice] = -1e6 / (adjacencymatrix[index, indice])
            # update to call resistance of vessel itself

    startlengthoutlets = len(patient.Topology.InletNodes) + len(patient.Topology.BifurcationNodes)

    # remove the wk nodes themselves, we already know the pressure there.
    resm = resistancematrix[:startlengthoutlets, :startlengthoutlets]

    inputvector = numpy.zeros((startlengthoutlets, 1))
    inputflowrate = patient.PatientData["StrokeVolume"] * (patient.PatientData["HeartRate"] / 60)
    inputvector[0, 0] = inputflowrate  # for all inlets

    for index, node in enumerate(patient.Topology.OutletNodes):
        for index2, node2 in enumerate(patient.Topology.BifurcationNodes):
            outlet = nodes.index(node)
            bif = nodes.index(node2)
            res = adjacencymatrix[outlet, bif]
            if res > 0:
                inputvector[bif, 0] += 1e6 * outletpressure[index] / res

    solution = scipy.linalg.solve(resm, inputvector)
    solvedsystem = numpy.allclose(numpy.dot(resm, solution), inputvector)

    # for index, node in enumerate(patient.Topology.OutletNodes):
    #     node.Pressure = outletpressure[index]

    # set pressure for each node
    for index in range(0, startlengthoutlets):
        nodes[index].Pressure = solution[index, 0]

    for bif in patient.Topology.BifurcationNodes:
        for node in bif.Connections:
            node.Pressure = bif.Pressure

    # determine flowrate for each vessel
    # flow in a vessel, Q=G*dP
    # we have R and P
    for vessel in patient.Topology.Vessels:
        vessel.FlowRate = 1e6 * (vessel.Nodes[0].Pressure - vessel.Nodes[-1].Pressure) / vessel.Resistance
        print(vessel.Name)
        print(vessel.FlowRate)

    # totalvolume = 0
    # for node in patient.Topology.OutletNodes:
    #     for vessel in patient.Topology.Vessels:
    #         if node in vessel.Nodes:
    #             totalvolume += vessel.FlowRate

    # plt.imshow(resistancematrix, cmap="gray")
    # plt.show()


class ResistanceNetwork:
    def __init__(self):
        self.AdjacencyMatrix = []
        self.ConductanceMatrix = []
        self.inputvector = []

    def SystemMatrix(self):
        diagonalelements = numpy.sum(self.ConductanceMatrix)
        systemmat = -1 * self.ConductanceMatrix
        for index, element in enumerate(diagonalelements):
            systemmat[index, index] = element
