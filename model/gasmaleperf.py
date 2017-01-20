import numpy as np
from numpy import pi
import gpkitmodels.aircraft.GP_submodels
from gpkitmodels.aircraft.GP_submodels.breguet_endurance import BreguetEndurance
from gpkitmodels.aircraft.GP_submodels.gas_engine import Engine
from gpkitmodels.aircraft.GP_submodels.wing import Wing
from gpkitmodels.aircraft.GP_submodels.cylindrical_fuselage import Fuselage
from gpkitmodels.aircraft.GP_submodels.empennage import Empennage
from gpkitmodels.aircraft.GP_submodels.tail_boom import TailBoomState
from gpkitmodels.aircraft.GP_submodels.tail_boom_flex import TailBoomFlexibility
from helpers import summing_vars
from gpkit import Model, Variable, Vectorize, units

# pylint: disable=invalid-name

class Aircraft(Model):
    "the JHO vehicle"
    def setup(self, Wfueltot, DF70=False):

        self.fuselage = Fuselage(Wfueltot)
        self.wing = Wing()
        self.engine = Engine(DF70)
        self.empennage = Empennage()
        self.pylon = Pylon()

        components = [self.fuselage, self.wing, self.engine, self.empennage,
                      self.pylon]
        self.smeared_loads = [self.fuselage, self.engine, self.pylon]

        Wzfw = Variable("W_{zfw}", "lbf", "zero fuel weight")
        Wpay = Variable("W_{pay}", 10, "lbf", "payload weight")
        Ppay = Variable("P_{pay}", 10, "W", "payload power")
        Wavn = Variable("W_{avn}", 8, "lbf", "avionics weight")
        lantenna = Variable("l_{antenna}", 13.4, "in", "antenna length")
        wantenna = Variable("w_{antenna}", 10.2, "in", "antenna width")
        # propr = Variable("r", 11, "in", "propellor radius")
        Volpay = Variable("\\mathcal{V}_{pay}", 1.0, "ft**3", "payload volume")
        Volavn = Variable("\\mathcal{V}_{avn}", 0.125, "ft**3",
                          "avionics volume")

        constraints = [
            Wzfw >= sum(summing_vars(components, "W")) + Wpay + Wavn,
            self.empennage.horizontaltail["V_h"] <= (
                self.empennage.horizontaltail["S"]
                * self.empennage.horizontaltail["l_h"]/self.wing["S"]**2
                * self.wing["b"]),
            self.empennage.verticaltail["V_v"] <= (
                self.empennage.verticaltail["S"]
                * self.empennage.verticaltail["l_v"]/self.wing["S"]
                / self.wing["b"]),
            self.wing["C_{L_{max}}"]/self.wing["m_w"] <= (
                self.empennage.horizontaltail["C_{L_{max}}"]
                / self.empennage.horizontaltail["m_h"]),
            # enforce antenna sticking on the tail
            self.empennage.verticaltail["c_{t_v}"] >= wantenna,
            self.empennage.verticaltail["b"] >= lantenna,
            # enforce a cruciform with the htail infront of vertical tail
            self.empennage.tailboom["l"] >= (
                self.empennage.horizontaltail["l_h"]
                + self.empennage.horizontaltail["c_{r_h}"]),
            4./6*pi*self.fuselage["k_{nose}"]*self.fuselage["R"]**3 >= Volpay,
            self.fuselage["\\mathcal{V}_{body}"] >= (
                self.fuselage.fueltank["\\mathcal{V}"] + Volavn),
            ]

        if DF70:
            constraints.extend([self.engine["h"] <= 2*self.fuselage["R"]])

        return components, constraints

    def flight_model(self, state):
        return AircraftPerf(self, state)

    def loading(self, Wcent):
        return AircraftLoading(self, Wcent)

class Pylon(Model):
    "attachment from fuselage to pylon"
    def setup(self):

        h = Variable("h", 7, "in", "pylon height")
        l = Variable("l", 32.8, "in", "pylon length")
        S = Variable("S", "ft**2", "pylon surface area")
        W = Variable("W", 1.42, "lbf", "pylon weight")

        constraints = [S >= 2*l*h]

        return constraints

    def flight_model(self, state):
        return PylonAero(self, state)

class PylonAero(Model):
    "pylon drag model"
    def setup(self, static, state):

        Cf = Variable("C_f", "-", "fuselage skin friction coefficient")
        Re = Variable("Re", "-", "fuselage reynolds number")

        constraints = [
            Re == state["V"]*state["\\rho"]*static["l"]/state["\\mu"],
            Cf >= 0.455/Re**0.3,
            ]

        return constraints

class AircraftLoading(Model):
    "aircraft loading model"
    def setup(self, aircraft, Wcent):

        loading = [aircraft.wing.loading(Wcent)]
        loading.append(aircraft.empennage.loading())
        loading.append(aircraft.fuselage.loading(Wcent))

        tbstate = TailBoomState()
        loading.append(TailBoomFlexibility(aircraft.empennage.horizontaltail,
                                           aircraft.empennage.tailboom,
                                           aircraft.wing, tbstate))

        return loading

class AircraftPerf(Model):
    "performance model for aircraft"
    def setup(self, static, state, **kwargs):

        self.wing = static.wing.flight_model(state)
        self.fuselage = static.fuselage.flight_model(state)
        self.engine = static.engine.flight_model(state)
        self.htail = static.empennage.horizontaltail.flight_model(state)
        self.vtail = static.empennage.verticaltail.flight_model(state)
        self.tailboom = static.empennage.tailboom.flight_model(state)
        self.pylon = static.pylon.flight_model(state)

        self.dynamicmodels = [self.wing, self.fuselage, self.engine,
                              self.htail, self.vtail, self.tailboom, self.pylon]
        areadragmodel = [self.fuselage, self.htail, self.vtail, self.tailboom,
                         self.pylon]
        areadragcomps = [static.fuselage, static.empennage.horizontaltail,
                         static.empennage.verticaltail,
                         static.empennage.tailboom, static.pylon]

        Wend = Variable("W_{end}", "lbf", "vector-end weight")
        Wstart = Variable("W_{start}", "lbf", "vector-begin weight")
        CD = Variable("C_D", "-", "drag coefficient")
        CDA = Variable("CDA", "-", "area drag coefficient")
        mfac = Variable("m_{fac}", 2.1, "-", "drag margin factor")

        dvars = []
        for dc, dm in zip(areadragcomps, areadragmodel):
            if "C_d" in dm.varkeys:
                dvars.append(dm["C_d"]*dc["S"]/static.wing["S"])
            elif "C_f" in dm.varkeys:
                dvars.append(dm["C_f"]*dc["S"]/static.wing["S"])

        constraints = [CDA/mfac >= sum(dvars),
                       CD >= CDA + self.wing["C_d"]]

        return self.dynamicmodels, constraints

class FlightState(Model):
    "define environment state during a portion of an aircraft mission"
    def setup(self, alt, wind, **kwargs):

        rho = Variable("\\rho", "kg/m^3", "air density")
        h = Variable("h", alt, "ft", "altitude")
        href = Variable("h_{ref}", 15000, "ft", "Reference altitude")
        psl = Variable("p_{sl}", 101325, "Pa", "Pressure at sea level")
        Latm = Variable("L_{atm}", 0.0065, "K/m", "Temperature lapse rate")
        Tsl = Variable("T_{sl}", 288.15, "K", "Temperature at sea level")
        temp = [(t.value - l.value*v.value).magnitude
                for t, v, l in zip(Tsl, h, Latm)]
        Tatm = Variable("t_{atm}", temp, "K", "Air temperature")
        mu = Variable("\\mu", "N*s/m^2", "Dynamic viscosity")
        musl = Variable("\\mu_{sl}", 1.789*10**-5, "N*s/m^2",
                        "Dynamic viscosity at sea level")
        Rspec = Variable("R_{spec}", 287.058, "J/kg/K",
                         "Specific gas constant of air")

        # Atmospheric variation with altitude (valid from 0-7km of altitude)
        constraints = [rho == psl*Tatm**(5.257-1)/Rspec/(Tsl**5.257),
                       (mu/musl)**0.1 == 0.991*(h/href)**(-0.00529),
                       Latm == Latm]

        V = Variable("V", "m/s", "true airspeed")

        if wind:

            V_wind = Variable("V_{wind}", 25, "m/s", "Wind speed")
            constraints.extend([V >= V_wind])

        else:

            V_wind = Variable("V_{wind}", "m/s", "Wind speed")
            V_ref = Variable("V_{ref}", 25, "m/s", "Reference wind speed")

            constraints.extend([(V_wind/V_ref) >= 0.6462*(h/href) + 0.3538,
                                V >= V_wind])
        return constraints

class FlightSegment(Model):
    "creates flight segment for aircraft"
    def setup(self, N, aircraft, alt=15000, wind=False, etap=0.7):

        self.aircraft = aircraft

        with Vectorize(N):
            self.fs = FlightState(alt, wind)
            self.aircraftPerf = self.aircraft.flight_model(self.fs)
            self.slf = SteadyLevelFlight(self.fs, self.aircraft,
                                         self.aircraftPerf, etap)
            self.be = BreguetEndurance(self.aircraftPerf)

        self.submodels = [self.fs, self.aircraftPerf, self.slf, self.be]

        Wfuelfs = Variable("W_{fuel-fs}", "lbf", "flight segment fuel weight")

        self.constraints = [Wfuelfs >= self.be["W_{fuel}"].sum()]

        if N > 1:
            self.constraints.extend([self.aircraftPerf["W_{end}"][:-1] >=
                                     self.aircraftPerf["W_{start}"][1:]])

        return self.aircraft, self.submodels, self.constraints

class Loiter(Model):
    "make a loiter flight segment"
    def setup(self, N, aircraft, alt=15000, wind=False, etap=0.7):
        fs = FlightSegment(N, aircraft, alt, wind, etap)

        t = Variable("t", "days", "time loitering")
        constraints = [fs.be["t"] >= t/N]

        return constraints, fs

class Cruise(Model):
    "make a cruise flight segment"
    def setup(self, N, aircraft, alt=15000, wind=False, etap=0.7, R=200):
        fs = FlightSegment(N, aircraft, alt, wind, etap)

        R = Variable("R", R, "nautical_miles", "Range to station")
        constraints = [R/N <= fs["V"]*fs.be["t"]]

        return fs, constraints

class Climb(Model):
    "make a climb flight segment"
    def setup(self, N, aircraft, alt=15000, wind=False, etap=0.7, dh=15000):
        fs = FlightSegment(N, aircraft, alt, wind, etap)

        with Vectorize(N):
            hdot = Variable("\\dot{h}", "ft/min", "Climb rate")

        deltah = Variable("\\Delta h", dh, "ft", "altitude difference")
        hdotmin = Variable("\\dot{h}_{min}", 100, "ft/min",
                           "minimum climb rate")

        constraints = [
            hdot*fs.be["t"] >= deltah/N,
            hdot >= hdotmin,
            fs.slf["T"] >= (0.5*fs["\\rho"]*fs["V"]**2*fs["C_D"]
                            * fs.aircraft.wing["S"] + fs["W_{start}"]*hdot
                            / fs["V"]),
            ]

        return fs, constraints

class SLFMaxSpeed(Model):
    "steady level flight model"
    def setup(self, state, aircraft, perf, etap):

        T = Variable("T", "N", "thrust")
        etaprop = Variable("\\eta_{prop}", etap, "-", "propulsive efficiency")

        constraints = [
            (perf["W_{end}"]*perf["W_{start}"])**0.5 <= (
                0.5*state["\\rho"]*state["V_{max}"]**2*perf["C_L"]
                * aircraft.wing["S"]),
            T >= (0.5*state["\\rho"]*state["V_{max}"]**2*perf["C_D"]
                  *aircraft.wing["S"]),
            perf["P_{shaft-max}"] >= T*state["V_{max}"]/etaprop]

        return constraints

class SteadyLevelFlight(Model):
    "steady level flight model"
    def setup(self, state, aircraft, perf, etap):

        T = Variable("T", "N", "thrust")
        etaprop = Variable("\\eta_{prop}", etap, "-", "propulsive efficiency")

        constraints = [
            (perf["W_{end}"]*perf["W_{start}"])**0.5 <= (
                0.5*state["\\rho"]*state["V"]**2*perf["C_L"]
                * aircraft.wing["S"]),
            T >= (0.5*state["\\rho"]*state["V"]**2*perf["C_D"]
                  *aircraft.wing["S"]),
            perf["P_{shaft}"] >= T*state["V"]/etaprop]

        return constraints

class Mission(Model):
    "creates flight profile"
    def setup(self, DF70=False, wind=False):

        mtow = Variable("MTOW", "lbf", "max-take off weight")
        Wcent = Variable("W_{cent}", "lbf", "center aircraft weight")
        Wfueltot = Variable("W_{fuel-tot}", "lbf", "total aircraft fuel weight")
        LS = Variable("(W/S)", "lbf/ft**2", "wing loading")

        JHO = Aircraft(Wfueltot, DF70)
        loading = JHO.loading(Wcent)

        climb1 = Climb(10, JHO, alt=np.linspace(0, 15000, 11)[1:], etap=0.508, wind=wind)
        cruise1 = Cruise(1, JHO, etap=0.684, R=180, wind=wind)
        loiter1 = Loiter(5, JHO, etap=0.647, wind=wind)
        cruise2 = Cruise(1, JHO, etap=0.684, wind=wind)
        mission = [climb1, cruise1, loiter1, cruise2]

        constraints = [
            mtow >= climb1["W_{start}"][0],
            Wfueltot >= sum(fs["W_{fuel-fs}"] for fs in mission),
            mission[-1]["W_{end}"][-1] >= JHO["W_{zfw}"],
            Wcent >= Wfueltot + sum(summing_vars(JHO.smeared_loads, "W")),
            LS == mtow/JHO.wing["S"],
            loiter1["P_{total}"] >= (loiter1["P_{shaft}"] + (
                loiter1["P_{avn}"] + JHO["P_{pay}"])
                                     / loiter1["\\eta_{alternator}"])
            ]

        for i, fs in enumerate(mission[1:]):
            constraints.extend([
                mission[i]["W_{end}"][-1] == fs["W_{start}"][0]
                ])

        return JHO, mission, loading, constraints

    def process_solution(self, sol):
        print sol("MTOW")

def test():
    M = Mission(DF70=True)
    M.cost = 1/M["t_Mission, Loiter"]
    subs = {"b_Mission, Aircraft, Wing": 24,
            "l_Mission, Aircraft, Empennage, TailBoom": 7.0,
            "AR_v": 1.5, "AR": 24, "SM_{corr}": 0.5, "AR_h": 4, "k": 0.0,
            "(1-k/2)": 1, "d_0": 1}
    M.substitutions.update(subs)
    for p in M.varkeys["P_{avn}"]:
        M.substitutions.update({p: 65})
    for t in M.varkeys["\\theta_{max}"]:
        M.substitutions.update({t: 65})
    M.localsolve()

if __name__ == "__main__":
    M = Mission(DF70=True)
    M.cost = 1/M["t_Mission, Loiter"]
    subs = {"b_Mission, Aircraft, Wing": 24,
            "l_Mission, Aircraft, Empennage, TailBoom": 7.0,
            "AR_v": 1.5, "AR": 24, "SM_{corr}": 0.5, "AR_h": 4, "k": 0.0,
            "(1-k/2)": 1, "d_0": 1}
    M.substitutions.update(subs)
    for p in M.varkeys["P_{avn}"]:
        M.substitutions.update({p: 65})
    for t in M.varkeys["\\theta_{max}"]:
        M.substitutions.update({t: 65})
    sol = M.localsolve("mosek")

    LD = sol("C_L_Mission, Loiter, FlightSegment, AircraftPerf, WingAero")/sol("C_D_Mission, Loiter, FlightSegment, AircraftPerf")
