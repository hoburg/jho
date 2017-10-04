" print performance numbers for jho "
import numpy as np
from jho import Mission

def jho_subs(model):
    """get solution for as-built Jungle Hawk Owl"""
    model.cost = 1/model["t_Mission/Loiter"]
    subs = {"b_Mission/Aircraft/Wing": 24,
            "l_Mission/Aircraft/Empennage/TailBoom": 7.0,
            "AR_v": 1.5, "c_{root}": 15./12, "SM_{corr}": 0.5, "AR_h": 4,
            "k": 0.0, "(1-k/2)": 1, "d_0": 1,
            "R_Mission/Aircraft/Fuselage": 7./12,
            "\\tau_Mission/Aircraft/Wing": 0.113661, "k_{nose}": 2.4055,
            "k_{bulk}": 4.3601, "k_{body}": 3.6518,
            "W_Mission/Aircraft/Empennage": 4.096,
            "W_Mission/Aircraft/Wing": 14.979,
            "W_Mission/Aircraft/Fuselage": 9.615}
    model.substitutions.update(subs)
    for p in model.varkeys["P_{avn}"]:
        model.substitutions.update({p: 65})
    for t in model.varkeys["\\theta_{max}"]:
        model.substitutions.update({t: 65})
    model.substitutions.update({"w_{lim}": 1})
    for vk in model.varkeys["w"]:
        model.substitutions.update({vk: 2})

    del model.substitutions["m_{fac}_Mission/Aircraft/Empennage"]
    del model.substitutions["m_{fac}_Mission/Aircraft/Wing"]
    del model.substitutions["m_{fac}_Mission/Aircraft/Fuselage"]
    model.cost = (model.cost/model["m_{fac}_Mission/Aircraft/Empennage"]
                  / model["m_{fac}_Mission/Aircraft/Wing"]
                  / model["m_{fac}_Mission/Aircraft/Fuselage"])
    sol = model.localsolve("mosek", verbosity=0)

    subs = {"m_{fac}_Mission/Aircraft/Wing":
            sol("m_{fac}_Mission/Aircraft/Wing"),
            "m_{fac}_Mission/Aircraft/Empennage":
            sol("m_{fac}_Mission/Aircraft/Empennage"),
            "m_{fac}_Mission/Aircraft/Fuselage":
            sol("m_{fac}_Mission/Aircraft/Fuselage")}
    model.substitutions.update(subs)

    del model.substitutions["W_Mission/Aircraft/Empennage"]
    del model.substitutions["W_Mission/Aircraft/Wing"]
    del model.substitutions["W_Mission/Aircraft/Fuselage"]

def perf_solve(model):
    "solve "
    model.cost = 1/model["t_Mission/Loiter"]
    sol = model.localsolve("mosek", verbosity=0)

    cmac = sol("c_{MAC}").magnitude
    print "mean aerodynamic chord [ft] = %.4f" % cmac

    Vy = sol("V_Mission/Climb/FlightSegment/FlightState")[0]
    print "speed for best rate of climb [m/s]: Vy = %.3f" % Vy.magnitude

    Vytop = sol("V_Mission/Climb/FlightSegment/FlightState")[-1]
    print "speed at top of climb [m/s] = %.3f" % Vytop.magnitude

    vloiter = np.average(
        sol("V_Mission/Loiter/FlightSegment/FlightState").magnitude)
    print "design loiter speed [m/s] = %.3f" % vloiter

    rho = sol("\\rho_{sl}").items()[0][1]
    S = sol("S_Mission/Aircraft/Wing")
    w55 = sol("W_{zfw}")*(sol("W_{zfw}").magnitude + 5)/sol("W_{zfw}").magnitude

    Vrot55 = ((2*w55/rho/S/1.39)**0.5).to("m/s")*1.5
    Vrot150 = ((2*sol("MTOW")/rho/S/1.39)**0.5).to("m/s")*1.5

    print "rotation speed at 55 lbs [m/s] = %.3f" % Vrot55.magnitude
    print "rotation speed at 150 lbs [m/s] = %.3f" % Vrot150.magnitude

    return sol

def max_speed(model):
    " find maximum speed at altitude "
    oldcost = model.cost
    model.cost = 1./np.prod(model["V_Mission/Loiter/FlightSegment/FlightState"])
    model.substitutions.update({"t_Mission/Loiter": 0.02})
    sol = model.localsolve("mosek")
    vmax = max(sol("V_Mission/Loiter/FlightSegment/FlightState"))
    print "Max Speed [m/s]: %.2f" % vmax.magnitude
    model.cost = oldcost
    return vmax

def optimum_speeds(model):
    " find optimum speeds "
    for v in model.varkeys["m_{fac}"]:
        mods = v.models
        if "Climb" in mods or "Loiter" in mods or "Cruise" in mods:
            if "FlightState" in mods:
                model.substitutions.update({v: 0.001})

    sol = model.localsolve("mosek", verbosity=0)

    vmins = sol("V_Mission/Loiter/FlightSegment/FlightState")[0].magnitude
    print ("optimum loiter speed for min power, "
           "start of loiter [m/s] = %.3f" % vmins)

    vmine = sol("V_Mission/Loiter/FlightSegment/FlightState")[-1].magnitude
    print ("optimum loiter speed for min power, "
           "end of loiter [m/s] = %.3f" % vmine)

    vstr = "V_Mission/Cruise/FlightSegment/FlightState"
    vcrin = sol(vstr).items()[0][1].magnitude
    print "optimum cruise speed, inbound [m/s] = %.3f" % vcrin

    vcrout = sol(vstr).items()[1][1].magnitude
    print "optimum cruise speed, outbound [m/s] = %.3f" % vcrout

    for v in model.varkeys["m_{fac}"]:
        mods = v.models
        if "Climb" in mods or "Loiter" in mods or "Cruise" in mods:
            if "FlightState" in mods:
                model.substitutions.update({v: 1})

def max_payload(model):
    " solve for maximum allowable payload "
    oldcost = model.cost
    model.cost = 1./model["W_{pay}"]
    oldsubw = model.substitutions["W_{pay}"]
    model.substitutions.update({"t_Mission/Loiter": 5.5})
    oldsubhdot = model.substitutions["\\dot{h}_{min}"]
    model.substitutions.update({"\\dot{h}_{min}": 10})
    del model.substitutions["W_{pay}"]
    sol = model.localsolve("mosek")
    wtot = sol("W_{pay}").magnitude
    wpay = (wtot + 14.0/3.0)/(7.0/5.0)
    mtow = sol("MTOW").magnitude
    print "Max payload weight [lbf] = %.3f" % wpay
    print "Max take off weight [lbf] = %.3f" % mtow
    model.substitutions.update({"W_{pay}": oldsubw})
    model.substitutions.update({"\\dot{h}_{min}": oldsubhdot})
    model.cost = oldcost

if __name__ == "__main__":
    M = Mission(DF70=True)
    jho_subs(M)
    Sol = perf_solve(M)
    _ = max_speed(M)
    max_payload(M)
    optimum_speeds(M)
