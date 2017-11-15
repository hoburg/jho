" print performance numbers for jho "
import numpy as np
import matplotlib.pyplot as plt
from jho import Mission
plt.rcParams.update({'font.size':15})

def jho_subs(model):
    """get solution for as-built Jungle Hawk Owl"""
    model.cost = 1/model["t_Mission/Loiter"]
    subs = {model.JHO.wing.planform.b: 24,
            model.JHO.emp.tailboom.l: 7.0,
            model.JHO.emp.vtail.lv: 7.0,
            model.JHO.emp.vtail.planform.AR: 15,
            model.JHO.wing.planform.croot: 15./12,
            "SM_{corr}": 0.5,
            model.JHO.emp.htail.planform.AR: 4,
            model.JHO.emp.tailboom.k: 0.0,
            model.JHO.emp.tailboom.d0: 1,
            "R_Mission/Aircraft/Fuselage": 7./12,
            model.JHO.wing.planform.tau: 0.113661,
            "k_{nose}": 2.4055,
            "k_{bulk}": 4.3601, "k_{body}": 3.6518,
            model.JHO.emp.W: 4.096,
            model.JHO.wing.W: 14.979,
            "W_Mission/Aircraft/Fuselage": 9.615}
    model.substitutions.update(subs)
    for p in model.varkeys["P_{avn}"]:
        model.substitutions.update({p: 65})
    for t in model.varkeys["\\theta_{max}"]:
        model.substitutions.update({t: 65})
    model.substitutions.update({model.JHO.wing.spar.wlim: 1})
    for vk in model.varkeys["w"]:
        model.substitutions.update({vk: 2})

    del model.substitutions[model.JHO.emp.mfac]
    del model.substitutions[model.JHO.wing.mfac]
    del model.substitutions["m_{fac}_Mission/Aircraft/Fuselage"]
    model.cost = (model.cost/model[model.JHO.emp.mfac]
                  / model[model.JHO.wing.mfac]
                  / model["m_{fac}_Mission/Aircraft/Fuselage"])
    sol = model.localsolve("mosek", verbosity=0)

    subs = {model.JHO.wing.mfac: sol(model.JHO.wing.mfac),
            model.JHO.emp.mfac: sol(model.JHO.emp.mfac),
            "m_{fac}_Mission/Aircraft/Fuselage":
            sol("m_{fac}_Mission/Aircraft/Fuselage")}
    model.substitutions.update(subs)

    del model.substitutions[model.JHO.emp.W]
    del model.substitutions[model.JHO.wing.W]
    del model.substitutions["W_Mission/Aircraft/Fuselage"]

def perf_solve(model):
    "solve "
    del model.substitutions["t_Mission/Loiter"]
    model.cost = 1/model["t_Mission/Loiter"]
    sol = model.localsolve("mosek", verbosity=0)

    mtow = sol("MTOW").magnitude
    print "MTOW [lbs] = %.2f" % mtow

    wzfw = sol("W_{zfw}").magnitude
    print "Zero fuel weight [lbs] = %.2f" % wzfw

    b = sol(model.JHO.wing.planform.b).magnitude
    print "Wing span [ft] = %.2f" % b

    lfuse = sol("l_Mission/Aircraft/Fuselage")
    ltail = sol(model.JHO.emp.tailboom.l)
    ljho = (lfuse + ltail).to("ft").magnitude
    print "Aicraft length [ft] = %.2f" % ljho

    AR = sol(model.JHO.wing.planform.AR)
    print "Aspect ratio = %.2f" % AR

    cmac = sol(model.JHO.wing.planform.cmac).magnitude
    print "mean aerodynamic chord [ft] = %.4f" % cmac

    croot = sol(model.JHO.wing.planform.croot).magnitude
    print "root chord [ft] = %.3f" % croot

    Vy = sol("V_Mission/Climb/FlightSegment/FlightState")[0]
    print "speed for best rate of climb [m/s]: Vy = %.3f" % Vy.magnitude

    Vytop = sol("V_Mission/Climb/FlightSegment/FlightState")[-1]
    print "speed at top of climb [m/s] = %.3f" % Vytop.magnitude

    vloiter = np.average(
        sol("V_Mission/Loiter/FlightSegment/FlightState").magnitude)
    print "design loiter speed [m/s] = %.3f" % vloiter

    rho = sol("rhosl").items()[0][1]
    S = sol(model.JHO.wing.planform.S)
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
    vmax = max(sol("V_Mission/Loiter/FlightSegment/FlightState")).magnitude
    rho = sol("\\rho_Mission/Loiter/FlightSegment/FlightState")[0].magnitude
    rhosl = 1.225
    print "Max Speed [m/s]: %.2f" % (vmax*rhosl/rho)
    model.cost = oldcost
    return vmax

def optimum_speeds(model):
    " find optimum speeds "
    for v in model.varkeys["m_{fac}"]:
        mods = v.models
        if "Climb" in mods or "Loiter" in mods or "Cruise" in mods:
            if "FlightState" in mods:
                model.substitutions.update({v: 0.001})

    model.cost = 1/model["t_Mission/Loiter"]
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

def plot_climbrate(result):
    fig, ax = plt.subplots()
    ax.plot(result("h_Mission/Climb/FlightSegment/FlightState").magnitude - 1500, result("\\dot{h}"), "k")
    ax.set_xlabel("Altitude [ft]")
    ax.set_ylabel("Climb Rate [ft/min]")
    ax.set_xlim([0, 10000])
    ax.grid()
    return fig, ax

def plot_glide(result):
    LoD = np.mean(np.hstack([result(l)/result(d) for l, d in
                             zip(result("CL"), result("C_D"))]))
    h = (result("h_Mission/Climb/FlightSegment/FlightState")
         - result("h_{ref}_Mission/Climb/FlightSegment/FlightState")/10)
    R = (LoD*h).to("nmi")
    fig, ax = plt.subplots()
    ax.plot(h, R, "k")
    ax.set_xlabel("Altidue [ft]")
    ax.set_ylabel("Glide Range [nmi]")
    ax.grid()
    ax.set_xlim([0, 10000])
    return fig, ax

def test():
    M = Mission(DF70=True)
    jho_subs(M)
    M.substitutions["t_Mission/Loiter"] = 5
    Sol = perf_solve(M)
    optimum_speeds(M)
    _ = max_speed(M)
    max_payload(M)

if __name__ == "__main__":
    M = Mission(DF70=True)
    jho_subs(M)
    M.substitutions["t_Mission/Loiter"] = 5
    Sol = perf_solve(M)
    optimum_speeds(M)
    _ = max_speed(M)
    max_payload(M)
    f, a = plot_climbrate(Sol)
    f.savefig("crateh.jpg")
    f, a = plot_glide(Sol)
    f.savefig("gliderange.jpg")
