from jho import Mission
from print_perf import perf_solve, max_speed, optimum_speeds
import numpy as np
import matplotlib.pyplot as plt

def fix_vars(model, result):

    fixvars = ["b_Mission/Aircraft/Wing",
               "l_Mission/Aircraft/Empennage/TailBoom", "AR_v", "c_{root}",
               "SM_{corr}", "AR_h", "k", "(1-k/2)", "d_0",
               "R_Mission/Aircraft/Fuselage", "\\tau_Mission/Aircraft/Wing",
               "k_{nose}", "k_{bulk}", "k_{body}",
               "W_Mission/Aircraft/Empennage", "W_Mission/Aircraft/Wing",
               "W_Mission/Aircraft/Fuselage", "W_{eng}"]

    for v in fixvars:
        val = result(v)
        model.substitutions.update({v: val})

def highpay_subs(model):
    model.substitutions["t_Mission/Loiter"] = 5
    model.cost = model["MTOW"]

    for mf in model.varkeys["m_{fac}"]:
        if len(mf.models) == 4:
            if "AircraftPerf" in mf.models:
                model.substitutions.update({mf: 1.05})

    model.substitutions.update({"W_{pay}": 100, "P_{pay}": 650})

if __name__ == "__main__":
    M = Mission(DF70=False)
    highpay_subs(M)
    sol = M.localsolve("mosek")
    mtow5 = sol("MTOW").magnitude

    fix_vars(M, sol)
    perf_solve(M)

    weng = sol("W_{eng}")
    print "Engine Weight [lbf] = %.3f" % weng.magnitude
    psl = sol("P_{sl-max}").magnitude
    print "Engine Power [hp] = %.3f" % psl

    optimum_speeds(M)
    _ = max_speed(M)

    M = Mission(DF70=False)
    highpay_subs(M)
    fig, ax = plt.subplots()
    style = [":", "-.", "--", "-", ":"]
    for t, st in zip([3, 4, 5, 6, 7], style):
        M.substitutions.update({"t_Mission/Loiter": t})
        M.substitutions.update({"W_{pay}": ("sweep", np.linspace(50, 300, 5))})
        sol = M.localsolve("mosek", skipsweepfailures=True)
        ax.plot(sol("W_{pay}"), sol("MTOW"), linestyle=st, color="k", linewidth=1.2, label="Endurance: %d [days]" % t)

    ax.set_xlabel("Payload Weight [lbf]")
    ax.set_ylabel("Max Takeoff Weight [lbf]")
    ax.legend(loc=2)
    # ax.plot([5, 5], [0, mtow5], "--r")
    # ax.plot([2, 5], [mtow5, mtow5], "--r")
    ax.set_ylim([0, 2000])
    ax.grid()
    fig.savefig("endurancetrade.pdf")
