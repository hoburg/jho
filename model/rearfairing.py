from gasmaleperf import Mission
import matplotlib.pyplot as plt
import numpy as np

M = Mission(DF70=True)
M.cost = 1/M["t_Mission, Loiter"]
subs = {"b_Mission, Aircraft, Wing": 24,
        "l_Mission, Aircraft, Empennage, TailBoom": 7.0,
        "AR_v": 1.5, "AR": 24, "SM_{corr}": 0.5, "AR_h": 4, "k": 0.0,
        "(1-k/2)": 1, "d_0": 1, "R_Mission, Aircraft, Fuselage": 7./12}
M.substitutions.update(subs)
for p in M.varkeys["P_{avn}"]:
    M.substitutions.update({p: 65})
for t in M.varkeys["\\theta_{max}"]:
    M.substitutions.update({t: 65})
sol = M.localsolve("mosek")

kbulk = sol("k_{bulk}")

M.substitutions.update({"k_{bulk}": ("sweep", np.linspace(2, kbulk, 15))})
sol = M.localsolve("mosek")

fig, ax = plt.subplots()
ax.plot((sol("k_{bulk}")*sol("R_Mission, Aircraft, Fuselage")).to("inches"), sol("t_Mission, Loiter"))
ax.set_ylabel("endurance")
ax.set_xlabel("rear fairing length")
fig.savefig("rearfairingtrade.jpg")
