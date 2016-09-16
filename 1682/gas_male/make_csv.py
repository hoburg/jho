import numpy as np
import pandas as pd
from gasmale import GasMALE
from gpkit.small_scripts import unitstr

def output_csv(path, M, sol, varnames, margins):
    """
    This ouputs variables relevant accross a mission
    """
    sens = sol["sensitivities"]["constants"]

    fseg = {}
    for subm in M.submodels:
        if subm.__class__.__name__ == "Mission":
            for fs in subm.submodels:
                fseg[fs.name] = {"index": [], "shape": [], "start": []}

    start = [1]
    for subm in M.submodels:
        if subm.__class__.__name__ == "Mission":
            for i, fs in enumerate(subm.submodels):
                fseg[fs.name]["index"].append(fs.num)
                fseg[fs.name]["shape"].append(fs.N)
                start.append(start[i] + fs.N)
                fseg[fs.name]["start"].append(start[i])

    colnames = ["Units"]
    for fs in fseg:
        for i in fseg[fs]["shape"]:
            colnames += [fs]*i
    colnames.append("Label")

    data = {}
    for vname in varnames:
        data[vname] = [0]*(start[-1] + 1)
        if vname in sens:
            data[vname + " sensitivity"] = [""] + [0]*(start[-1]-1) + [""]

    i = 0
    for vname in varnames:
        for sv in sol(vname):
            for fs in fseg:
                if fs not in sv.models:
                    continue
                ind = sv.models.index(fs)
                ifs = fseg[fs]["index"].index(sv.modelnums[ind])
                data[vname][fseg[fs]["start"][ifs]:fseg[fs]["start"][ifs]
                            + sv.shape[0]] = sol(sv).magnitude[0:]
                if vname in sens:
                    data[vname + " sensitivity"][fseg[fs]["start"][ifs]:
                            fseg[fs]["start"][ifs] + sv.shape[0]] = sens[sv][0:]


    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = colnames
    df.to_csv("%soutput1.csv" % path)

def bd_csv_output(path, sol, varname):

    if varname in sol["sensitivities"]["constants"]:
        colnames = ["Value", "Units", "Sensitivitiy", "Label"]
    else:
        colnames = ["Value", "Units", "Label"]

    data = {}
    for sv in sol(varname):
        name = sv
        data[name] = [sol(sv).magnitude]
        data[name].append(unitstr(sv.units))
        if varname in sol["sensitivities"]["constants"]:
            data[name].append(sol["sensitivities"]["constants"][sv])
        data[name].append(sv.label)

    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = colnames
    df.to_csv("%s%s_breakdown.csv" %
              (path, varname.replace("{", "").replace("}", "")))

if __name__ == "__main__":
    M = GasMALE()
    M.substitutions.update({"t_{loiter}": 6})
    M.cost = M["MTOW"]
    Sol = M.solve("mosek")
    PATH = "/Users/mjburton11/Dropbox (MIT)/16.82GasMALE/GpkitReports/csvs/"

    Mission_vars = ["RPM", "BSFC", "V", "P_{shaft}", "P_{shaft-tot}",
                    "h_{dot}", "h", "T_{atm}", "\\mu", "\\rho", "W_{fuel}",
                    "W_{N}", "W_{N+1}", "C_D", "C_L", "\\eta_{prop}", "T",
                    "h_{loss}", "P_{shaft-max}", "t", "Re", "C_{f-fuse}",
                    "C_{D-fuse}", "c_{dp}", "V_{wind}"]
    Margins = ["BSFC", "c_{dp}"]
    output_csv(PATH, M, Sol, Mission_vars, Margins)
    bd_csv_output(PATH, Sol, "W")
    bd_csv_output(PATH, Sol, "m_{fac}")