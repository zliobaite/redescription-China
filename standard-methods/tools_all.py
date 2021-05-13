import re
import numpy

import matplotlib.pyplot as plt
from matplotlib.patches import BoxStyle
# params = {'text.latex.preamble' : [r'\usepackage{ams}']}
# plt.rcParams.update(params)

### MAP TOOLS
##########################
import tools_geomap

# projs_loc = ["merc", "stere", "mill"]
# projs_glob = ["moll", "sinu", "robin", "hammer"]
map_prefs = {"map_proj": "mill", "map_resolution": "crude", "sea_color": "#F9FCFF",
             "map_elems": ["seas", "coasts"]} #, "meridians", "parallels"]}
#map_corners = [[64, 126], [0, 45]]
map_corners = [[45, 126], [0, 45]]


### COLOR TOOLS
##########################
import matplotlib.colors

ALPHA_DOTS = 0.55
primary_colors = {
    "K": (0., 0., 0.), # black
    "W": (1., 1., 1.), # white
    "G": (.8, .8, .8), # gray
    "H": (.5, .5, .5), # half-gray
    "U": (68/255.,119/255.,170/255.), # blue
    "R": (204/255.,102/255.,119/255.), # red
    "Y": (221/255.,204/255.,119/255.), # yellow
    "V": (17/255.,119/255.,51/255.) # green
    }

TOL_COLORS = [(1 , (51/255.,34/255.,136/255.)),
    (11 , (102/255.,153/255.,204/255.)),
    (5 , (136/255.,204/255.,238/255.)),
    (7 , (68/255.,170/255.,153/255.)),
    (4 , primary_colors["V"]),
    (8 , (153/255.,153/255.,51/255.)),
    (3 , primary_colors["Y"]),
    (10 , (102/255.,17/255.,0/255.)),
    (2 , primary_colors["R"]),
    (12 , (170/255.,68/255.,102/255.)),
    (9 , (136/255.,34/255.,85/255.)),
    (6 , (170/255.,68/255.,153/255.))]
TOL_NAMES = ["DarkPurple", "DarkBlue", "LightBlue",
             "LightGreen", "DarkGreen",
             "DarkBrown", "LightBrown",
             "DarkRed", "LightRed",
             "LightPink", "DarkPink", "LightPurple"]

TOL_NAMED = dict([(n, TOL_COLORS[i][1]) for (i,n) in enumerate(TOL_NAMES)])

#### Sequential
SEQ_HTML = ["FFFFE5", "FFF7BC", "FEE391", "FEC44F", "FB9A29", "EC7014", "CC4C02", "993404", "662506"]
SEQ_COLORS = [(i, tuple([int(c[2*i:2*(i+1)], 16)/255. for i in range(3)])) for i,c in enumerate(SEQ_HTML)]

#### Divergent
DIV_HTML = ["3D52A1", "3A89C9", "77B7E5", "B4DDF7", "E6F5FE", "FFFAD2", "FFE3AA", "F9BD7E", "ED875E", "D24D3E", "AE1C3E"]
DIV_COLORS = [(i, tuple([int(c[2*i:2*(i+1)], 16)/255. for i in range(3)])) for i,c in enumerate(DIV_HTML)]


def make_colormap(cs, extrapolate=True, name=None):
    if not extrapolate:
        if name is None:
            name = 'C#%dd' % len(cs)
        cs = [(c[0], c[1], c[2]) for c in cs]
        return matplotlib.colors.ListedColormap(cs, name)
    cdict = {'red': [], "green": [], "blue": []}
    if name is None:
        name = 'C#%dc' % len(cs)
    if len(cs) == 1:
        cs.insert(0, primary_colors["G"])
    dvf = len(cs)-1.
    for ci, c in enumerate(cs):
        cdict["red"].append((ci/dvf, c[0], c[0]))
        cdict["green"].append((ci/dvf, c[1], c[1]))
        cdict["blue"].append((ci/dvf, c[2], c[2]))
    return matplotlib.colors.LinearSegmentedColormap(name, cdict)
cmA = make_colormap([primary_colors["R"], primary_colors["H"], primary_colors["U"]])
cmB = make_colormap([d[1] for d in DIV_COLORS], extrapolate=True)

def make_layeredcmaps(depth_layers):
    cmaps = {}
    mx = max(depth_layers)
    for ci, c in enumerate(depth_layers):
        fgreen = 1.-c/mx
        cdict = {'red':   [(0.0,  0.0, 0.0),
                   (1.0,  1.0, 1.0)],
                'green': [(0.0,  fgreen, fgreen),
                   (1.0,  fgreen, fgreen)],
                'blue':  [(0.0,  1.0, 1.0),
                   (1.0,  0.0, 0.0)]}
        cmaps[c] = matplotlib.colors.LinearSegmentedColormap("Layer%d" % c, cdict)
    return cmaps


### DATA TOOLS
##########################
### VARIABLE NAMES
map_vnames = {
"MEAN_HYP": "HYP",
"MEAN_LOP": "LOP",
"MEAN_HOD": "HOD",
"MEAN_AL": "AL",
"MEAN_OL": "OL",
"MEAN_SF": "SF",
"MEAN_OT": "OT",
"MEAN_CM": "CM",
"MEAN_OO": "OO",
"MEAN_ETH": "ETH",
"MEAN_LOPT": "LOPT",
"NB_SPC": "nbSpc",
##############
"bio1:TMeanY": "T^{\\sim}{Y}",
"bio2:TMeanRngD": "T^{\\sim}{RngD}",
"bio3:TIso": "T{Iso}",
"bio4:TSeason": "T{Season}",
"bio5:TMaxWarmM": "T^{+}{WarmM}",
"bio6:TMinColdM": "T^{-}{ColdM}",
"bio7:TRngY": "T{RngY}",
"bio8:TMeanWetQ": "T^{\\sim}{WetQ}",
"bio9:TMeanDryQ": "T^{\\sim}{DryQ}",
"bio10:TMeanWarmQuarter": "T^{\\sim}{WarmQ}",
"bio11:TMeanColdQ": "T^{\\sim}{ColdQ}",
"bio12:PTotY": "P{TotY}",
"bio13:PWetM": "P{WetM}",
"bio14:PDryM": "P{DryM}",
"bio15:PSeason": "P{Season}",
"bio16:PWetQ": "P{WetQ}",
"bio17:PDryQ": "P{DryQ}",
"bio18:PWarmQ": "P{WarmQ}",
"bio19:PColdQ": "P{ColdQ}",
"bioX0:PDiffWetDryM": "P{Ampl}",
"bioX1:NPP": "{NPP}"}

map_vtex = {
"MEAN_HYP": "\\vLHSa{}",
"MEAN_LOP": "\\vLHSb{}",
"MEAN_HOD": "\\vLHSc{}",
"MEAN_AL": "\\vLHSd{}",
"MEAN_OL": "\\vLHSe{}",
"MEAN_SF": "\\vLHSf{}",
"MEAN_OT": "\\vLHSg{}",
"MEAN_CM": "\\vLHSh{}",
"MEAN_OO": "\\vLHSi{}",
"MEAN_ETH": "\\vLHSj{}",
"MEAN_LOPT": "\\vLHSk{}",
"NB_SPC": "\\vLHSl{}",
##############
"bio1:TMeanY": "\\vRHSa{}",
"bio2:TMeanRngD": "\\vRHSb{}",
"bio3:TIso": "\\vRHSc{}",
"bio4:TSeason": "\\vRHSd{}",
"bio5:TMaxWarmM": "\\vRHSe{}",
"bio6:TMinColdM": "\\vRHSf{}",
"bio7:TRngY": "\\vRHSg{}",
"bio8:TMeanWetQ": "\\vRHSh{}",
"bio9:TMeanDryQ": "\\vRHSi{}",
"bio10:TMeanWarmQuarter": "\\vRHSj{}",
"bio11:TMeanColdQ": "\\vRHSk{}",
"bio12:PTotY": "\\vRHSl{}",
"bio13:PWetM": "\\vRHSm{}",
"bio14:PDryM": "\\vRHSn{}",
"bio15:PSeason": "\\vRHSo{}",
"bio16:PWetQ": "\\vRHSp{}",
"bio17:PDryQ": "\\vRHSq{}",
"bio18:PWarmQ": "\\vRHSr{}",
"bio19:PColdQ": "\\vRHSs{}",
"bioX0:PDiffWetDryM": "\\vRHSt{}",
"bioX1:NPP": "\\vRHSu{}",
"bioX2:NPPt": "\\vRHSv{}",
"bioX3:NPPp": "\\vRHSw{}",
"bioX4:NPPtSlack": "\\vRHSx{}",
"bioX5:NPPpSlack": "\\vRHSy{}",
"Z:C": "\\vRHSba{}",
"Z:AsiaFocus": "\\vRHSbb{}",
"Z:SouthChina": "\\vRHSbc{}",
"Z:SouthEastAsia": "\\vRHSbd{}",
"Z:TibetPlateau": "\\vRHSbe{}",
"Z:India": "\\vRHSbf{}"}

def get_vname(xi, head_x):
    return map_vnames.get(head_x[xi], head_x[xi])
def get_vtex(xi, head_x):
    return map_vtex.get(head_x[xi], head_x[xi])
def get_vcolor(xi, head_x):
    if head_x[xi] == "bioX1:NPP":
        color = "Y"
    elif re.search("bio[0-9]*:P", head_x[xi]):
        color = "U"
    elif re.search("bio[0-9]*:", head_x[xi]):
        color = "R"
    else:
        color = "V"
    return color

### LOADING DATA
CLUSTER_CODES = {}
def conv_clu(ss):
    s = ss.decode("utf-8")
    if s == "nan":
        return -99
    elif re.match("S[0-9]+$", s):
        return int(s[1:])
    else:
        parts = s.split(":")
        if len(parts) == 1:
            g = "Z"
            v = s
        else:
            g = parts[0]
            v = ":".join(parts[1:])
        if g not in CLUSTER_CODES:
            CLUSTER_CODES[g] = {}
        if v not in CLUSTER_CODES[g]:
            CLUSTER_CODES[g][v] = len(CLUSTER_CODES[g])
        return CLUSTER_CODES[g][v]
    return -1
    
def conv_member(s):
    if s.strip() == b"F:out":
        return 0
    elif s.strip() == b"F:in":
        return 1
    return -1
def conv_cont(s):
    if s.strip() == b"C:AF":
        return 1
    elif s.strip() == b"C:EU":
        return 2
    elif s.strip() == b"C:NA":
        return 3
    elif s.strip() == b"C:SA":
        return 4
    return -1    

def load_ABC(file_A, file_B, file_C):
    with open(file_A) as fp:
        l = fp.readline()
    head_A = l.strip("# ").strip().split(",")
    map_A = dict([(v,k) for (k,v) in enumerate(head_A)])
    with open(file_B) as fp:
        l = fp.readline()
    head_B = l.strip("# ").strip().split(",")
    map_B = dict([(v,k) for (k,v) in enumerate(head_B)])
    with open(file_C) as fp:
        l = fp.readline()
    head_C = l.strip("# ").strip().split(",")
    map_C = dict([(v,k) for (k,v) in enumerate(head_C)])
    
    converters_B = {map_B["Z:C"]: conv_cont,
        map_B["Z:AsiaFocus"]: conv_member,
        map_B["Z:SouthChina"]: conv_member,
        map_B["Z:SouthEastAsia"]: conv_member,
        map_B["Z:TibetPlateau"]: conv_member,
        map_B["Z:India"]: conv_member}
    
    converters_C = dict([(k, conv_clu) for v,k in map_C.items() if v not in ["id","long","lat"]])
    A = numpy.loadtxt(file_A, delimiter=",", skiprows=2)
    B = numpy.loadtxt(file_B, delimiter=",", skiprows=2, converters=converters_B)
    C = numpy.loadtxt(file_C, delimiter=",", skiprows=2, converters=converters_C)
    return {"data": A, "head": head_A, "map": map_A}, {"data": B, "head": head_B, "map": map_B}, {"data": C, "head": head_C, "map": map_C}


### ANALYSIS TOOLS
##########################
def withen(mat):
    tt = numpy.std(mat, 0)
    tt[numpy.where(tt == 0)] = 1
    return (mat - numpy.tile(numpy.mean(mat, 0), (mat.shape[0], 1)))/numpy.tile(tt, (mat.shape[0], 1))

def regres_score(r): ## get score
    return r.aic
def regres_log(r, xs_lbl=None, y_lbl=None): ## get log string
    return "AIC=%d MSE=%s F=%d LL=%d\t[%s]" % (r.aic, r.mse_resid, r.fvalue, r.llf, regres_meq(r, xs_lbl, y_lbl))
def regres_log_tex(r, xs_lbl=None, y_lbl=None): ## get log string
    return "$%s$ & $%d$ & $%d$ & $%d$ & $%.2f$ \\\\" % (regres_meq(r, xs_lbl, y_lbl), r.aic, r.fvalue, r.llf, r.mse_resid)

def regres_leg(r, xs_lbl=None, y_lbl=None): ## get log string
    return "$AIC=%d\quadF=%d$" % (r.aic, r.fvalue)
def regres_meq(r, xs_lbl=None, y_lbl=None): ## get model equation
    if xs_lbl is None:
        if r.params.shape[0] == 2:
            xs_lbl = ["x"]
        else:
            xs_lbl = ["x_%d" % (i+1) for i in range(r.params.shape[0]-1)]
    if y_lbl is None:
        y_lbl = "y"        
    # eq = " ".join(["%+.2f %s" % (coeff, xl) for coeff, xl in list(zip(r.params[1:], xs_lbl))+ [(r.params[0], "")]])
    eq = " ".join(["%+.2f\\, %s" % (coeff, xl) for coeff, xl in [(r.params[0], "")]+list(zip(r.params[1:], xs_lbl))])
    if eq[0] == "+":
        eq = eq[1:]
    if len(y_lbl) > 0:
        eq = "%s = %s" % (y_lbl, eq)
    return eq.strip() 


### PLOTTING TOOLS
##########################
def plot_corrs(D, C, vrs, head, axes):
    nz = len(vrs)
    for i in range(nz):
        if i > 0:
            axes[i-1, 0].set_ylabel("$%s$" % get_vname(i, head))
        if i < nz-1:
            axes[nz-2, i].set_xlabel("$%s$" % get_vname(i, head))
            
        for j in range(i):
            axes[i-1, j].plot(D[:, vrs[j]], D[:, vrs[i]], '.', ms=1, color=cmA((C[i,j]+1)/2), alpha=ALPHA_DOTS)
            bx,tx = axes[i-1, j].get_xlim()
            by,ty = axes[i-1, j].get_ylim()
            axes[i-1, j].text((bx+tx)/2, (by+ty)/2, "%.3f" % C[i,j],
                                  horizontalalignment='center', verticalalignment='center', fontsize="large",
                                  bbox=dict(color='white', alpha=0.9, boxstyle=BoxStyle("Round", pad=0.1)))
            if j != i-1:
                axes[j, i-1].axis("off")
            else:
                if i > 1:
                    axes[i-1, j].set_ylabel("$%s$" % get_vname(i, head))
                    axes[i-1, j].yaxis.set_label_position('right') 
                if i < nz-1:
                    axes[i-1, j].set_xlabel("$%s$" % get_vname(j, head))
                    axes[i-1, j].xaxis.set_label_position('top')


def plot_proj(Xproj, Uproj, head, dname, axes, di, lbl_offsets={}):
    axes[1,di].plot(Xproj[1::1,0], Xproj[1::1,1], ".", ms=1, color=primary_colors["G"], alpha=ALPHA_DOTS)
    for i, c in enumerate(Uproj):
        color = primary_colors[get_vcolor(i, head)]
        axes[0,di].plot([0,c[0]], [0,c[1]], "-", color=color)
        axes[0,di].text(c[0], c[1], "$%s$" % get_vname(i, head), color=color)
        axes[0,di].spines['right'].set_visible(False)
        axes[0,di].spines['top'].set_visible(False)        
    axes[0,di].set_title(f"{dname}")


def plot_reg_one(X, y, results, xi=None, head_x={}, yi=None, head_y={}, ax=None, i=None, j=None):
    if ax is None:
        fig, ax = plt.subplots()        
    elif i is not None:
        if j is not None:
            ny = ax.shape[0]
            ax = ax[i,j]
        else:
            ax = ax[i]
    ax.plot(X[:,0], y, ".", ms=1, color=primary_colors["G"], alpha=ALPHA_DOTS)                
    xmm = [numpy.min(X[:,0]), numpy.max(X[:,0])]
    ymm = [numpy.min(y), numpy.max(y)]
    vs = numpy.vstack([numpy.ones(2), xmm])
    c = primary_colors["K"]
    # if results.aic > 0:
    #     c = primary_colors["G"]
    if results.f_pvalue < 0.001:
        l = "-"
    elif results.f_pvalue < 0.01:
        l = "--"
    else:
        l = ":"
    ax.plot(xmm, numpy.dot(results.params, vs), color=c)
    ax.set_xlim(xmm)
    ax.set_ylim([ymm[0]-0.15*(ymm[1]-ymm[0]), ymm[1]+0.15*(ymm[1]-ymm[0])])
    ## results.aic, results.f_pvalue, results.mse_total                
    ax.text((xmm[0]+xmm[1])/2, ymm[1]+0.12*(ymm[1]-ymm[0]), "$"+regres_meq(results)+"$",
                       horizontalalignment='center', verticalalignment='top', fontsize="medium",
                       bbox=dict(color='white', alpha=0.9, boxstyle=BoxStyle("Round", pad=0.1)))
    ax.text((xmm[0]+xmm[1])/2, ymm[0]-0.12*(ymm[1]-ymm[0]), regres_leg(results),
                       horizontalalignment='center', verticalalignment='bottom', fontsize="medium",
                       bbox=dict(color='white', alpha=0.9, boxstyle=BoxStyle("Round", pad=0.1)))
    if j is not None and j == 0:
        lbly = "$%s$" % get_vname(yi, head_y)
        ax.set_ylabel(lbly, fontsize="large")
    if i is not None and i+1 == ny:
        lblx = "$%s$" % get_vname(xi, head_x)
        ax.set_xlabel(lblx, fontsize="large")

def plot_clusters(axes, xi, yi, lng_lat, nm, nbc, cs, clusters_sets, clusters_info, clustering_info, cmap):
    axe, bm, bm_args = tools_geomap.prepare_map(map_prefs, map_corners[0]+map_corners[1], axe=axes[yi,xi])
    axes[yi,xi].set_xlabel("%s %d [%.3f]" % (nm, nbc, clustering_info["silhouette"]))
    xs, ys = bm(lng_lat[:,0], lng_lat[:,1])
    colors = -numpy.ones(lng_lat.shape[0])
    labels = []
    for cii, ci in enumerate(cs):
        ck = ((nm, nbc), ci)

        colors[list(clusters_sets[clusters_info[ck]["basis"]])] = clusters_info[ck]["ccolor"]
        labels.append((clusters_info[ck]["ccolor"], clusters_info[ck]["label"] if clusters_info[ck]["label"] != "Ao" else "---", clusters_info[ck]["size"]))
    for li, ll in enumerate(sorted(labels, reverse=True)):
        cx, cy = bm(map_corners[0][0]+0.01*(map_corners[0][1]-map_corners[0][0]),
                    map_corners[1][0]+((li+.5)/len(labels))*(map_corners[1][1]-map_corners[1][0]))
        axes[yi,xi].text(cx, cy, "%s:%d" % (ll[1], ll[2]), color=cmap(ll[0]), zorder=15,
                        horizontalalignment='left', verticalalignment='center', fontsize="small",
                        bbox=dict(color='white', alpha=0.9, boxstyle=BoxStyle("Round", pad=0.05)))
        # cx, cy = bm(clusters_info[ck]["center"][0], clusters_info[ck]["center"][1])
        # axes[yi,xi].text(cx, cy, clusters_info[ck]["label"], color=cmap(clusters_info[ck]["ccolor"]), zorder=15,
        #                 horizontalalignment='center', verticalalignment='center', # fontsize="large",
        #                 bbox=dict(color='white', alpha=0.9, boxstyle=BoxStyle("Round", pad=0.05)))
    axes[yi,xi].scatter(xs, ys, c=colors, s=.9, edgecolors='none', vmin=0, vmax=1, zorder=10, alpha=.85, cmap=cmap)
