import re, string, numpy, codecs, itertools, os.path
from classQuery import SYM
from classSParts import SSetts, tool_ratio
import pdb

SIDE_CHARS = {0:"L", 1:"R", -1: "C"}
HAND_SIDE = {"LHS": 0, "RHS": 1, "0": 0, "1": 1, "COND": -1, "-1": -1}
NUM_CHARS = dict([(numpy.base_repr(ii, base=25), "%s" % chr(ii+ord("a"))) for ii in range(25)])

WIDTH_MID = .5


####   TOOLS FOR LATEX PRINTING
def digit_to_char(n, pad=None):
    if pad is None:
        tmp = "".join([NUM_CHARS[t] for t in numpy.base_repr(n, base=25)])
    else:
        tmp = ("z"*pad+"".join([NUM_CHARS[t] for t in numpy.base_repr(n, base=25)]))[-pad:]
    # print "%s -> %s" % (n, tmp)
    return tmp

def side_ltx_cmd(sid, padss=None):
    if sid in SIDE_CHARS:
        schar = SIDE_CHARS[sid]
    else:
        schar = digit_to_char(sid, padss)
    return schar

def var_ltx_cmd(sid, vid, padsv=None, padss=None):
    schar = side_ltx_cmd(sid, padss)
    vchar = digit_to_char(vid, padsv)
    return "\\v%sHS%s" % (schar.upper(), vchar.lower())

#### TODO update TeX commands
TEX_PCK = "\\usepackage{amsmath}\n"+ \
              "\\usepackage{amsfonts}\n"+ \
              "\\usepackage{amssymb}\n"+ \
              "\\usepackage{booktabs}\n"+ \
              "\\usepackage[mathletters]{ucs}\n"+ \
              "\\usepackage[utf8x]{inputenc}\n"
TEX_CLR = "\\usepackage{color}\n"+ \
              "\\definecolor{LHSclr}{rgb}{.855, .016, .078} %% medium red\n"+ \
              "% \\definecolor{LHSclr}{rgb}{.706, .012, .063} %% dark red\n"+ \
              "% \\definecolor{LHSclr}{rgb}{.988, .345, .392} %% light red\n"+ \
              "\\definecolor{RHSclr}{rgb}{.055, .365, .827} %% medium blue\n"+ \
              "% \\definecolor{RHSclr}{rgb}{.043, .298, .682} %% dark blue\n"+ \
              "% \\definecolor{RHSclr}{rgb}{.455, .659, .965} %% light blue\n"+ \
              "\\definecolor{LCclr}{rgb}{.50,.50,.50} %% medium gray\n"+ \
              "\\definecolor{RNclr}{rgb}{.40, .165, .553} %% medium purple\n"
TEX_CMD = "\\newcommand{\\iLHS}{\\mathbf{L}} % index for left hand side\n"+ \
              "\\newcommand{\\iRHS}{\\mathbf{R}} % index for right hand side\n"+ \
              "\\newcommand{\\iCOND}{\\mathbf{C}} % index for conditional\n"+ \
              "\\newcommand{\\SubCond}{C} % index for learn subset\n"+ \
              "\\newcommand{\\RSetLearn}{\\mathcal{O}} % index for learn subset\n"+ \
              "\\newcommand{\\RSetTest}{\\mathcal{I}} % index for test subset\n"+ \
              "\\newcommand{\\RSetAll}{} % index for whole\n"+ \
              "\\newcommand{\\RSetAA}{\\mathcal{I}} % index for whole non empty\n"+ \
              "\\newcommand{\\query}{q} % index for whole\n"+ \
              "\\newcommand{\\RName}[1]{\\textcolor{RNclr}{R#1}} \n"+ \
              "%% \\renewcommand{\\land}{\\text{\\textcolor{LCclr}{~AND~}}} \n"+ \
              "%% \\renewcommand{\\lor}{\\text{\\textcolor{LCclr}{~OR~}}} \n\n"+ \
              "\\newcommand{\\abs}[1]{\\vert#1\\vert} % absolute value\n"+ \
              "\\newcommand{\\perc}[1]{\\% #1} % perc \n"+ \
              "\\newcommand{\\ratio}[1]{\\div #1} % ratio \n"+ \
              "\\DeclareMathOperator*{\\pValue}{pV}\n"+ \
              "\\DeclareMathOperator*{\\jacc}{J}\n"+ \
              "\\DeclareMathOperator*{\\supp}{supp}\n"+ \
              "\\DeclareMathOperator*{\\pr}{p}\n"

def openTexDocument(names, standalone=False):
    ####### prepare variable names    
    names_alts = []
    names_commands = ""
    numvs_commands = ""
    macvs_commands = ""
    macfld_commands = "\\newcommand{\\PP}[3]{#1{#2}#3} \n"
    padss = len(digit_to_char(len(names)-1))
    for i, ns in enumerate(names):
        if ns is not None:
            names_alts.append([])
            padsv = len(digit_to_char(len(ns)-1))
            sltx = side_ltx_cmd(i, padss)
            macvs_commands += "\\newcommand{\\varname%s}[1]{\\text{\\textcolor{%sHSclr}{#1}}} \n" % (sltx, sltx)
            for ni, n in enumerate(ns):
                vlc = var_ltx_cmd(i, ni, padsv, padss)
                names_alts[i].append("$"+vlc+"{}$")
                tmp = re.sub("#", "\#", re.sub("_", "\\_", re.sub("\\\\", "\\\\textbackslash{}", n)))
                names_commands += "\\newcommand{%s}{\\varname%s{%s}}\n" % (vlc, sltx, tmp)
                numvs_commands += "%% \\newcommand{%s}{\\varname%s{v%d}}\n" % (vlc, sltx, ni)
        else:
            names_alts.append(None)
    if standalone:
        str_out = "\\documentclass{article}\n"+ \
              TEX_PCK + TEX_CLR + TEX_CMD + \
              macvs_commands + \
              names_commands+ \
              numvs_commands+ \
              macfld_commands+ \
              "\\begin{document}\n"+ \
              "\\begin{table}[h]\n"+ \
              "\\scriptsize\n"
    else:
        str_out = "\\scriptsize\n"
    return str_out, names_alts

def openTexTabular(list_fields, nblines=1):
    str_out = ""
    with_fname = False
    if nblines == 3:
        #### SPREAD ON THREE LINES
        str_out += "\\begin{tabular}{@{\\hspace*{1ex}}p{0.05\\textwidth}@{\\hspace*{3ex}}"+ ("p{%.3f\\textwidth}" % WIDTH_MID)
        for i in range(len(list_fields)-2):
            str_out += "@{\\hspace*{2ex}}l"
        str_out += "@{\\hspace*{1ex}}}\n"
        str_out += "\\toprule\n"
        with_fname=True
        
    elif nblines == 2:
        #### SPREAD ON TWO LINES 
        str_out += "\\begin{tabular}{@{\\hspace*{1ex}}r@{\\hspace*{1ex}}p{0.67\\textwidth}@{}r"
        for i in range(len(list_fields)-2):
            str_out += "@{\\hspace*{2ex}}r"
        str_out += "@{\\hspace*{1ex}}}\n\\toprule\n"

        str_out += " & " + RedProps.dispHeaderFields(list_fields, style="tex") + " \\\\\n"
        str_out += "%%% & " + RedProps.dispHeaderFields(list_fields) + " \\\\\n"
        str_out += "\\midrule\n"

    else:
        #### SINGLE LINE
        str_out += "\\begin{tabular}{@{\\hspace*{1ex}}r@{\\hspace*{1ex}}p{0.35\\textwidth}@{\\hspace*{1em}}p{0.35\\textwidth}"
        for i in range(len(list_fields)-2):
            str_out += "@{\\hspace*{2ex}}r"
        str_out += "@{\\hspace*{1ex}}}\n\\toprule\n"

        str_out += " & " + RedProps.dispHeaderFields(list_fields, style="tex") + " \\\\\n"
        str_out += "%%% & " + RedProps.dispHeaderFields(list_fields) + " \\\\\n"
        str_out += "\\midrule\n"
    return str_out, with_fname

def closeTexTabular():
    return "\\bottomrule\n\\end{tabular}\n"

def closeTexDocument(standalone):
    str_out = ""
    if standalone:
        str_out += "\\end{table}\n\\end{document}"
    ### auctex vars
    str_out += "\n%%%%%% Local Variables:\n"+ \
      "%%%%%% mode: latex\n"+ \
      "%%%%%% TeX-master: t\n"+ \
      "%%%%%% End:\n"
    return str_out
##############################################
####   TOOLS FOR FORMATTING AND PRINTING
def prepareValFmt(val, fmt={}, to_str=True, replace_none="-"):
    if val is None:
        if to_str:
            return replace_none
        return val
    fmtf = "%"+fmt.get("fmt", "s")
    if "rnd" in fmt:
        val = round(val, fmt["rnd"])
    elif fmt.get("supp_set", False):
        fmtf = "%s"
        val = fmt.get("sep", ",").join([("%"+fmt.get("fmt", "s")) % v for v in val])
    if to_str:
        return fmtf % val
    return val
    
def prepareFmtString(nblines, nbstats, last_one, tex=False):
    nbc = "%d" % (nbstats+1)
    if nblines == 3:
        #### SPREAD ON THREE LINES            
        if tex:
            frmts = "%(rid)s & & %(stats)s \\\\ [.2em]\n & \\multicolumn{"+ nbc +"}{p{.9\\textwidth}}{ %(q0)s } \\\\ \n" + \
                    " & \\multicolumn{"+ nbc +"}{p{.9\\textwidth}}{ %(q1)s } \\\\"
            if not last_one:
                frmts +=  " [.32em] \cline{2-"+ nbc +"} \\\\ [-.88em]"
        else:
            frmts = "%(rid)s%(stats)s\n%(q0)s\n%(q1)s"
    elif nblines == 2:
        if tex:
            frmts = "%(rid)s & %(q0)s & & %(stats)s \\\\ \n & \\multicolumn{2}{r}{ %(q1)s } \\\\"
            if not last_one:
                frmts +=  " [.3em]"
        else:
            frmts = "%(rid)s%(q0)s\t%(q1)s\n%(stats)s"
    else:
        if tex:
            frmts = "%(rid)s & %(all)s \\\\"
        else:
            frmts = "%(rid)s%(all)s"
    return frmts
##############################################
def cust_eval(exp, loc):
    try:
        return eval(exp, {}, loc)
    except ZeroDivisionError:
        return float("Inf")
    except TypeError:
        return None

def expand_val(val):
    v = val.strip()
    if v == "Ex*":
        return SSetts.labels[:4]
    elif v == "Em*":
        return SSetts.labels[4:]
    return [val.strip()]


class dummyC:
    def __init__(self):
        print "WARNING: Class not properly set up!"
    def initParsed(self):
        print "WARNING: Class not properly set up!"
    def parse(self):
        print "WARNING: Class not properly set up!"

class RedProps(object):

    pref_dir = os.path.dirname(os.path.abspath(__file__))
    def_file_basic = pref_dir+"/fields_defs_basic.txt"
    default_def_files = [def_file_basic]
    
    rset_match = "("+ "|".join(["all","learn","test","cond"] + HAND_SIDE.keys()) +")"
    what_match = "\w+"
    which_match = "\w+" 
    match_primitive = "(?P<prop>((?P<rset_id>"+rset_match+"))?:(?P<what>"+what_match+"):(?P<which>"+which_match+")?)"
        
    all_what_match = ""
    all_which_match = ""
    all_match_primitive = "(?P<prop>((?P<rset_id>"+rset_match+"))?:(?P<what>"+all_what_match+"):(?P<which>"+all_which_match+")?)"
    lbl_match = "(?P<what>"+all_what_match+")[_]?(?P<which>"+all_which_match+")[_]?(?P<rset_id>"+rset_match+")?"
    
    @classmethod
    def setupRProps(tcl, Qclass=dummyC, Rclass=dummyC, elems_typs=[]):
        tcl.Qclass = Qclass
        tcl.Rclass = Rclass
        tcl.elems_typs = elems_typs
        tcl.setupRPMatches()
    @classmethod
    def setupRPMatches(tcl):
        tcl.all_what_match = "("+ "|".join(["(?P<"+ preff +"what>"+cc.Pwhat_match+")" for (preff, cc) in tcl.elems_typs]) +")"
        tcl.all_which_match = "("+ "|".join(["(?P<"+ preff +"which>"+cc.Pwhich_match+")" for (preff, cc) in tcl.elems_typs]) +")"   
        tcl.all_match_primitive = "(?P<prop>((?P<rset_id>"+tcl.rset_match+"))?:(?P<what>"+tcl.all_what_match+"):(?P<which>"+tcl.all_which_match+")?)"
        tcl.lbl_match = "(?P<what>"+tcl.all_what_match+")[_]?(?P<which>"+tcl.all_which_match+")[_]?(?P<rset_id>"+tcl.rset_match+")?"


    substs = [("queryLHS", "query_LHS"), ("queryRHS", "query_RHS"), ("card_", "len"), ("alpha", "Exo"), ("beta", "Eox"), ("gamma", "Exx"), ("delta", "Eoo"), ("mua", "Exm"), ("mub", "Emx"), ("muaB", "Eom"), ("mubB", "Emo"), ("mud", "Emm"), ("status_disabled", "status_enabled")]
    
    # comm_tex = "\\newcommand{\\PP}[3]{#1{#2}#3}",
    map_lbls = {}
    map_lbls["txt"] = {"what": {}, "which": {}, "rset": {"all" : ""}}
    map_lbls["txt"]["specials"] = {}
    map_lbls["tex"] = {"what": {"query": "\\query", "acc": "\\jacc", "pval": "\\pValue",
                                "len": "\\abs", "card": "\\abs",
                                "ratio": "\\ratio", "perc": "\\perc"},
                       "which": {"I": "\\supp"},
                       "rset": {"all" : "\\RSetAll", 
                                "cond" : "\\SubCond",
                                "learn" : "_{\\RSetLearn}", "test" : "_{\\RSetTest}",
                                "ratioTL" : "_{\\RSetTest/\\RSetLearn}",
                                "ratioTA" : "_{\\RSetTest/\\RSetAA}", "ratioLA" : "_{\\RSetLearn/\\RSetAA}",
                                "LHS" : "_\\iLHS", "RHS" : "_\\iRHS", "COND" : "_\\iCOND",
                                "0" : "_\\iLHS", "1" : "_\\iRHS", "-1" : "_\\iCOND",}}
    map_lbls["tex"]["specials"] = {}
    map_lbls["gui"] = {"what": {"acc": "J", "pval": "pV", "perc": "%", "ratio": "/",
                                "len": ("|", "|"), "card": ("|", "|"),
                                "set": "", "supp": ""},
                       "which": {"I": "supp"},
                       "rset": {"all" : "",
                                "cond" : "C",
                                "learn" : SYM.SYM_LEARN, "test" : SYM.SYM_TEST,
                                "ratioTL" : "T/L",
                                "ratioTA" : "T/A", "ratioLA" : "L/A",
                                "0" : "LHS", "1" : "RHS", "-1" : "COND",}}
    map_lbls["gui"]["specials"] = {}
    lbl_patts = {"txt": "%(what)s_%(which)s_%(rset)s",
                 "tex": "$\\PP{%(what)s}{%(which)s}{%(rset)s}$",
                 "gui": "%(rset)s %(what)s%(which)s%(what_1)s"}
    def_elems = {"what": "", "which": "", "rset": "", "what_1": ""}


    modifiers_defaults = {"wsplits": False, "wmissing": False, "wcond": False}
    @classmethod
    def updateModifiers(tcl, red_list, modifiers={}):
        if ("wcond" not in modifiers) and any([red[1].hasCondition() for red in red_list]):
            modifiers["wcond"] = True
        if ("wmissing" not in modifiers) and any([red[1].hasMissing() for red in red_list]):
            modifiers["wmissing"] = True
        if ("wsplits" not in modifiers) and any([red[1].hasRSets() for red in red_list]):
            modifiers["wsplits"] = True
        return modifiers
    def getModifiersForData(tcl, data):
        if data is None:
            return {}
        return {"wsplits": data.hasLT(),
                "wcond": data.isConditional(),
                "wmissing": data.hasMissing()}
    @classmethod
    def hashModifiers(tcl, modifiers):
        return "%(wmissing)d%(wcond)d%(wsplits)d" % modifiers
    @classmethod
    def getStyles(tcl, style):
        if style in ["tex", "txt"]:
            return style, style
        elif style == "gui":
            return style, "U"
        else:                
            return "txt", style

    
    @classmethod
    def isPrimitive(tcl, exp):
        return re.match(tcl.match_primitive, exp) is not None
    @classmethod
    def getPrimitiveWs(tcl, exp):
        # all_mtch = re.match(tcl.all_match_primitive, exp)
        mtch = re.match(tcl.match_primitive, exp)
        if mtch is not None:
            return (mtch.group("what"), mtch.group("which") or "", mtch.group("rset_id") or "all")
        return (None, None, None)

    ### PRIMITIVE FORMATS AND LABELS
    @classmethod
    def getPrimitiveFmt(tcl, exp):
        what, which, rset = tcl.getPrimitiveWs(exp)
        if what is not None:
            if what in ["len", "card", "depth", "width"] or re.search("nb$", what):
                return {"fmt": "d"}
            if what in ["supp", "set"] or re.search("set$", what):
                return {"fmt": "s", "supp_set": True, "sep": ", "}
            if what in ["perc"]:
                return {"fmt": ".2f", "rnd": 2}
            if what in ["ratio"]:
                return {"fmt": ".4f", "rnd": 4}
            if what in ["acc", "pval"]:
                return {"fmt": ".3f", "rnd": 3}
        return {"fmt": "s"}

    @classmethod        
    def getPrimitiveLbl(tcl, exp, style="txt"):
        return tcl.prepareWsLbl(tcl.getPrimitiveWs(exp), style=style)
    @classmethod        
    def prepareWsLbl(tcl, ws, style="txt"):
        (what, which, rset) = ws
        if what is not None:
            if (what, which, rset) in tcl.map_lbls[style]["specials"]:
                return tcl.map_lbls[style]["specials"][(what, which, rset)]
            lbl_elems = dict(tcl.def_elems)
            for pk, pp in [("rset", rset), ("what", what), ("which", which)]:
                tt = tcl.map_lbls[style][pk].get(pp, pp)
                if type(tt) is tuple:
                    lbl_elems.update(dict([(pk+("_%d" % ti), tk) for ti, tk in enumerate(tt)]))
                    tt = tt[0]
                lbl_elems.update({pk: tt})
            if rset in HAND_SIDE and what != "query" and lbl_elems["which"] == "":
                lbl_elems["which"] = "q"                
            tmp = tcl.lbl_patts[style] % lbl_elems
            return re.sub("__+", "_", tmp.strip("_"))
        return "??"
    @classmethod
    def getPrimitiveNameForLbl(tcl, ll):
        if ll in tcl.map_lbls["txt"]["specials"].values():
            return dict([(v,k) for (k,v) in tcl.map_lbls["txt"]["specials"].items()])[ll]
        fld = ll        
        for (sf, st) in tcl.substs:
            fld = fld.replace(sf, st)
        parts = {"rset_id": "all", "what": None, "which": ""}
        ptyps = {"what": None, "which": None}
        mtch_all = re.match(tcl.lbl_match+"$", fld)
        if mtch_all is not None:
            if mtch_all.group("rset_id") is not None:
                parts["rset_id"] = mtch_all.group("rset_id")
            for part in ["which", "what"]:
                if mtch_all.group(part) is not None:
                    parts[part] = mtch_all.group(part)
                    ptyps[part] = [preff for preff, cc in tcl.elems_typs if mtch_all.group(preff+part) is not None]
        if parts["what"] is not None:
            return "%(rset_id)s:%(what)s:%(which)s" % parts
        return "??"

    @classmethod
    def getFieldLbl(tcl, ff, style="txt"):
        tmp = tcl.getPrimitiveLbl(ff, style)
        if tmp == "??":
            prts = ff.split("_")
            if len(prts) == 1:
                ws = (prts[0], "", "")
            elif len(prts) == 2:
                ws = (prts[0], "", prts[1])
            elif len(prts) == 3:
                ws = (prts[0], prts[1], prts[2])
            tmp = tcl.prepareWsLbl(ws, style=style)
            if tmp == "??":
                return ff
            return tmp
        return tmp
    @classmethod
    def getFieldKeyforLbl(tcl, ll):
        tmp = tcl.getPrimitiveKeyForLbl(ll)
        if tmp == "??":
            return ff
        return tmp 

    ### PRIMITIVE EVAL
    @classmethod
    def getPrimitiveVal(tcl, red, exp, details={}):
        what, which, rset = tcl.getPrimitiveWs(exp)
        if what is not None:
            return red.getProp(what, which, rset, details)        
    @classmethod
    def getPrimitiveValFormatted(tcl, red, exp, details={}, to_str=True):
        fmt = tcl.getPrimitiveFmt(exp)
        val = tcl.getPrimitiveVal(red, exp, details)
        return prepareValFmt(val, fmt, to_str)        

    ### DERIVATIVE EVAL
    @classmethod
    def compEVal(tcl, red, exp, details={}):
        texp = exp
        Rd = {}
        for mtch in list(re.finditer(tcl.match_primitive, exp))[::-1]:
            Rd[mtch.group("prop")] = red.getProp(mtch.group("what"), mtch.group("which"), mtch.group("rset_id"), details)
            texp = texp[:mtch.start()] + ('R["%s"]' % mtch.group("prop")) + texp[mtch.end():]
        return cust_eval(texp, loc={"R": Rd})
    @classmethod
    def compEVals(tcl, red, exps, details={}):
        props_collect = set()
        trans_exps = {}
        for eid, exp in exps.items():
            texp = exp
            for mtch in list(re.finditer(tcl.match_primitive, exp))[::-1]:
                props_collect.add(mtch)
                texp = texp[:mtch.start()] + ('R["%s"]' % mtch.group("prop")) + texp[mtch.end():]
            trans_exps[eid] = texp        
        Rd = {}
        for prop in props_collect:
            Rd[prop.group("prop")] = red.getProp(prop.group("what"), prop.group("which"), prop.group("rset_id"), details)
        props = {}
        for eid, exp in trans_exps.items():
            props[eid] = cust_eval(exp, loc={"R": Rd})
        return props

    @classmethod
    def refreshEVals(tcl, red, exps, details={}):
        red.setCacheEVals(tcl.compEVals(red, exps, details))
    @classmethod
    def getEVal(tcl, red, k, exp=None, fresh=False, default=None, details={}):
        if (exp is not None) and (not red.hasCacheEVal(k) or fresh):
            red.setCacheEVal(k, tcl.compEVal(red, exp, details))
        return red.getCacheEVal(k)
    @classmethod
    def getEVals(tcl, red, exps, fresh=False, details={}):
        if fresh:
            revals = exps
        else:
            revals = dict([(k,v) for (k,v) in exps.items() if not red.hasCacheEVal(k)])
        if len(revals) > 0:
            red.updateCacheEVals(tcl.compEVals(red, revals, details))
        return dict([(k,red.getCacheEVal(k)) for k in exps.keys()])

    def getEValF(self, red, k, fresh=False, default=None, details={}):
        exp = self.getFieldExp(k)
        return self.getEVal(red, k, exp=exp, fresh=fresh, default=default, details=details)
    def getEValGUI(self, red, details):
        if "k" in details:
            val = self.getEValF(red, details["k"], details=details)
            return self.formatVal(val, details["k"], to_str=details.get("to_str", False), replace_none=details.get("replace_none", "-"))
        return None
    
    def __init__(self, fields_fns=None):
        self.derivatives = {}
        self.setupFDefs(fields_fns)

    def setupFDefs(self, fields_fns=None):
        if fields_fns is None:
            fields_fns = self.default_def_files
        self.field_lists = {}
        for ff in fields_fns:
            current_list_name, current_list_fields = (None, None)
            try:
                with open(ff) as fp:            
                    current_list_name, current_list_fields = self.readFieldsFile(fp)
            except IOError:
                print "Cannot read fields defs from file %s!" % ff

            else:
                if current_list_name is not None:
                    self.setFieldsList(current_list_name, current_list_fields)


                
    def readFieldsFile(self, fields_fp):
        current_list_fields = []
        current_list_name = None        
        for line in fields_fp:
            prts = line.strip().split("\t")
            if prts[0] == "field" and len(prts) > 2:
                dets = {"exp": prts[2], "fmt": "s"}
                if prts[3] == "set":
                    dets.update({"supp_set": True, "sep": ", "})
                else:
                    dets["fmt"] = prts[3]
                    mtc = re.search("\.(?P<rnd>[0-9]+)f", prts[3])
                    if mtc is not None:
                        dets["rnd"] = int(mtc.group("rnd"))
                self.addDerivativeField(prts[1], dets)
            if prts[0] == "fieldlist" and len(prts) == 2:
                if current_list_name is not None:
                    self.setFieldsList(current_list_name, current_list_fields)
                current_list_name = prts[1]
                current_list_fields = []
            elif current_list_name is not None:
                current_list_fields.extend(self.parseFieldsLine(line))
        return current_list_name, current_list_fields
        
                
    modifiers_mtch = "((\[(?P<modify>[^\]]*)\])|((?P<plc_hld>\w)\=)|(\{(?P<vals>[^}]*)\})|(?P<spc> +))"
    def parseFieldsLine(self, line):
        #### HERE        
        prts = line.strip().split("\t")
        if len(prts) == 1:
            return [{"field": prts[0]}]
        elif prts[0] == "list" and len(prts) > 1 and self.hasFieldsList(prts[1]):
            return [dict(f) for f in self.getFieldsList(prts[1])]
        elif len(prts) == 2:
            combs = [[]]
            pile_mods = []
            prev_mod = False
            current_plc_hld = None
            plc_hlds = []
            for mtch in re.finditer(self.modifiers_mtch, prts[1]):
                elem = [m for m in mtch.groupdict().items() if m[1] is not None][0]
                if elem[0] == "modify":
                    prev_mod = True
                    pile_mods.append(elem[1])
                elif elem[0] == "plc_hld":
                    current_plc_hld = elem[1]
                elif elem[0] == "vals":
                    mods = [p for p in pile_mods]
                    if not prev_mod and len(mods) > 0:
                        mods.append(-1)
                        pile_mods.pop()
                    prev_mod = False
                    for val in elem[1].split(","):
                        for v in expand_val(val):
                            combs[-1].append((v, mods))
                elif elem[0] == "spc":
                    plc_hlds.append(current_plc_hld)
                    if prev_mod and len(combs[-1]) == 0:
                        combs[-1].append((None, tuple(pile_mods)))
                    combs.append([])
                    pile_mods = []
                    prev_mod = False
                    current_plc_hld = None
            plc_hlds.append(current_plc_hld)
            if prev_mod and len(combs[-1]) == 0:
                combs[-1].append((None, tuple(pile_mods)))
            fields = []
            for p in itertools.product(*combs):
                vals, mods = ([], [])
                for pp in p:
                    vals.append(pp[0])
                    if len(pp[1]) > 1 and pp[1][-1] == -1:
                        mods.append("not (%s)" % " and ".join(pp[1][:-1]))
                    else:
                        mods.extend(pp[1])
                fld = prts[0]
                for vi, v in enumerate(vals):
                    if v is not None:
                        fld = fld.replace(plc_hlds[vi], v)
                fields.append({"field": fld})
                if len(mods) > 0:
                    fields[-1]["modify"] = " and ".join(mods)
            return fields
        return []
    
    def addDerivativeField(self, fk, ff):
        self.derivatives[fk] = ff
    def hasDerivativeField(self, fk):
        return fk in self.derivatives
    def getDerivativeField(self, fk):
        return self.derivatives[fk]
    def getAllDerivativeFields(self):
        return self.derivatives.keys()
    def hasFieldsList(self, fk):
        return fk in self.field_lists
    def getFieldsList(self, fk):
        return self.field_lists.get(fk, [])
    def setFieldsList(self, fk, flist):
        self.field_lists[fk] = flist
    def getListsKeys(self):
        return self.field_lists.keys()


    def getFieldExp(self, ff):
        if self.hasDerivativeField(ff):
            return self.getDerivativeField(ff)["exp"]
        elif self.isPrimitive(ff):            
            return ff
        elif ff in self.Rclass.info_what:
            return ":%s:" % ff

    def getFieldFmt(self, ff):
        if self.hasDerivativeField(ff):
            return self.getDerivativeField(ff)
        else:
            return self.getPrimitiveFmt(ff)
    def getFieldDelim(self, ff, style="txt"):
        if style == "tex" and  re.search("[df]$", self.getFieldFmt(ff).get("fmt")):
            return "$"
        return ""
    def getFieldDets(self, ff):
        tmp = {}
        if self.hasDerivativeField(ff):
            return self.getDerivativeField(ff)["exp"]
        elif self.isPrimitive(ff):
            return ff

    def formatVal(self, val, fk, to_str=True, replace_none="-"):
        fmt = self.getFieldFmt(fk)
        return prepareValFmt(val, fmt, to_str, replace_none)


    #### PREPARING DERIVATIVE FIELDS
    def getListFields(self, lfname, modifiers={}):
        ddmod = dict(self.modifiers_defaults)
        ddmod.update(modifiers)
        fields = []
        for f in self.getFieldsList(lfname):
            if cust_eval(f.get("modify", "True"), ddmod):
                fields.append(f["field"])
        return fields

    def getCurrentListFields(self, ckey, modifiers={}):
        ddmod = dict(self.modifiers_defaults)
        ddmod.update(modifiers)
        hmod = self.hashModifiers(ddmod)
        options = ["%s-%s" % (ckey, hmod), "%s" % ckey, "basic"]
        for opt in options:
            if self.hasFieldsList(opt):
                return self.getListFields(opt, modifiers)
    def setCurrentListFields(self, ll, ckey, modifiers={}):
        ddmod = dict(self.modifiers_defaults)
        ddmod.update(modifiers)
        hmod = self.hashModifiers(ddmod)
        ck = "%s-%s" % (ckey, hmod)
        flist = [{"field": f} for f in ll]
        self.setFieldsList(ck, flist)

    def getAllFields(self, ckey, modifiers={}):
        ddmod = dict(self.modifiers_defaults)
        ddmod.update(modifiers)
        all_fields = {}
        for lfname in self.getListsKeys():
            for fi, f in enumerate(self.getFieldsList(lfname)):
                if cust_eval(f.get("modify", "True"), ddmod):
                    field = f["field"]
                    if field not in all_fields:
                        all_fields[field] = [0.,0]
                    all_fields[field][0] += fi
                    all_fields[field][1] += 1
        for field in self.getAllDerivativeFields():
            if field not in all_fields:
                all_fields[field] = [1,.1]

        afs = all_fields.keys()
        afs.sort(key=lambda x: all_fields[x][0]/all_fields[x][1])
        return afs
        
    def getExpDict(self, lfields):
        return dict([(f, self.getFieldExp(f)) for f in lfields])

    ##############################################
    ####        PRINTING
    ##############################################
    @classmethod
    def dispHeaderFields(tcl, list_fields, style="txt", sep=None):
        dstyle, qstyle = tcl.getStyles(style)
        if sep is None:
            if dstyle == "tex":
                sep = " & "
            else:
                sep = "\t"
        return sep.join([tcl.getFieldLbl(f, dstyle) for f in list_fields])

    def dispHeader(self, list_fields="basic", modifiers={}, style="txt", sep=None):
        if not type(list_fields) is list:
            list_fields = self.getListFields(list_fields, modifiers)
        return self.dispHeaderFields(list_fields, style, sep)


    def disp(self, red, names=[None, None], row_names=None, with_fname=False, rid="", nblines=1, delim="", last_one=False, list_fields="basic", modifiers={}, style="txt", sep=None):
        dstyle, qstyle = self.getStyles(style)
        if not type(list_fields) is list:
            list_fields = self.getListFields(list_fields, modifiers)
        if sep is None:
            if dstyle == "tex":
                sep = " & "
            else:
                sep = "\t"
        details = {"style": qstyle}
        if names[0] is not None or names[1] is not None:
            details["names"] = names
            details["named"] = True
        if row_names is not None:
            details["row_names"] = row_names
            details["named"] = True

        exp_dict = self.getExpDict(list_fields)
        evals_dict = self.compEVals(red, exp_dict, details)

        if type(with_fname) == list and len(with_fname) == len(exp_dict):
            lbls = with_fname
            with_fname = True
        elif with_fname:
            lbls = ["%s=" % self.getFieldLbl(f, dstyle) for f in list_fields]
            
        dts = {}
        lbl = ""
        entries = {"stats": [], "all": [], "q0": "", "q1": "", "rid": rid}
        for fid, field in enumerate(list_fields):
            entry = self.formatVal(evals_dict[field], field, to_str=True)                    
            delim = self.getFieldDelim(field, dstyle)
            if with_fname:
                lbl = lbls[fid]
            dts[field] = (lbl+delim+entry+delim)
            if re.search("LHS:query:", field):
                entries["q0"] = dts[field]
            elif re.search("RHS:query:", field):
                entries["q1"] = dts[field]
            else:
                entries["stats"].append(dts[field])
            entries["all"].append(dts[field])
        nbstats = len(entries["stats"])
        nball = len(entries["all"])
        entries["all"] = sep.join(entries["all"])
        entries["stats"] = sep.join(entries["stats"])
        frmts = prepareFmtString(nblines, nbstats, last_one, dstyle == "tex")
        return frmts % entries
    
    ##############################################
    ####        PRINTING
    ##############################################
    def printRedList(self, reds, names=[None, None], fields=None, full_supp=False, supp_names=None, nblines=1, modifiers={}, style="txt"):
        try:
            red_list = sorted(reds.items())
        except AttributeError:
            red_list = list(enumerate(reds))

        modifiers = self.updateModifiers(red_list, modifiers)
        ckey = "txt"
        if style == "tex": ckey = "tex"
        all_fields = self.getCurrentListFields(ckey, modifiers)
        if type(fields) is list and len(fields) > 0:
            if fields[0] == -1:
                all_fields.extend(fields[1:])
            else:
                all_fields = fields
        if full_supp:
            all_fields.extend(self.getListFields("supps", modifiers))

        if style == "tex":
            return self.printTexRedList(reds=reds, names=names, list_fields=all_fields, nblines=nblines)
        str_out = self.dispHeader(all_fields, sep="\t", style=style) 
        for ri, red in red_list:
            str_out += "\n" + self.disp(red, names=names, row_names=supp_names, nblines=nblines, list_fields=all_fields, sep="\t", style=style)
        return str_out

    def printTexRedList(self, reds, names=[None, None], list_fields=None, nblines=1, standalone=False, modifiers={}):
        standalone = True
        try:
            red_list = sorted(reds.items())
            last_ri = red_list[-1][0]
        except AttributeError:
            red_list = list(enumerate(reds))
            last_ri = len(reds)-1

        if list_fields is None:
            modifiers = self.updateModifiers(red_list, modifiers)
            list_fields = self.getCurrentListFields("tex", modifiers)

        str_odoc, names_alts = openTexDocument(names, standalone)
        str_oth, with_fname = openTexTabular(list_fields, nblines)
        str_reds = ""
        for ri, red in red_list:
            ridstr = "\RName{%s}" %ri
            str_reds += self.disp(red, names_alts, list_fields=list_fields, style="tex", sep=" & ", with_fname=with_fname, rid=ridstr, nblines=nblines, last_one=(ri == last_ri)) + "\n" 
        str_cth = closeTexTabular()
        str_cdoc = closeTexDocument(standalone)
        return str_odoc + str_oth + str_reds + str_cth + str_cdoc


    ##############################################
    ####        PARSING
    ##############################################
    def parseHeader(self, string, sep=None):
        default_queries_fields = self.getListFields("queries", {})
        if sep is None:
            seps = ["\t", ",", ";"]
        else:
            seps = [sep]
        for ss in seps:
            fields = [self.getPrimitiveNameForLbl(s.strip()) for s in string.strip().split(ss)]
            if all([any([re.match(h, f) is not None for f in fields]) for h in default_queries_fields]):
                return fields, ss
        return None, None

    def parseQueries(self, string, list_fields=None, sep="\t", names=[None, None]):
        if type(string) is str:
            string = codecs.decode(string, 'utf-8','replace')

        if list_fields is None:
            list_fields = self.getListFields("basic", {})
        default_queries_fields = self.getListFields("queries", {})
        poplist_fields = list(list_fields) ### to pop out the query fields...
        map_fields = dict([(v,k) for (k,v) in enumerate(list_fields)])

        queries = [None, None, None]
        lpartsList = {}
        parts = [s.strip() for s in string.rsplit(sep)]
        for side, fldu in enumerate(default_queries_fields):
            flds = [(fldu, False)]
            if names[side] is not None:
                flds.append((fldu, True))
            while len(flds) > 0:
                (fld, named) = flds.pop(0)
                poplist_fields[map_fields[fld]] = None
                if map_fields[fld] >= len(parts):
                    raise Warning("Did not find expected query field for side %d (field %s expected at %d, found only %d fields)!" % (side, fld, map_fields[fld], len(parts)))
                else:
                    try:
                        if named:
                            query = self.Qclass.parse(parts[map_fields[fld]], names[side])
                        else:
                            query = self.Qclass.parse(parts[map_fields[fld]])
                    except Exception as e:
                        query = None
                        if len(flds) == 0:
                            raise Warning("Failed parsing query for side %d (field %s, string %s)!\n\t%s" % (side, fld, parts[map_fields[fld]], e))
                    if query is not None:
                        queries[side] = query
                        flds = []


        for pi, fk in enumerate(poplist_fields):
            if fk is not None and pi < len(parts):
                vs = parts[pi]
                if vs == '-':
                    continue
                fmt = self.getFieldFmt(fk)
                if fmt.get("supp_set", False):
                    sep = fmt.get("sep", ",")
                    vstr = vs.strip().split(sep)
                    try:
                        vs = map(int, vstr)
                    except ValueError:
                        vs = vstr
                if re.search("f$", fmt.get("fmt", "")):
                    vs = float(vs)
                elif re.search("d$", fmt.get("fmt", "")):
                    vs = int(vs)
                lpartsList[fk] = vs

        for side in [0, 1]:
            if queries[side] is None:
                queries[side] =  self.Qclass()
        if queries[-1] is not None:
            lpartsList["queryCOND"] = queries[-1]
        return (queries[0], queries[1], lpartsList)

    def parseRedList(self, fp, data, reds=None):
        if data is not None and data.hasNames():
            names = data.getNames()
        else:
            names = [None, None]
        
        list_fields = None
        sep = None
        more = []
        if reds is None:
            reds = []
        lid = 0
        for line in fp:
            lid += 1
            if len(line.strip()) > 0 and not re.match("^[ \t]*#", line):
                if list_fields is None:
                    list_fields, sep = self.parseHeader(line)
                else:
                    (queryL, queryR, lpartsList) = self.parseQueries(line, list_fields, sep, names)
                    r = self.Rclass.initParsed(queryL, queryR, lpartsList, data)
                    if r is not None:
                        reds.append(r)
                        more.append(line)
        return reds, {"fields": list_fields, "sep": sep, "lines": more}

    
if __name__ == '__main__':
    # print Redescription.exp_details.keys()
    from classBatch import Batch
    from classData import Data
    from classRedescription import Redescription
    # from classRedescription import parseRedList
    import sys

    rep = "/home/egalbrun/TKTL/misc/ecometrics/compare_more/v3/"
    # rep = "/home/egalbrun/short/raja_small/"
    # data = Data([rep+"data_LHS.csv", rep+"data_RHS.csv", {}, ""], "csv")
    # filename = rep+"redescriptions.csv"
    data = Data([rep+"data/IUCN_EU_nbspc3+_focus_agg.csv", rep+"data/IUCN_EU_nbspc3+_focus_iav.csv", {}, ""], "csv")
    filename = rep+"xps/biotraits-iav_IUCN_EU_focus_nbspc3+_i.01o.3m.1.queries"

    rp = Redescription.getRP()
    filep = open(filename, mode='r')
    redsO = Batch([])    
    rp.parseRedList(filep, data, redsO)
    filep.close()
    print redsO
    
    exit()

    
    lfname = "basic"
    modifiers = {"wsplits": True}
    tt = rp.getListFields(lfname, modifiers)
    print rp.getExpDict(tt)
    exit()
    
    filep = open(filename, mode='r')
    reds = Batch([])
    rp.parseRedList(filep, data, reds)
    filep.close()
    rid = 4
    # print reds[rid].disp(with_fname=True, style="U", modifiers={"wmissing": True})
    # print reds[rid].disp(with_fname=True, style="gui")
    # print reds[rid].disp(with_fname=True, style="T")
    # print reds[rid].disp(with_fname=True, style="txt")
    # print reds[rid].disp(with_fname=True, style="tex")
    # # th = rp.dispHeader("basic", style="txt")
    # # ff = rp.parseHeader(th)
    # # print ff

    # # print rp.printTexRedList(reds, names=data.getNames())
    print rp.printRedList(reds, names=data.getNames(), style="T")
    
    # print rp.dispHeader("basic", style="tex")
    # print rp.disp(reds[1], names=data.getNames(), style="tex", with_fname=True)

    # rid = 2
    # details = {"names": data.getNames(), "row_names": data.getRNames(), "named": True, "style": "U"}
    # exps = dict(enumerate(["all:acc:", "all:len:I", ":set:Exo", "test:ratio:Exx", "LHS:query:", "RHS:query:", "LHS:Lnb:", "LHS:len:", "RHS:Tset:", "RHS:depth:", "RHS:containsOR:"]))
    # # for ei, exp in exps.items():
    # #     ll = RedProps.getPrimitiveLbl(exp)
    # #     nn = RedProps.getPrimitiveNameForLbl(ll)
    # #     print "PRIM", exp, ll, nn
    
    # print RedProps.getEVals(reds[rid], exps, fresh=False, details=details)


        # print "-----------------"
        # print RedProps.getPrimitiveLbl(prim, style="txt"), RedProps.getPrimitiveLbl(prim, style="tex"), RedProps.getPrimitiveLbl(prim, style="gui")
        # print RedProps.getPrimitiveVal(reds[rid], prim), RedProps.getPrimitiveValFormatted(reds[rid], prim, details=details, style="txt")


