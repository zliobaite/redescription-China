import re, string, numpy, codecs
from classQuery import  *
from classSParts import  SParts, tool_pValSupp, tool_pValOver, tool_ratio
from classBatch import Batch
from classRedProps import RedProps
import toolRead
import pdb

ACTIVE_RSET_ID = "S0"
HAND_SIDE = {"LHS": 0, "RHS": 1, "0": 0, "1": 1, "COND": -1, "-1": -1}

def mapSuppNames(supp, details={}):
    if details.get("named", False) and "row_names" in details:
        return [details["row_names"][t] for t in supp]
    return supp


class Redescription(object):
    diff_score = Query.diff_length + 1
    
    ### PROPS WHAT
    info_what_dets = {}
    # info_what_dets = {"queryLHS": "self.prepareQueryLHS",
    #                   "queryRHS": "self.prepareQueryRHS",
    #                   "queryCOND": "self.prepareQueryCOND"}
    info_what = {"track": "self.getTrack()", "status_enabled": "self.getStatus()"}
    Pwhat_match = "("+ "|".join(info_what.keys()+info_what_dets.keys()) +")"

    ### PROPS WHICH
    which_rids = "rids"
    Pwhich_match = "("+ "|".join([which_rids]) +")"    
    @classmethod
    def hasPropWhich(tcl, which):
        return re.match(tcl.Pwhich_match, which) is not None


    RP = None
    @classmethod
    def setupRP(tcl, fields_fns=None):
        elems_typs = [("q", Query), ("s", SParts), ("r", Redescription)]
        RedProps.setupRProps(Query, Redescription, elems_typs)
        tcl.RP = RedProps(fields_fns)

    @classmethod
    def getRP(tcl, rp=None):
        if rp is None:
            if tcl.RP is None:
                tcl.setupRP()
            return tcl.RP
        return rp
    
    def __init__(self, nqueryL=None, nqueryR=None, nsupps = None, nN = -1, nPrs = [-1,-1], ssetts=None):
        self.resetRestrictedSuppSets()
        self.queries = [nqueryL, nqueryR]
        if nsupps is not None:
            self.sParts = SParts(ssetts, nN, nsupps, nPrs)
            self.dict_supp_info = None
        else:
            self.sParts = None
            self.dict_supp_info = {}
        self.lAvailableCols = [None, None]
        self.vectorABCD = None
        self.status = 1
        self.track = []
        self.cache_evals = {}
        self.condition = None

    def fromInitialPair(initialPair, data, dt={}):
        if initialPair[0] is None and initialPair[1] is None:
            return None
        supps_miss = [set(), set(), set(), set()]
        queries = [None, None]
        for side in [0,1]:
            suppS, missS = (set(), set())
            if type(initialPair[side]) is Query:
                queries[side] = initialPair[side]
                suppS, missS = initialPair[side].recompute(side, data)
            else:
                queries[side] = Query()
                if type(initialPair[side]) is Literal:
                    queries[side].extend(None, initialPair[side])                
                    suppS, missS = data.literalSuppMiss(side, initialPair[side])
            supps_miss[side] = suppS
            supps_miss[side+2] = missS

        r = Redescription(queries[0], queries[1], supps_miss, data.nbRows(), [len(supps_miss[0])/float(data.nbRows()),len(supps_miss[1])/float(data.nbRows())], data.getSSetts())
        r.track = [(-1, -1)]

        if dt.get("litC") is not None:
            litC = dt["litC"]
            if type(litC) is list:
                # if len(litC) > 1: pdb.set_trace()
                qC = Query(OR=False, buk=litC) 
            else:                
                qC = Query(buk=[litC])                
            supp_cond, miss_cond = qC.recompute(-1, data)
            r.setCondition(qC, supp_cond)
        if data.hasLT():
            r.setRestrictedSupp(data)
        return r
    fromInitialPair = staticmethod(fromInitialPair)

    def fromQueriesPair(queries, data):
        r = Redescription(queries[0].copy(), queries[1].copy())
        r.recompute(data)        
        r.track = [tuple([0] + sorted(r.queries[0].invCols())), tuple([1] + sorted(r.queries[1].invCols()))]
        if len(queries) > 2 and queries[2] is not None:
            qC = queries[2]
            supp_cond, miss_cond = qC.recompute(-1, data)
            r.setCondition(qC, supp_cond)            
        return r
    fromQueriesPair = staticmethod(fromQueriesPair)

    def dropSupport(self):
        if self.sParts is not None:
            self.dict_supp_info.toDict()
            self.sParts = None

    def compare(self, y):
        if self.score() > y.score():
            return Redescription.diff_score
        elif self.score() == y.score():
            return Query.comparePair(self.queries[0], self.queries[1], y.queries[0], y.queries[1])
        else:
            return -Redescription.diff_score

    def interArea(self, redB, side):
        if redB is not None:
            return len(redB.supp(side) & self.supp(side))* len(redB.invColsSide(side) & self.invColsSide(side))
        return 0
    def unionArea(self, redB, side):
        if redB is not None:
            return len(redB.supp(side) | self.supp(side))* len(redB.invColsSide(side) | self.invColsSide(side))
        return 0
    def overlapAreaSide(self, redB, side):
        if len(redB.invColsSide(side) & self.invColsSide(side)) == 0:
            return 0
        areaU = self.unionArea(redB, side)
        return tool_ratio(self.interArea(redB, side), areaU)
    def overlapAreaTotal(self, redB):
        areaUL = self.unionArea(redB, 0)
        areaUR = self.unionArea(redB, 1)
        return tool_ratio(self.interArea(redB, 0) + self.interArea(redB, 1),areaUL+areaUR)
    def overlapAreaL(self, redB):
        return self.overlapAreaSide(redB, 0)
    def overlapAreaR(self, redB):
        return self.overlapAreaSide(redB, 1)
    def overlapAreaMax(self, redB):
        return max(self.overlapAreaSide(redB, 0), self.overlapAreaSide(redB, 1))

    def overlapRows(self, redB):
        if redB is not None:
            return tool_ratio(len(redB.getSuppI() & self.getSuppI()), min(redB.getLenI(), self.getLenI()))
        return 0
    
    def oneSideIdentical(self, redescription):
        return self.queries[0] == redescription.queries[0] or self.queries[1] == redescription.queries[1]
    def bothSidesIdentical(self, redescription):
        return self.queries[0] == redescription.queries[0] and self.queries[1] == redescription.queries[1]

    def equivalent(self, y):
       return abs(self.compare(y)) < Query.diff_balance
        
    # def __hash__(self):
    #      return int(hash(self.queries[0])+ hash(self.queries[1])*100*self.score())
        
    def __len__(self):
        return self.length(0) + self.length(1)

    def usesOr(self, side=None):
        if side is not None:
            return self.queries[side].usesOr()
        return self.queries[0].usesOr() or self.queries[1].usesOr()

    def supp(self, side):
        return self.supports().supp(side)

    def miss(self, side):
        return self.supports().miss(side)
            
    def score(self):
        return self.getAcc()

    def supports(self):
        return self.sParts
    
    def partsAll(self):
        return self.supports().sParts

    def partsFour(self):
        return [self.supports().suppA(), self.supports().suppB(), self.supports().suppI(), self.supports().suppO()]

    def partsThree(self):
        return [self.supports().suppA(), self.supports().suppB(), self.supports().suppI()]
    
    def partsNoMiss(self):
        return self.supports().sParts[:4]
    
    def query(self, side=None):
        if side == -1:
            return self.getQueryC()
        return self.queries[side]

    def getQueries(self):
        return self.queries

    def getQueryC(self):
        if self.condition is not None:
            return self.condition.get("q")
    def getSupportsC(self):
        if self.condition is not None:
            return self.condition.get("sparts")
    def getSuppC(self):
        if self.condition is not None:
            return self.condition.get("supp")

    def hasCondition(self):
        return self.condition is not None
    def setCondition(self, qC=None, supp_cond=None): ### here
        self.condition = None
        if qC is not None:
            if supp_cond is None:
                sparts = None
            else:
                sparts = self.supports().copy()
                sparts.update(0, False, supp_cond)
                sparts.update(1, False, supp_cond)
            self.condition = {"q": qC, "supp": supp_cond, "sparts": sparts}
    
    def probas(self):
        return self.supports().probas()

    def probasME(self, dbPrs, epsilon=0):
        return [self.queries[side].probaME(dbPrs, side, epsilon) for side in [0,1]]

    def surpriseME(self, dbPrs, epsilon=0):
        #return [-numpy.sum(numpy.log(numpy.absolute(SParts.suppVect(self.supports().nbRows(), self.supports().suppSide(side), 0) - self.queries[side].probaME(dbPrs, side)))) for side in [0,1]]
        return -numpy.sum(numpy.log(numpy.absolute(SParts.suppVect(self.supports().nbRows(), self.supports().suppI(), 0) - self.queries[0].probaME(dbPrs, 0)*self.queries[1].probaME(dbPrs, 1))))

    def exME(self, dbPrs, epsilon=0):
        prs = [self.queries[side].probaME(dbPrs, side, epsilon) for side in [0,1]]
        surprises = []
        tmp = [i for i in self.supports().suppI() if prs[0][i]*prs[1][i] == 0]
        surprises.append(-numpy.sum(numpy.log([prs[0][i]*prs[1][i] for i in self.supports().suppI()])))
        surprises.extend([-numpy.sum(numpy.log([prs[side][i] for i in self.supports().suppSide(side)])) for side in [0,1]])

        return surprises + [len(tmp) > 0]

        N = self.supports().nbRows()
        margsPr = [numpy.sum([prs[side][i] for i in self.supports().suppSide(side)]) for side in [0,1]]
        pvals = [tool_pValOver(self.supports().lenI(), N, int(margsPr[0]), int(margsPr[1])), tool_pValSupp(N, self.supports().lenI(), margsPr[0]*margsPr[1]/N**2)]
        return surprises, pvals
    
    def length(self, side):
        return len(self.queries[side])
        
    def availableColsSide(self, side, data=None, single_dataset=False):
        if self.lAvailableCols[side] is not None and self.length(1-side) != 0:
            tt = set(self.lAvailableCols[side])
            if single_dataset:
                tt &= set(self.lAvailableCols[1-side])
            if data is not None:
                for ss in [0,1]:
                    if data.hasGroups(ss):
                        for c in self.queries[ss].invCols():
                            tt = [t for t in tt if data.areGroupCompat(t, c, side, ss)]
            return tt
        return set() 
    def nbAvailableCols(self):
        if self.lAvailableCols[0] is not None and self.lAvailableCols[1] is not None:
            return len(self.lAvailableCols[0]) + len(self.lAvailableCols[1])
        return -1
    def updateAvailable(self, souvenirs):
        nb_extensions = 0
        for side in [0, 1]:
            if self.lAvailableCols[side] is None or len(self.lAvailableCols[side]) != 0:
                self.lAvailableCols[side] =  souvenirs.availableMo[side] - souvenirs.extOneStep(self, side)
                nb_extensions += len(souvenirs.availableMo[side]) - self.length(side) - len(self.lAvailableCols[side])
        return nb_extensions
    def removeAvailables(self):
        self.lAvailableCols = [set(),set()]
            
    def update(self, data=None, side= -1, opBool = None, literal= None, suppX=None, missX=None):
        if side == -1 :
            self.removeAvailables()
        else:
            op = Op(opBool)
            self.queries[side].extend(op, literal)
            self.supports().update(side, op.isOr(), suppX, missX)
            self.dict_supp_info = None
            if self.lAvailableCols[side] is not None:
                self.lAvailableCols[side].remove(literal.colId())
            self.track.append(((1-side) * 1-2*int(op.isOr()), literal.colId()))

    def setFull(self, max_var=None):
        if max_var is not None:
            for side in [0,1]:
                if self.length(side) >= max_var[side]:
                    self.lAvailableCols[side] = set()
                
    def kid(self, data, side= -1, op = None, literal= None, suppX= None, missX=None):
        kid = self.copy()        
        kid.update(data, side, op, literal, suppX, missX)
        return kid
            
    def copy(self):
        r = Redescription(self.queries[0].copy(), self.queries[1].copy(), \
                             self.supports().supparts(), self.supports().nbRows(), self.probas(), self.supports().getSSetts())
        for side in [0,1]:
            if self.lAvailableCols[side] is not None:
                r.lAvailableCols[side] = set(self.lAvailableCols[side])
        r.status = self.status
        r.track = list(self.track)
        r.restricted_sets = {}
        for sid, rst in self.restricted_sets.items():
            r.restricted_sets[sid] = {"sParts": rst["sParts"],
                                         "prs": [rst["prs"][0], rst["prs"][1]],
                                         "rids": set(rst["rids"])}

        return r

    def recomputeQuery(self, side, data= None, restrict=None):
        return self.queries[side].recompute(side, data, restrict)
    
    def invLiteralsSide(self, side):
        return self.queries[side].invLiterals()

    def invLiterals(self):
        return [self.invLiteralsSide(0), self.invLiteralsSide(1)]
    
    def invColsSide(self, side):
        return self.queries[side].invCols()

    def invCols(self):
        return [self.invColsSide(0), self.invColsSide(1)]
        
    def setRestrictedSupp(self, data):
        ### USED TO BE STORED IN: self.restrict_sub, self.restricted_sParts, self.restricted_prs = None, None, None
        self.setRestrictedSuppSets(data, supp_sets=None)
        
    def resetRestrictedSuppSets(self):
        self.restricted_sets = {}

    def setRestrictedSuppSets(self, data, supp_sets=None):
        self.dict_supp_info = None
        if supp_sets is None:
            if data.hasLT():
                supp_sets = data.getLT()
            else:
                supp_sets = {ACTIVE_RSET_ID: data.nonselectedRows()}
        for sid, sset in supp_sets.items():
            if len(sset) == 0:
                self.restricted_sets[sid] = {"sParts": None,
                                             "prs": None,
                                             "rids": set()}
            elif sid not in self.restricted_sets or self.restricted_sets[sid]["rids"] != sset:
                (nsuppL, missL) = self.recomputeQuery(0, data, sset)
                (nsuppR, missR) = self.recomputeQuery(1, data, sset)
                if len(missL) + len(missR) > 0:
                    rsParts = SParts(data.getSSetts(), sset, [nsuppL, nsuppR, missL, missR])
                else:
                    rsParts = SParts(data.getSSetts(), sset, [nsuppL, nsuppR])

                self.restricted_sets[sid] = {"sParts": rsParts,
                                             "prs": [self.queries[0].proba(0, data, sset),
                                                     self.queries[1].proba(1, data, sset)],
                                             "rids": set(sset)}
            
    def getNormalized(self, data=None, side=None):
        if side is not None:
            sides = [side]
        else:
            sides = [0,1]
        queries = [self.queries[side] for side in [0,1]]
        c = [False, False]
        for side in sides:
            queries[side], c[side] = self.queries[side].algNormalized()
        if c[0] or c[1]:
            red = Redescription.fromQueriesPair(queries, data)
            ### check that support is same
            # if self.supports() != red.supports():
            #     print "ERROR ! SUPPORT CHANGED WHEN NORMALIZING..."
            #     pdb.set_trace()
            return red, True            
        else:
            return self, False
        
    def recompute(self, data):
        (nsuppL, missL) = self.recomputeQuery(0, data)
        (nsuppR, missR) = self.recomputeQuery(1, data)
#        print self.disp()
#        print ' '.join(map(str, nsuppL)) + ' \t' + ' '.join(map(str, nsuppR))
        if len(missL) + len(missR) > 0:
            self.sParts = SParts(data.getSSetts(), data.nbRows(), [nsuppL, nsuppR, missL, missR])
        else:
            self.sParts = SParts(data.getSSetts(), data.nbRows(), [nsuppL, nsuppR]) #TODO: recompute
        self.prs = [self.queries[0].proba(0, data), self.queries[1].proba(1, data)]
        if data.hasLT():
            self.setRestrictedSupp(data)
        if self.hasCondition():
            qC = self.getQueryC()
            supp_cond, miss_cond = qC.recompute(-1, data)
            self.setCondition(qC, supp_cond)                        
        self.dict_supp_info = None

    def check(self, data):
        result = 0
        details = None
        if self.supports() is not None: #TODO: sparts
            (nsuppL, missL) = self.recomputeQuery(0, data)
            (nsuppR, missR) = self.recomputeQuery(1, data)
            
            details = ( len(nsuppL.symmetric_difference(self.supports().supp(0))) == 0, \
                     len(nsuppR.symmetric_difference(self.supports().supp(1))) == 0, \
                     len(missL.symmetric_difference(self.supports().miss(0))) == 0, \
                     len(missR.symmetric_difference(self.supports().miss(1))) == 0 )        
            result = 1
            for detail in details:
                result*=detail
        return (result, details)

    def hasMissing(self):
        return self.supports().hasMissing()

    def getStatus(self):
        return self.status
    def getEnabled(self):
        return self.status

    def flipEnabled(self):
        self.status = -self.status

    def setEnabled(self):
        self.status = 1
    def setDisabled(self):
        self.status = -1

    def setDiscarded(self):
        self.status = -2

    ##### GET FIELDS INFO INVOLVING ADDITIONAL DETAILS (PRIMARILY FOR SIREN)
    def getQueriesU(self, details=None):
        if details is not None and "names" in details:
            return self.queries[0].disp(details["names"][0], style="U") + "---" + self.queries[1].disp(details["names"][1], style="U")
        else:
            return self.queries[0].disp(style="U") + "---" + self.queries[1].disp(style="U")

    def getQueryLU(self, details=None):
        if details is not None and "names" in details:
            return self.queries[0].disp(details["names"][0], style="U") #, unicd=True)
        else:
            return self.queries[0].disp(style="U")

    def getQueryRU(self, details=None):
        if details is not None and "names" in details:
            return self.queries[1].disp(details["names"][1], style="U") #, unicd=True)
        else:
            return self.queries[1].disp(style="U")

    def getTrack(self, details=None):
        if details is not None and ( details.get("aim", None) == "list" or details.get("format", None) == "str"):
            return ";".join(["%s:%s" % (t[0], ",".join(map(str,t[1:]))) for t in self.track])
        else:
            return self.track

    def getSortAble(self, details=None):
        if details.get("aim") == "sort":
            return (self.status, details.get("id", "?"))
        return ""

    def getShortRid(self, details=None):
        return "R%s" % details.get("id", "?")

    def getEnabled(self, details=None):
        return 1*(self.status>0)

    def getTypeParts(self, details=None):
        return self.supports().getTypeParts()
    def getMethodPVal(self, details=None):
        return self.supports().getMethodPVal()    

    def hasRSets(self):
        return len(self.restricted_sets) > 0
    def getRSet(self, details=None):
        if type(details) is dict:
            rset_id = details.get("rset_id")
        else:
            rset_id = details
        if rset_id is not None:
            if rset_id == "all":
                return {"sParts": self.supports()}
            elif rset_id == "cond" and self.hasCondition():
                return {"sParts": self.getSupportsC()}
            elif rset_id in self.restricted_sets:
                return self.restricted_sets[rset_id]
            return None
        elif ACTIVE_RSET_ID in self.restricted_sets:
            return self.restricted_sets[ACTIVE_RSET_ID]
        else:            
            return {"sParts": self.supports()}
    def getRSetParts(self, details=None):
        rset = self.getRSet(details)
        if rset is not None:
            return rset.get("sParts")
    def getRSetIds(self, details=None):
        rset = self.getRSet(details)
        if rset is not None and "rids" in rset:
            return sorted(rset["rids"])
        return None
       
    def getRSetABCD(self, details=None):
        ssp = self.getRSetParts(details)
        if ssp is not None:
            return ssp.get("sParts").getVectorABCD(force_list=True, rest_ids=ssp.get("rids"))
        
    def getAccRatio(self, details=None):
        if details is not None and (details.get("rset_id_num") in self.restricted_sets \
               or details.get("rset_id_den") in self.restricted_sets):
            acc_num = self.getRSetParts(details.get("rset_id_num")).acc()
            acc_den = self.getRSetParts(details.get("rset_id_den")).acc()
            return tool_ratio(acc_num, acc_den)
        return 1.

    def getLenRatio(self, details=None):
        if details is not None and (details.get("rset_id_num") in self.restricted_sets \
               or details.get("rset_id_den") in self.restricted_sets):
            len_num = self.getRSetParts(details.get("rset_id_num")).lenI()
            len_den = self.getRSetParts(details.get("rset_id_den")).lenI()
            return tool_ratio(len_num, len_den)
        return 1.
    
    def getAcc(self, details=None):
        return self.getRSetParts(details).acc()
    def getPVal(self, details=None):
        return self.getRSetParts(details).pVal()

    def getLenP(self, details=None):
        if "part_id" in details:
            return self.getRSetParts(details).lenP(details["part_id"])
        return -1

    def getLenI(self, details=None):
        return self.getRSetParts(details).lenI()
    def getLenU(self, details=None):
        return self.getRSetParts(details).lenU()
    def getLenL(self, details=None):
        return self.getRSetParts(details).lenL()
    def getLenR(self, details=None):
        return self.getRSetParts(details).lenR()
    def getLenO(self, details=None):
        return self.getRSetParts(details).lenO()
    def getLenN(self, details=None):
        return self.getRSetParts(details).lenN()
    def getLenA(self, details=None):
        return self.getRSetParts(details).lenA()
    def getLenB(self, details=None):
        return self.getRSetParts(details).lenB()
    
    def getSuppI(self, details=None):
        return self.getRSetParts(details).suppI()
    def getSuppU(self, details=None):
        return self.getRSetParts(details).suppU()
    def getSuppL(self, details=None):
        return self.getRSetParts(details).suppL()
    def getSuppR(self, details=None):
        return self.getRSetParts(details).suppR()
    def getSuppO(self, details=None):
        return self.getRSetParts(details).suppO()
    def getSuppN(self, details=None):
        return self.getRSetParts(details).suppN()
    def getSuppA(self, details=None):
        return self.getRSetParts(details).suppA()
    def getSuppB(self, details=None):
        return self.getRSetParts(details).suppB()
    
    def getProp(self, what, which=None, rset_id=None, details=None):
        if Query.hasPropWhat(what) and rset_id in HAND_SIDE:
            q = self.query(HAND_SIDE[rset_id])            
            if q is not None:
                dts = {"side": HAND_SIDE[rset_id]}
                dts.update(details)
                return q.getProp(what, which, dts)
            return None

        if rset_id is not None and which == self.which_rids: ### ids details for split sets
            rset_ids = self.getRSetIds(rset_id)
            if rset_ids is None:
                return None
            if what == "len" or what == "card":
                return len(rset_ids)
            elif what == "supp" or what == "set":
                return mapSuppNames(rset_ids, details)
            elif what == "perc":
                return tool_ratio(100.*len(rset_ids), self.supports().nbRows())
            elif what == "ratio":
                return tool_ratio(len(rset_ids), self.supports().nbRows())

        if SParts.hasPropWhat(what): ### info from supp parts
            rset_parts = self.getRSetParts(rset_id)
            if rset_parts is None:
                return None
            prp = self.getRSetParts(rset_id).getProp(what, which)
            if what == "supp" or what == "set":
                return mapSuppNames(prp, details)            
            return prp
        elif what in Redescription.info_what_dets: ### other redescription info
            methode = eval(Redescription.info_what_dets[what])
            if callable(methode):
                return methode(details)
        elif what in Redescription.info_what: ### other redescription info
            return eval(Redescription.info_what[what])
         
    def setCacheEVals(self, cevs):
        self.cache_evals = cevs
    def updateCacheEVals(self, cevs):
        self.cache_evals.update(cevs)
    def setCacheEVal(self, ck, cev):
        self.cache_evals[ck] = cev
    def getCacheEVal(self, ck):
        return self.cache_evals.get(ck)
    def hasCacheEVal(self, ck):
        return ck in self.cache_evals

    def getEValGUI(self, details):
        val = None
        if "rp" in details:
            val = details["rp"].getEValGUI(self, details)
        if val is None:
            return details.get('replace_none')
        return val

    
##### PRINTING AND PARSING METHODS
    #### FROM HERE ALL PRINTING AND READING
    def __str__(self):
        str_av = ["?", "?"]
        for side in [0,1]:
            if self.availableColsSide(side) is not None:
                str_av[side] = "%d" % len(self.availableColsSide(side))
        tmp = ('%s + %s terms:' % tuple(str_av)) + ('\t (%i): %s\t%s\t%s\t%s' % (len(self), self.dispQueries(sep=" "), self.dispStats(sep=" "), self.dispSuppL(sep=" "), self.getTrack({"format":"str"})))
        return tmp

    def dispQueries(self, names=[None, None, None], sep='\t'):
        sides = [0, 1]
        if self.hasCondition():
            sides.append(-1)
        return sep.join(["q%s:%s" % (side, self.query(side).disp(names[side])) for side in sides])
    def dispStats(self, sep="\t"):
        if self.hasCondition():
            return self.supports().dispStats(sep)+sep+"COND:"+self.getSupportsC().dispStats(sep)
        return self.supports().dispStats(sep)
    def dispSupp(self, sep="\t"):
        if self.hasCondition():
            return self.supports().dispSupp(sep)+sep+"COND:"+self.getSupportsC().dispSupp(sep)
        return self.supports().dispSupp(sep)        
    def dispSuppL(self, sep="\t"):
        if self.hasCondition():
            return self.supports().dispSuppL(sep)+sep+"COND:"+self.getSupportsC().dispSuppL(sep)
        return self.supports().dispSuppL(sep)
    
    def prepareQuery(self, side, details={}, named=False):
        style=details.get("style", "")
        q = self.query(side)
        if q is None: return ""
        if (named or details.get("named", False)) and "names" in details:
            return q.disp(names=details["names"][side], style=style)
        return q.disp(style=style)        
    def prepareQueryRHS(self, details):
        return self.prepareQuery(1, details)
    def prepareQueryLHS(self, details):
        return self.prepareQuery(0, details)
    def prepareQueryCOND(self, details):
        return self.prepareQuery(-1, details)    

    def disp(self, names=[None, None], row_names=None, with_fname=False, rid="", nblines=1, delim="", last_one=False, list_fields="basic", modifiers={}, style="txt", sep=None, rp=None):
        return self.getRP(rp).disp(self, names, row_names, with_fname, rid, nblines, delim, last_one, list_fields, modifiers, style, sep)
    @classmethod   
    def parse(tcl, stringQueries, stringSupp=None, data=None, rp=None):
        if data is not None and data.hasNames():
            names = data.getNames()
        else:
            names = [None, None]
        (queryL, queryR, lpartsList) = tcl.getRP(rp).parseQueries(stringQueries, names=names)
        supportsS = None
        if data is not None and stringSupp is not None and type(stringSupp) == str and re.search('\t', stringSupp):
            supportsS = SParts.parseSupport(stringSupp, data.nbRows(), data.getSSetts())
        return tcl.initParsed(queryL, queryR, lpartsList, data, supportsS)
    @classmethod   
    def initParsed(tcl, queryL, queryR, lpartsList={}, data = None, supportsS=None):
        status_enabled = None
        if "status_enabled" in lpartsList:
            status_enabled = int(lpartsList.pop("status_enabled"))
        r = None
        if supportsS is not None:
            r = tcl(queryL, queryR, supportsS.supparts(), data.nbRows(), [set(),set()], [ queryL.proba(0, data), queryR.proba(1, data)], data.getSSetts())

            for key, v in lpartsList.items():
                tv = RedProps.getEVal(r, key)
                if tv != v:
                    raise Warning("Something wrong in the supports ! (%s: %s ~ %s)\n" \
                                    % (key, v, tv))

        if r is None:
            r = tcl(queryL, queryR)
            if data is not None:
                r.recompute(data)
            else:
                r.cache_evals = lpartsList

        if "queryCOND" in lpartsList:
            qC = lpartsList["queryCOND"]
            supp_cond = None
            if data is not None:
                supp_cond, miss_cond = qC.recompute(-1, data)
            r.setCondition(qC, supp_cond)
        if r is not None and status_enabled is not None:
            r.status = status_enabled
        return r
    
        
if __name__ == '__main__':
    # print Redescription.exp_details.keys()
    from classData import Data
    from classQuery import Query
    import sys

    rep = "/home/egalbrun/short/vaalikone_FILES/"
    rep = "/home/egalbrun/short/raja_time/"
    data = Data([rep+"data_LHS_lngtid.csv", rep+"data_RHS_lngtid.csv", {}, "NA"], "csv")

    rp = Redescription.getRP()
    # filename = rep+"redescriptions.csv"
    filename = rep+"tid_test.queries"
    filep = open(filename, mode='r')
    reds = Batch([])
    rp.parseRedList(filep, data, reds)
    with open("/home/egalbrun/short/tmp_queries.txt", mode='w') as fp:
        # pdb.set_trace()
        fp.write(rp.printRedList(reds, missing=True))
        ## fp.write(rp.printTexRedList(reds, names = [data.getNames(0), data.getNames(1)], nblines=3, standalone=True))
        ## fields=[-1, "CUST:XX=q0:containsC:0", "Lnb_queryLHS", "Lset_queryLHS", "Lnb_queryRHS", "Lset_queryRHS", "containsAND_queryRHS", "containsOR_queryRHS"]
    exit()
    # rep = "/home/galbrun/TKTL/redescriptors/data/vaalikone/"
    # data = Data([rep+"vaalikone_profiles_test.csv", rep+"vaalikone_questions_test.csv", {}, "NA"], "csv")

    # reds = []
    # with codecs.open("../../bazar/queries.txt", encoding='utf-8', mode='r') as f:
    #     for line in f:
    #         if len(line.strip().split("\t")) >= 2:
    #             try:
    #                 tmpLHS = Query.parse(line.strip().split("\t")[0], data.getNames(0))
    #                 tmpRHS = Query.parse(line.strip().split("\t")[1], data.getNames(1))
    #             except:
    #                 continue
    #             r = Redescription.fromQueriesPair([tmpLHS, tmpRHS], data)
    #             reds.append(r)

    # with codecs.open("../../bazar/queries_list2.txt", encoding='utf-8', mode='w') as f:
    #     f.write(rp.printRedList(reds))

    # with codecs.open("../../bazar/queries_list2.txt", encoding='utf-8', mode='r') as f:
    #     reds, _ = rp.parseRedList(f, data)

    # for red in reds:
    #     print red.disp()
