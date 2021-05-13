from classQuery import Op, Query
from classRedescription import Redescription
import pdb

class ExtensionError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Extension(object):

    def __init__(self, ssetts, adv=None, clp=None, sol=None):
        ### self.adv is a tuple: acc, varBlue, varRed, contrib, fixBlue, fixRed
        self.ssetts = ssetts
        self.condition = None
        if adv is not None and len(adv) == 3 and len(adv[2]) == 4 and clp==None and sol==None:
            self.adv = adv[0]
            self.clp = adv[1]
            self.side, self.op, tmp_neg, self.literal = adv[2]
        else:
            self.adv = adv
            self.clp = clp
            if sol is not None:
                self.side, self.op, tmp_neg, self.literal = sol
            else:
                self.side, self.op, self.literal = None, None, None

    def setClp(self, clp, neg=False):
        if clp is None:
            self.clp = clp
        else:
            if len(clp) == 2:
                lin = clp[0]
                lparts = clp[1]
                lout = [lparts[i] - lin[i] for i in range(len(lparts))]
                if neg:
                    self.clp = [lout, lin, lparts]
                else:
                    self.clp = [lin, lout, lparts]

    def setCondition(self, condition):
        self.condition = condition

    def getCondition(self):
        return self.condition
    def hasCondition(self):
        return self.condition is not None

    def dispLiteral(self):
        lit = self.getLiteral()
        if lit is None:
            return "-"
        elif type(lit) is list:
            return " AND ".join(["%s" % ll for ll in lit])
        return "%s" % lit

    def getLiteral(self):
        if self.isValid():
            return self.literal
    def setLiteral(self, literal):
        self.literal = literal
        
    def getOp(self):
        if self.isValid():
            return self.op

    def getSide(self):
        if self.isValid():
            return self.side

    def getAcc(self):
        if self.isValid():
            return self.adv[0]

    def getVarBlue(self):
        if self.isValid():
            return self.adv[1]

    def getVarRed(self):
        if self.isValid():
            return self.adv[2]

    def getCLP(self):
        if self.isValid():
            return self.clp
    
    def kid(self, red, data):
        supp = data.supp(self.getSide(), self.getLiteral())
        miss = data.miss(self.getSide(), self.getLiteral())
        tmp = red.kid(data, self.getSide(), self.getOp(), self.getLiteral(), supp, miss)
        if self.hasCondition():
            litC = self.getCondition().getLiteral()
            if type(litC) is list:
                # if len(litC) > 1: pdb.set_trace()
                qC = Query(OR=False, buk=litC) 
            else:
                qC = Query(buk=[litC])
            supp_cond, miss_cond = qC.recompute(-1, data)
            tmp.setCondition(qC, supp_cond)
        return tmp

    def isValid(self):
        return self.adv is not None and len(self.adv) > 2

    def isNeg(self):
        if self.isValid():
            if type(self.getLiteral()) is list:
                return False
            return self.getLiteral().isNeg()

    def __str__(self):
        tmp = "Empty extension"
        if self.isValid():
            tmp = ("Extension:\t (%d, %s, %s) -> %f CLP:%s ADV:%s" % (self.getSide(), Op(self.getOp()), self.dispLiteral(), self.getAcc(), str(self.clp), str(self.adv)))
            if self.hasCondition():
                tmp += "\nCOND_%s" % self.getCondition()
        return tmp 

    def disp(self, base_acc=None, N=0, prs=None, coeffs=None):
        strPieces = ["", "", "", ""]
        score_of = self
        if self.isValid():
            strPieces[self.getSide()] = "%s %s" % (Op(self.getOp()), self.dispLiteral()) 
            if self.hasCondition():
                # score_of = self.getCondition()
                strPieces[2] = self.getCondition().dispLiteral()

            if base_acc is None:
                strPieces[-1] = '----\t%1.7f\t----\t----\t% 5i\t% 5i' \
                                % (score_of.getAcc(), score_of.getVarBlue(), score_of.getVarRed())
            else:
                strPieces[-1] = '\t\t%+1.7f \t%1.7f \t%1.7f \t%1.7f\t% 5i\t% 5i' \
                                % (score_of.score(base_acc, N, prs, coeffs), score_of.getAcc(), \
                                   score_of.pValQuery(N, prs), score_of.pValRed(N, prs) , score_of.getVarBlue(), score_of.getVarRed())

        return '* %20s <==> * %20s %20s %s' % tuple(strPieces) # + "\n\tCLP:%s" % str(self.clp)
            
    def score(self, base_acc, N, prs, coeffs):
        ### HERE: HOW TO SCORE WITH CONDITION?
        if self.isValid():
            sc = coeffs["impacc"]*self.impacc(base_acc) \
                   + coeffs["rel_impacc"]*self.relImpacc(base_acc) \
                   + self.pValRedScore(N, prs, coeffs) \
                   + self.pValQueryScore(N, prs, coeffs)
            if False: #self.hasCondition():
                sc += self.getCondition().score(base_acc, N, prs, coeffs)
            return sc

                   
    def relImpacc(self, base_acc=0):
        if self.isValid():
            if base_acc != 0:
                return (self.adv[0] - base_acc)/base_acc
            else:
                return self.adv[0]
        
    def impacc(self, base_acc=0):
        if self.isValid():
            return (self.adv[0] - base_acc)
        
    def pValQueryScore(self, N, prs, coeffs=None):
        if self.isValid():
            if coeffs is None or coeffs["pval_query"] < 0:
                return coeffs["pval_query"] * self.ssetts.pValQuery(N, prs)
            elif coeffs["pval_query"] > 0:
                return -coeffs["pval_fact"]*(coeffs["pval_query"] < self.pValQuery(N, prs))
            else:
                return 0

    def pValRedScore(self, N, prs, coeffs=None):
        if self.isValid():
            if coeffs is None or coeffs["pval_red"] < 0:
                return coeffs["pval_red"] * self.pValRed(N, prs)
            elif coeffs["pval_red"] > 0:
                return -coeffs["pval_fact"]*(coeffs["pval_red"] < self.pValRed(N, prs))
            else:
                return 0

    def pValQuery(self, N=0, prs=None):
        if self.isValid():
            return self.ssetts.pValQueryCand(self.side, self.op, self.isNeg(), self.clp, N, prs)

    def pValRed(self, N=0, prs=None):
        if self.isValid():
            try:
                return self.ssetts.pValRedCand(self.side, self.op, self.isNeg(), self.clp, N, prs)
            except IndexError:                
                pdb.set_trace()
                return self.ssetts.pValRedCand(self.side, self.op, self.isNeg(), self.clp, N, prs)
                print self

    def __cmp__(self, other):
        return self.compare(other)
    
    def compare(self, other):
        tmp = self.compareAdv(other)
        if tmp == 0:
            return cmp(self.getLiteral(), other.getLiteral())

    def compareAdv(self, other):
        if other is None:
            return 1
        if not self.isValid():
            return -1

        if type(other) in [tuple, list]:
            other_adv = other
        else:
            other_adv = other.adv

        return cmp(self.adv, other_adv)
        
    
class ExtensionsBatch(object):
    def __init__(self, N=0, coeffs=None, current=None):
        self.current = current
        self.base_acc = self.current.getAcc()
        self.N = N
        self.prs = self.current.probas()
        self.coeffs = coeffs
        self.bests = {}
        self.tmpsco = {}
    def scoreCand(self, cand):
        if cand is not None:
            return cand.score(self.base_acc, self.N, self.prs, self.coeffs)

    def pos(self, cand):
        if cand.isValid():
            return (cand.getSide(), cand.getOp())

    def get(self, pos):
        if pos in self.bests:
            return self.bests[pos]
        else:
            return None
        
    def update(self, cands):
        for cand in cands:
            pos = self.pos(cand)
            self.scoreCand(cand)
            if pos is not None and ( pos not in self.bests or self.scoreCand(cand) > self.scoreCand(self.bests[pos])):
                self.bests[pos] = cand

    
    def updateDL(self, cands, rm, data):
        for cand in cands:
            kid = cand.kid(self.current, data)
            top = rm.getTopDeltaRed(kid, data)
            pos = self.pos(cand)
            self.scoreCand(cand)
            if pos is not None and ( pos not in self.bests or -top[0] > self.tmpsco[pos]):
                self.bests[pos] = cand
                self.tmpsco[pos] = -top[0]

    def improving(self, min_impr=0):
        return dict([(pos, cand)  for (pos, cand) in self.bests.items() \
                     if self.scoreCand(cand) >= min_impr])

    def improvingKids(self, data, min_impr=0, max_var=[-1,-1]):
        kids = []
        for (pos, cand) in self.bests.items():
            if self.scoreCand(cand) >= min_impr:
                kid = cand.kid(self.current, data)
                kid.setFull(max_var)
                if kid.getAcc() != cand.getAcc():
                    raise ExtensionError("[in Extension.improvingKids]\n%s\n\t%s\n\t~> %s" % (self.current, cand, kid))
                if kid.hasCondition() and kid.getAcc("cond") != cand.getCondition().getAcc():
                    raise ExtensionError("[in Extension.improvingKids COND]\n%s\n\t%s\n\t~> %s" % (self.current, cand, kid))
                kids.append(kid)
        return kids
    

    def improvingKidsDL(self, data, min_impr=0, max_var=[-1,-1], rm=None):
        tc = rm.getTopDeltaRed(self.current, data)
        min_impr = -tc[0]
        # print "DL impr---", min_impr, self.tmpsco
        kids = []
        for (pos, cand) in self.bests.items():
            if self.tmpsco[pos] >= min_impr:
                kid = cand.kid(self.current, data)
                kid.setFull(max_var)
                if kid.getAcc() != cand.getAcc():
                    raise ExtensionError("[in Extension.improvingKidsDL]\n%s\n\t%s\n\t~> %s" % (self.current, cand, kid))
            
                kids.append(kid)
        return kids

        
    def __str__(self):
        dsp  = 'Extensions Batch:\n' 
        dsp += 'Redescription: %s' % self.current
        dsp += '\n\t  %20s        %20s        %20s' \
                  % ('LHS extension', 'RHS extension', 'Condition')
            
        dsp += '\t\t%10s \t%9s \t%9s \t%9s\t% 5s\t% 5s' \
                      % ('score', 'Accuracy',  'Query pV','Red pV', 'toBlue', 'toRed')
            
        for k,cand in self.bests.items(): ## Do not print the last: current redescription
            dsp += '\n\t%s' % cand.disp(self.base_acc, self.N, self.prs, self.coeffs)
        return dsp
