from classData import CatColM
from classConstraints import Constraints
from classCharbon import CharbonGreedy
from classExtension import Extension
from classSParts import SParts, tool_ratio
from classQuery import  *
import numpy
import pdb

class CharbonGMiss(CharbonGreedy):

    name = "GreedyMiss"
    def handlesMiss(self):
        return True

    def getCandidates(self, side, col, supports, red, colsC=None):
        currentRStatus = Constraints.getStatusRed(red, side)
        method_string = 'self.getCandidates%i' % col.typeId()
        try:
            method_compute =  eval(method_string)
        except AttributeError:
              raise Exception('Oups No candidates method for this type of data (%i)!'  % col.typeId())
        cands = method_compute(side, col, supports, currentRStatus)
        # print "=======MISS=========="
        for cand in cands:            
            supp = col.suppLiteral(cand.getLiteral())
            lparts = supports.lparts()
            lin = supports.lpartsInterX(supp)
            # print cand.getLiteral(), cand.isNeg()
            # print lparts, lin, len(supp)
            # print "Before", cand
            # # print "CLP \tB: %s\tA:%s" % (cand.clp, [lin, lparts])
            cand.setClp([lin, lparts], cand.isNeg())
            if colsC is not None and self.constraints.getCstr("add_condition"):
                ss = supports.copy()
                ss.update(side, cand.getOp(), supp)
                cand.setCondition(self.getCondition(colsC, ss))
            # print cand
        # print "================="            
        return cands
    
    def getCandidates1(self, side, col, supports, currentRStatus=0):
        cands = []
        lparts = supports.lparts()
        lmiss = supports.lpartsInterX(col.miss())
        lin = supports.lpartsInterX(col.supp())

        for op in self.constraints.getCstr("allw_ops", side=side, currentRStatus=currentRStatus):
            for neg in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus):
                adv, clp = self.getAC(side, op, neg, lparts, lmiss, lin, self.isCond(currentRStatus))
                if adv is not None :
                    cands.append(Extension(self.constraints.getSSetts(), adv, clp, (side, op, neg, Literal(neg, BoolTerm(col.getId())))))
        return cands

    def getCandidates2(self, side, col, supports, currentRStatus=0):
        return self.getCandidatesNonBool(side, col, supports, currentRStatus)

    def getCandidates3(self, side, col, supports, currentRStatus=0):
        return self.getCandidatesNonBool(side, col, supports, currentRStatus)

    def getCandidatesNonBool(self, side, col, supports, currentRStatus=0):
        cands = []
        lparts = supports.lparts()
        lmiss = supports.lpartsInterX(col.miss())

        for cand in self.findCover(side, col, lparts, lmiss, supports, currentRStatus):
            cands.append(cand[1])
        return cands

    ##################################################
    ###### CONDITIONAL
    def isCond(self, currentRStatus=0):
        return self.constraints.isStatusCond(currentRStatus)
    
    def getCondition(self, colsC, supports):
        sI = supports.part_union(self.constraints.getSSetts().IDS_inter)
        sD = supports.part_union(self.constraints.getSSetts().IDS_diff)

        cond_sparts = SParts(self.constraints.getSSetts(), supports.nbRows(), [sI, sD.union(sI)])
        lparts = cond_sparts.lparts()
        cond_cand = self.getConditionCand(colsC, cond_sparts, lparts)
        if cond_cand is not None:
            supp = self.getCCandSupp(colsC, cond_cand)
            lin = cond_sparts.lpartsInterX(supp)
            cond_cand.setClp([lin, lparts], False)
            return cond_cand

    def getConditionCand(self, colsC, cond_sparts, lparts):
        cis = range(len(colsC))
        prev = None
        best = ([], 0)
        while len(cis) > 0:
            current = []
            for ci in cis:
                lmiss = cond_sparts.lpartsInterX(colsC[ci].miss())
                cands = self.findCover(1, colsC[ci], lparts, lmiss, cond_sparts, Constraints.getStatusCond())
                if len(cands) == 1:
                    cand = cands[0][1]
                    if best[1] < cand.getAcc():
                        best = ([len(current)], cand.getAcc())
                    elif best[1] == cand.getAcc():
                        best[0].append(len(current))
                    current.append(cand)
            if len(best[0]) == 0:
                cis = []
            else:
                basis = (None, None, cond_sparts.nbRows(), 0.)
                for cc in best[0]:
                    cand = current[cc]
                    supp = colsC[cand.getLiteral().colId()].suppLiteral(cand.getLiteral())
                    if cand.getVarRed() > basis[-1] or (cand.getVarRed() == basis[-1] and len(supp) < basis[-2]):
                        basis = (cc, supp, len(supp), cand.getVarRed())
                cis = [ci for cii, ci in enumerate(cis) if basis[0] != cii]
                keep_cand, keep_supp = current[basis[0]], basis[1]
                if prev is None:
                    keep_cand.setLiteral([keep_cand.getLiteral()])
                else:
                    keep_cand.setLiteral([keep_cand.getLiteral()]+prev.getLiteral())
                prev = keep_cand
                cond_sparts.update(1, False, keep_supp)
                lparts = cond_sparts.lparts()
                best = ([], prev.getAcc())
        return prev
            
    def getCCandSupp(self, colsC, cond_cand):
        lits = cond_cand.getLiteral()
        if type(lits) is Literal:
            lits = [lits]
        return set.intersection(*[colsC[lit.colId()].suppLiteral(lit) for lit in lits])
    ##################################################

    
####################### COVER METHODS
    def findCover(self, side, col, lparts, lmiss, supports, currentRStatus=0):
        method_string = 'self.findCover%i' % col.typeId()
        try:
            method_compute =  eval(method_string)
        except AttributeError:
              raise Exception('Oups No covering method for this type of data (%i)!'  % col.typeId())
        return method_compute(side, col, lparts, lmiss, supports, currentRStatus)

    def findCover1(self, side, col, lparts, lmiss, supports, currentRStatus=0):
        cands = []
        lin = supports.lpartsInterX(col.supp())
        for op in self.constraints.getCstr("allw_ops", side=side, currentRStatus=currentRStatus):
            for neg in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus):
                tmp_adv, tmp_clp  = self.getAC(side, op, neg, lparts, lmiss, lin, self.isCond(currentRStatus))
                if tmp_adv is not None:
                    cands.append((False, Extension(self.constraints.getSSetts(), tmp_adv, tmp_clp, [side, op, neg, Literal(neg, BoolTerm(col.getId()))])))

                ### to negate the other side when looking for initial pairs
                if self.constraints.getCstr("neg_query_init", side=side, currentRStatus=currentRStatus):
                    tmp_adv, tmp_clp  = self.getAC(side, op, neg, \
                                                       self.constraints.getSSetts().negateParts(1-side, lparts), self.constraints.getSSetts().negateParts(1-side, lmiss), self.constraints.getSSetts().negateParts(1-side, lin), self.isCond(currentRStatus))
                    if tmp_adv is not None:
                        cands.append((True, Extension(self.constraints.getSSetts(), tmp_adv, tmp_clp, [side, op, neg, Literal(neg, BoolTerm(col.getId()))])))

        return cands

    def findCover2(self, side, col, lparts, lmiss, supports, currentRStatus=0):
        cands = []
        allw_neg = True in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus)
        negs = self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus)
        if self.constraints.getCstr("multi_cats"): ## redundant negation
            negs = [False]
            
        for op in self.constraints.getCstr("allw_ops", side=side, currentRStatus=currentRStatus):
            for neg in negs:
                best = (None, None, None)
                bestNeg = (None, None, None)
                collect_goods = []
                collect_goodsNeg = []
                for (cat, supp) in col.iter_cats():
                    lin = supports.lpartsInterX(supp)
                    ## best = self.updateACT(best, Literal(neg, CatTerm(col.getId(), cat)), side, op, neg, lparts, lmiss, lin)
                    ######################
                    tmp_adv = self.getAC(side, op, neg, lparts, lmiss, lin, self.isCond(currentRStatus))
                    if best[0] < tmp_adv[0]:
                        best = (tmp_adv[0], tmp_adv[1], [side, op, neg, Literal(neg, CatTerm(col.getId(), cat))])
                    if tmp_adv[0] is not None:
                        collect_goods.append((tmp_adv[0], cat, lparts, lmiss, lin))
                    ######################
                    
                    ### to negate the other side when looking for initial pairs
                    if self.constraints.getCstr("neg_query_init", side=side, currentRStatus=currentRStatus):
                        # bestNeg = self.updateACT(bestNeg, Literal(neg, CatTerm(col.getId(), cat)), side, op, neg, \
                        #                       self.constraints.getSSetts().negateParts(1-side, lparts), self.constraints.getSSetts().negateParts(1-side, lmiss), self.constraints.getSSetts().negateParts(1-side, lin))
                        ######################
                        Nlparts = self.constraints.getSSetts().negateParts(1-side, lparts)
                        Nlmiss = self.constraints.getSSetts().negateParts(1-side, lmiss)
                        Nlin = self.constraints.getSSetts().negateParts(1-side, lin)
                        tmp_adv = self.getAC(side, op, neg, Nlparts, Nlmiss, Nlin, self.isCond(currentRStatus))
                        if bestNeg[0] < tmp_adv[0]:
                            bestNeg = (tmp_adv[0], tmp_adv[1], [side, op, neg, Literal(neg, CatTerm(col.getId(), cat))])
                        if tmp_adv[0] is not None:
                            collect_goodsNeg.append((tmp_adv[0], cat, Nlparts, Nlmiss, Nlin))
                        ######################
                    
                #### HERE CATS
                if best[0] is not None:
                    bb = self.combCats(best, allw_neg, side, op, neg, col, collect_goods, currentRStatus=currentRStatus)
                    cands.append((False, Extension(self.constraints.getSSetts(), bb)))
                    # cands.append((False, Extension(self.constraints.getSSetts(), best)))
                if bestNeg[0] is not None:
                    bb = self.combCats(best, allw_neg, side, op, neg, col, collect_goodsNeg, currentRStatus=currentRStatus)
                    cands.append((True, Extension(self.constraints.getSSetts(), bb)))
                    # cands.append((True, Extension(self.constraints.getSSetts(), bestNeg)))
        return cands

    def additionsLParts(self, lparts, lmiss, lin, nextc, other_side=0):
        if other_side == 0:
            ccum_lparts = lparts
            ccum_lmiss = lmiss
            ccum_lin = self.constraints.getSSetts().addition(lin, nextc[-1])
        else:
            ccum_lparts = self.constraints.getSSetts().additionOtherSide(lparts, nextc[-3], other_side < 0)
            ccum_lmiss = self.constraints.getSSetts().additionOtherSide(lmiss, nextc[-2], other_side < 0)       
            ccum_lin = self.constraints.getSSetts().additionOtherSide(lin, nextc[-1], other_side < 0)
            
        # print "--- SUM", other_side
        # print lparts, lmiss, lin
        # print nextc[-3], nextc[-2], nextc[-1]
        # print ">>"
        # print ccum_lparts, ccum_lmiss, ccum_lin
        return ccum_lparts, ccum_lmiss, ccum_lin
            
    def combCats(self, best, allw_neg, side, op, neg, col, collected, other_side=0, currentRStatus=0):
        if not self.constraints.getCstr("multi_cats"):
            return best
        collected.sort()
        nextc = collected.pop()
        best_cat = [nextc[1]]
        best_score = (nextc[0], None)
        cum_lparts, cum_lmiss, cum_lin = (list(nextc[-3]), list(nextc[-2]), list(nextc[-1]))
        while len(collected) > 0:
            nextc = collected.pop()
            ccum_lparts, ccum_lmiss, ccum_lin = self.additionsLParts(cum_lparts, cum_lmiss, cum_lin, nextc, other_side)
            tmp_adv = self.getAC(side, op, neg, ccum_lparts, ccum_lmiss, ccum_lin, self.isCond(currentRStatus))
            if best_score[0] < tmp_adv[0]:
                best_cat.append(nextc[1])
                cum_lparts = ccum_lparts
                cum_lmiss = ccum_lmiss
                cum_lin = ccum_lin                    
                best_score = tmp_adv                
        if len(best_cat) > 1:
            if col is None:
                best = (best_score[0], best_score[1], [side, op, neg, set(best_cat)])
            else:
                lit = col.makeCatLit(best_cat, neg, allw_neg)
                best = (best_score[0], best_score[1], [side, op, lit.isNeg(), lit])
        return best
        
    def findCover3(self, side, col, lparts, lmiss, supports, currentRStatus=0):
        cands = []
        if self.constraints.isStatusCond(currentRStatus) or self.inSuppBounds(side, True, lparts) or self.inSuppBounds(side, False, lparts):  ### DOABLE
            segments = col.makeSegments(self.constraints.getSSetts(), side, supports, self.constraints.getCstr("allw_ops", side=side, currentRStatus=currentRStatus))
            for cand in self.findCoverSegments(side, col, segments, lparts, lmiss, currentRStatus):
                cands.append((False, cand))

            ### to negate the other side when looking for initial pairs
            if self.constraints.getCstr("neg_query_init", side=side, currentRStatus=currentRStatus):
                nlparts = self.constraints.getSSetts().negateParts(1-side, lparts)
                nlmiss = self.constraints.getSSetts().negateParts(1-side, lmiss)

                if self.inSuppBounds(side, True, nlparts): ### DOABLE
                    nsegments = col.makeSegments(self.constraints.getSSetts(), side, supports.negate(1-side), self.constraints.getCstr("allw_ops", side=side, currentRStatus=currentRStatus))
                    ##H pdb.set_trace()
                    for cand in self.findCoverSegments(side, col, nsegments, nlparts, nlmiss, currentRStatus):
                        cands.append((True, cand))
        return cands
        
    def findCoverSegments(self, side, col, segments, lparts, lmiss, currentRStatus=0):
        cands = []
        for op in self.constraints.getCstr("allw_ops", side=side, currentRStatus=currentRStatus):
            if len(segments[op]) < self.constraints.getCstr("max_seg"):
                cands.extend(self.findCoverFullSearch(side, op, col, segments, lparts, lmiss, currentRStatus))
            else:
                if (False in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus)):
                    cands.extend(self.findPositiveCover(side, op, col, segments, lparts, lmiss, currentRStatus))
                if (True in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus)):
                    cands.extend(self.findNegativeCover(side, op, col, segments, lparts, lmiss, currentRStatus))
        return cands

    def findCoverFullSearch(self, side, op, col, segments, lparts, lmiss, currentRStatus=0):
        cands = []
        bests = {False: (None, None, None), True: (None, None, None)}

        for seg_s in range(len(segments[op])):
            lin = self.constraints.getSSetts().makeLParts()
            for seg_e in range(seg_s,len(segments[op])):
                lin = self.constraints.getSSetts().addition(lin, segments[op][seg_e][2])
                for neg in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus):
                    bests[neg] = self.updateACT(bests[neg], (seg_s, seg_e), side, op, neg, lparts, lmiss, lin, self.isCond(currentRStatus))

        for neg in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus):
            if bests[neg][0]:
                bests[neg][-1][-1] = col.getLiteralSeg(neg, segments[op], bests[neg][-1][-1])
                if bests[neg][-1][-1] is not None:
                    cands.append(Extension(self.constraints.getSSetts(), bests[neg]))
        return cands

    def findNegativeCover(self, side, op, col, segments, lparts, lmiss, currentRStatus=0):
        is_cond = self.isCond(currentRStatus)
        cands = []
        lin_f = self.constraints.getSSetts().makeLParts()
        # bests_f = [(0, 0, self.constraints.getSSetts().makeLParts())]
        bests_f = [(self.advAcc(side, op, False, lparts, lmiss, lin_f, is_cond), 0, lin_f)] 
        best_track_f = [0]
        lin_b = self.constraints.getSSetts().makeLParts()
        # bests_b = [(0, 0, self.constraints.getSSetts().makeLParts())]
        bests_b = [(self.advAcc(side, op, False, lparts, lmiss, lin_b, is_cond), 0, lin_b)]
        best_track_b = [0]

        for  i in range(len(segments[op])):
            # FORWARD
            lin_f = self.constraints.getSSetts().addition(lin_f, segments[op][i][2])
            if  self.advRatioVar(side, op, lin_f, is_cond) > bests_f[-1][0]:
                lin_f = self.constraints.getSSetts().addition(lin_f, bests_f[-1][2])
                bests_f.append((self.advAcc(side, op, False, lparts, lmiss, lin_f, is_cond), i+1, lin_f))
                lin_f = self.constraints.getSSetts().makeLParts()
            best_track_f.append(len(bests_f)-1)

            # BACKWARD
            lin_b = self.constraints.getSSetts().addition(lin_b, segments[op][-(i+1)][2])
            if  self.advRatioVar(side, op, lin_b, is_cond) > bests_b[-1][0]:
                lin_b = self.constraints.getSSetts().addition(lin_b, bests_b[-1][2])
                bests_b.append((self.advAcc(side, op, False, lparts, lmiss, lin_b, is_cond), i+1, lin_b))
                lin_b = self.constraints.getSSetts().makeLParts()
            best_track_b.append(len(bests_b)-1)

        #pdb.set_trace()
        best_t = (None, None, None)
        for b in bests_b:
            if b[1] == len(segments[op]):
                f = bests_f[0]
            else:
                f = bests_f[best_track_f[len(segments[op])-(b[1]+1)]]
            if self.advRatioVar(side, op, f[2], is_cond) > b[0]:
                best_t = self.updateACT(best_t, (f[1], len(segments[op]) - (b[1]+1)), side, op, False, lparts, lmiss, self.constraints.getSSetts().addition(f[2], b[2]), is_cond)
            else:
                best_t = self.updateACT(best_t, (0, len(segments[op]) - (b[1]+1)), side, op, False, lparts, lmiss, b[2], is_cond)

        for f in bests_f:
            if f[1] == len(segments[op]):
                b = bests_b[0]
            else:
                b = bests_b[best_track_b[len(segments[op])-(f[1]+1)]]
            if self.advRatioVar(side, op, b[2], is_cond) > f[0]: 
                best_t = self.updateACT(best_t, (f[1], len(segments[op]) - (b[1]+1)), side, op, False, lparts, lmiss, self.constraints.getSSetts().addition(f[2], b[2]), is_cond)
            else:
                best_t = self.updateACT(best_t, (f[1], len(segments[op])-1), side, op, False, lparts, lmiss, f[2], is_cond)

        if best_t[0] is not None:
            tmp = best_t[-1][-1]
            best_t[-1][-1] = col.getLiteralSeg(True, segments[op], best_t[-1][-1])
            if best_t[-1][-1] is not None:
                cands.append(Extension(self.constraints.getSSetts(), best_t))
        return cands

    def findPositiveCover(self, side, op, col, segments, lparts, lmiss, currentRStatus=0):
        is_cond = self.isCond(currentRStatus)
        cands = []
        lin_f = self.constraints.getSSetts().makeLParts()
        nb_seg_f = 0
        best_f = (None, None, None)
        lin_b = self.constraints.getSSetts().makeLParts()
        nb_seg_b = 0
        best_b = (None, None, None)

        for i in range(len(segments[op])-1):
            # FORWARD
            if i > 0 and self.advAcc(side, op, False, lparts, lmiss, segments[op][i][2], is_cond) < self.advRatioVar(side, op, lin_f, is_cond):
                lin_f = self.constraints.getSSetts().addition(lin_f, segments[op][i][2])
                nb_seg_f += 1
            else: 
                lin_f = segments[op][i][2]
                nb_seg_f = 0
            best_f = self.updateACT(best_f, (i - nb_seg_f, i), side, op, False, lparts, lmiss, lin_f, is_cond)

            # BACKWARD
            if i > 0 and self.advAcc(side, op, False, lparts, lmiss, segments[op][-(i+1)][2], is_cond) < self.advRatioVar(side, op, lin_b, is_cond):
                lin_b = self.constraints.getSSetts().addition(lin_b, segments[op][-(i+1)][2])
                nb_seg_b += 1
            else:
                lin_b = segments[op][-(i+1)][2]
                nb_seg_b = 0
            best_b = self.updateACT(best_b, (len(segments[op])-(1+i), len(segments[op])-(1+i) + nb_seg_b), \
                                    side, op, False, lparts, lmiss, lin_b, is_cond)

        if best_b[0] is not None and best_f[0] is not None:
            bests = [best_b, best_f]

            if best_b[-1][-1][0] > best_f[-1][-1][0] and best_b[-1][-1][1] > best_f[-1][-1][1] and best_b[-1][-1][0] <= best_f[-1][-1][1]:
                lin_m = self.constraints.getSSetts().makeLParts()
                for seg in segments[op][best_b[-1][-1][0]:best_f[-1][-1][1]+1]:
                    lin_m = self.constraints.getSSetts().addition(lin_m, seg[2])
                tmp_adv_m, tmp_clp_m  = self.getAC(side, op, False, lparts, lmiss, lin_m, is_cond)
                if tmp_adv_m is not None:
                    bests.append((tmp_adv_m, tmp_clp_m, [side, op, False, (best_b[-1][-1][0], best_f[-1][-1][1])]))

            bests.sort()
            best = bests[-1]
            
        elif not best_f[0] is not None:
            best = best_f
        else:
            best = best_b

        if best[0] is not None:
            best[-1][-1] = col.getLiteralSeg(False, segments[op], best[-1][-1])
            if best[-1][-1] is not None:
                cands.append(Extension(self.constraints.getSSetts(), best))
        return cands

################################################################### PAIRS METHODS
###################################################################

    def computePair(self, colL, colR, colsC=None):
        min_type = min(colL.typeId(), colR.typeId())
        max_type = max(colL.typeId(), colR.typeId())
        method_string = 'self.do%i%i' % (min_type, max_type)
        try:
            method_compute =  eval(method_string)
        except AttributeError:
              raise Exception('Oups this combination does not exist (%i %i)!'  % (min_type, max_type))
        if colL.typeId() == min_type:
            (scores, literalsL, literalsR) = method_compute(colL, colR, 1)
        else:
            (scores, literalsR, literalsL) =  method_compute(colL, colR, 0)
        pairs = []
        for i in range(len(scores)):
            pair = {"litL": literalsL[i], "litR": literalsR[i], "score": scores[i]}

            #### compute additional condition
            if colsC is not None and self.constraints.getCstr("add_condition"):
                rsparts = SParts(self.constraints.getSSetts(), colL.nbRows(),
                                 [colL.suppLiteral(pair["litL"]), colR.suppLiteral(pair["litR"]), colL.miss(), colR.miss()])
                
                sI = rsparts.part_union(self.constraints.getSSetts().IDS_inter)
                sD = rsparts.part_union(self.constraints.getSSetts().IDS_diff)
                
                cond_sparts = SParts(self.constraints.getSSetts(), colL.nbRows(), [sI, sD.union(sI)])
                lparts = cond_sparts.lparts()        
                cond_cand = self.getConditionCand(colsC, cond_sparts, lparts)
                if cond_cand is not None:
                    pair["litC"] = cond_cand.getLiteral()
                    pair["scoreC"] = cond_cand.getAcc()
            pairs.append(pair)
        return pairs

    def doBoolStar(self, colL, colR, side):
        if side == 1:
            (supports, fixTerm, extCol) = (SParts(self.constraints.getSSetts(), colL.nbRows(), [colL.supp(), set(), colL.miss(), set()]), BoolTerm(colL.getId()), colR)
        else:
            (supports, fixTerm, extCol) = (SParts(self.constraints.getSSetts(), colL.nbRows(), [set(), colR.supp(), set(), colR.miss()]), BoolTerm(colR.getId()), colL)

        return self.fit(extCol, supports, side, fixTerm)

    def fit(self, col, supports, side, termX):
        (scores, literalsFix, literalsExt) = ([], [], [])   
        lparts = supports.lparts()
        lmiss = supports.lpartsInterX(col.miss())
        currentRStatus = Constraints.getStatusPair(col, side, termX)
        cands = self.findCover(side, col, lparts, lmiss, supports, currentRStatus=currentRStatus)
        for cand in cands:
            scores.append(cand[1].getAcc())
            literalsFix.append(Literal(cand[0], termX))
            literalsExt.append(cand[1].getLiteral())
        return (scores, literalsFix, literalsExt)

    def do11(self, colL, colR, side):
        return self.doBoolStar(colL, colR, side)

    def do12(self, colL, colR, side):
        return self.doBoolStar(colL, colR, side)
        
    def do13(self, colL, colR, side):
        return self.doBoolStar(colL, colR, side)
        
    def do22(self, colL, colR, side):
        return self.subdo22Full(colL, colR, side)

    def do23(self, colL, colR, side):
        return self.subdo23Full(colL, colR, side)
    
    def do33(self, colL, colR, side):
#         if len(init_info[1-side]) == 3: # fit FULL
#             if len(init_info[1-side][0]) > len(init_info[side][0]): 
#                 (scores, literalsB, literalsA)= self.subdo33Full(colL, colR, side)
#                 return (scores, literalsA, literalsB)
#             else:
#                 return self.subdo33Full(colL, colR, side)
#         else:
#             return self.subdo33Heur(colL, colR, side)
        if len(colL.interNonMode(colR.nonModeSupp())) >= self.constraints.getCstr("min_itm_in") :
            return self.subdo33Full(colL, colR, side)
        else:
            return ([], [], [])
    
    def subdo33Heur(self, colL, colR, side):
        ### Suitable for sparse data
        bestScore = None
        if True: ### DOABLE
            ## FIT LHS then RHS
            supports = SParts(self.constraints.getSSetts(), colL.nbRows(), [set(), colR.nonModeSupp(), set(), colR.miss()])
            (scoresL, literalsFixL, literalsExtL) = self.fit(colL, supports, 0, idR)
            for tL in literalsExtL:
                suppL = colL.suppLiteral(tL)
                supports = SParts(self.constraints.getSSetts(), colL.nbRows(), [suppL, set(), colL.miss(), set()])
                (scoresR, literalsFixR, literalsExtR) = self.fit(colR, supports, 1, tL)
                for i in range(len(scoresR)):
                    if scoresR[i] > bestScore:
                        (scores, literalsL, literalsR) = ([scoresR[i]], [literalsFixR[i]], [literalsExtR[i]])
                        bestScore = scoresR[i]
                        
            ## FIT RHS then LHS
            supports = SParts(self.constraints.getSSetts(), colL.nbRows(), [colL.nonModeSupp(), set(), colL.miss(), set()])
            (scoresR, literalsFixR, literalsExtR) = self.fit(colR, supports, 1, idL)
            for tR in literalsExtR:
                suppR = colR.suppLiteral(tR)
                supports = SParts(self.constraints.getSSetts(), colL.nbRows(), [set(), suppR, set(), colR.miss()])
                (scoresL, literalsFixL, literalsExtL) = self.fit(colL, supports, 0, tR)
                for i in range(len(scoresL)):
                    if scoresL[i] > bestScore:
                        (scores, literalsL, literalsR) = ([scoresR[i]], [literalsExtL[i]], [literalsFixL[i]])
                        bestScore = scoresL[i]
                        
#             if len(scores) > 0:
#                print "%f: %s <-> %s" % (scores[0], literalsA[0], literalsB[0])
        return (scores, literalsL, literalsR)

    def subdo33Full(self, colL, colR, side):
        org_side = side
        best = []
        bUpE=1
        bUpF=1
        interMat = []
        bucketsL = colL.buckets()
        bucketsR = colR.buckets()

        if len(bucketsL[0]) > len(bucketsR[0]):
            bucketsF = bucketsR; colF = colR; bucketsE = bucketsL; colE = colL; side = 1-side; flip_side = True
        else:
            bucketsF = bucketsL; colF = colL; bucketsE = bucketsR; colE = colR; flip_side = False
            
        (scores, literalsF, literalsE) = ([], [], [])
        ## DOABLE
        
        # print "Nb buckets: %i x %i"% (len(bucketsF[1]), len(bucketsE[1]))
        # if ( len(bucketsF[1]) * len(bucketsE[1]) > self.constraints.getCstr("max_prodbuckets") ): 
        nbb = self.constraints.getCstr("max_prodbuckets") / float(len(bucketsF[1]))
        if len(bucketsE[1]) > nbb: ## self.constraints.getCstr("max_sidebuckets"):

            if len(bucketsE[1])/nbb < self.constraints.getCstr("max_agg"):
                ### collapsing buckets on the largest side is enough to get within the reasonable size
                bucketsE = colE.collapsedBuckets(self.constraints.getCstr("max_agg"), nbb)
                bUpE=3 ## in case of collapsed bucket the threshold is different

            else:
                ### collapsing buckets on the largest side is NOT enough to get within the reasonable size
                bucketsE = None

                #### try cats
                bbs = [dict([(bi, es) for (bi, es) in enumerate(bucketsL[0]) \
                                 if ( len(es) > self.constraints.getCstr("min_itm_in") and \
                                          colL.nbRows() - len(es) > self.constraints.getCstr("min_itm_out"))]),
                       dict([(bi, es) for (bi, es) in enumerate(bucketsR[0]) \
                                 if ( len(es) > self.constraints.getCstr("min_itm_in") and \
                                          colR.nbRows() - len(es) > self.constraints.getCstr("min_itm_out"))])]

                ## if len(bbs[0]) > 0 and ( len(bbs[1]) == 0 or len(bbs[0])/float(len(bucketsL[0])) < len(bbs[1])/float(len(bucketsR[0]))):
                
                nbes = [float(max(sum([len(v) for (k,v) in bbs[s].items()]), .5)) for s in [0,1]]
                side = None
                if len(bbs[0]) > 0 and ( len(bbs[1]) == 0 or nbes[0]/len(bbs[0]) > nbes[1]/len(bbs[1]) ):
                    ccL, ccR, side = (CatColM(bbs[0], colL.nbRows(), colL.miss()), colR, 1) 
                elif len(bbs[1]) > 0:
                    ccL, ccR, side = (colL, CatColM(bbs[1], colR.nbRows(), colR.miss()), 0) 

                if side is not None:
                    #### working with on variable as categories is workable
                    ## print "Trying cats...", len(bucketsL[0]), len(bucketsR[0]), len(bbs[0]), len(bbs[1])
                    (scores, literalsFix, literalsExt) = self.subdo23Full(ccL, ccR, side, try_comb=False)
                    if side == 1:
                        literalsL = []
                        literalsR = literalsExt
                        for ltc in literalsFix:
                            c = ltc.getTerm().getCat()
                            if type(c) is set and len(c) > 0:
                                c = sorted(c)[0]
                            val = bucketsL[1][c]
                            literalsL.append( Literal(ltc.isNeg(), NumTerm(colL.getId(), val, val)) )
                    else:
                        literalsL = literalsExt                    
                        literalsR = []
                        for ltc in literalsFix:
                            c = ltc.getTerm().getCat()
                            if type(c) is set and len(c) > 0:
                                c = sorted(c)[0]
                            val = bucketsR[1][c]
                            literalsR.append( Literal(ltc.isNeg(), NumTerm(colR.getId(), val, val)) )

                    return (scores, literalsL, literalsR)

                else:
                    #### working with on variable as categories is NOT workable
                    ### the only remaining solution is aggressive collapse of buckets on both sides
                    nbb = numpy.sqrt(self.constraints.getCstr("max_prodbuckets"))
                    bucketsE = colE.collapsedBuckets(self.constraints.getCstr("max_agg"), nbb)
                    bUpE=3 ## in case of collapsed bucket the threshold is different
                    bucketsF = colF.collapsedBuckets(self.constraints.getCstr("max_agg"), nbb)
                    bUpF=3 ## in case of collapsed bucket the threshold is different
                    side = org_side
                    ## print "Last resort solution...", nbb, len(bucketsL[0]), len(bucketsR[0])                    

        if bucketsE is not None and ( len(bucketsF[1]) * len(bucketsE[1]) < self.constraints.getCstr("max_prodbuckets") ):
            partsMubB = len(colF.miss())
            missMubB = len(colF.miss() & colE.miss())
            totInt = colE.nbRows() - len(colF.miss()) - len(colE.miss()) + missMubB
            #margE = [len(intE) for intE in bucketsE[0]]
            
            lmissFinE = [len(colF.miss() & bukE) for bukE in bucketsE[0]]
            lmissEinF = [len(colE.miss() & bukF) for bukF in bucketsF[0]]
            margF = [len(bucketsF[0][i]) - lmissEinF[i] for i in range(len(bucketsF[0]))]
            totMissE = len(colE.miss())
            totMissEinF = sum(lmissEinF)
            
            for bukF in bucketsF[0]: 
                interMat.append([len(bukF & bukE) for bukE in bucketsE[0]])
            
            if bucketsF[2] is not None :
                margF[bucketsF[2]] += colF.lenMode()
                for bukEId in range(len(bucketsE[0])):
                    interMat[bucketsF[2]][bukEId] += len(colF.interMode(bucketsE[0][bukEId])) 

            if bucketsE[2] is not None :
                #margE[bucketsE[2]] += colE.lenMode()
                for bukFId in range(len(bucketsF[0])):
                    interMat[bukFId][bucketsE[2]] += len(colE.interMode(bucketsF[0][bukFId]))        

            if bucketsF[2] is not None and bucketsE[2] is not None:
                interMat[bucketsF[2]][bucketsE[2]] += len(colE.interMode(colF.modeSupp()))

#             ### check marginals
#             totF = 0
#             for iF in range(len(bucketsF[0])):
#                 sF = sum(interMat[iF])
#                 if sF != margF[iF]:
#                     raise Error('Error in computing the marginals (1)')
#                 totF += sF

#             totE = 0
#             for iE in range(len(bucketsE[0])):
#                 sE = sum([intF[iE] for intF in interMat])
#                 if sE != margE[iE]:
#                     raise Error('Error in computing the marginals (2)')
#                 totE += sE

#             if totE != totF or totE != colE.nbRows():
#                 raise Error('Error in computing the marginals (3)')


            belowF = 0
            lowF = 0
            while lowF < len(interMat) and totInt - belowF >= self.constraints.getCstr("min_itm_in"):

                aboveF = 0
                upF = len(interMat)-1
                while upF >= lowF and totInt - belowF - aboveF >= self.constraints.getCstr("min_itm_in"):
                    if belowF + aboveF  >= self.constraints.getCstr("min_itm_out"):
                        EinF = [sum([interMat[iF][iE] for iF in range(lowF,upF+1)]) for iE in range(len(interMat[lowF]))]
                        EoutF = [sum([interMat[iF][iE] for iF in range(0,lowF)+range(upF+1,len(interMat))]) for iE in range(len(interMat[lowF]))]
                        lmissE = sum(lmissEinF[lowF:upF+1])
                        #totEinF = sum(EinF)
                        lparts = self.constraints.getSSetts().makeLParts([(self.constraints.getSSetts().partId(self.constraints.getSSetts().Exo, 1-side), totInt - aboveF - belowF + lmissE), \
                                                    (self.constraints.getSSetts().partId(self.constraints.getSSetts().Emo, 1-side), partsMubB ), \
                                                    (self.constraints.getSSetts().partId(self.constraints.getSSetts().Eoo, 1-side), aboveF + belowF + totMissEinF - lmissE)], 0)
                        lmiss  = self.constraints.getSSetts().makeLParts([(self.constraints.getSSetts().partId(self.constraints.getSSetts().Exo, 1-side), lmissE ), \
                                                    (self.constraints.getSSetts().partId(self.constraints.getSSetts().Emo, 1-side), missMubB ), \
                                                    (self.constraints.getSSetts().partId(self.constraints.getSSetts().Eoo, 1-side), totMissEinF - lmissE )], 0)

                        belowEF = 0
                        outBelowEF = 0
                        lowE = 0
                        while lowE < len(interMat[lowF]) and totInt - belowF - aboveF - belowEF >= self.constraints.getCstr("min_itm_in"):
                            aboveEF = 0
                            outAboveEF = 0
                            upE = len(interMat[lowF])-1
                            while upE >= lowE and totInt - belowF - aboveF - belowEF - aboveEF >= self.constraints.getCstr("min_itm_in"):
                                
                                lmissF = sum(lmissFinE[lowE:upE+1])
                                lin = self.constraints.getSSetts().makeLParts([(self.constraints.getSSetts().partId(self.constraints.getSSetts().Exo, 1-side), totInt - belowF - aboveF - belowEF - aboveEF), \
                                                         (self.constraints.getSSetts().partId(self.constraints.getSSetts().Emo, 1-side), lmissF ), \
                                                         (self.constraints.getSSetts().partId(self.constraints.getSSetts().Eoo, 1-side), belowF + aboveF - outAboveEF - outBelowEF)], 0)

                                best = self.updateACTP33(best, (lowF, upF, lowE, upE), side, True, False, lparts, lmiss, lin)
                                aboveEF+=EinF[upE]
                                outAboveEF+=EoutF[upE]
                                upE-=1
                            belowEF+=EinF[lowE]
                            outBelowEF+=EoutF[lowE]
                            lowE+=1
                    aboveF+=margF[upF]
                    upF-=1
                belowF+=margF[lowF]
                lowF+=1

        for b in best:
            tF = colF.getLiteralBuk(False, bucketsF[1], b[-1][-1][0:2], bucketsF[bUpF])
            tE = colE.getLiteralBuk(False, bucketsE[1], b[-1][-1][2:], bucketsE[bUpE])
            if tF is not None and tE is not None:
                literalsF.append(tF)
                literalsE.append(tE)
                scores.append(b[0][0])

        if flip_side:
            return (scores, literalsE, literalsF)
        else:
            return (scores, literalsF, literalsE)

    def subdo22Full(self, colL, colR, side):
        configs = [(0, False, False), (1, False, True), (2, True, False), (3, True, True)]
        allw_neg = True
        if True not in self.constraints.getCstr("neg_query", side=side, type_id=2):
            configs = configs[:1]
            allw_neg = False
        best = [[] for c in configs]

        # print "--------------------------------------"
        # print "\t".join(["", ""]+[catR for catR in colR.cats()])
        # print "\t".join(["", ""]+["%d" % len(colR.suppCat(catR)) for catR in colR.cats()])
        # for catL in colL.cats():
        #     print "\t".join([catL, "%d" % len(colL.suppCat(catL))]+["%d" % len(colL.suppCat(catL).intersection(colR.suppCat(catR))) for catR in colR.cats()])
        # print "--------------------------------------"
        
        for catL in colL.cats():
            ### TODO DOABLE
            supports = SParts(self.constraints.getSSetts(), colL.nbRows(), [colL.suppCat(catL), set(), colL.miss(), set()])
            lparts = supports.lparts()
            lmiss = supports.lpartsInterX(colR.miss())
            
            for catR in colR.cats():
                lin = supports.lpartsInterX(colR.suppCat(catR))
                for (i, nL, nR) in configs:
                    if nL:
                        tmp_lparts = self.constraints.getSSetts().negateParts(0, lparts)
                        tmp_lmiss = self.constraints.getSSetts().negateParts(0, lmiss)
                        tmp_lin = self.constraints.getSSetts().negateParts(0, lin)
                    else:
                        tmp_lparts = lparts
                        tmp_lmiss = lmiss
                        tmp_lin = lin

                    # print "--", catL, catR, nL, nR
                    # print tmp_lparts, tmp_lmiss, tmp_lin
                    best[i] = self.updateACTP22(best[i], (catL, catR), side, True, nR, tmp_lparts, tmp_lmiss, tmp_lin)

        (scores, literalsFix, literalsExt) = ([], [], [])
        if self.constraints.getCstr("multi_cats"):
            # print "PAIR ----", colL, colR
            (scores, literalsFix, literalsExt) = self.combPairsCats(best, [colL, colR], configs, allw_neg)
            # if len(scores) > 0:
            #     print "---- Multi cats:"
            #     for ii in range(len(scores)):
            #         print scores[ii], literalsFix[ii], literalsExt[ii]
            
        for (i, nL, nR) in configs:
            for b in best[i]:
                scores.append(b[0][0])
                literalsFix.append(Literal(nL, CatTerm(colL.getId(), b[-1][-1][0])))
                literalsExt.append(Literal(nR, CatTerm(colR.getId(), b[-1][-1][1])))
        return (scores, literalsFix, literalsExt)

    def combPairsCats(self, best, cols, configs, allw_neg):
        (scores, literalsFix, literalsExt) = ([], [], [])
        for (i, nL, nR) in configs:
            map_cat = [{}, {}]
            for b in best[i]:                
                for ss in [0,1]:
                    if b[-1][-1][ss] not in map_cat[ss]:
                        map_cat[ss][b[-1][-1][ss]] = []
                    map_cat[ss][b[-1][-1][ss]].append((b[0][0], b[-1][-1][1-ss], b[1][2], b[1][3], b[1][0]))
                        
            for ss in [0,1]:
                cats_loc = [None, None]
                kk = map_cat[ss].keys()
                for k in kk:
                    tmp_cats = [] #(c[0], c[1]) for c in map_cat[ss][k]] 
                    if len(map_cat[ss][k]) > 1:
                        other_side = 0
                        if (ss==1):
                            other_side = 1
                            if nL: other_side = -1
                        # print "==============="
                        # print "GO", k, nR, nL, other_side, "\n*".join([""]+[str(mm) for mm in map_cat[ss][k]])
                        bb = self.combCats(None, allw_neg, 1, True, nR, None, map_cat[ss][k], other_side=other_side)
                        if bb is not None:
                            tmp_cats = [tmp_cat for tmp_cat in tmp_cats if tmp_cat[1] not in bb[-1][-1]]
                            tmp_cats.append((bb[0][0], bb[-1][-1]))
                        for tmp_cat in tmp_cats:
                            cats_loc[ss] = k
                            cats_loc[1-ss] = tmp_cat[1]
                            scores.append(tmp_cat[0])
                            literalsFix.append(cols[0].makeCatLit(cats_loc[0], nL, allw_neg))
                            literalsExt.append(cols[1].makeCatLit(cats_loc[1], nR, allw_neg))                        
        return (scores, literalsFix, literalsExt)

    

    def subdo23Full(self, colL, colR, side, try_comb=True):
        if side == 0:
            (colF, colE) = (colR, colL)
        else:
            (colF, colE) = (colL, colR)

        configs = [(0, False, False), (1, False, True), (2, True, False), (3, True, True)]
        allw_neg = True
        if True not in self.constraints.getCstr("neg_query", side=side, type_id=2) or True not in self.constraints.getCstr("neg_query", side=side, type_id=3):
            configs = configs[:1]
            allw_neg = False
        best = [[] for c in configs]

        buckets = colE.buckets()
        bUp = 1
        nbb = self.constraints.getCstr("max_prodbuckets") / float(len(colF.cats()))
        if len(buckets[1]) > nbb: ## self.constraints.getCstr("max_sidebuckets"):
             bUp=3 ## in case of collapsed bucket the threshold is different
             buckets = colE.collapsedBuckets(self.constraints.getCstr("max_agg"), nbb)
             #pdb.set_trace()

        ### TODO DOABLE
        if buckets is not None and ( len(buckets[1]) * len(colF.cats()) <= self.constraints.getCstr("max_prodbuckets")): 
            partsMubB = len(colF.miss())
            missMubB = len(colF.miss() & colE.miss())
            
            missMat = [len(colF.miss() & buk) for buk in buckets[0]]
            totMiss = sum(missMat)

            marg = [len(buk) for buk in buckets[0]]
            if buckets[2] is not None :
                marg[buckets[2]] += colE.lenMode()

            for cat in colF.cats():
                lparts = self.constraints.getSSetts().makeLParts([(self.constraints.getSSetts().Exo, len(colF.suppCat(cat)) ), (self.constraints.getSSetts().Emo, partsMubB ), (self.constraints.getSSetts().Eoo, - colF.nbRows())], 1-side)
                lmiss  = self.constraints.getSSetts().makeLParts([(self.constraints.getSSetts().Exo, len(colF.suppCat(cat) & colE.miss()) ), (self.constraints.getSSetts().Emo, missMubB ), (self.constraints.getSSetts().Eoo, -len(colE.miss()) )], 1-side)

                interMat = [len(colF.suppCat(cat) & buk) for buk in buckets[0]]
                if buckets[2] is not None :
                    interMat[buckets[2]] += len(colE.interMode(colF.suppCat(cat)))        

                totIn = sum(interMat) 
                below = 0
                missBelow = 0
                low = 0
                while low < len(interMat) and \
                          (totIn - below >= self.constraints.getCstr("min_itm_in") or totIn - below >= self.constraints.getCstr("min_itm_out")):
                    above = 0
                    missAbove = 0
                    up = len(interMat)-1
                    while up >= low and \
                          (totIn - below - above >= self.constraints.getCstr("min_itm_in") or totIn - below - above >= self.constraints.getCstr("min_itm_out")):
                        lin = self.constraints.getSSetts().makeLParts([(self.constraints.getSSetts().Exo, totIn - below - above), (self.constraints.getSSetts().Emo, totMiss - missBelow - missAbove ), (self.constraints.getSSetts().Eoo, -sum(marg[low:up+1]))], 1-side)
                        for (i, nF, nE) in configs:
                            if nF:
                                tmp_lparts = self.constraints.getSSetts().negateParts(1-side, lparts)
                                tmp_lmiss = self.constraints.getSSetts().negateParts(1-side, lmiss)
                                tmp_lin = self.constraints.getSSetts().negateParts(1-side, lin)
                            else:
                                tmp_lparts = lparts
                                tmp_lmiss = lmiss
                                tmp_lin = lin

                            best[i] = self.updateACTP23(best[i], (cat, low, up), side, True, nE, tmp_lparts, tmp_lmiss, tmp_lin)
                            
                        above+=interMat[up]
                        missAbove+=missMat[up]
                        up-=1
                    below+=interMat[low]
                    missBelow+=missMat[low]
                    low+=1

        
        (scores, literalsFix, literalsExt) = ([], [], [])
        if try_comb and self.constraints.getCstr("multi_cats"):
            (scores, literalsFix, literalsExt) = self.combNumCats(best, [colF, colE], configs, allw_neg, side, buckets, bUp)
            # if len(scores) > 0:
            #     print "---- Multi cats:"
            #     for ii in range(len(scores)):
            #         print scores[ii], literalsFix[ii], literalsExt[ii]

        for (i, nF, nE) in configs:
            for b in best[i]:
                tE = colE.getLiteralBuk(nE, buckets[1], b[-1][-1][1:], buckets[bUp])
                #tE = colE.getLiteralBuk(nE, buckets[1], idE, b[-1][-1][1:])
                if tE is not None:
                    literalsExt.append(tE)
                    literalsFix.append(Literal(nF, CatTerm(colF.getId(), b[-1][-1][0])))
                    scores.append(b[0][0])
        return (scores, literalsFix, literalsExt)


    def combNumCats(self, best, cols, configs, allw_neg, side, buckets, bUp):
        (scores, literalsFix, literalsExt) = ([], [], [])
        for (i, nF, nE) in configs:
            map_cat = {}
            for b in best[i]:
                buk = b[-1][-1][1:]
                if buk not in map_cat:
                    map_cat[buk] = []
                map_cat[buk].append((b[0][0], b[-1][-1][0], b[1][2], b[1][3], b[1][0]))

            kk = map_cat.keys()
            for k in kk:
                tmp_cats = [] #(c[0], c[1]) for c in map_cat[k]] 
                if len(map_cat[k]) > 1:
                    other_side = 1
                    if nE: other_side = -1
                    bb = self.combCats(None, allw_neg, side, True, nF, None, map_cat[k], other_side=other_side)
                    if bb is not None:
                        tmp_cats = [tmp_cat for tmp_cat in tmp_cats if tmp_cat[1] not in bb[-1][-1]]
                        tmp_cats.append((bb[0][0], bb[-1][-1]))
                    for tmp_cat in tmp_cats:
                        tE = cols[1].getLiteralBuk(nE, buckets[1], k, buckets[bUp])
                        if tE is not None:
                            literalsExt.append(tE)
                            literalsFix.append(cols[0].makeCatLit(tmp_cat[1], nF, allw_neg))
                            scores.append(tmp_cat[0])                                                    
        return (scores, literalsFix, literalsExt)

    
##### TOOLS METHODS
    # compute the advance resulting of appending X on given side with given operator and negation
    # from intersections of X with parts (clp)            
    def getExtAC(self, ext):
        return self.getAC(ext.side, ext.op, ext.literal.isNeg(), ext.clp[2], ext.clp[3], ext.clp[0])

    def advRatioVar(self, side, op, lin_f, is_cond=False):
        if is_cond:
            den = self.constraints.getSSetts().sumPartsId(side, self.constraints.getSSetts().IDS_varden[False], lin_f)
            num = self.constraints.getSSetts().sumPartsId(side, self.constraints.getSSetts().IDS_varnum[False], lin_f)
            return tool_ratio(num, den+num)
        return self.constraints.getSSetts().advRatioVar(side, op, lin_f)
    def advAcc(self, side, op, neg, lparts, lmiss, lin, is_cond=False):
        if is_cond:
            return self.getAC_cond(side, op, neg, lparts, lmiss, lin, no_const=True)
        return self.constraints.getSSetts().advAcc(side, op, neg, lparts, lmiss, lin)
        
    def getAC_cond(self, side, op, neg, lparts, lmiss, lin, no_const=False):
        lout = [lparts[i] - lmiss[i] - lin[i] for i in range(len(lparts))]
        clp = (lin, lout, lparts, lmiss)

        contri = self.constraints.getSSetts().sumPartsIdInOut(side, neg, self.constraints.getSSetts().IDS_varnum[False], clp)
        # contriAND = self.constraints.getSSetts().sumPartsIdInOut(side, neg, self.constraints.getSSetts().IDS_cont[False], clp)
        # contriOR = self.constraints.getSSetts().sumPartsIdInOut(side, neg, self.constraints.getSSetts().IDS_cont[True], clp)
        # print "Contri COND=%d AND=%d OR=%d" % (contri, contriAND, contriOR)
        if no_const or (contri >= self.constraints.getCstr("min_itm_in")):
            # pdb.set_trace()
            fixBlue, fixRed = (0, 0)
            varBlue = self.constraints.getSSetts().sumPartsIdInOut(side, neg, self.constraints.getSSetts().IDS_varnum[False], clp)
            varRed = varBlue + self.constraints.getSSetts().sumPartsIdInOut(side, neg, self.constraints.getSSetts().IDS_varden[False], clp)
            acc = tool_ratio(varBlue, varRed)
            # print "PIECES", varBlue, varRed, contri, fixBlue, fixRed
            return (acc, varBlue, varRed, contri, fixBlue, fixRed), clp
        return None, clp

    
    def getAC(self, side, op, neg, lparts, lmiss, lin, is_cond=False, no_const=False):
        if is_cond:
            return self.getAC_cond(side, op, neg, lparts, lmiss, lin, no_const)

        lout = [lparts[i] - lmiss[i] - lin[i] for i in range(len(lparts))]
        clp = (lin, lout, lparts, lmiss)
        
        contri = self.constraints.getSSetts().sumPartsIdInOut(side, neg, self.constraints.getSSetts().IDS_cont[op], clp)        
        if no_const or (contri >= self.constraints.getCstr("min_itm_c")):
            varBlue = self.constraints.getSSetts().sumPartsIdInOut(side, neg, self.constraints.getSSetts().IDS_varnum[op], clp)
            fixBlue = self.constraints.getSSetts().sumPartsIdInOut(side, neg, self.constraints.getSSetts().IDS_fixnum[op], clp)
            if no_const or (varBlue+fixBlue >= self.constraints.getCstr("min_itm_in")):
                sout = self.constraints.getSSetts().sumPartsIdInOut(side, neg, self.constraints.getSSetts().IDS_out[op], clp)
                if no_const or (sout >= self.constraints.getCstr("min_itm_out")):
                    varRed = self.constraints.getSSetts().sumPartsIdInOut(side, neg, self.constraints.getSSetts().IDS_varden[op], clp)
                    fixRed = self.constraints.getSSetts().sumPartsIdInOut(side, neg, self.constraints.getSSetts().IDS_fixden[op], clp)
                    acc = tool_ratio(varBlue + fixBlue, varRed + fixRed)
                    # print "PIECES", sout, varBlue, varRed, contri, fixBlue, fixRed
                    return (acc, varBlue, varRed, contri, fixBlue, fixRed), clp
        return None, clp

    def updateACT(self, best, lit, side, op, neg, lparts, lmiss, lin, is_cond=False):
        tmp_adv = self.getAC(side, op, neg, lparts, lmiss, lin, is_cond)
        if best[0] < tmp_adv[0]:
            return tmp_adv[0], tmp_adv[1], [side, op, neg, lit]
        else:
            return best
        ### EX: best = self.updateACT(best, Literal(neg, BoolTerm(col.getId())), side, op, neg, lparts, lmiss, lin, col.nbRows())

    def updateACTP(self, best, lit, side, op, neg, lparts, lmiss, lin, conflictF):
        tmp_adv = self.getAC(side, op, neg, lparts, lmiss, lin)
        if tmp_adv[0] is None:
            return best
        inserted = False
        i = 0
        while i < len(best):
            if best[i][0] > tmp_adv[0]:
                if conflictF(best[i][-1][-1], lit):  ## found conflicting of better quality 
                    return best
            else:
                if not inserted: 
                    best.insert(i,(tmp_adv[0], tmp_adv[1], [side, op, neg, lit]))
                    inserted = True
                elif conflictF(best[i][-1][-1], lit): ## found conflicting of lesser quality, remove
                    best.pop(i)
                    i -=1
            i+=1
        if not inserted:
            best.append((tmp_adv[0], tmp_adv[1], [side, op, neg, lit]))
        return best


    def conflictP22(self, litA, litB):
        # return True
        if self.constraints.getCstr("multi_cats"):
            return (litA[0] == litB[0]) and (litA[1] == litB[1])
        else:
            return (litA[0] == litB[0]) or (litA[1] == litB[1])
    def conflictP23(self, litA, litB):
        # return True
        if self.constraints.getCstr("multi_cats"):
            return (litA[0] == litB[0]) and not (litA[1] > litB[2] or litB[1] > litA[2])
        else:
            return (litA[0] == litB[0]) or not (litA[1] > litB[2] or litB[1] > litA[2])
    def conflictP33(self, litA, litB):
        # return True
        return not ((litA[0] > litB[1] or litB[0] > litA[1]) and (litA[2] > litB[3] or litB[2] > litA[3]))

    def updateACTP22(self, best, lit, side, op, neg, lparts, lmiss, lin):
        return self.updateACTP(best, lit, side, op, neg, lparts, lmiss, lin, self.conflictP22)
    def updateACTP23(self, best, lit, side, op, neg, lparts, lmiss, lin):
        return self.updateACTP(best, lit, side, op, neg, lparts, lmiss, lin, self.conflictP23)
    def updateACTP33(self, best, lit, side, op, neg, lparts, lmiss, lin):
        return self.updateACTP(best, lit, side, op, neg, lparts, lmiss, lin, self.conflictP33)


    def inSuppBounds(self, side, op, lparts):
        return self.constraints.getSSetts().sumPartsId(side, self.constraints.getSSetts().IDS_varnum[op] + self.constraints.getSSetts().IDS_fixnum[op], lparts) >= self.constraints.getCstr("min_itm_in") \
               and self.constraints.getSSetts().sumPartsId(side, self.constraints.getSSetts().IDS_cont[op], lparts) >= self.constraints.getCstr("min_itm_c")
