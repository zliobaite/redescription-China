from classData import CatColM
from classConstraints import Constraints
from classCharbon import CharbonGreedy
from classExtension import Extension
from classSParts import SParts, tool_ratio
from classQuery import  *
import numpy
import pdb

class CharbonGStd(CharbonGreedy):

    name = "GreedyStd"
    
    def getCandidates(self, side, col, supports, red, colsC=None):
        currentRStatus = Constraints.getStatusRed(red, side)
        method_string = 'self.getCandidates%i' % col.typeId()
        try:
            method_compute =  eval(method_string)
        except AttributeError:
              raise Exception('Oups No candidates method for this type of data (%i)!'  % col.typeId())
        cands = method_compute(side, col, supports, currentRStatus)
        # print "======STD==========="
        for cand in cands:
            supp = col.suppLiteral(cand.getLiteral())
            # pdb.set_trace()
            lparts = supports.lparts()
            lin = supports.lpartsInterX(supp)
            # print cand.getLiteral(), cand.isNeg()
            # print lparts, lin, len(supp)
            cand.setClp([lin, lparts], cand.isNeg())
            # print cand
            if colsC is not None and self.constraints.getCstr("add_condition"):
                ss = supports.copy()
                ss.update(side, cand.getOp(), supp)
                cand.setCondition(self.getCondition(colsC, ss))
        # print "================="
        return cands
        
    def getCandidates1(self, side, col, supports, currentRStatus=0):
        cands = []

        lparts = supports.lparts()
        lin = supports.lpartsInterX(col.supp())
        if side:
            fixed_colors = [[lparts[2], lparts[1]], [lparts[0], lparts[3]]]
            var_colors = [[lin[2], lin[1]], [lin[0], lin[3]]]
        else:
            fixed_colors = [[lparts[2], lparts[0]], [lparts[1], lparts[3]]]
            var_colors = [[lin[2], lin[0]], [lin[1], lin[3]]]

        for op in self.constraints.getCstr("allw_ops", side=side, currentRStatus=currentRStatus):
            for neg in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus):
                try:
                    adv = self.getAdv(side, op, neg, fixed_colors, var_colors[op], self.isCond(currentRStatus))
                except TypeError:
                    pdb.set_trace()
                    print side, op, neg, fixed_colors, var_colors[op]
                if adv is not None :
                    cands.append(Extension(self.constraints.getSSetts(), adv, None, (side, op, neg, Literal(neg, BoolTerm(col.getId())))))
        return cands

    def getCandidates2(self, side, col, supports, currentRStatus=0):
        return self.getCandidatesNonBool(side, col, supports, currentRStatus)

    def getCandidates3(self, side, col, supports, currentRStatus=0):
        return self.getCandidatesNonBool(side, col, supports, currentRStatus)

    def getCandidatesNonBool(self, side, col, supports, currentRStatus=0):
        cands = []
        lparts = supports.lparts()

        for cand in self.findCover(side, col, lparts, supports, currentRStatus):
            cands.append(cand[1])
        return cands

    ##################################################
    ###### CONDITIONAL
    def isCond(self, currentRStatus=0):
        return self.constraints.isStatusCond(currentRStatus)
    
    def getCondition(self, colsC, supports):
        cond_sparts = SParts(self.constraints.getSSetts(), supports.nbRows(), [supports.suppI(), supports.suppU()])
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
                cands = self.findCover(1, colsC[ci], lparts, cond_sparts, Constraints.getStatusCond())
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
    def findCover(self, side, col, lparts, supports, currentRStatus=0):
        method_string = 'self.findCover%i' % col.typeId()
        try:
            method_compute =  eval(method_string)
        except AttributeError:
              raise Exception('Oups No covering method for this type of data (%i)!'  % col.typeId())
        return method_compute(side, col, lparts, supports, currentRStatus)

    def findCover1(self, side, col, lparts, supports, currentRStatus=0):
        cands = []
        lin = supports.lpartsInterX(col.supp())
        if side:
            fixed_colors = [[lparts[2], lparts[1]], [lparts[0], lparts[3]]]
            var_colors = [[lin[2], lin[1]], [lin[0], lin[3]]]
        else:
            fixed_colors = [[lparts[2], lparts[0]], [lparts[1], lparts[3]]]
            var_colors = [[lin[2], lin[0]], [lin[1], lin[3]]]

        if self.constraints.getCstr("neg_query_init", side=side, currentRStatus=currentRStatus):
            if side:
                nfixed_colors = [[lparts[1], lparts[2]], [lparts[3], lparts[0]]]
                nvar_colors = [[lin[1], lin[2]], [lin[3], lin[0]]]
            else:
                nfixed_colors = [[lparts[0], lparts[2]], [lparts[3], lparts[1]]]
                nvar_colors = [[lin[0], lin[2]], [lin[3], lin[1]]]

        for op in self.constraints.getCstr("allw_ops", side=side, currentRStatus=currentRStatus):            
            for neg in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus):
                adv = self.getAdv(side, op, neg, fixed_colors, var_colors[op], self.isCond(currentRStatus))
                if adv is not None:
                    cands.append((False, Extension(self.constraints.getSSetts(), adv, None, (side, op, neg, Literal(neg, BoolTerm(col.getId()))))))

                ### to negate the other side when looking for initial pairs
                if self.constraints.getCstr("neg_query_init", side=side, currentRStatus=currentRStatus):
                    adv = self.getAdv(side, op, neg, nfixed_colors, nvar_colors[op], self.isCond(currentRStatus))
                    if adv is not None :
                        cand = Extension(self.constraints.getSSetts(), adv, None, (side, op, neg, Literal(neg, BoolTerm(col.getId()))))
                        cands.append((True, cand))  

        return cands

    def findCover2(self, side, col, lparts, supports, currentRStatus=0):
        cands = []
        allw_neg = True in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus)
        negs = self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus)
        if self.constraints.getCstr("multi_cats"): ## redundant negation
            negs = [False]

        if side:
            fixed_colors = [[lparts[2], lparts[1]], [lparts[0], lparts[3]]]
        else:
            fixed_colors = [[lparts[2], lparts[0]], [lparts[1], lparts[3]]]

        if self.constraints.getCstr("neg_query_init", side=side, currentRStatus=currentRStatus):
            if side:
                nfixed_colors = [[lparts[1], lparts[2]], [lparts[3], lparts[0]]]
            else:
                nfixed_colors = [[lparts[0], lparts[2]], [lparts[3], lparts[1]]]
                
        for op in self.constraints.getCstr("allw_ops", side=side, currentRStatus=currentRStatus):            
            for neg in negs:
                best = (None, None, None)
                bestNeg = (None, None, None)
                collect_goods = []
                collect_goodsNeg = []
                for (cat, supp) in col.iter_cats():
                    lin = supports.lpartsInterX(supp)
                    if side:
                        var_colors = [[lin[2], lin[1]], [lin[0], lin[3]]]
                    else:
                        var_colors = [[lin[2], lin[0]], [lin[1], lin[3]]]

                    # best = self.updateACTColors(best, , side, op, neg, fixed_colors, var_colors[op], self.isCond(currentRStatus))
                    ###################### 
                    tmp_adv = self.getAdv(side, op, neg, fixed_colors, var_colors[op], self.isCond(currentRStatus))
                    if best[0] < tmp_adv:
                        best = (tmp_adv, None, [side, op, neg, Literal(neg, CatTerm(col.getId(), cat))]) ## [fixed_colors, tuple(var_colors)]
                    if tmp_adv is not None:
                        collect_goods.append((tmp_adv, cat, fixed_colors, var_colors[op]))
                    ######################
                    
                    ### to negate the other side when looking for initial pairs
                    if self.constraints.getCstr("neg_query_init", side=side, currentRStatus=currentRStatus):
                        if side:
                            nvar_colors = [[lin[1], lin[2]], [lin[3], lin[0]]]
                        else:
                            nvar_colors = [[lin[0], lin[2]], [lin[3], lin[1]]]

                        # bestNeg = self.updateACTColors(bestNeg, Literal(neg, CatTerm(col.getId(), cat)), side, op, neg, nfixed_colors, nvar_colors[op])
                        ######################
                        tmp_adv = self.getAdv(side, op, neg, nfixed_colors, nvar_colors[op], self.isCond(currentRStatus))
                        if bestNeg[0] < tmp_adv:
                            bestNeg = (tmp_adv, None, [side, op, neg, Literal(neg, CatTerm(col.getId(), cat))]) ## [fixed_colors, tuple(var_colors)]
                        if tmp_adv is not None:
                            collect_goodsNeg.append((tmp_adv, cat, nfixed_colors, nvar_colors[op]))
                        ######################
                        
                if best[0] is not None:
                    bb = self.combCats(best, allw_neg, side, op, neg, col, collect_goods, currentRStatus=currentRStatus)
                    cands.append((False, Extension(self.constraints.getSSetts(), bb)))
                    ## cands.append((False, Extension(self.constraints.getSSetts(), best)))
                if bestNeg[0] is not None:
                    bb = self.combCats(best, allw_neg, side, op, neg, col, collect_goodsNeg, currentRStatus=currentRStatus)
                    cands.append((True, Extension(self.constraints.getSSetts(), bb)))
                    # cands.append((True, Extension(self.constraints.getSSetts(), bestNeg)))
        return cands

    def additionsLParts(self, lparts, lin, nextc, other_side=0):
        if other_side ==  0:
            ccum_lparts = lparts
            ccum_lin = [lin[0]+nextc[-1][0], lin[1]+nextc[-1][1]]
        elif other_side == 1:
            ccum_lparts = [[lparts[0][0]+nextc[-2][0][0], lparts[0][1]+nextc[-2][0][1]],
                           [lparts[1][0]+nextc[-2][1][0], lparts[1][1]-(nextc[-2][0][0]+nextc[-2][0][1]+nextc[-2][1][0])]]
            ccum_lin = [lin[0]+nextc[-1][0], lin[1]-nextc[-1][0]]
        else:
            ccum_lparts = [[lparts[0][0]+nextc[-2][0][0], lparts[0][1]+nextc[-2][0][1]],
                           [lparts[1][0]-(nextc[-2][0][0]+nextc[-2][0][1]+nextc[-2][1][1]), lparts[1][1]+nextc[-2][1][1]]]
            ccum_lin = [lin[0]-nextc[-1][1], lin[1]+nextc[-1][1]]
            
        # print "--- SUM", other_side
        # print lparts, lin
        # print nextc[-2], nextc[-1]
        # print ">>"
        # print ccum_lparts, ccum_lin
        return ccum_lparts, ccum_lin
    
    def combCats(self, best, allw_neg, side, op, neg, col, collected, other_side=0, currentRStatus=0):
        if not self.constraints.getCstr("multi_cats"):
            return best
        collected.sort()
        nextc = collected.pop()
        best_cat = [nextc[1]]
        best_score = nextc[0]
        cum_lparts, cum_lin = (list(nextc[-2]), list(nextc[-1]))
        while len(collected) > 0:
            nextc = collected.pop()
            ccum_lparts, ccum_lin = self.additionsLParts(cum_lparts, cum_lin, nextc, other_side)
            tmp_adv = self.getAdv(side, op, neg, ccum_lparts, ccum_lin, self.isCond(currentRStatus))
            if best_score < tmp_adv:
                best_cat.append(nextc[1])
                cum_lparts = ccum_lparts
                cum_lin = ccum_lin
                best_score = tmp_adv
        if len(best_cat) > 1:
            if col is None:
                best = (best_score, None, [side, op, neg, set(best_cat)])
            else:
                lit = col.makeCatLit(best_cat, neg, allw_neg)
                best = (best_score, None, [side, op, lit.isNeg(), lit])
        return best

    
    def findCover3(self, side, col, lparts, supports, currentRStatus=0):
        cands = []
        if self.constraints.isStatusCond(currentRStatus) or self.inSuppBounds(side, True, lparts) or self.inSuppBounds(side, False, lparts):  ### DOABLE
            segments = col.makeSegmentsColors(self.constraints.getSSetts(), side, supports, self.constraints.getCstr("allw_ops", side=side, currentRStatus=currentRStatus))
            if side:
                fixed_colors = [[lparts[2], lparts[1]], [lparts[0], lparts[3]]]
            else:
                fixed_colors = [[lparts[2], lparts[0]], [lparts[1], lparts[3]]]
 
            for cand in self.findCoverSegments(side, col, segments, fixed_colors, currentRStatus):
                cands.append((False, cand))

            ### to negate the other side when looking for initial pairs
            if self.constraints.getCstr("neg_query_init", side=side, currentRStatus=currentRStatus):
                if side:
                    nfixed_colors = [[lparts[1], lparts[2]], [lparts[3], lparts[0]]]
                else:
                    nfixed_colors = [[lparts[0], lparts[2]], [lparts[3], lparts[1]]]
                
                if self.inSuppBounds(side, True, lparts[-1::-1]): ### DOABLE
                    nsegments = [[[i, j, [k[1], k[0]]] for (i,j,k) in segments[0]], [[i, j, [k[1], k[0]]] for (i,j,k) in segments[1]]]
                    ##pdb.set_trace()
                    for cand in self.findCoverSegments(side, col, nsegments, nfixed_colors, currentRStatus):
                        cands.append((True, cand))
        return cands
        
    def findCoverSegments(self, side, col, segments, fixed_colors, currentRStatus=0):
        cands = []
        for op in self.constraints.getCstr("allw_ops", side=side, currentRStatus=currentRStatus): #TODO: limit this? Maybe or just add a different parameter/new method even
            if len(segments[op]) < self.constraints.getCstr("max_seg"):
                cands.extend(self.findCoverFullSearch(side, op, col, segments, fixed_colors, currentRStatus))
            else:
                if (False in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus)):
                    cands.extend(self.findPositiveCover(side, op, col, segments, fixed_colors, currentRStatus))
                if (True in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus)):
                    cands.extend(self.findNegativeCover(side, op, col, segments, fixed_colors, currentRStatus))
        return cands

    def findCoverFullSearch(self, side, op, col, segments, fixed_colors, currentRStatus=0):
        cands = []
        bests = {False: (None, None, None), True: (None, None, None)}

        for seg_s in range(len(segments[op])):
            var_colors = [0, 0] 
            for seg_e in range(seg_s,len(segments[op])):
                var_colors[0] += segments[op][seg_e][2][0]
                var_colors[1] += segments[op][seg_e][2][1]
                for neg in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus):
                    bests[neg] = self.updateACTColors(bests[neg], (seg_s, seg_e), side, op, neg, fixed_colors, var_colors, self.isCond(currentRStatus))

        for neg in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus):
            if bests[neg][0]:
                bests[neg][-1][-1] = col.getLiteralSeg(neg, segments[op], bests[neg][-1][-1])
                if bests[neg][-1][-1] is not None:
                    cands.append(Extension(self.constraints.getSSetts(), bests[neg]))
        return cands

    def findNegativeCover(self, side, op, col, segments, fixed_colors, currentRStatus=0):
        is_cond = self.isCond(currentRStatus)
        
        cands = []
        var_colors_f = [0, 0]
        # bests_f = [(0, 0, [0, 0])]
        bests_f = [(self.getAdv(side, op, False, fixed_colors, var_colors_f, is_cond=is_cond, no_const=True)[0], 0, [0, 0])]
        best_track_f = [0]
        var_colors_b = [0, 0]
        # bests_b = [(0, 0, [0, 0])]
        bests_b = [(self.getAdv(side, op, False, fixed_colors, var_colors_b, is_cond=is_cond, no_const=True)[0], 0, [0, 0])]
        best_track_b = [0]

        for i in range(len(segments[op])):
            # FORWARD
            var_colors_f[0] += segments[op][i][2][0]
            var_colors_f[1] += segments[op][i][2][1]
            if self.advRatioVar(var_colors_f, is_cond) > bests_f[-1][0]:
                var_colors_f[0] += bests_f[-1][2][0]
                var_colors_f[1] += bests_f[-1][2][1]

                bests_f.append((self.getAdv(side, op, False, fixed_colors, var_colors_f, is_cond=is_cond, no_const=True)[0], i+1, var_colors_f))
                var_colors_f = [0, 0]
            best_track_f.append(len(bests_f)-1)

            # BACKWARD
            var_colors_b[0] += segments[op][-(i+1)][2][0]
            var_colors_b[1] += segments[op][-(i+1)][2][1]

            if self.advRatioVar(var_colors_b, is_cond) > bests_b[-1][0]:
                var_colors_b[0] += bests_b[-1][2][0]
                var_colors_b[1] += bests_b[-1][2][1]

                bests_b.append((self.getAdv(side, op, False, fixed_colors, var_colors_b, is_cond=is_cond, no_const=True)[0], i+1, var_colors_b))
                var_colors_b = [0, 0]
            best_track_b.append(len(bests_b)-1)

        ##pdb.set_trace()
        best_t = (None, None, None)
        for b in bests_b:
            if b[1] == len(segments[op]):
                f = bests_f[0]
            else:
                f = bests_f[best_track_f[len(segments[op])-(b[1]+1)]]
            if self.advRatioVar(f[2], is_cond) > b[0]:
                best_t = self.updateACTColors(best_t, (f[1], len(segments[op]) - (b[1]+1)), side, op, False, fixed_colors, [f[2][0] + b[2][0], f[2][1] + b[2][1]], self.isCond(currentRStatus))
            else:
                best_t = self.updateACTColors(best_t, (0, len(segments[op]) - (b[1]+1)), side, op, False, fixed_colors, b[2], self.isCond(currentRStatus))

        for f in bests_f:
            if f[1] == len(segments[op]):
                b = bests_b[0]
            else:
                b = bests_b[best_track_b[len(segments[op])-(f[1]+1)]]
            if self.advRatioVar(b[2], is_cond) > f[0]: 
                best_t = self.updateACTColors(best_t, (f[1], len(segments[op]) - (b[1]+1)), side, op, False, fixed_colors, [f[2][0] + b[2][0], f[2][1] + b[2][1]], self.isCond(currentRStatus))
            else:
                best_t = self.updateACTColors(best_t, (f[1], len(segments[op])-1), side, op, False, fixed_colors, f[2], self.isCond(currentRStatus))

        if best_t[0] is not None:
            tmp = best_t[-1][-1] 
            best_t[-1][-1] = col.getLiteralSeg(True, segments[op], best_t[-1][-1])
            if best_t[-1][-1] is not None:
                cands.append(Extension(self.constraints.getSSetts(), best_t))
        return cands

    def findPositiveCover(self, side, op, col, segments, fixed_colors, currentRStatus=0):
        is_cond = self.isCond(currentRStatus)
        
        cands = []
        var_colors_f = [0.0, 0]
        nb_seg_f = 0
        best_f = (None, None, None)
        var_colors_b = [0.0, 0]
        nb_seg_b = 0
        best_b = (None, None, None)

        for i in range(len(segments[op])-1):
            # FORWARD
            if i > 0 and self.getAdv(side, op, False, fixed_colors, segments[op][i][2], is_cond=is_cond, no_const=True)[0] < self.advRatioVar(var_colors_f, is_cond):
                var_colors_f[0] += segments[op][i][2][0]
                var_colors_f[1] += segments[op][i][2][1]
                nb_seg_f += 1
            else:
                var_colors_f[0] = segments[op][i][2][0]
                var_colors_f[1] = segments[op][i][2][1]
                nb_seg_f = 0
                
            best_f = self.updateACTColors(best_f, (i - nb_seg_f, i), side, op, False, fixed_colors, var_colors_f, is_cond)

            # BACKWARD
            if i > 0 and self.getAdv(side, op, False, fixed_colors, segments[op][-(i+1)][2], is_cond=is_cond, no_const=True)[0] < self.advRatioVar(var_colors_b, is_cond):
                var_colors_b[0] += segments[op][-(i+1)][2][0]
                var_colors_b[1] += segments[op][-(i+1)][2][1]
                nb_seg_b += 1
            else:
                var_colors_b[0] = segments[op][-(i+1)][2][0]
                var_colors_b[1] = segments[op][-(i+1)][2][1]
                nb_seg_b = 0

            best_b = self.updateACTColors(best_b, (len(segments[op])-(1+i), len(segments[op])-(1+i) + nb_seg_b), \
                                    side, op, False, fixed_colors, var_colors_b, is_cond)

        if best_b[0] is not None and best_f[0] is not None:
            bests = [best_b, best_f]

            if best_b[-1][-1][0] > best_f[-1][-1][0] and best_b[-1][-1][1] > best_f[-1][-1][1] and best_b[-1][-1][0] <= best_f[-1][-1][1]:
                var_colors_m = [0,0]
                for seg in segments[op][best_b[-1][-1][0]:best_f[-1][-1][1]+1]:
                    var_colors_m[0] += seg[2][0]
                    var_colors_m[1] += seg[2][1]
                tmp_adv_m  = self.getAdv(side, op, False, fixed_colors, var_colors_m, is_cond)
                if tmp_adv_m is not None:
                    bests.append((tmp_adv_m, [fixed_colors, tuple(var_colors_m)], \
                                  [side, op, False, (best_b[-1][-1][0], best_f[-1][-1][1])]))

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
                suppL = colL.suppLiteral(pair["litL"])
                suppR = colR.suppLiteral(pair["litR"])
                cond_sparts = SParts(self.constraints.getSSetts(), colL.nbRows(), [suppL.intersection(suppR), suppL.union(suppR)])
                lparts = cond_sparts.lparts()        
                cond_cand = self.getConditionCand(colsC, cond_sparts, lparts)
                if cond_cand is not None:
                    pair["litC"] = cond_cand.getLiteral()
                    pair["scoreC"] = cond_cand.getAcc()
            pairs.append(pair)
        return pairs

    def doBoolStar(self, colL, colR, side):
        if side == 1:
            (supports, fixTerm, extCol) = (SParts(self.constraints.getSSetts(), colL.nbRows(), [colL.supp(), set()]), BoolTerm(colL.getId()), colR)
        else:
            (supports, fixTerm, extCol) = (SParts(self.constraints.getSSetts(), colL.nbRows(), [set(), colR.supp()]), BoolTerm(colR.getId()), colL)

        return self.fit(extCol, supports, side, fixTerm)

    def fit(self, col, supports, side, termX):
        (scores, literalsFix, literalsExt) = ([], [], [])   
        lparts = supports.lparts()
        currentRStatus = Constraints.getStatusPair(col, side, termX)
        cands = self.findCover(side, col, lparts, supports, currentRStatus=currentRStatus)
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
            supports = SParts(self.constraints.getSSetts(), colL.nbRows(), [set(), colR.nonModeSupp()])
            (scoresL, literalsFixL, literalsExtL) = self.fit(colL, supports, 0, idR)
            for tL in literalsExtL:
                suppL = colL.suppLiteral(tL)
                supports = SParts(self.constraints.getSSetts(), colL.nbRows(), [suppL, set()])
                (scoresR, literalsFixR, literalsExtR) = self.fit(colR, supports, 1, tL)
                for i in range(len(scoresR)):
                    if scoresR[i] > bestScore:
                        (scores, literalsL, literalsR) = ([scoresR[i]], [literalsFixR[i]], [literalsExtR[i]])
                        bestScore = scoresR[i]
                        
            ## FIT RHS then LHS
            supports = SParts(self.constraints.getSSetts(), colL.nbRows(), [colL.nonModeSupp(), set()])
            (scoresR, literalsFixR, literalsExtR) = self.fit(colR, supports, 1, idL)
            for tR in literalsExtR:
                suppR = colR.suppLiteral(tR)
                supports = SParts(self.constraints.getSSetts(), colL.nbRows(), [set(), suppR])
                (scoresL, literalsFixL, literalsExtL) = self.fit(colL, supports, 0, tR)
                for i in range(len(scoresL)):
                    if scoresL[i] > bestScore:
                        (scores, literalsL, literalsR) = ([scoresR[i]], [literalsExtL[i]], [literalsFixL[i]])
                        bestScore = scoresL[i]
                        
#             if len(scores) > 0:
#                print "%f: %s <-> %s" % (scores[0], literalsA[0], literalsB[0])
        return (scores, literalsL, literalsR)

    def subdo33Full(self, colL, colR, side):
        ## print "SUBDO 33 FULL", colL.getId(), colR.getId()
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
                    
        ## print "buckets lengths\t(0,%d) %d\t(1,%d) %d\tcollapsed %d -- product %d" % (colL.id, len(bucketsL[1]), colR.id, len(bucketsR[1]), len(bucketsE[1]), len(bucketsF[1]) * len(bucketsE[1]))
        if bucketsE is not None and ( len(bucketsF[1]) * len(bucketsE[1]) < self.constraints.getCstr("max_prodbuckets") ):
            ## print "Trying buckets...", len(bucketsF[0]), len(bucketsE[0])
            totInt = colE.nbRows()
            #margE = [len(intE) for intE in bucketsE[0]]
            
            margF = [len(bucketsF[0][i]) for i in range(len(bucketsF[0]))]
            
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
                        #totEinF = sum(EinF)
                        
                        fixed_colors = [[0, 0],[totInt - aboveF - belowF, aboveF + belowF]]

                        belowEF = 0
                        outBelowEF = 0
                        lowE = 0
                        while lowE < len(interMat[lowF]) and totInt - belowF - aboveF - belowEF >= self.constraints.getCstr("min_itm_in"):
                            aboveEF = 0
                            outAboveEF = 0
                            upE = len(interMat[lowF])-1
                            while upE >= lowE and totInt - belowF - aboveF - belowEF - aboveEF >= self.constraints.getCstr("min_itm_in"):
                                var_colors = [totInt - belowF - aboveF - belowEF - aboveEF, belowF + aboveF - outAboveEF - outBelowEF]
                                best = self.updateACTColorsP33(best, (lowF, upF, lowE, upE), side, True, False, fixed_colors, var_colors)
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
        ##### THIS NEEDS CHANGE PARTS
        configs = [(0, False, False), (1, False, True), (2, True, False), (3, True, True)]
        allw_neg = True
        if True not in self.constraints.getCstr("neg_query", side=side, type_id=2):
            configs = configs[:1]
            allw_neg = False
        best = [[] for c in configs]
        tot = colL.nbRows()        
        for catL in colL.cats():
            totL = len(colL.suppCat(catL))
            lparts = [0, totL, 0, colL.nbRows()-totL]
            
            for catR in colR.cats():
                totR = len(colR.suppCat(catR))
                interLR = len(colL.suppCat(catL) & colR.suppCat(catR))
                lin = [0, interLR, 0, totR - interLR]
 
                for (i, nL, nR) in configs:
                    if nL:
                        fixed_colors = [[lparts[0], lparts[2]], [lparts[3], lparts[1]]]
                        var_colors = [lin[3], lin[1]]                            
                    else:
                        fixed_colors = [[lparts[2], lparts[0]], [lparts[1], lparts[3]]]
                        var_colors = [lin[1], lin[3]]                            
                    
                    best[i] = self.updateACTColorsP22(best[i], (catL, catR), 0, True, nR, fixed_colors, var_colors)

        (scores, literalsFix, literalsExt) = ([], [], [])
        if self.constraints.getCstr("multi_cats"):
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
                    map_cat[ss][b[-1][-1][ss]].append((b[0][0], b[-1][-1][1-ss], b[1][0], b[1][1]))

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
        ##### THIS NEEDS CHANGE PARTS
        if side == 0:
            (colF, colE) = (colR, colL)
        else:
            (colF, colE) = (colL, colR)

        configs = [(0, False, False), (1, False, True), (2, True, False), (3, True, True)]
        allw_neg = True
        if True not in self.constraints.getCstr("neg_query", side=side, type_id=3):
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

            marg = [len(buk) for buk in buckets[0]]
            if buckets[2] is not None :
                marg[buckets[2]] += colE.lenMode()

            for cat in colF.cats():
                totF = len(colF.suppCat(cat))
                lparts = [0, totF, 0, colF.nbRows() - totF]

                interMat = [len(colF.suppCat(cat) & buk) for buk in buckets[0]]
                if buckets[2] is not None :
                    interMat[buckets[2]] += len(colE.interMode(colF.suppCat(cat)))        

                totIn = sum(interMat)
                below = 0
                low = 0
                while low < len(interMat) and \
                          (totIn - below >= self.constraints.getCstr("min_itm_in") or totIn - below >= self.constraints.getCstr("min_itm_out")):
                    above = 0
                    up = len(interMat)-1

                    while up >= low and \
                              (totIn - below - above >= self.constraints.getCstr("min_itm_in") or totIn -below-above >= self.constraints.getCstr("min_itm_out")):
                        pin = totIn - below - above

                        lin = [0, pin, 0, sum(marg[low:up+1]) - pin]
                        for (i, nF, nE) in configs:
                            if nF:
                                fixed_colors = [[lparts[0], lparts[2]], [lparts[3], lparts[1]]]
                                var_colors = [lin[3], lin[1]]                            
                            else:
                                fixed_colors = [[lparts[2], lparts[0]], [lparts[1], lparts[3]]]
                                var_colors = [lin[1], lin[3]]                            

                        # for (i, nF, nE) in configs:
                        #     if nF and nE:
                        #         fixed_colors = [[lparts[3], lparts[1]], [lparts[0], lparts[2]]]
                        #         var_colors = [lin[0], lin[2]]                            
                        #     elif nE:
                        #         fixed_colors = [[lparts[1], lparts[3]], [lparts[2], lparts[0]]]
                        #         var_colors = [lin[2], lin[0]]
                        #     elif nF:
                        #         fixed_colors = [[lparts[0], lparts[2]], [lparts[3], lparts[1]]]
                        #         var_colors = [lin[3], lin[1]]                            
                        #     else:
                        #         fixed_colors = [[lparts[2], lparts[0]], [lparts[1], lparts[3]]]
                        #         var_colors = [lin[1], lin[3]]                            

                            best[i] = self.updateACTColorsP23(best[i], (cat, low, up), side, True, nE, fixed_colors, var_colors)
                            
                        above+=interMat[up]
                        up-=1
                    below+=interMat[low]
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
                if tE is not None:
                    literalsExt.append(tE)                    
                    literalsFix.append(colF.makeCatLit(b[-1][-1][0], nF))
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
                map_cat[buk].append((b[0][0], b[-1][-1][0], b[1][0], b[1][1]))

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

    def advRatioVar(self, var_colors, is_cond=False):
        if is_cond:
            return tool_ratio(var_colors[0], var_colors[0]+var_colors[1])
        return tool_ratio(var_colors[0], var_colors[1])
                
    def getAdv(self, side, op, neg, fixed_colors, var_colors, is_cond=False, no_const=False):
        fix_num = None
        if neg:
            tmp_var = [fixed_colors[op][i] - var_colors[i] for i in range(len(var_colors))]
        else:
            # tmp_var = [var_colors[i] for i in range(len(var_colors))]
            tmp_var = var_colors

        if is_cond:
            ### SPECIAL AND
            # if fixed_colors[op][1] - tmp_var[1] >= self.constraints.getCstr("min_itm_c") \
            #    and tmp_var[0] >= self.constraints.getCstr("min_itm_in") \
            #    and fixed_colors[1-op][1] + fixed_colors[op][1] - tmp_var[1] >= self.constraints.getCstr("min_itm_out"):
            if no_const or (tmp_var[0] >= self.constraints.getCstr("min_itm_in")):
                contri = tmp_var[0]
                fix_num, fix_den = (0, 0)
                # if self.constraints.getCstr("constraint_score") == "comp_acc":
                #     ### COMP ACC
                #     var_num = float(tmp_var[0])/(tmp_var[1] + tmp_var[0])
                #     tt = (fixed_colors[False][1]-tmp_var[1] + fixed_colors[False][0]-tmp_var[0])
                #     if tt > 0:
                #         var_den = float(fixed_colors[False][0]-tmp_var[0])/tt
                #     else:
                #         var_den = 0
                                            
                # else:
                if True:
                    #### ACCURACY
                    var_num, var_den = (tmp_var[0], tmp_var[1] + tmp_var[0])
                # print "PIECES", var_num, var_den, contri, fix_num, fix_den
                # print "%d/%d : %d/%d" % (tmp_var[0], tmp_var[1] + tmp_var[0], fixed_colors[False][0]-tmp_var[0], tt)
            
        elif op:
            #OR
            if no_const or (tmp_var[0] >= self.constraints.getCstr("min_itm_c") \
               and fixed_colors[op][1] - tmp_var[1] >= self.constraints.getCstr("min_itm_out") \
               and fixed_colors[1-op][0] + tmp_var[0] >= self.constraints.getCstr("min_itm_in")):
                contri = tmp_var[0]
                fix_num = fixed_colors[1-op][0]
                var_num = tmp_var[0]
                fix_den = fixed_colors[1-op][0] + fixed_colors[1-op][1] + fixed_colors[op][0]
                var_den = tmp_var[1]
                # sout = fixed_colors[op][1] - tmp_var[1]
                # print "PIECES", sout, var_num, var_den, contri, fix_num, fix_den
        else:
            # AND
            if no_const or (fixed_colors[op][1] - tmp_var[1] >= self.constraints.getCstr("min_itm_c") \
               and tmp_var[0] >= self.constraints.getCstr("min_itm_in") \
               and fixed_colors[1-op][1] + fixed_colors[op][1] - tmp_var[1] >= self.constraints.getCstr("min_itm_out")):
                contri = fixed_colors[op][1] - tmp_var[1]
                fix_num = 0
                var_num = tmp_var[0]
                fix_den = fixed_colors[op][0] + fixed_colors[1-op][0]
                var_den = tmp_var[1]
                # sout = fixed_colors[1-op][1] + fixed_colors[op][1] - tmp_var[1]
                # print "PIECES", sout, var_num, var_den, contri, fix_num, fix_den

        if fix_num is not None:
            acc = tool_ratio(fix_num + var_num, fix_den + var_den)
            return (acc, var_num, var_den, contri, fix_num, fix_den)
        return None
                
    # def getClp(self, side, op, neg, fixed_colors, var_colors):
    #     lparts = [fixed_colors[0][1], fixed_colors[1][0], fixed_colors[0][0], fixed_colors[1][1]]
    #     if op:
    #         lin = [0, var_colors[0], 0, var_colors[1]]
    #     else:
    #         lin = [var_colors[1], 0, var_colors[0], 0]

    #     if side == 1:
    #         lin[0], lin[1] = lin[1], lin[0]
    #         lparts[0], lparts[1] = lparts[1], lparts[0] 

    #     lout = [lparts[i] - lin[i] for i in range(len(lparts))]
    #     return [lin, lout, lparts]

    # def getBestClp(self, bestColors):
    #     return bestColors[0],  self.getClp(bestColors[2][0], bestColors[2][1], bestColors[2][2], bestColors[1][0], bestColors[1][1]), bestColors[2] 
        
    def updateACTColors(self, best, lit, side, op, neg, fixed_colors, var_colors, is_cond=False):
        tmp_adv = self.getAdv(side, op, neg, fixed_colors, var_colors, is_cond)
        if best[0] < tmp_adv:
            return tmp_adv, None, [side, op, neg, lit] ## [fixed_colors, tuple(var_colors)]
        else:
            return best

    def updateACTColorsP(self, best, lit, side, op, neg, fixed_colors, var_colors, conflictF):
        tmp_adv = self.getAdv(side, op, neg, fixed_colors, var_colors)
        if tmp_adv is None:
            return best
        inserted = False
        i = 0
        while i < len(best):
            if best[i][0] > tmp_adv:
                if conflictF(best[i][-1][-1], lit):  ## Best already contains conflicting of better quality 
                    return best
            else:
                if not inserted: 
                    #best.insert(i,(tmp_adv, None, [side, op, neg, lit]))
                    best.insert(i,(tmp_adv, (tuple(fixed_colors), tuple(var_colors)), [side, op, neg, lit]))
                    inserted = True
                elif conflictF(best[i][-1][-1], lit): ## Best contains conflicting of lesser quality than inserted, remove
                    best.pop(i)
                    i -=1
            i+=1
        if not inserted:
            best.append((tmp_adv, (tuple(fixed_colors), tuple(var_colors)), [side, op, neg, lit]))
            #best.append((tmp_adv, None, [side, op, neg, lit]))
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

    def updateACTColorsP22(self, best, lit, side, op, neg, fixed_colors, var_colors):
        return self.updateACTColorsP(best, lit, side, op, neg, fixed_colors, var_colors, self.conflictP22)            
    def updateACTColorsP23(self, best, lit, side, op, neg, fixed_colors, var_colors):
        return self.updateACTColorsP(best, lit, side, op, neg, fixed_colors, var_colors, self.conflictP23)
    def updateACTColorsP33(self, best, lit, side, op, neg, fixed_colors, var_colors):
        return self.updateACTColorsP(best, lit, side, op, neg, fixed_colors, var_colors, self.conflictP33)

### Reintegrate as well as linter
    def inSuppBounds(self, side, op, lparts):
        return self.constraints.getSSetts().sumPartsId(side, self.constraints.getSSetts().IDS_varnum[op] + self.constraints.getSSetts().IDS_fixnum[op], lparts) >= self.constraints.getCstr("min_itm_in") \
               and self.constraints.getSSetts().sumPartsId(side, self.constraints.getSSetts().IDS_cont[op], lparts) >= self.constraints.getCstr("min_itm_c")
