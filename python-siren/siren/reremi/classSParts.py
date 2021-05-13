##### RENAMED PARTS
## grep -E "gamma|beta|alpha|delta|mua|mub|muaB|mubB|mud" ../reremi/*.py
## sed -i -e 's/alpha/Exo/g' -e 's/beta/Eox/g' -e 's/gamma/Exx/g' -e 's/delta/Eoo/g' -e 's/mua/Exm/g' -e 's/mub/Emx/g' -e 's/muaB/Eom/g' -e 's/mubB/Emo/g' -e 's/mud/Emm/g' ../reremi/*.py
## grep -E "SSetts.gamma|SSetts.beta|SSetts.alpha|SSetts.delta|SSetts.mua|SSetts.mub|SSetts.muaB|SSetts.mubB|SSetts.mud" ../views/*.py
## sed -i -e 's/SSetts.alpha/SSetts.Exo/g' -e 's/SSetts.beta/SSetts.Eox/g' -e 's/SSetts.gamma/SSetts.Exx/g' -e 's/SSetts.delta/SSetts.Eoo/g' -e 's/SSetts.mua/SSetts.Exm/g' -e 's/SSetts.mub/SSetts.Emx/g' -e 's/SSetts.muaB/SSetts.Eom/g' -e 's/SSetts.mubB/SSetts.Emo/g' -e 's/SSetts.mud/SSetts.Emm/g' ../views/*.py
## sed -i -e 's/ExmB/Eom/g' -e 's/EmxB/Emo/g' ../*/*.py

from scipy.special import gammaln
from scipy.stats import binom
import numpy, random, re
import pdb

def tool_ratio(num, den):
    if num is None or den is None:
        return None
    if den == 0:
        if num > 0:
            return float("Inf")
        else:
            return 0.
    else:
        return float(num)/den

def tool_hypergeomPMF(k, M, n, N):
    tot, good = M, n
    bad = tot - good
    return numpy.exp(gammaln(good+1) - gammaln(good-k+1) - gammaln(k+1) + gammaln(bad+1) \
                              - gammaln(bad-N+k+1) - gammaln(N-k+1) - gammaln(tot+1) \
                              + gammaln(tot-N+1) + gammaln(N+1))
#same as the following but numerically more precise
#return comb(good,k) * comb(bad,N-k) / comb(tot,N)

def tool_pValOver(kInter, nbRows, suppL, suppR, lU=None):
    ## probability that two sets of these size have intersection equal or larger than kInter
    if lU is None:
        return sum([ tool_hypergeomPMF(k, nbRows, suppL, suppR) for k in range(kInter, min(suppL, suppR)+1)])
    else:
        return sum([ tool_hypergeomPMF(k, nbRows, suppL, suppR) for k in range(0, lU - kInter+1)])
        
def tool_pValSupp(nbRows, supp, pr, lU=None):
    ## probability that an itemset with marginal probability pr has support equal or larger than supp
    if lU is None:
        return 1-binom.cdf(supp-1, nbRows, pr)
    else:
        return binom.cdf(lU-supp, nbRows, pr)         


class SSetts(object):

    # labels = ['\t|  \n', '\t  |\n', '\t| |\n', '\t   \n', '\t| :\n', '\t: |\n', '\t  :\n', '\t:  \n', '\t: :\n' ]
    # labels = ['**', '__', '==', '  ', '*.', '"_', '..', '""', '::' ]

    status = [(True, False), (False, True), (True, True), (False, False),
              (True, None), (None, True), (False, None), (None, False), (None, None)]
    labels_status = {True:"1", False:"0", None:"?"}
    labelsu_status = {True:"x", False:"o", None:"?"}
    labelsm_status = {True:"x", False:"o", None:"m"}
    labels = ["E%s%s" % (labelsm_status[slhs], labelsm_status[srhs]) for (slhs, srhs) in status]
    # labelsu_status = {True:u"\u2081", False:u"\u2080", None:u"\u2098"}
    labels_sparts = ["E%s%s" % (labels_status[slhs], labels_status[srhs]) for (slhs, srhs) in status]
    ## WITH UNICODE
    sym_status = labelsu_status
    # ## WITHOUT UNICODE
    # sym_status = labels_status
    
    sym_sparts = [sym_status[slhs]+sym_status[srhs] for (slhs, srhs) in status]
    
    ### define the numerical values for the parts
    for i, l in enumerate(labels):
        exec("%s = %d" % (l,i))

    io_labels = ["into", "out", "tot", "imiss"]
    for i, l in enumerate(io_labels):
        exec("%s = %d" % (l,i))

    map_label_part = dict([(s,p) for (p,s) in enumerate(labels)])
    map_status_part = dict([(s,p) for (p,s) in enumerate(status)])
    @classmethod
    def mapStatusToSPart(tcl, status):
        return tcl.map_status_part.get(status, -1)

    @classmethod
    def mapStatusToSPart(tcl, status):
        return tcl.map_status_part.get(status, -1)

        
    # # indexes of the parts
    # (Exo, Eox, Exx, Eoo, Exm, Emx, Eom, Emo, Emm) = range(9)
    # (into, out, tot, imiss) = range(4)

    # (self.Exo, self.Eox, self.Exx, self.Eoo, self.Exm, self.Emx, self.Eom, self.Emo, self.Emm) = range(9)
    # (self.into, self.out, self.tot, self.imiss) = range(4)


    ## TRUTH TABLE:
    ## A B    OR    AND
    ## T T    T     T
    ## T F    T     F
    ## T M    T     M
    ## F F    F     F
    ## F M    M     F
    ## M M    M     M

    ## PARTS:
    ##        A  |  B
    ## ----------------
    ## Exo   T  |  F
    ## Eox   F  |  T
    ## Exx   T  |  T
    ## Eoo   F  |  F
    ## Exm   T  |  M
    ## Emx   M  |  T
    ## Eom   F  |  M
    ## Emo   M  |  F
    ## Emm   M  |  M

    ## EXTENDING: A op X
    ##       |     op = OR    |     op = AND   |
    ##       |   x    o    m  |   x    o    m  |
    ## -----------------------------------------
    ## Exo   |  ---- xo ----  |  xo   oo   mo  |
    ## Eox   |  xx   ox   mx  |  ---- ox ----  |
    ## Exx   |  ---- xx ----  |  xx   ox   mx  |
    ## Eoo   |  xo   oo   mo  |  ---- oo ----  |
    ## Exm   |  ---- xm ----  |  xm   om   mm  |
    ## Emx   |  xx   mx   mx  |  mx   ox   mx  |
    ## Eom   |  xm   om   mm  |  ---- om ----  |
    ## Emo   |  xo   mo   mo  |  mo   oo   mo  |
    ## Emm   |  xm   mm   mm  |  mm   om   mm  |

    

    def __init__(self, type_parts="none", methodpVal="Marg"):
        self.type_parts = None
        self.resetPartsIds(type_parts)
        self.setMethodPVal(methodpVal)

    def getTypeParts(self):
        return self.type_parts
    def getMethodPVal(self):
        return self.methodpVal

    def reset(self, type_parts=None, methodpVal=None):
        if type_parts is not None:
            self.resetPartsIds(type_parts)
        if methodpVal is not None:
            self.setMethodPVal(methodpVal)
        
    def resetPartsIds(self, type_parts):
        if type_parts is True:
            type_parts = "rejective"
        if type_parts is False:
            type_parts = "none"
            
        if self.type_parts == type_parts:
            return
        elif self.type_parts == "none" and type_parts != "exclu":
            return

        self.type_parts = type_parts
        
        # indexes from the parts when looking from the right (A=L, B=R) or the left (A=R,B=L) 
        self.side_index = [[0,1,2,3,4,5,6,7,8], [1,0,2,3,5,4,7,6,8]]

        # indexes for the intersections with parts
        # (into: part inter X_True, out: part inter X_False, miss: part inter X_Missing, tot: total part = into + out + miss)
        # indexed for the intersections with parts when considering positive or negative X
        self.neg_index = [[0, 1, 2, 3], [1, 0, 2, 3]]

        ## and same for initializing relative parts
        for i, l in enumerate(self.labels):
            exec("%s = %d" % (l,i))
        for i, l in enumerate(self.io_labels):
            exec("%s = %d" % (l,i))

        # (Exo, Eox, Exx, Eoo, Exm, Emx, Eom, Emo, Emm) = range(9)
        self.last_nonmiss = Eoo

        if type_parts == "none":

            #####################################################################################
            ###############             WITHOUT MISSING VALUES                      #############
            #####################################################################################

            self.bottom = Exo
            self.top = Eoo
            ##############################################################
            ####         J=       |Exx| /
            ####              |Exx|+|Exo|+|Eox|
            ##############################################################

            ##### TO COMPUTE ADVANCE while building, INDEXED BY OPERATOR (0: AND, 1: OR)
            # Parts in numerator (BLUE), independent of X 
            self.IDS_fixnum = [[], [(tot, Exx)]]
            # Parts in numerator (BLUE), dependent of X: toBlue
            self.IDS_varnum = [[(into, Exx)] ,[(into, Eox)]]
            # Parts in denominator (RED), independent of X
            self.IDS_fixden = [[(tot, Eox), (tot, Exx)], [(tot, Exx), (tot, Exo), (tot, Eox)]]
            # Parts in denominator (RED), dependent of X: toRed
            self.IDS_varden = [[(into, Exo)], [(into, Eoo)]]
            # Parts left uncovered (OUT), (always dependent of X)
            self.IDS_out = [[(out, Exo), (tot, Eoo)], [(out, Eoo)]]
            # Parts in contribution (CONT), (always dependent of X)
            # Contribution: AND entities removed from Exo, OR: entities added to Exx
            self.IDS_cont = [[(out, Exo)], [(into, Eox)]]
            # Parts in the new support of the extended query
            self.IDS_nsupp = [[(into, Exo), (into, Exx)], [(tot, Exo), (tot, Exx), (into, Eox), (into, Eoo)]]

            #### TO COMPUTE ACCURACY after building
            self.IDS_dL = [Exo]
            self.IDS_dR = [Eox]
            self.IDS_diff = list(self.IDS_dL + self.IDS_dR)
            self.IDS_inter = [Exx]
            self.IDS_uncovered = [Eoo]

            ##############################################################

            #### TO COMPUTE SUPPORTS, no index
            self.IDS_supp = (Exx, Exo)
            self.IDS_miss = ()
            # indexes swaping when negating one side (0: negating A, 1: negating B)
            self.IDS_negated = [(Eoo, Exx, Eox, Exo), \
                                (Exx, Eoo, Exo, Eox)]

            #### END NO MISSING VALUES
            #####################################################################################


        else:

            #####################################################################################
            ###############                WITH MISSING VALUES                      #############
            #####################################################################################
            self.bottom = Exo
            self.top = Emm
            ##############################################################
            #### REJECTIVE  J=       |Exx| /
            ####              |Exx|+|Exo|+|Eox|
            ##############################################################

            if type_parts == "rejective" or type_parts == "grounded":

                ##### TO COMPUTE ADVANCE while building, INDEXED BY OPERATOR (0: AND, 1: OR)
                # Parts in numerator (BLUE), independent of X 
                self.IDS_fixnum = [[],
                              [(tot, Exx)]]
                # Parts in numerator (BLUE), dependent of X: toBlue
                self.IDS_varnum = [[(into, Exx)] ,
                              [(into, Eox), (into, Emx)]]
                # Parts in denominator (RED), independent of X
                self.IDS_fixden = [[(tot, Eox)],
                              [(tot, Exx), (tot, Exo)]]
                # Parts in denominator (RED), dependent of X: toRed
                self.IDS_varden = [[(into, Exo), (into, Exx), (out, Exx), (out, Emx)],
                              [(into, Emx), (into, Emo), (into, Eoo), (into, Eox), (out, Eox)]]

                # Parts left uncovered (OUT), (always dependent of X)
                self.IDS_out = [[(out, Exo), (out, Emo), (tot, Eoo)],
                           [(out, Eoo)]]
                # Parts in contribution (CONT), (always dependent of X)
                # Contribution: AND entities removed from Exo, OR: entities added to Exx
                self.IDS_cont = [[(out, Exo)],
                            [(into, Eox), (into, Emx)]]
                # Parts in the new support of the extended query
                self.IDS_nsupp = [[(into, Exo), (into, Exx)],
                             [(tot, Exo), (tot, Exx), (tot, Exm), (into, Emx),
                                  (into, Eox), (into, Eoo), (into, Emo), (into, Eom), (into, Emm)]]

                #### TO COMPUTE ACCURACY after building
                self.IDS_dL = [Exo]
                self.IDS_dR = [Eox]
                self.IDS_diff = list(self.IDS_dL + self.IDS_dR)                
                self.IDS_inter = [Exx]
                self.IDS_uncovered = [Eoo]


            ##############################################################
            #### OPTIMISTIC    J= |Exx|+|Exm|+|Emx|+|Emm| /
            ####            |Exo|+|Eox|+|Exx|+|Exm|+|Emx|+|Emm|
            ##############################################################

            elif type_parts == "optimistic":

                ##### TO COMPUTE ADVANCE while building, INDEXED BY OPERATOR (0: AND, 1: OR)
                self.IDS_fixnum = [[(imiss, Exm), (imiss, Emx), (imiss, Emm), (imiss, Exx)],
                                     [(tot, Exm), (tot, Exx), (tot, Emm), (tot, Emx), (imiss, Eom), (imiss, Eox)]]
                self.IDS_varnum = [[(into, Exm), (into, Exx), (into, Emx), (into, Emm)],
                                     [(into, Eom), (into, Eox), (imiss, Eom), (imiss, Eox)]]
                self.IDS_fixden = [[(tot, Exx), (tot, Emx), (tot, Eox)],
                                     [(tot, Exm), (tot, Exx), (tot, Emx), (tot,  Emm), (imiss, Eom), (tot, Eox), (tot, Exo)]]
                self.IDS_varden = [[(into, Exo), (into, Exm), (into, Emm), (imiss, Exm), (imiss, Emm)],
                                     [(into, Eoo), (into, Emo), (into, Eom)]]

                self.IDS_out = [[(out, Exo), (imiss, Exo), (out, Exm), (tot, Eoo), (tot, Emo), (out, Emm), (tot, Eom)],
                                  [(out, Eoo), (imiss, Eoo), (out, Emo), (imiss,Emo), (out,Eom)]]
                self.IDS_cont = [[(out, Exo)],
                                   [(into, Eox), (into, Emx)]]
                self.IDS_nsupp = [[(into, Exo), (into, Exm), (into, Exx), (into, Emx), (into, Emm)],
                                    [(tot, Exo), (tot, Exx), (tot, Exm), (into, Emx), (into, Eox), (into, Eoo), (into, Emo), (into, Emm), (into, Eom)]]

                #### TO COMPUTE ACCURACY after building
                self.IDS_dL = [Exo]
                self.IDS_dR = [Eox]
                self.IDS_diff = list(self.IDS_dL + self.IDS_dR)                
                self.IDS_inter = [Exx, Emx, Exm, Emm]
                self.IDS_uncovered = [Eoo, Emo, Eom]


            ##############################################################
            #### PESSIMISTIC    J=          |Exx|   /
            ####            |Exo|+|Eox|+|Exx|+|Exm|+|Emx|+|Eom|+|Emo|+|Emm|
            ##############################################################
            #--- corrected (oct.17) from    J=          |Exx|   /
            #---            |Exo|+|Eox|+|Exx|+|Eom|+|Emo|+|Emm|
            #--- the varibles below were right, the formula above wasn't, and didn't match...

            elif type_parts == "pessimistic":

                ##### TO COMPUTE ADVANCE while building, INDEXED BY OPERATOR (0: AND, 1: OR)
                self.IDS_fixnum = [[],
                                     [(tot, Exx)]]
                self.IDS_varnum = [[(into, Exx)] ,
                                     [(into, Eox), (into, Emx)]]
                self.IDS_fixden = [[(imiss, Exo), (tot, Exm), (tot, Exx), (tot, Emx), (tot, Eox), (imiss, Emo), (tot, Emm), (tot, Eom)],
                                     [(tot, Exx), (tot, Exo), (tot, Eox), (tot, Exm), (tot, Eom), (tot, Emx), (tot, Emo), (tot, Emm)]]
                self.IDS_varden = [[(into, Exo), (into, Emo)],
                                     [(into, Eoo), (imiss, Eoo)]]

                self.IDS_out = [[(out, Exo), (out, Emo), (tot, Eoo)],
                                  [(out, Eoo)]]
                self.IDS_cont = [[(out, Exo)],
                                   [(into, Eox), (into, Emx)]]
                self.IDS_nsupp = [[(into, Exo), (into, Exx), (into, Exm)],
                                    [(tot, Exo), (tot, Exx), (tot, Exm), (into, Emx), (into, Eox), (into, Eoo), (into, Emo), (into, Emm)]]

                #### TO COMPUTE ACCURACY after building
                self.IDS_dL = [Exo, Exm, Emo]
                self.IDS_dR = [Eox, Emx, Eom, Emm]
                self.IDS_diff = list(self.IDS_dL + self.IDS_dR)
                self.IDS_inter =  [Exx]
                self.IDS_uncovered = [Eoo]

            ##############################################################
            #### BASIC     J=              |Exx|  /
            ####               |Exo|+|Eox|+|Exx|+|Exm|+|Emx|
            ##############################################################

            elif type_parts == "basic":

                ##### TO COMPUTE ADVANCE while building, INDEXED BY OPERATOR (0: AND, 1: OR)
                self.IDS_fixnum = [[], [(tot, Exx)]]
                self.IDS_varnum = [[(into, Exx)] ,
                                   [(into, Eox), (into, Emx)]]
                self.IDS_fixden = [[(tot, Exx), (tot, Emx), (tot, Eox)],
                                   [(tot, Exo), (tot, Eox), (tot, Exx), (tot, Exm), (tot, Emx)]]
                self.IDS_varden = [[(into, Exo), (into, Exm)],
                                   [(into, Eoo), (into, Emo), (into, Eom)]]

                self.IDS_out = [[(out, Exo), (imiss, Exo), (tot, Eoo), (out, Exm), (imiss, Exm), (tot, Eom), (tot, Emo), (out, Emm)],
                                [(out, Eoo), (imiss, Eoo), (out, Emo), (imiss, Emo), (out, Eom), (imiss, Eom), (out, Emm), (imiss, Emm)]]
                self.IDS_cont = [[(out, Exo)],
                                 [(into, Eox), (into, Emx)]]
                self.IDS_nsupp = [[(into, Exo), (into, Exx), (into, Exm)],
                                  [(tot, Exo), (tot, Exx), (tot, Exm), (into, Eox), (into, Eoo), (into, Emx), (into, Emo), (into, Eom), (into, Emm)]]

                #### TO COMPUTE ACCURACY after building
                self.IDS_dL = [Exo, Exm]
                self.IDS_dR = [Eox, Emx]
                self.IDS_diff = list(self.IDS_dL + self.IDS_dR)
                self.IDS_inter = [Exx]
                self.IDS_uncovered = [Eoo, Emo, Eom, Emm]

            ##############################################################
            #### TEST     J=      |Exo|+|Eox|  /
            ####               |Exo|+|Eox|+|Exx|
            ##############################################################

            elif type_parts == "exclu":

                ##### TO COMPUTE ADVANCE while building, INDEXED BY OPERATOR (0: AND, 1: OR)
                self.IDS_fixnum = [[(tot, Eox)],
                                   [(tot, Exo)]]
                self.IDS_varnum = [[(into, Exo), (out, Exx), (out, Emx)],
                                   [(into, Eoo), (out, Eox), (into, Emo)]]
                
                self.IDS_fixden = [[(tot, Eox)],
                                   [(tot, Exx), (tot, Exo)]]
                self.IDS_varden = [[(into, Exo), (into, Exx), (out, Exx), (out, Emx)],
                                   [(into, Eox), (out, Eox), (into, Eoo), (into, Emx), (into, Emo)]]

                self.IDS_out = [[(out, Exo), (out, Emo), (tot, Eoo)],
                           [(out, Eoo)]]
                self.IDS_cont = [[(out, Exx)],
                            [(out, Eox), (into, Emo)]]
                self.IDS_nsupp = [[(into, Exo), (into, Exx)],
                             [(tot, Exo), (tot, Exx), (tot, Exm), (into, Emx),
                                  (into, Eox), (into, Eoo), (into, Emo), (into, Eom), (into, Emm)]]

                #### TO COMPUTE ACCURACY after building
                self.IDS_dL = [Exo]
                self.IDS_dR = [Eox]
                self.IDS_diff = [Exx]
                self.IDS_inter = [Exo, Eox]
                self.IDS_uncovered = [Eoo]

                

            ##############################################################

            ### TO COMPUTE SUPPORTS, no index
            self.IDS_supp = (Exx, Exo, Exm)
            self.IDS_miss = (Emx, Emo, Emm)
            ### indexes swaping when negating one side (0: negating A, 1: negating B)
            self.IDS_negated = [(Eoo, Exx, Eox, Exo, Eom, Emx, Exm, Emo, Emm), \
                           (Exx, Eoo, Exo, Eox, Exm, Emo, Eom, Emx, Emm)]

         #### END WITH MISSING VALUES
         ############################################################################

    ### return part label
    def getLabels(self):
        return self.labels
    
    def getLabel(self, id):
        return self.labels[id]

    # return the index corresponding to part_id when looking from given side 
    def partId(self, part_id, side=0):
        return self.side_index[side][part_id]

    # return the index corresponding to part_id when negating given side 
    def negatedPartId(self, part_id, side=0):
        return self.IDS_negated[side][part_id]

    
    # return the index corresponding to inout and possible negation
    def inOutId(self, inout_id, neg=0):
        return self.neg_index[neg][inout_id]


    # sums the values in parts that correspond to part_id indexes given in parts_id
    ## parts_id can be
    ##  * a list of pairs (inout, part_id), inout are then ignored 
    ##  * a list of values part_id
    def sumPartsId(self, side, parts_id, parts):
        if type(parts) == list  and len(parts_id) > 0:
            if type(parts_id[0]) == int:
                 ids = parts_id    
            elif len(parts_id[0]) == 2:
                (inout, ids) = zip(*parts_id)
            else:
                ids = []
            return sum([parts[self.partId(part_id, side)] for part_id in set(ids)])
        elif type(parts) == int :
            return 1*(parts in [self.partId(part_id[1], side) for part_id in parts_id])
        return 0

    def suppPartRange(self):
        return range(self.bottom, self.top+1)
    def suppPartRangeNoMiss(self):
        return range(self.bottom, self.last_nonmiss+1)

    
    # sums the values in parts that correspond to inout and part_id indexes given in parts_id
    ## parts_id must be
    ##  * a list of pairs (inout, part_id)
    def sumPartsIdInOut(self, side, neg, parts_id, parts):
        return sum([parts[self.inOutId(part_id[0], neg)][self.partId(part_id[1], side)] for part_id in parts_id])


    # return parts reordered to match the new indexes of parts corresponding to negation of given side
    def negateParts(self, side, parts):
        return [parts[self.negatedPartId(p, side)] for p in range(len(parts))]

    # compute the ratio of BLUE/RED parts depending on intersection with X
    def advRatioVar(self, side, op, parts):
        den = self.sumPartsId(side, self.IDS_varden[op], parts)
        num = self.sumPartsId(side, self.IDS_varnum[op], parts)
        return tool_ratio(num, den)

    # compute the accuracy resulting of appending X on given side with given operator and negation
    # from intersections of X with parts (clp)
    def advAcc(self, side, op, neg, lparts, lmiss, lin):
        lout = [lparts[i] - lmiss[i] - lin[i] for i in range(len(lparts))]
        clp = (lin, lout, lparts, lmiss)
        return tool_ratio(self.sumPartsIdInOut(side, neg, self.IDS_varnum[op] + self.IDS_fixnum[op], clp),
                          self.sumPartsIdInOut(side, neg, self.IDS_varden[op] + self.IDS_fixden[op], clp))

    # sets the method to compute p-values
    def setMethodPVal(self, methodpVal='Marg'):
        try:
            self.methodpVal = methodpVal.capitalize()
            eval('self.pVal%sQueryCand' % (self.methodpVal))
            eval('self.pVal%sRedCand' % (self.methodpVal))
            # self.pValQueryCand = eval('self.pVal%sQueryCand' % (self.methodpVal))
            # self.pValRedCand = eval('self.pVal%sRedCand' % (self.methodpVal))
        except AttributeError:
            raise Exception('Oups method to compute the p-value does not exist !')


    def pValRedCand(self, side, op, neg, lParts, N, prs = None, method=""):
        meth = eval('self.pVal%sRedCand' % (self.methodpVal))
        return meth(side, op, neg, lParts, N, prs)

    def pValQueryCand(self, side, op, neg, lParts, N, prs = None):
        meth = eval('self.pVal%sQueryCand' % (self.methodpVal))
        return meth(side, op, neg, lParts, N, prs)
        # return 0 # self.pValSuppQueryCand(side, op, neg, lParts, N, prs)
    
    # query p-value using support probabilities (binomial), for candidates
    def pValSuppQueryCand(self, side, op, neg, lParts, N, prs = None):
        if prs is None:
            return 0
        else:
            lInter = self.sumPartsId(side, self.IDS_supp, lParts[self.inOutId(self.into, neg)])
            lX = float(sum(lParts[self.inOutId(self.into, neg)]))     
            if op:
                return 1-tool_pValSupp(N, lInter, prs[side] + lX/N - prs[side]*lX/N)
            else: 
                return tool_pValSupp(N, lInter, prs[side]*lX/N)

    # query p-value using marginals (binomial), for candidates
    def pValMargQueryCand(self, side, op, neg, lParts, N, prs = None):
        if prs is None:
            return 0
        else:
            lInter = self.sumPartsId(side, self.IDS_supp, lParts[self.inOutId(self.into, neg)])
            lsupp = self.sumPartsId(side, self.IDS_supp, lParts[self.inOutId(self.tot, neg)])
            lX = float(sum(lParts[self.inOutId(self.into, neg)]))     
            if op:
                return 1-tool_pValSupp(N, lInter, lsupp*lX/(N*N))
            else: 
                return tool_pValSupp(N, lInter, lsupp*lX/(N*N))

    # query p-value using support sizes (hypergeom), for candidates
    def pValOverQueryCand(self, side, op, neg, lParts, N, prs = None):
        if prs is None:
            return 0
        else:
            lInter = self.sumPartsId(side, self.IDS_supp, lParts[self.inOutId(self.into, neg)])
            lsupp = self.sumPartsId(side, self.IDS_supp, lParts[self.inOutId(self.tot, neg)])
            lX = sum(lParts[self.inOutId(self.into, neg)])
            if op:
                return 1-tool_pValOver(lInter, N, lsupp, lX)
            else: 
                return tool_pValOver(lInter, N, lsupp, lX)

            
    # redescription p-value using support probabilities (binomial), for candidates
    def pValSuppRedCand(self, side, op, neg, lParts, N, prs = None):
        lInter = self.sumPartsIdInOut(side, neg, self.IDS_fixnum[op] + self.IDS_varnum[op], lParts)
        lUnion = None
        if self.type_parts == "exclu":
            lUnion = self.sumPartsIdInOut(side, neg, self.IDS_fixden[op] + self.IDS_varden[op], lParts)
                
        lX = float(sum(lParts[self.inOutId(self.into, neg)]))
        # if self.pValOut: pdb.set_trace()
        if prs is None :
            lO = self.sumPartsId(1-side, self.IDS_supp, lParts[self.inOutId(self.tot, neg)])
            return tool_pValSupp(N, lInter, float(lO*lX)/(N*N), lU=lUnion)
        elif op:
            return tool_pValSupp(N, lInter, prs[1-side]*(prs[side] + lX/N - prs[side]*lX/N), lU=lUnion)
        else: 
            return tool_pValSupp(N, lInter, prs[1-side]*(prs[side] * lX/N), lU=lUnion)


    # redescription p-value using marginals (binomial), for candidates
    def pValMargRedCand(self, side, op, neg, lParts, N, prs = None):
        lInter = self.sumPartsIdInOut(side, neg, self.IDS_fixnum[op] + self.IDS_varnum[op], lParts)
        lUnion = None
        if self.type_parts == "exclu":
            lUnion = self.sumPartsIdInOut(side, neg, self.IDS_fixden[op] + self.IDS_varden[op], lParts)

        lO = self.sumPartsId(1-side, self.IDS_supp, lParts[self.inOutId(self.tot, neg)])
        lS = self.sumPartsIdInOut(side, neg, self.IDS_nsupp[op], lParts)
        return tool_pValSupp(N, lInter, float(lO*lS)/(N*N), lU=lUnion)
    
    # redescription p-value using support sizes (hypergeom), for candidates
    def pValOverRedCand(self, side, op, neg, lParts, N, prs = None):
        lInter = self.sumPartsIdInOut(side, neg, self.IDS_fixnum[op] + self.IDS_varnum[op], lParts)
        lUnion = None
        if self.type_parts == "exclu":
            lUnion = self.sumPartsIdInOut(side, neg, self.IDS_fixden[op] + self.IDS_varden[op], lParts)

        lO = self.sumPartsId(1-side, self.IDS_supp, lParts[self.inOutId(self.tot, neg)])
        lS = self.sumPartsIdInOut(side, neg, self.IDS_nsupp[op], lParts)
        return tool_pValOver(lInter, N, lO, lS, lU=lUnion)


    # initialize parts counts
    # default count for every part is zero
    # pairs contains a list of (part_id, value)
    # if value is non negative, the count of part_id is set to that value
    # if value is negative, the count of part_id is set to - value - sum of the other parts set so far
    def makeLParts(self, pairs=[], side=0):
        lp = [0 for i in range(self.top+1)]
        for (part_id, val) in pairs:
            if self.partId(part_id, side) < len(lp):
                if val < 0:
                    tmp = sum(lp)
                    lp[self.partId(part_id, side)] = -val- tmp
                else:
                    lp[self.partId(part_id, side)] = val
            else:
                if val > 0:
                    raise Exception("Some missing data where there should not be any!")
        return lp
    

    # adds to parts counts
    # lpartsY can be a part_id in wich case the result of the addition
    # is lpartsX where that part in incremented by one
    def addition(self, lpartsX, lpartsY):
        if type(lpartsY) == list:
            lp = [lpartsX[i]+lpartsY[i] for i in range(len(lpartsX))]    
        else:
            lp = list(lpartsX)
            if type(lpartsY) == int :
                lp[lpartsY] += 1
        return lp

    def additionOtherSide(self, lpartsX, lpartsY, neg=False):
        lp = [lpartsX[i] for i in range(len(lpartsX))]
        if neg:
            XX = self.Exo
        else:
            XX = self.Eoo
        for x in range(self.last_nonmiss+1):
            if x != XX:
                lp[x] += lpartsY[x]
                lp[XX] -= lpartsY[x]
        return lp

    
class SParts(object):

    ### PROPS WHAT
    info_what = {"acc": "self.acc()", "pval": "self.pVal()"}
    props_what = ["len", "card", "supp", "set", "perc", "ratio"]    
    Pwhat_match = "("+ "|".join(info_what.keys()+ props_what) +")"
    @classmethod
    def hasPropWhat(tcl, what):
        return re.match(tcl.Pwhat_match, what) is not None

    ### PROPS WHICH
    sets_letters = "PIULROABN"        
    Pwhich_match = "("+ "|".join(["["+sets_letters+"]"] + SSetts.map_label_part.keys()) +")"    
    @classmethod
    def hasPropWhich(tcl, which):
        return re.match(tcl.Pwhich_match, which) is not None

    props_stats = [("acc", None), ("len", "I"), ("pval", None)]
    
    def __init__(self, ssetts, N, supports, prs = [1,1]):
        #### init from dict_info
        self.ssetts = ssetts
        if type(N) == dict:
            sdict = N
            self.missing = False
            self.sParts = [set() for i in range(len(self.ssetts.getLabels()))]
            self.prs = [-1, -1]
            self.N = 0
            supp_keys = sdict.keys()
            for i, supp_key in enumerate(self.ssetts.getLabels()):
                if supp_key in sdict:
                    if i > 3 and len(sdict[supp_key]) > 0:
                        self.missing = True
                    self.sParts[i] = set(sdict.pop(supp_key))

            if 'pr_0' in sdict:
                self.prs[0] = sdict.pop('pr_0')
            if 'pr_1' in sdict:
                self.prs[1] = sdict.pop('pr_1')
            if 'N' in sdict:
                self.N = sdict.pop('N')
            if not self.missing:
                del self.sParts[4:]
        else:
            if type(N) is set:
                self.N = len(N)
                bk = N
            else:
                self.N = N
                bk = None
            self.prs = prs
            self.vect = None
            ### if include all empty missing parts, remove 
            if type(supports) == list and len(supports) == 4 and len(supports[2]) + len(supports[3]) == 0 :
                supports = supports[0:2]
            elif type(supports) == list and len(supports) == 9 and len(supports[8]) + len(supports[7]) + len(supports[6]) + len(supports[5]) + len(supports[4]) == 0 :
                supports = supports[0:3]

            ### sParts is a partition of the rows (Eoo is not explicitely stored when there are no missing values)
            ## two supports: interpreted as (suppL, suppR)
            if type(supports) == list and len(supports) == 2 :
                (suppL, suppR) = supports
                self.missing = False
                self.sParts = [ set(suppL - suppR), \
                           set(suppR - suppL), \
                           set(suppL & suppR)]
            ## three supports: interpreted as (Exo, Eox, Exx)
            elif type(supports) == list and len(supports) == 3:
                self.missing = False
                self.sParts = [ set(supports[0]), set(supports[1]), set(supports[2])]
            ## four supports: interpreted as (suppL, suppR, missL, missR)
            elif type(supports) == list and len(supports) == 4:
                self.missing = True
                (suppL, suppR, missL, missR) = supports
                self.sParts = [ set(suppL - suppR - missR), \
                           set(suppR - suppL - missL), \
                           set(suppL & suppR), \
                           set(range(self.N)) - suppL -suppR - missL - missR, \
                           set(suppL & missR), \
                           set(suppR & missL), \
                           set(missR - suppL - missL), \
                           set(missL - suppR - missR), \
                           set(missL & missR) ]
            ## nine supports: interpreted as (Exo, Eox, Exx, Eoo, Exm, Emx, Eom, Emo, Emm)
            elif type(supports) == list and len(supports) == 9:
                self.missing = True
                self.sParts = [set(support) for support in supports]
            ## else: set all empty
            else:
                self.missing = False
                self.sParts = [set(), set(), set(), set(), set(), set(), set(), set(), set()]
                bk = None
            if bk is not None:
                if len(self.sParts) == 3:
                    self.sParts.append(set(bk))
                else:
                    self.sParts[self.ssetts.Eoo] = set(bk)
                for si, sp in enumerate(self.sParts):
                    if si != self.ssetts.Eoo:
                        self.sParts[self.ssetts.Eoo] -= sp
                        
    def copy(self):        
        return SParts(self.ssetts, self.N, self.sParts, prs = list(self.prs))
        
    # def __eq__(self, other):
    #     print "Calling EQ"
    #     if isinstance(other, SParts) and len(other.sParts) == len(self.sParts):
    #         for i, p in enumerate(self.sParts):
    #             if other.sParts[i] != p:                    
    #                 return False
    #         return True
    #     return False

    def __cmp__(self, other):
        if isinstance(other, SParts) and len(other.sParts) == len(self.sParts):
            lps = [len(p) for p in self.sParts]
            lpo = [len(p) for p in other.sParts]
            if lps == lpo:
                for i, p in enumerate(self.sParts):
                    if other.sParts[i] != p:                    
                        return cmp(p, other.sParts[i])
                return 0
            return cmp(lps, lpo)
        return -1

    def getTypeParts(self):
        return self.ssetts.getTypeParts()
    def getMethodPVal(self):
        return self.ssetts.getMethodPVal()

    def proba(self, side):
        return self.prs[side]

    def pVal(self):
        try:
            return eval('self.pVal%s()' % self.getMethodPVal())
        except AttributeError:
            raise Exception('Oups method to compute the p-value does not exist !')

    def getSSetts(self):
        return self.ssetts

    def nbRows(self):
        return self.N

    def toDict(self, with_Eoo=False):
        sdict = {}
        for i in range(len(self.sParts)):
                 sdict[self.ssetts.getLabel(i)] = self.part(i)
                 sdict["card_" + self.ssetts.getLabel(i)] = self.lpart(i)
                 sdict["perc_" + self.ssetts.getLabel(i)] = self.lpart(i) * 100. / self.N
        if with_Eoo:
                 sdict[self.ssetts.getLabel(SSetts.Eoo)] = self.part(SSetts.Eoo)
                 sdict["card_" + self.ssetts.getLabel(SSetts.Eoo)] = self.lpart(SSetts.Eoo)
                 sdict["perc_" + self.ssetts.getLabel(SSetts.Eoo)] = self.lpart(SSetts.Eoo) * 100. / self.N
        for side in [0, 1]:
                 if self.prs[side] != -1:
                     sdict["pr_" + str(side)] = self.prs[side]
        sdict["N"] = self.N
        for info_key, info_meth in SParts.info_what.items():
            sdict[info_key] = eval(info_meth)
        return sdict
            
    # contains missing values
    def hasMissing(self):
        return self.missing

    # return copy of the probas
    def probas(self):
        return list(self.prs)

    # return support (used to create new instance of SParts)
    def supparts(self):
        return self.sParts

    # return new instance of SParts corresponding to negating given side
    def negate(self, side=0):
        if self.missing:
            return SParts(self.ssetts, self.N, self.ssetts.negateParts(side, self.sParts))
        else:
            self.sParts.append(self.part(self.ssetts.Eoo))
            n = self.ssetts.negateParts(side, self.sParts)
            return SParts(self.ssetts, self.N, n[0:-1])

    def part(self, part_id, side=0):
        pid = self.ssetts.partId(part_id, side)
        if pid < len(self.sParts):
            return self.sParts[pid]
        elif part_id == self.ssetts.Eoo:
            return set(range(self.N)) - self.sParts[0] - self.sParts[1] - self.sParts[2]
        else:
            return set()
        
    def lpart(self, part_id, side=0):
        pid = self.ssetts.partId(part_id, side)
        if pid < len(self.sParts):
            return len(self.sParts[pid])
        elif part_id == self.ssetts.Eoo:
            return self.N - len(self.sParts[0]) - len(self.sParts[1]) - len(self.sParts[2])
        else:
            return 0

    def parts(self, side=0):
        return [self.part(i, side) for i in range(self.ssetts.top+1)]

    def parts4M(self, side=0):
        if self.missing:
            return [self.part(i, side) for i in range(self.ssetts.Eoo+1)]+[set().union(*[self.part(i, side) for i in range(self.ssetts.Eoo+1, self.ssetts.top+1)])]
        else:
            return self.parts(side)
            
    def lparts(self, side=0):
        return [self.lpart(i, side) for i in range(self.ssetts.top+1)]
    
    def partInterX(self, suppX, part_id, side=0):
        pid = self.ssetts.partId(part_id, side)
        if pid < len(self.sParts):
            return set(suppX & self.sParts[pid])
        elif part_id == self.ssetts.Eoo:
            return set(suppX - self.sParts[0] - self.sParts[1] - self.sParts[2])
        else:
            return set()
        
    def lpartInterX(self, suppX, part_id, side=0):
        pid = self.ssetts.partId(part_id, side)
        if pid < len(self.sParts):
            return len(suppX & self.sParts[pid])
        elif part_id == self.ssetts.Eoo:
            return len(suppX - self.sParts[0] - self.sParts[1] - self.sParts[2])
        else:
            return 0

    def partsInterX(self, suppX, side=0):
        return [self.partInterX(suppX, i, side) for i in range(self.ssetts.top+1)]
    
    def lpartsInterX(self, suppX, side=0):
        if self.missing:
            return [self.lpartInterX(suppX, i, side) for i in range(self.ssetts.top+1)]
        else:
            la = self.lpartInterX(suppX, self.ssetts.Exo, side)
            lb = self.lpartInterX(suppX, self.ssetts.Eox, side)
            lc = self.lpartInterX(suppX, self.ssetts.Exx, side)
            tmp = [la, lb, lc, len(suppX) - la - lb - lc]
            for i in range(len(tmp), self.ssetts.top+1):
                tmp.append(0)
            return tmp

    def nbParts(self):
        return self.ssetts.top+1
        
    def lparts_union(self, ids, side=0):
        return sum([self.lpart(i, side) for i in ids])

    def part_union(self, ids, side=0):
        union = set()
        for i in ids:
            union |= self.part(i, side)
        return union

    def supp(self, side=0):
        return self.part_union(self.ssetts.IDS_supp, side)
    def nonSupp(self, side=0):
        if not self.missing:
            return set(range(self.N)) - self.supp(side)
        else:
            return self.part_union(set(range(self.ssetts.top+1)) - set(self.ssetts.IDS_supp + self.ssetts.IDS_miss), side)
    def miss(self, side=0):
        if not self.missing:
            return set()
        else:
            return self.part_union(self.ssetts.IDS_miss, side)

    def lenSupp(self, side=0):
        return self.lparts_union(self.ssetts.IDS_supp, side)
    def lenNonSupp(self, side=0):
        return self.N - self.lenSupp(side) - self.lenMiss(side)
    def lenMiss(self, side=0):
        if not self.missing:
            return 0
        else:
            return self.lparts_union(self.ssetts.miss_ids, side)

    ### SUPPORTS
    def suppSide(self, side):
        if side == 0:
            return self.part_union(self.ssetts.IDS_dL+self.ssetts.IDS_inter, 0)
        else:
            return self.part_union(self.ssetts.IDS_dR+self.ssetts.IDS_inter, 0)
    def suppP(self, i, side=0):
        return self.part(i, side)
    def suppD(self, side=0):
        return self.part_union(self.ssetts.IDS_diff, side)    
    def suppI(self, side=0):
        return self.part_union(self.ssetts.IDS_inter, side)
    def suppU(self, side=0):
        return self.part_union(self.ssetts.IDS_inter+self.ssetts.IDS_diff, side)
    def suppL(self, side=0):
        return self.suppSide(0)
    def suppR(self, side=0):
        return self.suppSide(1)
    def suppO(self, side=0):
        return self.part_union(self.ssetts.IDS_uncovered, side)
    def suppA(self, side=0):
        return self.part_union(self.ssetts.IDS_dL+self.ssetts.IDS_inter, side)
    def suppB(self, side=0):
        return self.part_union(self.ssetts.IDS_dR+self.ssetts.IDS_inter, side)
    def suppN(self, side=0):
        if len(self.sParts) == 4:
            return self.part_union(range(4), side)
        else:
            return set(range(self.N))

    ### LENGHTS
    ## corresponding lengths
    def lenSide(self, side):
        if side == 0:
            return self.lparts_union(self.ssetts.IDS_dL+self.ssetts.IDS_inter, 0)
        else:
            return self.lparts_union(self.ssetts.IDS_dR+self.ssetts.IDS_inter, 0)
    def lenP(self, i, side=0):
        return self.lpart(i, side)
    def lenD(self, side=0):
        return self.lparts_union(self.ssetts.IDS_diff, side)
    def lenI(self, side=0):
        return self.lparts_union(self.ssetts.IDS_inter, side)
    def lenU(self, side=0):
        return self.lenD(side)+self.lenI(side)
    def lenL(self, side=0):
        return self.lparts_union(self.ssetts.IDS_dL+self.ssetts.IDS_inter, side)
    def lenR(self, side=0):
        return self.lparts_union(self.ssetts.IDS_dR+self.ssetts.IDS_inter, side)
    def lenO(self, side=0):
        return self.lparts_union(self.ssetts.IDS_uncovered, side)    
    def lenA(self, side=0):
        return self.lparts_union(self.ssetts.IDS_dL+self.ssetts.IDS_inter, side)
    def lenB(self, side=0):
        return self.lparts_union(self.ssetts.IDS_dR+self.ssetts.IDS_inter, side)
    def lenN(self, side=0):
        if len(self.sParts) == 4:
            return self.lparts_union(range(4), side)
        else:
            return self.N
    
    # def lenD(self, side=0):
    #     return self.lparts_union(self.ssetts.IDS_diff, side)    
    # def lenI(self, side=0):
    #     return self.lparts_union(self.ssetts.IDS_inter, side)
    # def lenU(self, side=0):
    #     return self.lparts_union(self.ssetts.IDS_inter+self.ssetts.IDS_diff, side)
    #     return self.suppI(side) | self.suppD(side)
    # def lenL(self, side=0):
    #     return self.lenSide(0)
    # def lenR(self, side=0):
    #     return self.lenSide(1)
    # def lenO(self, side=0):
    #     return self.lparts_union(self.ssetts.IDS_uncovered, side)

    def getProp(self, what, which=None):
        if what in SParts.info_what:
            return eval(SParts.info_what[what])
        wt = what
        if what == "card":
            wt = "len"
        elif what == "supp":
            wt = "set"
        methode = eval("self.%s" % wt)
        if callable(methode):
            return methode(which)


    def len(self, which="I"):
        if which in SSetts.map_label_part:
            return self.lenP(SSetts.map_label_part[which])
        elif which in SParts.sets_letters:
            return eval("self.len%s()" % which)
    def set(self, which="I"):
        if which in SSetts.map_label_part:
            return self.suppP(SSetts.map_label_part[which])
        elif which in SParts.sets_letters:
            return eval("self.supp%s()" % which)   
    def ratio(self, which="I"):
        return tool_ratio(self.len(which), self.nbRows())
    def perc(self, which="I"):
        return tool_ratio(self.len(which), self.nbRows()/100.)
    
    # accuracy
    def acc(self, side=0):
        lenI = self.lenI(side)
        return tool_ratio(lenI, lenI+self.lenD(side))

        # redescription p-value using support probabilities (binomial), for redescriptions
    def pValSupp(self):
        if self.prs == [-1,-1] or self.N == -1:
            return -1
        elif self.lenSupp(0)*self.lenSupp(1) == 0:
            return 0
        else:
            lUnion = self.lenU() if self.getTypeParts() == "exclu" else None                
            return tool_pValSupp(self.N, self.lenI(), self.prs[0]*self.prs[1], lU=lUnion)

    # redescription p-value using marginals (binomial), for redescriptions
    def pValMarg(self):
        if self.N == -1:
            return -1
        elif self.lenSupp(0)*self.lenSupp(1) == 0:
            return 0
        else:
            lUnion = self.lenU() if self.getTypeParts() == "exclu" else None                
            return tool_pValSupp(self.N, self.lenI(), float(self.lenSupp(0)*self.lenSupp(1))/(self.N*self.N), lU=lUnion)

    # redescription p-value using support sizes (hypergeom), for redescriptions
    def pValOver(self):
        if self.N == -1:
            return -1
        elif self.lenSupp(0)*self.lenSupp(1) == 0:
            return 0
        else:
            lUnion = self.lenI() if self.getTypeParts() == "exclu" else None
            return tool_pValOver(self.lenI(), self.N, self.lenSupp(0) ,self.lenSupp(1), lU=lUnion)
    
    # moves the instersection of supp with part with index id_from to part with index id_to
    def moveInter(self, side, id_from, id_to, supp):
        self.sParts[self.ssetts.partId(id_to, side)] |= (self.sParts[self.ssetts.partId(id_from,side)] & supp)
        self.sParts[self.ssetts.partId(id_from,side)] -= supp

    # update support probabilities
    def updateProba(prA, prB, OR):
        if type(prA) == int and prA == -1:
            return prB
        elif OR :
            return prA + prB - prA*prB
        else :
            return prA*prB
    updateProba = staticmethod(updateProba)

    # update support probabilities
    def updateProbaMass(prs, OR):
        if len(prs) == 1:
            return prs[0]
        elif OR :
            return reduce(lambda x, y: x+y-x*y, prs)
        else :
            return numpy.prod(prs)
    updateProbaMass = staticmethod(updateProbaMass)

    # update supports and probabilities resulting from appending X to given side with given operator
    def update(self, side, OR, suppX, missX=None):
        self.vect = None
        union = None
        self.prs[side] = SParts.updateProba(self.prs[side], len(suppX)/float(self.N), OR)
            
        if not self.missing and (type(missX) == set and len(missX) > 0):
            self.missing = True
            if len(self.sParts) == 3:
                self.sParts.append(set(range(self.N)) - self.sParts[0] - self.sParts[1] -self.sParts[2])
            else:
                union = set(self.sParts[0] | self.sParts[1] | self.sParts[2] | self.sParts[3])
            self.sParts.extend( [set(), set(), set(), set(), set() ])
            
        if self.missing and self.ssetts.top > self.ssetts.Eoo:
            if OR : ## OR
                ids_from_to_supp = [(self.ssetts.Eox, self.ssetts.Exx ), (self.ssetts.Eoo, self.ssetts.Exo ),
                                    (self.ssetts.Emx, self.ssetts.Exx ), (self.ssetts.Emo, self.ssetts.Exo ),
                                    (self.ssetts.Eom, self.ssetts.Exm ), (self.ssetts.Emm, self.ssetts.Exm )]
                for (id_from, id_to) in ids_from_to_supp:
                    self.moveInter(side, id_from, id_to, suppX)

                if (type(missX) == set and len(missX) > 0):
                    ids_from_to_miss = [(self.ssetts.Eox, self.ssetts.Emx ), (self.ssetts.Eoo, self.ssetts.Emo ),
                                        (self.ssetts.Eom, self.ssetts.Emm )]
                    for (id_from, id_to) in ids_from_to_miss:
                        self.moveInter(side, id_from, id_to, missX)
            
            else: ## AND
                if (type(missX) == set and len(missX) > 0):
                    suppXB  = set(range(self.N)) - suppX - missX
                else:
                    suppXB  = set(range(self.N)) - suppX
                ids_from_to_suppB = [(self.ssetts.Exo, self.ssetts.Eoo ), (self.ssetts.Exx, self.ssetts.Eox ),
                                    (self.ssetts.Exm, self.ssetts.Eom ), (self.ssetts.Emx, self.ssetts.Eox ),
                                    (self.ssetts.Emo, self.ssetts.Eoo ), (self.ssetts.Emm, self.ssetts.Eom )]
                for (id_from, id_to) in ids_from_to_suppB:
                    self.moveInter(side, id_from, id_to, suppXB)
                
                if (type(missX) == set and len(missX) > 0):
                    ids_from_to_miss = [(self.ssetts.Exo, self.ssetts.Emo ), (self.ssetts.Exx, self.ssetts.Emx ),
                                        (self.ssetts.Exm, self.ssetts.Emm )]
                    for (id_from, id_to) in ids_from_to_miss:
                        self.moveInter(side, id_from, id_to, missX)
                
        else :
            if OR : ## OR
                self.sParts[self.ssetts.partId(self.ssetts.Exo,side)] |= (suppX
                                                                       - self.sParts[self.ssetts.partId(self.ssetts.Eox, side)]
                                                                       - self.sParts[self.ssetts.partId(self.ssetts.Exx, side)])
                self.sParts[self.ssetts.partId(self.ssetts.Exx,side)] |= (suppX
                                                                       & self.sParts[self.ssetts.partId(self.ssetts.Eox, side)])
                self.sParts[self.ssetts.partId(self.ssetts.Eox,side)] -= suppX
            
            else: ## AND
                self.sParts[self.ssetts.partId(self.ssetts.Eox,side)] |= (self.sParts[self.ssetts.partId(self.ssetts.Exx, side)]
                                                                       - suppX )
                self.sParts[self.ssetts.partId(self.ssetts.Exx,side)] &= suppX
                self.sParts[self.ssetts.partId(self.ssetts.Exo,side)] &= suppX
        if union is not None:
            self.sParts[self.ssetts.Eoo] = union - self.sParts[self.ssetts.Exx] - self.sParts[self.ssetts.Eox] - self.sParts[self.ssetts.Exo]
        
    # computes vector ABCD (vector containg for each row the index of the part it belongs to)
    def makeVectorABCD(self, force_list=False, rest_ids=None):
        if self.vect is None or (force_list and type(self.vect) is not list):
            if len(self.sParts) == 4 and not force_list:
                # svect = {}
                self.vect = {}
                for partId in range(len(self.sParts)):
                    for i in self.sParts[partId]:
                        self.vect[i] = partId
            else:
                self.vect = [self.ssetts.Eoo for i in range(self.N)]
                map_rest = {}
                if rest_ids is not None:
                    map_rest = dict([(vvv, vvi) for (vvi,vvv) in enumerate(sorted(rest_ids))])
                for partId in range(len(self.sParts)):
                    for i in self.sParts[partId]:
                        self.vect[map_rest.get(i, i)] = partId
                        
                        
    def getVectorABCD(self, force_list=False, rest_ids=None):
        self.makeVectorABCD(force_list, rest_ids)
        if type(self.vect) is dict:
            return None
        return list(self.vect)

    # returns the index of the part the given row belongs to, vectorABCD need to have been computed 
    def partRow(self, row):
        return self.vect[row]

    # return the index of the part the given row belongs to
    # or the intersection of the mode of X with the different parts if row == -1, vectorABCD need to have been computed 
    def lpartsRow(self, row, X=None):
        lp = None
        if row == -1 and X is not None :
            if self.missing:
                lp = [len(X.interMode(self.sParts[i])) for i in range(self.ssetts.top+1)]
            else:
                lp = [0 for i in range(self.nbParts())]
                lp[0] = len(X.interMode(self.sParts[0]))
                lp[1] = len(X.interMode(self.sParts[1]))
                lp[2] = len(X.interMode(self.sParts[2]))
                lp[3] = X.lenMode() - lp[0] - lp[1] - lp[2]
        elif row is not None:
            lp = self.vect[row]
        return lp

############## PRINTING
##############
    # def __str__(self):
#         s = '|'
#         r = '||\n|'
#         if self.missing: up_to = self.ssetts.Emm
#         else: up_to = self.ssetts.Eoo
#         for i in range(up_to+1):
#             s += '|%s' % (3*self.ssetts.getLabel(i))
#             r += '| % 4i ' % self.lpart(i,0)
#         return s+r+'||'

    def __str__(self):
        return "SUPPORT:" + self.dispSuppL(sep=" ")
    
    def dispSuppL(self, sep="\t"):
        return sep.join(["card_" + self.ssetts.getLabel(i)+":" + str(len(self.sParts[i]))         for i in range(len(self.sParts))])
    
    def dispSupp(self, sep="\t"):
        supportStr = ''
        for i in sorted(self.supp(0)): supportStr += "%i "%i
        supportStr += sep
        for i in sorted(self.supp(1)): supportStr += "%i "%i
        if self.missing:
            supportStr += sep
            for i in sorted(self.miss(0)): supportStr += "%i "%i
            supportStr += sep
            for i in sorted(self.miss(1)): supportStr += "%i "%i
        return supportStr

    def dispStats(self, sep="\t"):
        return sep.join(["%s%s:%s" % (what, which or "", self.getProp(what, which)) for (what, which) in self.props_stats])

    # compute the resulting support and missing when combining X and Y with given operator
    def partsSuppMiss(OR, XSuppMiss, YSuppMiss):
        if XSuppMiss is None:
            return YSuppMiss
        elif YSuppMiss is None:
            return XSuppMiss
        elif OR:
            supp = set(XSuppMiss[0] | YSuppMiss[0])
            miss = set(XSuppMiss[1] | YSuppMiss[1]) - supp
        else:
            miss = set((XSuppMiss[1] & YSuppMiss[1]) | (XSuppMiss[1] & YSuppMiss[0]) | (YSuppMiss[1] & XSuppMiss[0]))
            supp = set(XSuppMiss[0] & YSuppMiss[0])
        return (supp, miss)
    partsSuppMiss = staticmethod(partsSuppMiss)

    def partsSuppMissMass(OR, SuppMisses):
        if len(SuppMisses) == 1:
            return SuppMisses[0]
        elif len(SuppMisses) > 1:
            if OR:
                supp = reduce(set.union, [X[0] for X in SuppMisses])
                miss = reduce(set.union, [X[1] for X in SuppMisses]) - supp
            else:
                supp = reduce(set.intersection, [X[0] for X in SuppMisses])
                miss = reduce(set.intersection, [X[0].union(X[1]) for X in SuppMisses]) - supp
            return (supp, miss)
    partsSuppMissMass = staticmethod(partsSuppMissMass)

        # Make binary out of supp set
    def suppVect(self, N, supp, val=1):
        vect = None
        if 2*len(supp) < N:
            st = supp
            v = val
            if val == 1:
                vect = numpy.zeros(N)
            else:
                vect = numpy.ones(N)
        else:
            st = set(range(N)) - supp
            v = 1-val
            if val == 0:
                vect = numpy.zeros(N)
            else:
                vect = numpy.ones(N)
        for i in st:
            vect[i] = v
        return vect
    suppVect = staticmethod(suppVect)

    def parseSupport(stringSupp, N, ssetts):
        partsSupp = stringSupp.rsplit('\t')
        if len(partsSupp) == 2:
            return SParts(ssetts, N, [SParts.parseSupportPart(partsSupp[0]), SParts.parseSupportPart(partsSupp[1])])
        elif len(partsSupp) == 4:
            return SParts(ssetts, N, [SParts.parseSupportPart(partsSupp[0]), SParts.parseSupportPart(partsSupp[1]), \
                          SParts.parseSupportPart(partsSupp[2]), SParts.parseSupportPart(partsSupp[3])])
        return None
    parseSupport = staticmethod(parseSupport)

    def parseSupportPart(string):
        nsupp = set()
        for i in string.strip().rsplit():
            try:
                nsupp.add(int(i))
            except TypeError, detail:
                raise Exception('Unexpected element in the support: %s\n' %i)
        return nsupp
    parseSupportPart = staticmethod(parseSupportPart)

