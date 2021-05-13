from __future__ import unicode_literals
import wx
### from wx import ALIGN_CENTER, ALL, EXPAND, HORIZONTAL
### from wx import FONTFAMILY_DEFAULT, FONTSTYLE_NORMAL, FONTWEIGHT_NORMAL
### from wx import BoxSizer, Button, Font, NewId
### from wx import EVT_BUTTON, EVT_LEFT_DCLICK

import numpy
# The recommended way to use wx with mpl is with the WXAgg
# backend. 
# import matplotlib
# matplotlib.use('WXAgg')

from ..reremi.classQuery import QTree
from ..reremi.classSParts import SSetts
from ..reremi.classRedescription import Redescription
from classDrawerBasis import DrawerEntitiesTD

import pdb

class DrawerRedTree(DrawerEntitiesTD):
    
    all_width = 1.
    height_inter = [2., 3.] ### starting at zero creates troubles with supp drawing, since it's masking non zero values..
    maj_space = 0.02
    min_space = 0.01
    flat_space = 0.03
    margins_sides = 0.5
    margins_tb = 0.1
    margin_hov = min_space/2.
    missing_yy = -1./6

    ann_xy = (10,0)
    
    def __init__(self, view):
        self.view = view
        self.trees = None
        self.store_supp = None
        self.elements = {"active_info": False, "act_butt": [1]}
        self.initPlot()
        self.plot_void()
        ## self.draw()
        
    def getCanvasConnections(self):
        return [('pick_event', self.OnPick),
                ('button_release_event', self.on_click),
                ('motion_notify_event', self.on_motion)]

        
    def plotTrees(self, trees, rids=None):        
        draw_settings = self.getDrawSettings()
        self.plotTreesT(trees, draw_settings, rids)
        # self.plotTreesBasic(trees, draw_settings)
        for side in [0,1]:
            self.plotTree(side, trees[side], None, draw_settings)
            if self.hasMissingPoint(side):
                b = trees[side].getBottomX()
                self.axe.plot(b, self.height_inter[0]+self.missing_yy, 'o', mec='k', mfc=draw_settings[-1]["color_e"], zorder=10)

    def plotTreesBasic(self, trees, draw_settings):
        ##### DRAW each line
        keys = []
        for side in [0,1]:
            for k in trees[side].getLeaves():
                keys.append((side, k, trees[side].getNodeXY(k)[1]))

        mat = numpy.zeros((self.view.getParentData().nbRows(), len(keys)+1))
        for ki, (side, k, pos) in enumerate(keys):
            for si, supp_part in enumerate(trees[side].getNodeSuppSets(k)):
                mat[list(supp_part), ki] = pos
                mat[list(supp_part), -1] = si

        mask = mat > 0
        parts = list(mat[:,-1])
        oids = numpy.argsort((mat[:,:-1].sum(axis=1)+(mat[:,-1]+1)*mat.max()*mat.shape[1])/(mask.sum(axis=1)))
        mat[oids,-1] = numpy.linspace(self.height_inter[0], self.height_inter[1], mat.shape[0])
        for li, k in zip(*numpy.where(mat[:, :-1]> 0)):
            x,y = trees[keys[k][0]].getNodeXY(keys[k][1])
            b = trees[keys[k][0]].getBottomX()
            self.axe.plot((b, 0), (y, mat[li,-1]), color=draw_settings[parts[li]]["color_e"])
            ## self.axe.plot((b, 0), (y, mat[li,-1]), color=draw_settings[parts[li]]["color_e"])
            
    def plotTreesT(self, trees, draw_settings, rids=None):
        ##### DRAW polygons
        map_ids = None
        if rids is None:
            nbrows = self.view.getParentData().nbRows()
        else:
            nbrows = len(rids)
            map_ids = dict([(v,k) for (k,v) in enumerate(rids)])
        mat = numpy.zeros(nbrows, dtype=int)
        tdata = []
        for side in [0,1]:
            for leaf in trees[side].getLeaves():
                for pid, supp_part in enumerate(trees[side].getNodeSuppSets(leaf)):
                    if len(supp_part) > 0:
                        tdata.append((side, leaf, pid))
                        #### ASSIGN ROWS TO LEAVES GROUPS
                        if map_ids is None: 
                            mat[list(supp_part)] += 2**(len(tdata))        
                        else:
                            mat[[map_ids[ii] for ii in supp_part]] += 2**(len(tdata))
                       
                            
        ### gather data for the blocks
        blocks = {}
        map_v_to_bid = {}
        scores_pos = []
        has_miss_points = [False, False]
        for v in numpy.unique(mat):
            sidesb = [-1,-1]
            leaves = []
            pstatus = [None, None]
            vert_pos = 0.
            idxs = tuple([i-1 for (i,vb) in enumerate(bin(v)[::-1]) if vb == '1'])
            for idx in idxs:                
                (side, leaf, pid) = tdata[idx]
                sidesb[side] = idx
                if pstatus[side] is not None:
                    pstatus[side] |= trees[side].isLeafInNode(leaf)
                else:
                    pstatus[side] = trees[side].isLeafInNode(leaf)
                y_leaf = trees[side].getNodeXY(leaf)[1]
                leaves.append({"y": y_leaf, "id": leaf, "part": pid, "side": side})
                
                vert_pos += y_leaf                

            for side in [0,1]:
                if sidesb[side] == -1:
                    has_miss_points[side] = True
                    leaves.append({"y": self.height_inter[0]+self.missing_yy, "id": -1, "part": -1, "side": side})
                
            eids = numpy.where((mat == v))[0]
            bid = tuple(sidesb)
            scores_pos.append((vert_pos/len(leaves), bid))
            map_v_to_bid[v] = bid

            part = SSetts.mapStatusToSPart(tuple(pstatus))
            blocks[bid] = {"eids": eids, "part": part,
                           "pstatus": tuple(pstatus), "leaves": leaves}
        
        #### compute block sizes and positions
        scores_pos.sort()
        bids = [s[1] for s in scores_pos]
        sizes = [len(blocks[bid]["eids"]) for bid in bids]
        tot_spaces = (len(sizes)-1)*self.maj_space
        scale_h = numpy.sum(sizes)/(1-tot_spaces)
        tot_height = self.height_inter[1]-self.height_inter[0]
        for i in range(len(sizes)):
            blocks[bids[i]]["y_bot"] = (i*self.maj_space+numpy.sum(sizes[:i])/scale_h)*tot_height+self.height_inter[0]
            blocks[bids[i]]["y_top"] = (i*self.maj_space+numpy.sum(sizes[:i+1])/scale_h)*tot_height+self.height_inter[0]
            blocks[bids[i]]["x_mid"] = 0.
            for side in [0,1]:
                blocks[bids[i]]["x_leaf%d" %side] = trees[side].getBottomX()
                blocks[bids[i]]["x_flt%d" %side] = -numpy.sign(0.5-side)*self.flat_space

        ### actual plotting
        for bi, bid in enumerate(bids):
            block = blocks[bid]
            # print "BLOCK", bid,block["part"], block["y_leaf0"], block["y_leaf1"]
            for li, leaf in enumerate(block["leaves"]):
                coords_poly = ((block["x_mid"], block["x_flt%d" % leaf["side"]], block["x_leaf%d" % leaf["side"]],
                                block["x_flt%d" % leaf["side"]], block["x_mid"], block["x_mid"]),
                               (block["y_bot"], block["y_bot"], leaf["y"],
                                block["y_top"], block["y_top"], block["y_bot"]))
                # (r,g,b,a) = draw_settings[block["part"]]["color_f"]
                self.axe.fill(coords_poly[0], coords_poly[1],
                               color=draw_settings[block["part"]]["color_f"], linewidth=0) #, linewidth=nb/nbtot)
                # coords_poly[0][:-1], coords_poly[1][:-1], color=draw_settings[block["part"]]["color_e"], linewidth=1
                # if None in block["pstatus"]: ### if group involves missing values 
                #     self.axe.fill(coords_poly[0], coords_poly[1],
                #                color=draw_settings[block["part"]]["color_e"], linewidth=0, alpha=.5) #, linewidth=nb/nbtot)
        self.store_supp = {"has_miss_points": has_miss_points, "blocks": blocks, "mat": mat,  "map_v_to_bid": map_v_to_bid, "map_ids": map_ids, "rids": rids}
        
    def plotTree(self, side, tree, node, ds=None):
        
        color_dot = {0: ds[0]["color_e"], 1: ds[1]["color_e"]}
        line_style = {QTree.branchY: "-", QTree.branchN: "--"}
        # rsym = {0: ">", 1: "<"}
        x, y = tree.getNodeXY(node)

        if tree.isLeafNode(node):
            b = tree.getBottomX()
            self.axe.plot((x, b), (y, y), 'k:')
            self.axe.plot(x, y, 'k.')
            if node is None:
                pass
            elif tree.isLeafInNode(node):
                self.axe.plot(b, y, 'ko', picker=5, gid="%d:%d:-1.T" % (side, node), zorder=10)
            else:
                self.axe.plot(b, y, 'wo', picker=5, gid="%d:%d:+1.T" % (side, node), zorder=10)
                
        else:
            if tree.isParentNode(node):
                for ynb in [0,1]:                        
                    for ci, child in enumerate(tree.getNodeChildren(node, ynb)):
                        xc, yc = tree.getNodeXY(child)
                        self.axe.plot((x, xc), (y, yc), 'k'+line_style[ynb], linewidth=1.5)
                        self.plotTree(side, tree, child, ds=ds)

            if tree.isSplitNode(node):
                self.axe.plot(x, y, color=color_dot[side], marker='s')
                ant = tree.getNodeSplit(node).disp(names=self.getParentData().getNames(side))
                ant = ant.replace("<", "$\leq$").replace(">", "$\geq$")
                self.axe.annotate(ant, xy =(x, y), xytext =(x, y+0.02),
                                  horizontalalignment='center', color=color_dot[side],
                                  bbox=dict(boxstyle="round", fc="w", ec="none", alpha=0.7),**self.view.getFontProps()
                                  )
                self.axe.annotate(ant, xy =(x, y), xytext =(x, y+0.02),
                                  horizontalalignment='center', color=color_dot[side],
                                  bbox=dict(boxstyle="round", fc=color_dot[side], ec="none", alpha=0.3),**self.view.getFontProps()
                                  )

                
    def okTrees(self):
        return self.trees is not None \
          and self.trees[0] is not None and self.trees[1] is not None \
          and not self.trees[0].isBroken() and not self.trees[1].isBroken()

    def getVecAndDets(self, inter_params=None):
        vec = self.getPltDtH().getSuppABCD()
        vec_dets = self.getPltDtH().getVecDets()
        return vec, vec_dets

    def update(self, update_trees=True):
        if self.view.wasKilled():
            return

        if self.isReadyPlot():
            red = self.getPltDtH().getRed()

            if update_trees:
                red.queries[0].side = 0
                red.queries[1].side = 1
                self.trees = [red.queries[0].toTree(self.getParentData()), red.queries[1].toTree(self.getParentData())]

            for butt in self.getElement("buttons"):
                if not self.okTrees():
                    butt["element"].Disable()
                else:
                    butt["element"].Enable()
                    
            if self.okTrees():
                #### WARNING: CONFUSION RESTRICTED SUPPORT <> MISSING VALUES
                ## rsupp = red.supports().parts4M()
                rids = red.getRSetIds(self.getPltDtH().getDetailsSplit())                

                rsupp = red.getRSetParts(self.getPltDtH().getDetailsSplit()).parts4M()
                for side in [0,1]:
                    self.trees[side].computeSupps(side, self.getParentData(), rsupp)
                    if update_trees:
                        self.trees[side].positionTree(side, all_width=self.all_width, height_inter=self.height_inter)

                self.clearPlot()

                inter_params = self.getParamsInter()
                vec, vec_dets = self.getVecAndDets(inter_params)
                draw_settings = self.getDrawSettings()
                ### for highlights
                self.dots_draws = self.prepareEntitiesDots(vec, vec_dets, draw_settings)
                # print "========= LEFT =======\n%s\n%s\n" % (red.queries[0], self.trees[0])
                # print "========= RIGHT =======\n%s\n%s\n" % ( red.queries[1], self.trees[1])
                self.plotTrees(self.trees, rids)

                self.axe.set_xlim([-self.all_width-self.margins_sides, self.all_width+self.margins_sides])
                if self.hasMissingPoints():
                    self.axe.set_ylim([self.height_inter[0]+self.missing_yy-self.margins_tb,
                                       self.height_inter[1]+self.margins_tb])
                else:
                    self.axe.set_ylim([self.height_inter[0]-self.margins_tb, self.height_inter[1]+self.margins_tb])

                self.axe.set_xticks([])
                self.axe.set_yticks([])

                self.draw()
                self.setFocus()
            else:
                self.plot_void()
                
    def hasMissingPoints(self):
        return self.hasMissingPoint(0) or self.hasMissingPoint(1)

    def hasMissingPoint(self, side):
        if self.store_supp is None:
            return False
        return self.store_supp["has_miss_points"][side]

    def makeAdditionalElements(self, panel=None):
        if panel is None:
            panel = self.getLayH().getPanel()
        flags = wx.ALIGN_CENTER | wx.ALL # | wx.EXPAND

        buttons = []
        buttons.extend([{"element": wx.Button(panel, wx.NewId(), size=(self.getLayH().butt_w,-1), label="Expand"),
                         "function": self.view.OnExpandSimp},
                        {"element": wx.Button(panel, wx.NewId(), size=(self.getLayH().butt_w,-1), label="Simplify LHS"),
                         "function": self.OnSimplifyLHS},
                        {"element": wx.Button(panel, wx.NewId(), size=(self.getLayH().butt_w,-1), label="Simplify RHS"),
                         "function": self.OnSimplifyRHS}])
        for i in range(len(buttons)):
            buttons[i]["element"].SetFont(wx.Font(8, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
            
        ##############################################
        add_boxB = wx.BoxSizer(wx.HORIZONTAL)
        add_boxB.AddSpacer((self.getLayH().getSpacerWn()/2.,-1))
        add_boxB.Add(buttons[1]["element"], 0, border=1, flag=flags)
        add_boxB.AddSpacer((self.getLayH().getSpacerWn(),-1))
        add_boxB.Add(buttons[2]["element"], 0, border=1, flag=flags)
        add_boxB.AddSpacer((self.getLayH().getSpacerWn(),-1))
        add_boxB.Add(buttons[0]["element"], 0, border=1, flag=flags)
        add_boxB.AddSpacer((self.getLayH().getSpacerWn()/2.,-1))
        
        self.setElement("buttons", buttons)
        self.setElement("inter_elems", {})
        return [add_boxB]


    def OnPick(self, event):
        if event.mouseevent.button in self.getElement("act_butt"):
            gid_parts = event.artist.get_gid().split(".")
    #def sendOtherPick(self, gid_parts):
            if gid_parts[-1] == "T":
                pp = gid_parts[0].split(":")
                if len(pp) == 3:
                    side, node, onoff = map(int, pp)
                    if self.trees[side] is not None and self.trees[side].hasNode(node):
                        if onoff == -1:
                            self.removeBranchQ(side, node)
                        else:
                            self.addBranchQ(side, node)
            elif gid_parts[-1] == "S":
                self.simplify(int(gid_parts[0]))

    def removeBranchQ(self, side, node):
        if self.trees[side].isLeafInNode(node):
            bidd = self.trees[side].getNodeLeaf(node)
            qu = self.getPltDtH().getQuery(side).copy()
            qu.buk.pop(bidd)                        
            if qu != self.getPltDtH().getQuery(side):
                for n in self.trees[side].getLeaves():
                    if self.trees[side].getNodeLeaf(n) > bidd:
                        self.trees[side].setNodeLeaf(n, self.trees[side].getNodeLeaf(n)-1)
                self.trees[side].setNodeLeaf(node, -1)
                self.getPltDtH().updateQuery(side, query=qu, update_trees=False)

    def addBranchQ(self, side, node):
        if self.trees[side].isLeafOutNode(node):
            buk = self.trees[side].getBranchQuery(node)
            qu = self.getPltDtH().getQuery(side).copy()
            bid = qu.appendBuk(buk)
            if qu != self.getPltDtH().getQuery(side):
                self.trees[side].setNodeLeaf(node, bid)
                self.getPltDtH().updateQuery(side, query=qu, update_trees=False)


    def OnSimplifyLHS(self, event):
        if self.okTrees():
            self.simplify(0)
    def OnSimplifyRHS(self, event):
        if self.okTrees():
            self.simplify(1)
    def simplify(self, side):
        qu = self.trees[side].getSimpleQuery()
        self.getPltDtH().updateQuery(side, query=qu, force=True)
    
    def getCoordsXY(self, idp):
        return []            
    def getCoordsXYA(self, idp):
        center, bid = self.getCenterForId(idp)
        return (self.all_width+self.margins_sides, center)
    
    def drawEntity(self, idp, fc, ec, sz=1, zo=5, dsetts={}):
        ## print "DRAW", idp
        center, bid = self.getCenterForId(idp)
        block = self.store_supp["blocks"][bid]
        lines = []
        for leaf in block["leaves"]:
            coords_poly = ((block["x_mid"], block["x_flt%d" % leaf["side"]], block["x_leaf%d" % leaf["side"]]),
                           (center, center, leaf["y"]))
            lines.extend(self.axe.plot(coords_poly[0], coords_poly[1], color=fc, linewidth=1, zorder=zo))
        return lines

    def hasDotsReady(self):
        return self.store_supp is not None

    
    def drawAnnotation(self, xy, ec, tag, xytext=None):
        if xytext is None:
            xytext = self.getAnnXY()
        dsetts = self.getDrawSettings()
        lines = []
        lines.extend(self.axe.plot((self.flat_space, xy[0]), (xy[1], xy[1]), color=ec, linewidth=1, alpha=0.5))
        lines.extend(self.axe.plot((self.flat_space, xy[0]), (xy[1], xy[1]), color=dsetts["colhigh"], linewidth=1, alpha=0.3))

        lines.extend(DrawerEntitiesTD.drawAnnotation(self, xy, ec, tag, xytext))
        return lines

    def inCapture(self, event):
        return self.okTrees() and event.inaxes == self.getAxe() and numpy.abs(event.xdata) < self.flat_space and event.ydata > self.height_inter[0] and event.ydata < self.height_inter[1]
               
    def getLidAt(self, x, y):
        bid, rp = self.getRPPoint(y)
        if bid is not None:
            if self.store_supp.get("rids") is None:
                return self.store_supp["blocks"][bid]["eids"][rp]
            else:
                return self.store_supp["rids"][self.store_supp["blocks"][bid]["eids"][rp]]


    def getRPPoint(self, y):
        ### Get the code part and line rank for position y
        if self.store_supp is not None and "blocks" in self.store_supp:
            bounds = sorted([(bb["y_bot"], bb["y_top"], bid) for bid, bb in self.store_supp["blocks"].items()])
            i = 0
            while i < len(bounds) and y > bounds[i][1]:
                i += 1
            if i < len(bounds) and y > bounds[i][0]:
                bid = bounds[i][-1]
                nbr = len(self.store_supp["blocks"][bid]["eids"])
                return bid, min(int(nbr*(y - bounds[i][0])/(bounds[i][1] - bounds[i][0])), nbr-1)
            # else:
            #     if i < len(pmap)-1 and y + self.margin_hov > pmap[i][0] and pmap[i+1][-1] == 1:
            #         return pmap[i+1][1], 0

            #     elif i > 0 and y - self.margin_hov <= pmap[i-1][0] and pmap[i-1][-1] == 1:
            #         pp = pmap[i-1][1]
            #         return pp, self.store_supp["pos"][pp][2]-1
        return None, 0

    def getCenterForId(self, idp):
        ## pdb.set_trace()
        idpm = idp
        if self.store_supp.get("map_ids") is not None:
            idpm =  self.store_supp["map_ids"][idp]
        v = self.store_supp["mat"][idpm]
        bid = self.store_supp["map_v_to_bid"][v]
        block = self.store_supp["blocks"][bid]
        pls = numpy.where(block["eids"] == idpm)[0]
        if len(pls) == 1:
            pl = pls[0]
            center = block["y_bot"] + (block["y_top"]-block["y_bot"])*(pl+.5)/(len(block["eids"])+1.)
            return center, bid
        return block["y_bot"], None

