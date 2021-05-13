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
import matplotlib
matplotlib.use('WXAgg')

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

from ..reremi.classSParts import SSetts
from ..reremi.classRedescription import Redescription
from classDrawerBasis import DrawerEntitiesTD

import pdb

class DrawerRedCorrel(DrawerEntitiesTD):
    
    # all_width = 1.
    # height_inter = [2., 3.] ### starting at zero creates troubles with supp drawing, since it's masking non zero values..
    # maj_space = 0.02
    # min_space = 0.01
    # flat_space = 0.03
    # margins_sides = 0.5
    # margins_tb = 0.1
    # margin_hov = min_space/2.
    # missing_yy = -1./6

    # ann_xy = (10,0)
    
    def __init__(self, view):
        self.view = view
        self.store_supp = None
        self.elements = {"active_info": False, "act_butt": [1]}
        self.parts_in = {SSetts.Exx: True, SSetts.Exo: False, SSetts.Eox: False, SSetts.Eoo: False}
        self.parts_out = {SSetts.Exx: True, SSetts.Exo: True, SSetts.Eox: True, SSetts.Eoo: True} 
        self.initPlot()
        self.plot_void()
        ## self.draw()
        
    def getCanvasConnections(self):
        return []

    def getVecAndDets(self, inter_params=None):
        vec_dets = self.getPltDtH().getVecDets()
        return vec, vec_dets
        

    def OnRedraw(self, event=None):
        self.update()
    
    def update(self, update_trees=True):
        if self.view.wasKilled():
            return

        if self.isReadyPlot():

            elems = self.getElement("inter_elems")
            ks = self.parts_in.keys()
            if elems is not None:
                for k in ks:
                    lbl_in = "in_%s" % SSetts.labels[k]
                    lbl_out = "out_%s" % SSetts.labels[k]
                    if lbl_in in elems:
                        if elems[lbl_in].IsChecked():
                            self.parts_in[k] = True
                        else:
                            self.parts_in[k] = False
                    if lbl_out in elems:
                        if elems[lbl_out].IsChecked():
                            self.parts_out[k] = True
                        else:
                            self.parts_out[k] = False
                            
            vec = self.getPltDtH().getSuppABCD()
            mat, details, mcols = self.getParentData().getMatrix(nans=numpy.nan, only_able=True)
            self.clearPlot()

            pos_in = numpy.zeros(vec.shape, dtype=bool)
            pos_out = numpy.zeros(vec.shape, dtype=bool)
            for part, v in self.parts_in.items():
                if v:
                    pos_in |= (vec == part)
            for part, v in self.parts_out.items():
                if v:
                    pos_out |= (vec == part)

            Rall = numpy.corrcoef(mat)
            if numpy.sum(pos_in) > 0:
                Rin = numpy.corrcoef(mat[:, pos_in])
            else:
                Rin = Rall
            if numpy.sum(pos_out) > 0:
                Rout = numpy.corrcoef(mat[:, pos_out])
            else:
                Rout = Rout

            cmap = cm.get_cmap('PuOr')
            xs, ys = numpy.meshgrid(numpy.arange(Rout.shape[0]), numpy.arange(Rout.shape[0]))
            flt_xs, flt_ys, flt_Rout, flt_Rin = numpy.ravel(xs), numpy.ravel(ys), numpy.ravel(Rout), numpy.ravel(Rin)
            ids_under = (flt_ys == 0) | (flt_xs == Rout.shape[0]-1)
            ids_under = flt_ys > flt_xs
            flt_xs, flt_ys, flt_Rout, flt_Rin = flt_xs[ids_under], flt_ys[ids_under], flt_Rout[ids_under], flt_Rin[ids_under]
            angle = numpy.pi*1/4
            rot = numpy.array([[numpy.cos(angle), -numpy.sin(angle)],[numpy.sin(angle), numpy.cos(angle)]])
            rot_xys = numpy.dot(numpy.vstack([flt_xs, flt_ys]).T, rot)           

            labels = [d["name"] for d in details]
            width_band = .45
            margband_bot = width_band
            margband_top = 1.
            
            for i, lbl in enumerate(labels):
                cline = "#888888" if i %2 == 0 else "#CCCCCC"
                if i > 0:
                    xys = numpy.dot(numpy.array([[-1, i], [-margband_top, i-width_band], [-margband_top, i+width_band], [(i-1)+margband_bot, i+width_band], [(i-1)+margband_bot, i-width_band]]), rot)
                    self.axe.fill(xys[1:, 0], xys[1:, 1], color=cline, alpha=.7, zorder=1)
                    self.axe.text(xys[0, 0], xys[0, 1], lbl, rotation=-numpy.degrees(angle), ha="right", va="bottom")

                if i < len(labels)-1:
                    xys = numpy.dot(numpy.array([[i, Rout.shape[0]], [i-width_band, (Rout.shape[0]-1)+margband_top], [i+width_band, (Rout.shape[0]-1)+margband_top], [i+width_band, (i+1)-margband_bot], [i-width_band, (i+1)-margband_bot]]), rot)
                    self.axe.fill(xys[1:, 0], xys[1:, 1], color=cline, alpha=.7, zorder=1)
                    self.axe.text(xys[0, 0], xys[0, 1], lbl, rotation=numpy.degrees(angle), ha="left", va="bottom")

            patches = [Circle((rot_xys[i,0], rot_xys[i,1]), radius=.48*numpy.abs(flt_Rin[i])) for i in range(flt_Rin.shape[0])]
            ## .1+.4*numpy.abs(flt_Rin[i]-flt_Rout[i])
            fcolors = [cmap(.5*(flt_Rin[i]+1)) for i in range(flt_Rin.shape[0])]
            ecolors = [cmap(.5*(flt_Rout[i]+1)) for i in range(flt_Rin.shape[0])]
            # ecolors = [ecmap(.5*numpy.abs(flt_Rin[i]-flt_Rout[i])) for i in range(flt_Rin.shape[0])]
            p = PatchCollection(patches, alpha=1., zorder=10, facecolors = fcolors, edgecolors = ecolors, linewidths=(2.,))
            self.axe.add_collection(p)

            diag_size = numpy.sqrt(2)*(Rout.shape[0]-1)
            nb_bins = 100
            width_bin = diag_size/nb_bins
            for i in range(nb_bins):
                self.axe.fill([i*width_bin, (i+1)*width_bin, (i+1)*width_bin, i*width_bin], [-.75, -.75, -1.25, -1.25], color=cmap(i/(nb_bins-1.)))

            self.axe.text(-1, -1., "-1", ha="center", va="center")
            self.axe.text(diag_size+1, -1., "+1", ha="center", va="center")
                    
            self.axe.set_xlim([0-.2*Rout.shape[0], diag_size+.2*Rout.shape[0]])
            self.axe.set_ylim([-2, diag_size/2+.2*Rout.shape[0]])
            self.axe.set_xticks([])
            self.axe.set_yticks([])

            self.draw()
            self.setFocus()
        else:
            self.plot_void()
                

    # def makeAdditionalElements(self, panel=None):
    #     self.setElement("buttons", [])
    #     self.setElement("inter_elems", {})
    #     return []

    def makeAdditionalElements(self, panel=None):
        if panel is None:
            panel = self.getLayH().getPanel()
        flags = wx.ALIGN_CENTER | wx.ALL # | wx.EXPAND

        buttons = []
        buttons.extend([{"element": wx.Button(panel, size=(self.getLayH().butt_w,-1), label="Expand"),
                         "function": self.view.OnExpandSimp},
                        {"element": wx.Button(panel, size=(self.getLayH().butt_w,-1), label="Redraw"),
                         "function": self.OnRedraw}])
        

        for i in range(len(buttons)):
            buttons[i]["element"].SetFont(wx.Font(8, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))

        inter_elems = {}

        ##############################################
        add_boxB = wx.BoxSizer(wx.HORIZONTAL)
        add_boxB.AddSpacer((self.getLayH().getSpacerWn()/2.,-1))

        for part, v in self.parts_in.items():
            sub = SSetts.labels[part] 
            inter_elems["in_"+sub] = wx.CheckBox(panel, wx.NewId(), "", style=wx.ALIGN_RIGHT)
            inter_elems["in_"+sub].SetValue(v)
            inter_elems["out_"+sub] = wx.CheckBox(panel, wx.NewId(), "", style=wx.ALIGN_RIGHT)
            inter_elems["out_"+sub].SetValue(self.parts_out[part])

            v_box = wx.BoxSizer(wx.VERTICAL)
            label = wx.StaticText(panel, wx.ID_ANY, sub)
            label.SetFont(wx.Font(8, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
            v_box.Add(label, 0, border=1, flag=flags) #, userData={"where": "*"})
            v_box.Add(inter_elems["in_"+sub], 0, border=1, flag=flags) #, userData={"where":"*"})
            v_box.Add(inter_elems["out_"+sub], 0, border=1, flag=flags) #, userData={"where":"*"})
            add_boxB.Add(v_box, 0, border=1, flag=flags)

        add_boxB.AddSpacer((self.getLayH().getSpacerWn(),-1))
        add_boxB.Add(buttons[0]["element"], 0, border=1, flag=flags)
        add_boxB.AddSpacer((self.getLayH().getSpacerWn(),-1))
        add_boxB.Add(buttons[1]["element"], 0, border=1, flag=flags)

        add_boxB.AddSpacer((self.getLayH().getSpacerWn()/2.,-1))

        self.setElement("buttons", buttons)
        self.setElement("inter_elems", inter_elems)
        return [add_boxB]


    # def OnPick(self, event):

    # def hasDotsReady(self):
    #     return self.store_supp is not None
    
    # def drawAnnotation(self, xy, ec, tag, xytext=None):

    # def inCapture(self, event):
               
    # def getLidAt(self, x, y):
