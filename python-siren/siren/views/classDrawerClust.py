import wx, numpy, re
# The recommended way to use wx with mpl is with the WXAgg backend. 
import matplotlib
matplotlib.use('WXAgg')

import matplotlib.pyplot as plt
import matplotlib.colors

from classDrawerBasis import DrawerBasis, DrawerEntitiesTD
import pdb

    
class DrawerClustTD(DrawerEntitiesTD):
        
    @classmethod
    def getCMap(tcl, ltid):
        return plt.get_cmap("rainbow")
    def drawPoly(self):
        return False
    
    def getVecAndDets(self, inter_params=None):
        vec, vec_dets = self.getPltDtH().getVecAndDets(inter_params.get("choice_nbc"))
        self.setElement("vec", vec)
        self.setElement("vec_dets", vec_dets)
        return vec, vec_dets

    def getAxisLims(self):
        return self.getPltDtH().getParentCoordsExtrema()

    def makeFinish(self, xylims=(0,1,0,1), xybs=(.1,.1)):
        self.axe.axis([xylims[0], xylims[1], xylims[2], xylims[3]])

    
    #### SEC: ACTIONS
    ######################################
    def makeAdditionalElements(self, panel=None):
        if panel is None:
            panel = self.getLayH().getPanel()
        flags = wx.ALIGN_CENTER | wx.ALL # | wx.EXPAND

        buttons = []

        inter_elems = {}
        inter_elems["slide_opac"] = wx.Slider(panel, -1, 10, 0, 100, wx.DefaultPosition, (self.getLayH().sld_w, -1), wx.SL_HORIZONTAL)
        inter_elems["choice_nbc"] = wx.Choice(panel, -1)
        inter_elems["choice_nbc"].SetItems(["1"])
        inter_elems["choice_nbc"].SetSelection(0)


        ##############################################
        add_boxB = wx.BoxSizer(wx.HORIZONTAL)
        add_boxB.AddSpacer((self.getLayH().getSpacerWn()/2.,-1))

        v_box = wx.BoxSizer(wx.VERTICAL)
        label = wx.StaticText(panel, wx.ID_ANY,u"- opac. disabled +")
        label.SetFont(wx.Font(8, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
        v_box.Add(label, 0, border=1, flag=flags) #, userData={"where": "*"})
        v_box.Add(inter_elems["slide_opac"], 0, border=1, flag=flags) #, userData={"where":"*"})
        add_boxB.Add(v_box, 0, border=1, flag=flags)

        add_boxB.AddSpacer((self.getLayH().getSpacerWn(),-1))
        v_box = wx.BoxSizer(wx.VERTICAL)
        label = wx.StaticText(panel, wx.ID_ANY, "dist. inter c")
        label.SetFont(wx.Font(8, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
        add_boxB.Add(label, 0, border=1, flag=flags)
        add_boxB.Add(inter_elems["choice_nbc"], 0, border=1, flag=flags)   

        add_boxB.AddSpacer((self.getLayH().getSpacerWn()/2.,-1))

        self.setElement("buttons", buttons)
        self.setElement("inter_elems", inter_elems)        
        return [add_boxB]


    def plotMapperHist(self, axe, vec, vec_dets, mapper, nb_bins, corners, draw_settings):
        norm = matplotlib.colors.Normalize(vmin=0, vmax=1, clip=True)
        mappers = [matplotlib.cm.ScalarMappable(norm=norm, cmap="Purples"),
                   matplotlib.cm.ScalarMappable(norm=norm, cmap="binary")]

        x0, x1, y0, y1, bx, by = corners
        fracts = [.25, .05] ## ratio bars occ/fixed
        nbc = len(vec_dets["binLbls"])        
        bins_ticks = numpy.arange(nbc)
        tmpb = [b-0.5 for b in bins_ticks]
        tmpb.append(tmpb[-1]+1)

        # norm_bins_ticks = [(bi-tmpb[0])/float(tmpb[-1]-tmpb[0]) * 0.95*float(y1-y0) + y0 + 0.025*float(y1-y0) for bi in bins_ticks]
        # norm_bins = [(bi-tmpb[0])/float(tmpb[-1]-tmpb[0]) * 0.95*float(y1-y0) + y0 + 0.025*float(y1-y0) for bi in tmpb]
        norm_bins_ticks = [(bi-tmpb[0])/float(tmpb[-1]-tmpb[0]) *float(y1-y0) + y0 for bi in bins_ticks]
        norm_bins = [(bi-tmpb[0])/float(tmpb[-1]-tmpb[0]) *float(y1-y0) + y0 for bi in tmpb]
        left = [norm_bins[i] for i in range(nbc)]
        width = [norm_bins[i+1]-norm_bins[i] for i in range(nbc)]


        nbr = vec_dets["more"][0]["occ_cnt"].shape[0]
        h_occ = (fracts[0]*(x1-x0))/nbr
        h_hist = fracts[1]*(x1-x0)+2*bx
        bottom_occ = x1
        bottom_hist = bottom_occ+nbr*h_occ
        top_hist = bottom_hist+h_hist
        btms = [bottom_occ+i*h_occ for i in range(nbr)]
        
        bckc = "white"        
        bins_lbl = vec_dets["binLbls"]
        #vvmax = int(numpy.max(vec))
        colors = [mapper.to_rgba(i) for i in vec_dets["binVals"]]        
        # colors[-1] = draw_settings["default"]["color_f"]
        
        axe.barh(y0, nbr*h_occ+h_hist, y1-y0, x1, color=bckc, edgecolor=bckc)
        # axe.plot([bottom_occ, bottom_occ], [y0, y1-y0], color="blue")
        # axe.plot([bottom_hist, bottom_hist], [y0, y1-y0], color="red")
        # axe.plot([bottom+nbr*h, bottom+nbr*h], [y0, y1-y0], color="red")
        axe.barh(left, numpy.ones(nbc)*h_hist, width, numpy.ones(nbc)*bottom_hist, color=colors, edgecolor=bckc, linewidth=2)
        axe.plot([bottom_hist, bottom_hist], [norm_bins[0], norm_bins[-1]], color="black", linewidth=.2)
        axe.plot([bottom_occ, bottom_occ], [norm_bins[0], norm_bins[-1]], color="black", linewidth=.2)
        

        for pi, i in enumerate(vec_dets["more"]["orids"]):
            clrs = [mappers[int(vec_dets["more"][i]["occ_cnt"][j])].to_rgba(vec_dets["more"][i]["occ_avg"][j]) for j,v in enumerate(vec_dets["more"][i]["occ_avg"])]
            axe.barh(numpy.ones(nbr)*left[pi], numpy.ones(nbr)*h_occ, numpy.ones(nbr)*width[pi], btms, color=clrs, edgecolor=bckc, linewidth=0)
        
        x1 += nbr*h_occ+h_hist #(fracts[0]+fracts[1])*(x1-x0)+2*bx

        self.hist_click_info = {"left_edge_map": x0, "right_edge_map": bottom_occ, "right_edge_occ": bottom_hist, "right_edge_hist": x1,
                                 "hedges_hist": norm_bins, "vedges_occ": btms}
        
        axe.set_yticks(norm_bins_ticks)
        axe.set_yticklabels(bins_lbl, **self.view.getFontProps())
        # self.axe.yaxis.tick_right()
        axe.tick_params(direction="inout", left="off", right="on",
                            labelleft="off", labelright="on")
        return (x0, x1, y0, y1, bx, by)
        

    def on_click(self, event):
        # print "Event location:", event.xdata, event.ydata
        if self.clickActive() and self.inCapture(event):
            if event.xdata > self.hist_click_info['right_edge_occ'] and event.xdata < self.hist_click_info['right_edge_hist'] and \
              event.ydata > self.hist_click_info['hedges_hist'][0] and event.ydata < self.hist_click_info['hedges_hist'][-1]:
                self.on_click_hist(event)
            elif event.xdata > self.hist_click_info['right_edge_map'] and event.xdata < self.hist_click_info['right_edge_occ'] and \
              event.ydata > self.hist_click_info['hedges_hist'][0] and event.ydata < self.hist_click_info['hedges_hist'][-1]:
                self.on_click_occ(event)
            elif event.xdata > self.hist_click_info['left_edge_map'] and event.xdata < self.hist_click_info['right_edge_map'] and \
              event.ydata > self.hist_click_info['hedges_hist'][0] and event.ydata < self.hist_click_info['hedges_hist'][-1]:
                lid = self.getLidAt(event.xdata, event.ydata)
                if lid is not None:
                    self.sendEmphasize([lid])

    def on_click_hist(self, event):
        bini = 0
        while event.ydata > self.hist_click_info['hedges_hist'][bini]:
            bini += 1
        bval = self.getElement("vec_dets")["binVals"][bini-1]
        lids = numpy.where(self.getElement("vec") == bval)[0]
        if len(lids) > 0:
            self.sendEmphasize(lids)

    def on_click_occ(self, event):
        bini = 0
        while bini < len(self.hist_click_info['hedges_hist']) and event.ydata > self.hist_click_info['hedges_hist'][bini]:
            bini += 1
        ri = 0
        while ri < len(self.hist_click_info['vedges_occ']) and event.xdata > self.hist_click_info['vedges_occ'][ri]:
            ri += 1
        # status = 1
        # if event.ydata < (self.hist_click_info['hedges_hist'][bini]+self.hist_click_info['hedges_hist'][bini-1])/2.:
        #     status = 0
        bval = self.getElement("vec_dets")["binVals"][bini-1]
        etor = self.getPltDtH().getEtoR()
        lids = numpy.where(etor[:,ri-1] & (self.getElement("vec") == bval))[0]
        if len(lids) > 0:
            self.sendEmphasize(lids)
            
    def makeEmphTag(self, lid):
        tag = "%s" % self.getParentData().getRName(lid)
        if self.getElement("vec") is not None:
            c = self.getElement("vec")[lid]
            if c >= 0:
                tag += ": c%s" % c
        return tag
