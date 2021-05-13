import wx
### from wx import ALIGN_BOTTOM, ALIGN_CENTER, ALIGN_LEFT, ALIGN_RIGHT, ALL, HORIZONTAL, VERTICAL, ID_ANY, EXPAND, RAISED_BORDER, SL_HORIZONTAL
### from wx import EVT_BUTTON, EVT_SCROLL_THUMBRELEASE, FONTFAMILY_DEFAULT, FONTSTYLE_NORMAL, FONTWEIGHT_NORMAL
### from wx import BoxSizer, Button, CallLater, CheckBox, Choice, DefaultPosition, Font, NewId, Panel,  Slider, StaticText, TextCtrl

import numpy
# The recommended way to use wx with mpl is with the WXAgg backend. 
# import matplotlib
# matplotlib.use('WXAgg')
from classDrawerBasis import DrawerEntitiesTD, DrawerBasis
from classDrawerClust import DrawerClustTD

import pdb

class DrawerProj(DrawerBasis):

    #info_band_height = 240
    margin_hov = 0.01

    def makeAdditionalElements(self, panel=None):
        if panel is None:
            panel = self.getLayH().getPanel()
        flags = wx.ALIGN_CENTER | wx.ALL # | wx.EXPAND

        buttons = []
        buttons.extend([{"element": wx.Button(panel, size=(self.getLayH().butt_w,-1), label="Expand"),
                         "function": self.view.OnExpandSimp},
                        {"element": wx.Button(panel, size=(self.getLayH().butt_w,-1), label="Reproject"),
                         "function": self.view.OnReproject}])

        for i in range(len(buttons)):
            buttons[i]["element"].SetFont(wx.Font(8, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))

        inter_elems = {}
        inter_elems["slide_opac"] = wx.Slider(panel, -1, 10, 0, 100, wx.DefaultPosition, (self.getLayH().sld_w, -1), wx.SL_HORIZONTAL)

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
        add_boxB.Add(buttons[0]["element"], 0, border=1, flag=flags)
        add_boxB.AddSpacer((self.getLayH().getSpacerWn(),-1))
        add_boxB.Add(buttons[1]["element"], 0, border=1, flag=flags)

        add_boxB.AddSpacer((self.getLayH().getSpacerWn()/2.,-1))

        self.setElement("buttons", buttons)
        self.setElement("inter_elems", inter_elems)
        self.setElement("rep_butt", buttons[-1]["element"])
        return [add_boxB]

    def getProj(self):
        return self.view.getProj()
                        
    def makeFinish(self, xylims, xybs):
        if self.getProj().getCoords() is not None:
            if self.getProj().getAxisLabel(0) is not None:
                self.axe.set_xlabel(self.getProj().getAxisLabel(0),fontsize=12)
            if self.getProj().getAxisLabel(1) is not None:
                self.axe.set_ylabel(self.getProj().getAxisLabel(1),fontsize=12)
            self.axe.axis([xylims[0]-xybs[0], xylims[1]+xybs[0], xylims[2]-xybs[1], xylims[3]+xybs[1]])

    def isReadyPlot(self):
        return self.getProj() is not None    
    def getAxisLims(self):
        return self.getProj().getAxisLims()
    def drawPoly(self):
        return False

    def getCoordsXY(self, id):
        if self.getProj() is None:
            return (0,0)
        else:
            return (self.getProj().getCoords(0, ids=id), self.getProj().getCoords(1, ids=id))
    def getCoords(self, axi=None, ids=None):
        if self.getProj() is None:
            return None
        else:
            return self.getProj().getCoords(axi, ids)
    def getCoordsXYA(self, idp):
        return self.getCoordsXY(idp)

class DrawerEntitiesProj(DrawerProj, DrawerEntitiesTD): pass

class DrawerClustProj(DrawerProj, DrawerClustTD):
    
    def makeAdditionalElements(self, panel=None):
        if panel is None:
            panel = self.getLayH().getPanel()
        flags = wx.ALIGN_CENTER | wx.ALL # | wx.EXPAND

        buttons = []
        buttons.extend([{"element": wx.Button(panel, size=(self.getLayH().butt_w,-1), label="Reproject"),
                         "function": self.view.OnReproject}])

        for i in range(len(buttons)):
            buttons[i]["element"].SetFont(wx.Font(8, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))

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
        add_boxB.Add(buttons[0]["element"], 0, border=1, flag=flags)

        add_boxB.AddSpacer((self.getLayH().getSpacerWn(),-1))
        v_box = wx.BoxSizer(wx.VERTICAL)
        label = wx.StaticText(panel, wx.ID_ANY, "dist. inter c")
        label.SetFont(wx.Font(8, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
        add_boxB.Add(label, 0, border=1, flag=flags)
        add_boxB.Add(inter_elems["choice_nbc"], 0, border=1, flag=flags)
        
        add_boxB.AddSpacer((self.getLayH().getSpacerWn()/2.,-1))

        self.setElement("buttons", buttons)
        self.setElement("inter_elems", inter_elems)
        self.setElement("rep_butt", buttons[-1]["element"])
        return [add_boxB]
