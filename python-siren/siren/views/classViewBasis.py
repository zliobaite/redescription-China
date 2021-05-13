import wx
import numpy
 
from ..reremi.classSParts import SSetts

from classLayoutHandler import LayoutHandlerBasis, LayoutHandlerQueries

from classDrawerBasis import DrawerBasis, DrawerEntitiesTD
from classDrawerPara import DrawerRedPara
from classDrawerTree import DrawerRedTree
from classDrawerMap import DrawerEntitiesMap, DrawerClustMap
from classDrawerMappoly import DrawerEntitiesMappoly
from classDrawerProj import DrawerEntitiesProj, DrawerClustProj

from classDrawerCorrel import DrawerRedCorrel

from classPltDtHandler import PltDtHandlerBasis, PltDtHandlerRed, PltDtHandlerRedWithCoords
from classPltDtHList import PltDtHandlerListClust
from classProj import ProjFactory

import pdb

class ViewBasis(object):
    """
    The parent class of all visualizations.
    """

    colors_ord = ["color_l", "color_r", "color_i", "color_o"]
    colors_def = {"color_l": (255,0,0), "color_r": (0,0,255), "color_i": (160,32,240), "color_o": (153, 153, 153),
                  "grey_basic": (127,127,127), "grey_light": (153,153,153), "grey_dark": (85,85,85),
                  "color_h": (255, 255, 0), -1: (127, 127, 127)}
    DOT_ALPHA = 0.6
    ## 153 -> 99, 237 -> ed
    DOT_SHAPE = 's'
    DOT_SIZE = 3

    DELTA_ON = False
    DEF_ZORD = 3
    
    TID = "-"
    SDESC = "-"
    ordN = 0
    title_str = "Basis View"
    geo = False
    typesI = ""

    subcl_layh = LayoutHandlerBasis
    subcl_drawer = DrawerBasis
    subcl_pltdt = PltDtHandlerBasis
    
    def __init__(self, parent, vid, more=None):
        self.parent = parent
        self.vid = vid
        self.data = {}
        self.layH = self.subcl_layh(self)
        self.pltdtH = self.subcl_pltdt(self)
        self.drawer = self.subcl_drawer(self)
        self.layH.setToolbarDrawer(self.drawer)
        boxes = self.drawer.prepareInteractive(self.layH.getPanel())
        self.layH.finalize_init(boxes)
        
    def getLayH(self):
        return self.layH
    def getDrawer(self):
        return self.drawer
    def getPltDtH(self):
        return self.pltdtH
    def getParent(self):
        return self.parent

    
    #### SEC: VIEW IDENTIFICATION
    ###########################################

    @classmethod
    def getViewsDetails(tcl):
        if tcl.TID is not None:
            return {tcl.TID: {"title": tcl.title_str, "class": tcl, "more": None, "ord": tcl.ordN}}
        return {}
    
    @classmethod
    def suitableView(tcl, geo=False, what=None, tabT=None):
        return (tabT is None or tabT in tcl.typesI) and (not tcl.geo or geo)
    
    def getItemId(self):
        if self.hasParent():
            return self.getParentViewsm().getItemId(self.getId())
        return self.vid
    
    def getShortDesc(self):
        return "%s %s" % (self.getItemId(), self.SDESC)

    def getTitleDesc(self):
        return "%s %s" % (self.getItemId(), self.title_str)

    def getId(self):
        return (self.TID, self.vid)

    def getVId(self):
        return self.vid

    
    #### SEC: PARENT ACCESS
    ###########################################

    def hasParent(self):
        return self.parent is not None
    def getParentVizm(self):
        if self.hasParent():
            return self.parent.getVizm()
    def getParentViewsm(self):
        if self.hasParent():
            return self.parent.getViewsm()
    def getParentTab(self, which):
        if self.hasParent():
            return self.parent.getTab(which)
    def getParentData(self):
        if self.hasParent():
            return self.parent.dw.getData()
    def getParentPreferences(self):
        if self.hasParent():
            return self.parent.dw.getPreferences()
        return {}
    def getParentIcon(self, key=None):
        if self.hasParent():
            return self.parent.icons.get(key)
    def getParentTitlePref(self):
        if self.hasParent():
            return self.parent.titlePref
        return "Standalone "

    def isIntab(self):
        return self.getLayH().isInTab
    def toTop(self):
        self.getLayH().toTop()
    def _SetSize(self):
        self.getLayH()._SetSize()
    def updateTitle(self):
        self.getLayH().updateTitle()
    def getGPos(self):
        return self.getLayH().getGPos()
    def popSizer(self):
        return self.getLayH().popSizer()
    def destroy(self):
        self.getLayH().destroy()
    def wasKilled(self):
        return self.getLayH().wasKilled()
    def emphasizeOnOff(self, turn_on=set(), turn_off=set(), hover=False, review=True):
        self.getDrawer().emphasizeOnOff(turn_on, turn_off, hover, review)
    def updateRSets(self, new_rsets=None):
        self.getPltDtH().updateRSets(new_rsets)
        
        
    def lastStepInit(self, blocking=False):
        pass
    def OnQuit(self, event=None, upMenu=True, freeing=True):
        if self.hasParent():
            self.getParentViewsm().deleteView(self.getId(), freeing)
            self.getParentViewsm().unregisterView(vkey=self.getId(), upMenu=upMenu)
        else:
            self.getLayH().Destroy()
            
    #### SEC: MENU
    ######################################
    def makeMenu(self, frame=None):
        """
        Prepare the menu for this view.

        @type  frame: wx.Frame
        @param frame: The frame in which the menu resides
        """
        
        if self.isIntab():
            return
        
        if frame is None:
            frame = self.getLayH().frame

        menuBar = wx.MenuBar()
        if self.hasParent():
            menuBar.Append(self.parent.makeFileMenu(frame), "&File")
        menuBar.Append(self.getDrawer().makeActionsMenu(frame), "&Edit")
        menuBar.Append(self.makeVizMenu(frame), "&View")
        menuBar.Append(self.makeProcessMenu(frame), "&Process")
        
        if self.hasParent():
            menuBar.Append(self.parent.makeViewsMenu(frame), "&Windows")
            menuBar.Append(self.parent.makeHelpMenu(frame), "&Help")
        frame.SetMenuBar(menuBar)
        frame.Layout()

    def enumerateVizItems(self):
        if self.hasParent():
            return self.getParentViewsm().getViewsItems(vkey=self.getId())
        return []
    def makeVizMenu(self, frame, menuViz=None):
        """
        Prepare the visualization sub-menu for this view.

        @type  frame: wx.Frame
        @param frame: The frame in which the menu resides
        @type  menuViz: wx.Menu
        @param menuViz: Existing menu, if any, where entries will be appended
        @rtype:   wx.Menu
        @return:  the sub-menu, menuViz extended
        """

        self.ids_viewT = {}
        if menuViz is None:
            menuViz = wx.Menu()
        for item in self.enumerateVizItems():
            ID_NEWV = wx.NewId()
            m_newv = menuViz.Append(ID_NEWV, "%s" % item["title"],
                                    "Plot %s." % item["title"])
            if not item["suitable"]:
                m_newv.Enable(False)

            frame.Bind(wx.EVT_MENU, self.OnOtherV, m_newv)
            self.ids_viewT[ID_NEWV] = item["viewT"]
        if menuViz.GetMenuItemCount() == 0:
            self.getParent().appendEmptyMenuEntry(menuViz, "No Views", "There are no other views.")
        return menuViz       
    def OnOtherV(self, event):
        if self.hasParent():
            self.getParentViewsm().viewOther(viewT=self.ids_viewT[event.GetId()], vkey=self.getId())

    def makeProcessMenu(self, frame, menuPro=None):
        self.menu_map_pro = {}
        if menuPro is None:
            menuPro = wx.Menu()

        for process, details in self.getLayH().getProcesses():
            ID_PRO = wx.NewId()
            m_pro = menuPro.Append(ID_PRO, details["label"], details["legend"])
            if self.q_expand(details["more"]):
                frame.Bind(wx.EVT_MENU, self.OnExpandAdv, m_pro)
                self.menu_map_pro[ID_PRO] = process
            else:
                menuPro.Enable(ID_PRO, False)
        ct = menuPro.GetMenuItemCount()
        if self.hasParent():
            menuPro = self.parent.makeStoppersMenu(frame, menuPro)
        if ct < menuPro.GetMenuItemCount():
            menuPro.InsertSeparator(ct)
        return menuPro

    def OnExpandAdv(self, event):
        if self.getPltDtH().hasQueries():
            params = {"red": self.getPltDtH().getCopyRed()}
            if event.GetId() in self.menu_map_pro:
                params = self.getLayH().getProcessesParams(self.menu_map_pro[event.GetId()], params)
            self.getParent().expandFromView(params)

    def OnExpandSimp(self, event):
        if self.getPltDtH().hasQueries():
            params = {"red": self.getPltDtH().getCopyRed()}
            self.getParent().expandFromView(params)
    def getWeightCover(self, params):
        params["area"] = self.getDrawer().getHighlightedIds()
        return params
    def q_expand(self, more):
        if not self.getPltDtH().hasQueries():
            return False
        if more is None:
            return True
        res = True
        if "side" in more:
            res &= len(self.getPltDtH().getQuery(1-more["side"])) > 0
        if "in_weight" in more or "out_weight" in more:
            res &= self.getDrawer().q_has_selected()
        return res
            
    #### SEC: HANDLING SETTINGS
    ###########################################
    def getSettV(self, key, default=False):
        t = self.getParentPreferences()
        try:
            v = t[key]["data"]
        except:            
            v = default
        return v

    def getFontProps(self):
        return {"size": self.getSettV("plot_fontsize")}
    def getColorKey1(self, key, dsetts=None):
        if dsetts is None:
            dsetts = self.getParentPreferences()
        if key in dsetts:
            tc = dsetts[key]["data"]
        elif key in self.colors_def:
            tc = self.colors_def[key]
        else:
            tc = self.colors_def[-1]
        return [i/255.0 for i in tc]+[1.]
    def getColorKey255(self, key, dsetts=None):
        if dsetts is None:
            dsetts = self.getParentPreferences()
        if key in dsetts:
            tc = dsetts[key]["data"]
        elif key in self.colors_def:
            tc = self.colors_def[key]
        else:
            tc = self.colors_def[-1]
        return tc

    def getAlpha(self, alpha=None, color=None):
        if self.alpha_off:
            alpha = 1.
        else:
            if alpha is None:
                 alpha = self.DOT_ALPHA
            elif alpha < -1 or alpha > 1:
                alpha = numpy.sign(alpha)*(numpy.abs(alpha)%1)*self.DOT_ALPHA
            if alpha < 0:
                alpha = -color[3]*alpha
        return alpha
    
    def getColorA(self, color, alpha=None):
        alpha = self.getAlpha(alpha, color)
        return tuple([color[0],color[1], color[2], alpha])
    
    def getColorHigh(self):
        return self.getColorA(self.getColorKey1("color_h"))

    def getColors255(self):
        return  [ self.getColorKey255(color_k) for color_k in self.colors_ord ]

    def getColors1(self):
        return  [ self.getColorKey1(color_k) for color_k in self.colors_ord ]
    
    def getDrawSettDef(self):
        t = self.getParentPreferences()
        try:
            dot_shape = t["dot_shape"]["data"]
            dot_size = t["dot_size"]["data"]
        except:
            dot_shape = self.DOT_SHAPE
            dot_size = self.DOT_SIZE

        return {"color_f": self.getColorA(self.getColorKey1("grey_basic")),
                "color_e": self.getColorA(self.getColorKey1("grey_basic"), 1.),
                "shape": dot_shape, "size": dot_size, "zord": self.DEF_ZORD}

    def setAlphaOnOff(self):
        t = self.getParentPreferences()
        if t["alpha_off"]["data"] == 'yes':
            self.alpha_off = True
        else:
            self.alpha_off = False
            
    def getDrawSettings(self):
        self.setAlphaOnOff()
        colors = self.getColors1()
        colhigh = self.getColorHigh()
        fontprops = self.getFontProps()
        defaults = self.getDrawSettDef()
        if self.getSettV('miss_details'):
            zord_miss = self.DEF_ZORD
        else:
            zord_miss = -1       
        draw_pord = dict([(v,p) for (p,v) in enumerate([SSetts.Emm, SSetts.Exm, SSetts.Emx,
                                                        SSetts.Eom, SSetts.Emo,
                                                        SSetts.Eoo, SSetts.Eox,
                                                        SSetts.Exo, SSetts.Exx])])
            
        dd = numpy.nan*numpy.ones(numpy.max(draw_pord.keys())+1)
        for (p,v) in enumerate([SSetts.Eoo, SSetts.Eox, SSetts.Exo, SSetts.Exx]):
            dd[v] = p

        css = {"fontprops": fontprops, "draw_pord": draw_pord, "draw_ppos": dd, "shape": defaults["shape"], "colhigh": colhigh,
               "delta_on": self.getSettV('draw_delta', self.DELTA_ON)}
        for (p, iid) in enumerate([SSetts.Exo, SSetts.Eox, SSetts.Exx, SSetts.Eoo]):
            css[iid] = {"color_f": self.getColorA(colors[p]),
                        "color_e": self.getColorA(colors[p], 1.),
                        "shape": defaults["shape"], "size": defaults["size"],
                        "zord": self.DEF_ZORD}
        for (p, iid) in enumerate([SSetts.Exm, SSetts.Emx]):
            css[iid] = {"color_f": self.getColorA(colors[SSetts.Eoo], -.9),
                        "color_e": self.getColorA(colors[p], .9),
                        "shape": defaults["shape"], "size": defaults["size"]-1,
                        "zord": zord_miss}
        for (p, iid) in enumerate([SSetts.Emo, SSetts.Eom]):
            css[iid] = {"color_f": self.getColorA(colors[p], -.9),
                        "color_e": self.getColorA(colors[SSetts.Eoo], .9),
                        ## "color_e": self.getColorA(defaults["color_e"], .9),
                        "shape": defaults["shape"], "size": defaults["size"]-1,
                        "zord": zord_miss}
        css[SSetts.Emm] = {"color_f": self.getColorA(colors[SSetts.Eoo], -.9),
                           "color_e": self.getColorA(colors[SSetts.Eoo], .9),
                           "shape": defaults["shape"], "size": defaults["size"]-1,
                           "zord": zord_miss}
        # css[SSetts.Eoo] = {"color_f": self.getColorA(defaults["color_f"]),
        #                      "color_e": self.getColorA(defaults["color_e"], 1.),
        #                      "color_l": self.getColorA(defaults["color_l"]),
        #                      "shape": defaults["shape"], "size": defaults["size"]-1,
        #                      "zord": self.DEF_ZORD}
        css[-1] = {"color_f": self.getColorA(defaults["color_f"], .5),
                   "color_e": self.getColorA(defaults["color_e"], .5),
                   "shape": defaults["shape"], "size": defaults["size"]-1,
                   "zord": self.DEF_ZORD}
        css["default"] = defaults
        css[SSetts.Exo]["zord"] += 1
        css[SSetts.Eox]["zord"] += 1
        css[SSetts.Exx]["zord"] += 2
        css[SSetts.Eoo]["zord"] -= 1
        # print "---- COLOR SETTINGS"
        # for k,v in css.items():
        #     if type(k) is int:
        #         print "* %s" % k
        #         for kk,vv in v.items():
        #             print "\t%s\t%s" % (kk,vv)
        return css

    #### SEC: DATA HANDLING
    ###########################################   
    def setCurrent(self, qr=None):
        return self.getPltDtH().setCurrent(qr)

    def isSingleVar(self):
        return False
        
    def refresh(self):
        self.getLayH().autoShowSplitsBoxes()
        if self.isIntab():
            self.getLayH()._SetSize()

    def addStamp(self, pref=""):
        pass
            
class ViewEntitiesProj(ViewBasis):
    
    TID = "-"
    SDESC = "-"
    ordN = 0
    what = "entities"
    title_str = "Entities Projection"
    typesI = "evr"
    defaultViewT = ProjFactory.defaultView.PID + "_" + what

    subcl_drawer = DrawerEntitiesProj
    
    @classmethod
    def getViewsDetails(tcl):
        return ProjFactory.getViewsDetails(tcl, what=tcl.what)
    
    def __init__(self, parent, vid, more=None):
        self.parent = parent
        self.vid = vid
        self.data = {}
        self.initProject(more)
        self.layH = self.subcl_layh(self)
        self.pltdtH = self.subcl_pltdt(self)
        self.drawer = self.subcl_drawer(self)
        self.layH.setToolbarDrawer(self.drawer)
        boxes = self.additionalElements()
        boxes.extend(self.drawer.prepareInteractive(self.layH.getPanel()))
        self.layH.finalize_init(boxes)

    def getShortDesc(self):
        return "%s %s" % (self.getItemId(), self.getProj().SDESC)

    def getTitleDesc(self):
        return "%s %s" % (self.getItemId(), self.getProj().getTitle())

    def getId(self):
        return (self.getProj().PID, self.vid)
            
    def lastStepInit(self, blocking=False):
        if not self.wasKilled():
            if self.getProj().getCoords() is None:
                self.runProject(blocking)
            else:
                self.readyProj(self.proj)

    def getProj(self):
        return self.proj
            
    def OnReproject(self, rid=None, blocking=False):
        self.getProj().initParameters(self.boxes)
        # self.getProj().addParamsRandrep()
        # tmpp_id = self.projkeyf.GetValue().strip(":, ")
        # if (self.proj is None and len(tmpp_id) > 0) or tmpp_id != self.proj.getCode():
        #     self.initProject(tmpp_id)
        # else:
        #     self.initProject()
        self.runProject(blocking)

    def initProject(self, rid=None):
        ### print ProjFactory.dispProjsInfo()
        self.proj = ProjFactory.getProj(self.getParentData(), rid)
        
    def runProject(self, blocking=False):
        self.drawer.init_wait()
        if self.drawer.hasElement("rep_butt"):
            self.drawer.getElement("rep_butt").Disable()
            self.drawer.getElement("rep_butt").SetLabel("Wait...")
        if self.getPltDtH().hasQueries():
            self.getProj().addParamsRandrep({"vids": self.getPltDtH().getQCols()})
        if blocking:
            try:
                self.proj.do()
            except ValueError as e: #Exception as e:
                self.proj.clearCoords()
            self.readyProj(self.proj)
        else:
            self.parent.project(self.getProj(), self.getId())
        
    def readyProj(self, proj):
        if proj is not None:
            self.proj = proj
        elif self.proj is not None:
            self.proj.clearCoords()
        self.drawer.kill_wait()
        self.drawer.update()
        if self.drawer.hasElement("rep_butt"):
            self.drawer.getElement("rep_butt").Enable()
            self.drawer.getElement("rep_butt").SetLabel("Reproject")
            
    def makeBoxes(self, frame, proj):
        boxes = []
        for kp in proj.getTunableParamsK():
            label = wx.StaticText(frame, wx.ID_ANY, kp.replace("_", " ").capitalize()+":")
            ctrls = []
            value = proj.getParameter(kp)
            if type(value) in [int, float, str]:
                type_ctrl = "text"
                ctrls.append(wx.TextCtrl(frame, wx.NewId(), str(value)))
            elif type(value) is bool:
                type_ctrl = "checkbox" 
                ctrls.append(wx.CheckBox(frame, wx.NewId(), "", style=wx.ALIGN_RIGHT))
                ctrls[-1].SetValue(value)
            elif type(value) is list and kp in proj.options_parameters:
                type_ctrl = "checkbox"
                for k,v in proj.options_parameters[kp]:
                    ctrls.append(wx.CheckBox(frame, wx.NewId(), k, style=wx.ALIGN_RIGHT))
                    ctrls[-1].SetValue(v in value)
            elif kp in proj.options_parameters:
                type_ctrl = "choice" 
                ctrls.append(wx.Choice(frame, wx.NewId()))
                strs = [k for k,v in proj.options_parameters[kp]]
                ctrls[-1].AppendItems(strings=strs)
                try:
                    ind = strs.index(value)
                    ctrls[-1].SetSelection(ind)
                except ValueError:
                    pass
            boxes.append({"key": kp, "label": label, "type_ctrl": type_ctrl, "ctrls":ctrls, "value":value})
        return boxes
    
    def additionalElements(self):
        setts_boxes = []
        max_w = self.getLayH().getFWidth()-50
        current_w = 1000
        flags = wx.ALIGN_CENTER | wx.ALL

        self.boxes = self.makeBoxes(self.getLayH().getPanel(), self.getProj())
        # self.boxes = self.getProj().makeBoxes(self.panel)
        self.boxes.sort(key=lambda x : x["type_ctrl"])
        for box in self.boxes:
            block_w = box["label"].GetBestSize()[0] + sum([c.GetBestSize()[0] for c in box["ctrls"]])
            if current_w + block_w + 10 > max_w:
                setts_boxes.append(wx.BoxSizer(wx.HORIZONTAL))
                setts_boxes[-1].AddSpacer((10,-1))
                current_w = 10
            current_w += block_w + 10
            box["label"].SetFont(wx.Font(8, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
            setts_boxes[-1].Add(box["label"], 0, border=0, flag=flags | wx.ALIGN_RIGHT)
            for c in box["ctrls"]:
                c.SetFont(wx.Font(8, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
                setts_boxes[-1].Add(c, 0, border=0, flag=flags | wx.ALIGN_BOTTOM | wx.ALIGN_LEFT)
            setts_boxes[-1].AddSpacer((10,-1))

        self.nbadd_boxes = len(setts_boxes) 
        return setts_boxes
            
class ViewRed(ViewBasis):

    TID = None
    SDESC = "-"
    title_str = "2D"
    ordN = 1
    geo = False #True
    typesI = "vr"

    subcl_layh = LayoutHandlerQueries
    subcl_drawer = DrawerEntitiesTD
    subcl_pltdt = PltDtHandlerRed
    
    def updateQuery(self, sd=None, query=None):
        return self.getPltDtH().updateQuery(sd, query)
            
    def isSingleVar(self):
        return (len(self.data["queries"][0]) == 0 and self.data["queries"][1].isBasis(1, self.getParentData())) or \
          (len(self.data["queries"][1]) == 0 and self.data["queries"][0].isBasis(0, self.getParentData()))

    def addStamp(self, pref=""):
        self.getDrawer().addStamp(pref)

class ViewRedMap(ViewRed):

    TID = "MAP"
    SDESC = "Map"
    title_str = "Map"
    ordN = 1
    geo = True
    typesI = "vr"

    subcl_layh = LayoutHandlerQueries
    subcl_drawer = DrawerEntitiesMap
    subcl_pltdt = PltDtHandlerRedWithCoords

class ViewRedMappoly(ViewRed):

    TID = "MPP"
    SDESC = "MPoly"
    title_str = "MPoly"
    ordN = 1
    geo = True
    typesI = "r"

    subcl_layh = LayoutHandlerQueries
    subcl_drawer = DrawerEntitiesMappoly
    subcl_pltdt = PltDtHandlerRedWithCoords

    
class ViewRedPara(ViewRed):
    
    TID = "PC"
    SDESC = "Pa.Co."
    ordN = 2
    title_str = "Parallel Coordinates"
    typesI = "vr"
    geo = False

    subcl_drawer = DrawerRedPara

class ViewRedCorrel(ViewRed):
    
    TID = "CC"
    SDESC = "Correl"
    ordN = 6
    title_str = "Variable Correlations"
    typesI = "r"
    geo = False

    subcl_drawer = DrawerRedCorrel
    
class ViewRedTree(ViewRed):

    TID = "TR"
    SDESC = "Tree"
    ordN = 5
    title_str = "Decision Tree"
    typesI = "vr"
    geo = False
    
    subcl_drawer = DrawerRedTree
    subcl_pltdt = PltDtHandlerRed
    
    @classmethod
    def suitableView(tcl, geo=False, what=None, tabT=None):
        return (tabT is None or tabT in tcl.typesI) and (not tcl.geo or geo) and \
               ( what is None or (what[0].isTreeCompatible() and what[1].isTreeCompatible()))


class ViewRedProj(ViewEntitiesProj, ViewRed):

    TID = "EPJ"
    SDESC = "E.Proj."
    what = "entities"
    title_str = "Entities Projection"
    ordN = 10

    subcl_layh = LayoutHandlerQueries
    subcl_pltdt = PltDtHandlerRedWithCoords
    subcl_drawer = DrawerEntitiesProj
    
               
class ViewList(ViewBasis):
    
    TID = "L"
    SDESC = "LViz"
    ordN = 0
    title_str = "List View"
    geo = False
    typesI = "r"

    
    @classmethod
    def suitableView(tcl, geo=False, what=None, tabT=None):
        return tabT is None or tabT in tcl.typesI


class ViewClustMap(ViewList):
    
    TID = "CLM"
    SDESC = "CluMapLViz"
    ordN = 0
    title_str = "Map"
    typesI = "r"
    geo = True
    
    subcl_drawer = DrawerClustMap
    subcl_pltdt = PltDtHandlerListClust
    subcl_layh = LayoutHandlerBasis


class ViewClustProj(ViewEntitiesProj, ViewList):


    TID = "CLP"
    SDESC = "CluProjLViz"
    ordN = 10
    what = "cluster"
    title_str = "Cluster Proj View"
    typesI = "r"
    geo = False

    subcl_drawer = DrawerClustProj
    subcl_pltdt = PltDtHandlerListClust
    subcl_layh = LayoutHandlerBasis

