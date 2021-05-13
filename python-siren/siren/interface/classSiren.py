import os, os.path, random, re, numpy, glob
import time, math
import sys
import pickle

import wx
# import wx.lib.agw.pybusyinfo as PBI

### from wx import AboutBox, AboutDialogInfo, Bitmap, BoxSizer, BusyInfo, Button, CallLater, DefaultPosition, DisplaySize, FileDialog, Frame, Gauge, GridBagSizer, GridSizer, Icon
### from wx import ALIGN_CENTER, BITMAP_TYPE_PNG, CANCEL, CHANGE_DIR, DEFAULT_FRAME_STYLE, VERTICAL, VSCROLL, YES_NO, EXPAND, GA_HORIZONTAL, GA_SMOOTH, HSCROLL
### from wx import NB_TOP,  NO_BORDER, NO_DEFAULT, OK, OPEN, SAVE, STAY_ON_TOP, TE_MULTILINE, TE_READONLY, TE_RICH, TOP
### from wx import EVT_BUTTON, EVT_CLOSE, EVT_ENTER_WINDOW, EVT_LEFT_UP, EVT_MENU, EVT_NOTEBOOK_PAGE_CHANGED, EVT_SIZE, EVT_SPLITTER_UNSPLIT
### from wx import ICON_EXCLAMATION, ICON_INFORMATION
### from wx import ID_ABOUT, ID_COPY, ID_CUT, ID_EXIT, ID_HELP, ID_NO, ID_OK, ID_OPEN, ID_PASTE, ID_PREFERENCES, ID_SAVE, ID_SAVEAS
### from wx import Menu, MenuBar, MessageDialog, NewId, Notebook, NullBitmap, Panel, ScrolledWindow, SplitterWindow, StaticBitmap, TextCtrl, ToggleButton

import wx.lib.dialogs

from ..reremi.toolLog import Log
from ..reremi.classRedescription import Redescription
from ..reremi.classData import Data, DataError
from ..reremi.classConstraints import Constraints
from ..reremi.classBatch import Batch
from ..reremi.toolICList import ICList

from DataWrapper import DataWrapper, findFile
from classGridTable import VarTable, RowTable
from classCtrlTable import RedsManager, VarsManager
from classPreferencesDialog import PreferencesDialog
from classConnectionDialog import ConnectionDialog
from classSplitDialog import SplitDialog
from miscDialogs import ImportDataCSVDialog, ExportFigsDialog, FindDialog, MultiSelectorDialog, ChoiceElement
from ..views.factView import ViewFactory
from ..views.classVizManager import VizManager
from ..views.classViewsManager import ViewsManager
from ..work.toolWP import WorkPlant
from ..work.classWorkClient import WorkClient
from ..common_details import common_variables

import pdb

# try:
#     import test_code
# except ImportError:
#     test_code = None

  
def getRandomColor():
    return (random.randint(0,255), random.randint(0,255), random.randint(0,255))

class ERCache():

    def __init__(self, parent):
        self.parent = parent
        self.rids = []
        self.etor = None
        self.ddER = None

    def getRids(self):
        return self.rids

    def gatherReds(self):
        reds_map = self.parent.getRTab().getAllReds()
        self.setReds(reds_map)

    def needsRecompute(self):
        if self.etor is None:
            return True
        reds_map = self.parent.getRTab().getAllReds()
        rids = [rid for (rid, red) in reds_map]
        return rids != self.getRids()
        
        
    def setReds(self, reds_map):
        self.etor = None
        self.ddER = None
        reds = dict(reds_map)
        self.rids = [rid for (rid, red) in reds_map]
        
        nbE = 0
        if len(self.rids) > 0:
            nbE = reds[self.rids[0]].sParts.nbRows()
        self.etor = numpy.zeros((nbE, len(self.rids)), dtype=bool)
        for r, rid in enumerate(self.rids):
            self.etor[list(reds[rid].getSuppI()), r] = True

    def getRPos(self, rids):
        r_to_p = dict([(r,p) for (p,r) in enumerate(self.rids)])
        return [r_to_p[rid] for rid in rids if rid in r_to_p]
            
    def getEtoR(self, rids=None, eids=None):
        if self.needsRecompute():
            self.gatherReds()
            
        sub_etor = self.etor
        if eids is not None:
            sub_etor = sub_etor[eids,:]
        if rids is not None:
            ps = self.getRPos(rids)
            return sub_etor[:,ps]
        return sub_etor
    
    def computeDeduplicateER(self, etor=None):
        if etor is None:
            if self.etor is None:
                self.gatherReds()
            etor = self.etor
            
        #### DE-DUPLICATE REDS
        rb = {}
        for x,y in zip(*numpy.where(numpy.triu(numpy.dot(1.-etor.T, 1.-etor) + numpy.dot(etor.T*1., etor*1.),1)>=etor.shape[0])):
            rb[x] = y
        keep_rs = numpy.array([i for i in range(etor.shape[1]) if i not in rb], dtype=int)

        r_to_rep = -numpy.ones(etor.shape[1], dtype=int)
        r_to_rep[keep_rs] = numpy.arange(len(keep_rs))
        if len(rb) > 0:
            rfrm, rtt = zip(*rb.items())
            r_to_rep[numpy.array(rfrm)] = r_to_rep[numpy.array(rtt)]

        #### DE-DUPLICATE ENTITIES
        etor = etor[:,keep_rs]
        c01 = numpy.dot(1.-etor, 1.-etor.T) + numpy.dot(etor*1., etor.T*1.)
        nbr = etor.shape[1]
        eb = {}
        for x,y in zip(*numpy.where(numpy.triu(c01,1)>=nbr)):
            eb[x] = y            
        keep_es = numpy.array([i for i in range(etor.shape[0]) if i not in eb and numpy.sum(etor[i,:])>0], dtype=int)

        e_to_rep = -numpy.ones(etor.shape[0], dtype=int)
        e_to_rep[keep_es] = numpy.arange(len(keep_es))
        if len(eb) > 0:
            efrm, ett = zip(*eb.items())
            e_to_rep[numpy.array(efrm)] = e_to_rep[numpy.array(ett)]

        dists = nbr - c01[keep_es,:][:,keep_es]
        return {"e_rprt": keep_es, "r_rprt": keep_rs, "r_to_rep": r_to_rep, "e_to_rep": e_to_rep, "dists": dists,
                    "has_dup_r": len(rb) > 0, "has_dup_e": len(eb) > 0}

    def getDeduplicateER(self, rids=None, eids=None):
        if self.etor is None:
            self.gatherReds()

        if eids is None:
            eids = numpy.arange(self.etor.shape[0])
        if rids is None:
            rids = numpy.arange(self.etor.shape[1])
        ps = self.getRPos(rids)
        return self.computeDeduplicateER(self.etor[eids,:][:,ps])


class ProjCache():

    def __init__(self, capac=10):
        self.cache = {}
        ### for now, unused
        self.capac = capac

    def queryPC(self, proj, vid):
        phsh = "%s:%s" % (vid[0], proj.getParamsHash())
        if phsh in self.cache:
            if self.cache[phsh]["coords"] is not None:
                proj.setCoords(self.cache[phsh]["coords"])
                self.cache[phsh]["served"].append((vid, proj)) 
                return 0
            else:
                self.cache[phsh]["waiting"].append((vid, proj)) 
                return 1
        elif "random_state" not in proj.getParameters():
            self.cache[phsh] = {"coords": None, "waiting": [], "served": []}
        return -1

    def incomingPC(self, proj, vid):
        phsh = "%s:%s" % (vid[0], proj.getParamsHash())
        if phsh in self.cache:
            self.cache[phsh]["coords"] = proj.getCoords()
            while len(self.cache[phsh]["waiting"]) > 0:
                tmp =  self.cache[phsh]["waiting"].pop()
                tmp[1].setCoords(self.cache[phsh]["coords"])
                self.cache[phsh]["served"].append(tmp)
            return self.cache[phsh]["served"] 
        return []
 
class Siren():
    """ The main frame of the application
    """
    
    siren_srcdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if sys.platform == 'win32' and siren_srcdir.find("library.zip") > -1:
        siren_srcdir = re.sub("library.zip.*", "", siren_srcdir)
    if 'SIRENDATA' in os.environ:
        siren_datadir = os.environ['SIRENDATA']
        if not os.path.isdir(siren_datadir):
            raise RuntimeError('Path in environment SIRENDATA not a directory')
    else:
        siren_datadir = os.path.join(siren_srcdir, 'data')
        if not os.path.isdir(siren_datadir):
            siren_datadir = siren_srcdir    


        
    @classmethod
    def searchData(tcl, filen, folder=None, path=[], pref='data/'):
        ff = []
        for fldr in path:
            ff += ['../'+fldr, tcl.siren_srcdir+'/'+fldr]
        if folder is not None:
            ff += ['../'+pref+folder, tcl.siren_datadir+'/'+folder]
        return findFile(filen, path+ff)
        # tmp = findFile(filen, ff)
        # if tmp is None:
        #    raise Exception("findFile %s %s\t-->\t%s" % (filen, path+ff, tmp))
        #    print "findFile %s %s\t-->\t%s" % (filen, ff, tmp)
        # return tmp
    
    @classmethod
    def initIcons(tcl, icons_setts, path=[]):
        icons = {}
        for icon_name, icon_file in icons_setts.items():
            tmp = tcl.searchData(icon_file+".png", 'icons', path)
            if tmp is not None:    
                icons[icon_name] = wx.Bitmap(tmp)
            else:
                icons[icon_name] = wx.NullBitmap
        return icons    
    @classmethod
    def initConfs(tcl, cfiles):
        conf_defs = []
        for (filen, folder) in cfiles:
            cf = tcl.searchData(filen, "conf", path=[folder])
            if cf is not None:
                conf_defs.append(cf)
        return conf_defs
    @classmethod
    def initRedProps(tcl, ffiles):
        fields_fns = []
        for filen, folder in ffiles:
            cf = tcl.searchData(filen, "conf", path=[folder])
            if cf is not None:                
                fields_fns.append(cf)
                for mf in glob.glob(re.sub(filen, tcl.ff_mtch, cf)):
                    if mf != cf:
                        fields_fns.append(mf)
        Redescription.setupRP(fields_fns)
    
    titleTool = common_variables["PROJECT_NAME"]+' :: tools'
    titlePref = common_variables["PROJECT_NAME"]+' :: '
    titleHelp = common_variables["PROJECT_NAME"]+' :: help'
    helpInternetURL = common_variables["PROJECT_URL"]+'help'
    
    # For About dialog
    name = common_variables["PROJECT_NAME"]    
    programURL = common_variables["PROJECT_URL"]
    version = common_variables["VERSION"]
    cpyright = '(c) '+common_variables["COPYRIGHT_YEAR_FROM"]+'-' \
               +common_variables["COPYRIGHT_YEAR_TO"]+' ' \
               +common_variables["PROJECT_AUTHORS"]
    about_text = common_variables["PROJECT_DESCRIPTION_LINE"]+"\n"

    icons_setts = {"split_frame": "split",
                   "unsplit_frame": "unsplit",
                   "learn_act": "learn_act", 
                   "test_act": "test_act",
                   "learn_dis": "learn_dis",
                   "test_dis": "test_dis",
                   "kil": "cross",
                   "stamp": "stamp",
                   "inout": "up_right",
                   "outin": "down_right",
                   "save": "savefig"}

    main_tabs_ids = {"r": "reds", "e": "rows", "t": "log", "z": "viz", "v0": 0, "v1": 1, "v": "vars"}

    external_licenses = ['basemap', 'matplotlib', 'python', 'wx', 'grako']
    ffiles = [('fields_defs_basic.txt', 'reremi')]
    ff_mtch = 'fields_defs_*.txt'
    cfiles = [('miner_confdef.xml', 'reremi'), ('views_confdef.xml', 'views'), ('ui_confdef.xml', 'interface')]
    cfiles_io = [('inout_confdef.xml', 'reremi')]
    
    results_delay = 1000
         
    def __init__(self):
        self.initialized = True
        self.busyDlg = None
        self.findDlg = None
        self.proj_cache = ProjCache()
        self.er_cache = ERCache(self)
        self.dw = None
        self.vizm = None
        self.plant = WorkPlant()
        self.viewsm = ViewsManager(self)

        self.conf_defs = Siren.initConfs(self.cfiles)
        self.conf_defs_io = Siren.initConfs(self.cfiles_io)
        self.icon_file = Siren.searchData('siren_icon32x32.png', 'icons')
        self.license_file = Siren.searchData('LICENSE', 'licenses')
        self.helpURL = Siren.searchData('index.html', 'help')

                    # {"id": self.getDefaultTabId("v0"), "title":"LHS Variables",
                    #  "short": "LHS", "type":"v", "hide":False, "style":None},
                    # {"id": self.getDefaultTabId("v1"), "title":"RHS Variables",
                    #  "short": "RHS", "type":"v", "hide":False, "style":None},
        
        tmp_tabs = [{"id": self.getDefaultTabId("e"), "title":"Entities",
                     "short": "Ent", "type":"e", "style":None},
                    {"id": self.getDefaultTabId("v"), "title":"Variables",
                     "short": "Vars", "type":"v", "style":None},
                    {"id": self.getDefaultTabId("r"), "title":"Redescriptions",
                     "short": "Reds", "type":"r", "style":None},
                    {"id": self.getDefaultTabId("z"), "title":"Visualizations",
                     "short": "Viz", "type":"z", "style":None},
                    {"id": self.getDefaultTabId("t"), "title":"Log",
                     "short": "Log", "type":"t", "style": wx.TE_READONLY|wx.TE_MULTILINE},
            ]
        
        self.tabs = dict([(p["id"], p) for p in tmp_tabs])
        self.tabs_keys = [p["id"] for p in tmp_tabs]
        stn = self.tabs.keys()[0]
        if self.getDefaultTabId("v") in self.tabs:
            stn = self.getDefaultTabId("v")
        elif self.getDefaultTabId("e") in self.tabs:
            stn = self.getDefaultTabId("e")
        self.selectedTab = self.tabs[stn]

        self.logger = Log()
        self.icons = Siren.initIcons(self.icons_setts)
        Siren.initRedProps(self.ffiles)
        tmp = wx.DisplaySize()
        self.toolFrame = wx.Frame(None, -1, self.titleTool, pos = wx.DefaultPosition,
                                  size=(tmp[0]*0.66,tmp[1]*0.9), style = wx.DEFAULT_FRAME_STYLE)

        self.toolFrame.Bind(wx.EVT_CLOSE, self.OnQuit)
        self.toolFrame.Bind(wx.EVT_SIZE, self.OnSize)
        self.toolFrame.SetIcon(wx.Icon(self.icon_file, wx.BITMAP_TYPE_PNG))

        self.buffer_copy = None
        
        self.call_check = None

        self.create_tool_panel()
        self.changePage(stn)
        
        self.dw = DataWrapper(self.logger, conf_defs=self.conf_defs)

        ### About dialog
        self.info =  wx.AboutDialogInfo()
        self.info.SetName(self.name)
        self.info.SetWebSite(self.programURL)
        self.info.SetCopyright(self.cpyright)
        self.info.SetVersion(self.version)
        self.info.SetIcon(wx.Icon(self.icon_file, wx.BITMAP_TYPE_PNG))
        self.info.SetDescription(self.about_text)
        #with open(self.licence_file) as f:
        #    self.info.SetLicence(f.read())

        self.helpFrame = None
        
        ## Register file reading message functions to DataWrapper
        self.dw.registerStartReadingFileCallback(self.startFileActionMsg)
        self.dw.registerStopReadingFileCallback(self.stopFileActionMsg)
        self.readingFileDlg = None

        ### INITIALISATION OF DATA
        self.toolFrame.Show()

        self.resetConstraints()
        self.resetCoordinates()
        self.reloadAll()

        ### W/O THIS DW THINK IT'S CHANGED!
        self.dw.isChanged = False
        self.plant.setUpCall([self.doUpdates, self.resetLogger])
        self.resetLogger()

        self.initialized = True

    def getTab(self, which):
        if which in self.tabs:
            return self.tabs[which]["tab"]
    def getRTab(self):
        stn = self.getDefaultTabId("r")
        return self.tabs[stn]["tab"]

    def sysTLin(self):
        return sys.platform not in ["darwin", 'win32']

    def isInitialized(self):
        return self.initialized

    def hasDataLoaded(self):
        if self.dw is not None:
            return self.dw.getData() is not None
        return False

    def getReds(self):
        if self.dw is not None:
            return self.dw.getReds()
        return []
    def getNbReds(self):
        if self.dw is not None:
            return self.dw.getNbReds()
        return -1
    def getSaveListInfo(self):
        if self.matchTabType("vr"):
            return self.selectedTab["tab"].getSaveListInfo()
    def getToExportReds(self):
        if self.getDefaultTabId("r") in self.tabs and "tab" in self.tabs[self.getDefaultTabId("r")]:
            return self.tabs[self.getDefaultTabId("r")]["tab"].getItemsToExport(inc_hist=False)
        else:
            return []
    def getNbToExportReds(self):
        if self.getDefaultTabId("r") in self.tabs and "tab" in self.tabs[self.getDefaultTabId("r")]:
            return self.tabs[self.getDefaultTabId("r")]["tab"].getNbItemsToExport(inc_hist=False)
        else:
            return 0
    def getData(self):
        if self.dw is not None:
            return self.dw.getData()

    def updatePackReds(self):
        if self.getDefaultTabId("r") in self.tabs and "tab" in self.tabs[self.getDefaultTabId("r")]:
            rr = self.tabs[self.getDefaultTabId("r")]["tab"].getListsToPack(empty_add_electible=True)
            self.dw.resetRedescriptions(rr)
    def markPackRedsWritten(self):
        if self.getDefaultTabId("r") in self.tabs and "tab" in self.tabs[self.getDefaultTabId("r")]:
            self.tabs[self.getDefaultTabId("r")]["tab"].markSavedPack()

    def getPreferences(self):
        if self.dw is not None:
            return self.dw.getPreferences()
    def getLogger(self):
        return self.logger
    def getVizm(self):
        return self.vizm
    def getViewsm(self):
        return self.viewsm

        
######################################################################
###########     TOOL PANEL
######################################################################
## main panel, contains GridTables for the variables and redescriptions plus settings, log, etc.
        
    def create_tool_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas 
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.makeStatus(self.toolFrame)
        self.doUpdates()
        if self.hasSplit():
            self.splitter = wx.SplitterWindow(self.toolFrame)
            self.tabbed = wx.Notebook(self.splitter, -1, style=(wx.NB_TOP)) #, size=(3600, 1200))
        else:
            self.tabbed = wx.Notebook(self.toolFrame, -1, style=(wx.NB_TOP)) #, size=(3600, 1200))
        # self.tabbed.Bind(wx.EVT_LEFT_DOWN, self.testLeftD)

        tmp_keys = list(self.tabs_keys)
        #### Draw tabs
        for tab_id in tmp_keys:
            if self.tabs[tab_id]["type"] == "r":
                self.tabs[tab_id]["tab"] = RedsManager(self, tab_id, self.tabbed, self.tabs[tab_id]["short"])
                self.tabbed.AddPage(self.tabs[tab_id]["tab"].getSW(), self.tabs[tab_id]["title"])

            elif self.tabs[tab_id]["type"] == "v":
                self.tabs[tab_id]["tab"] = VarsManager(self, tab_id, self.tabbed, self.tabs[tab_id]["short"])
                self.tabbed.AddPage(self.tabs[tab_id]["tab"].getSW(), self.tabs[tab_id]["title"])

            # elif self.tabs[tab_id]["type"] == "v":
            #     self.tabs[tab_id]["tab"] = VarTable(self, tab_id, self.tabbed, self.tabs[tab_id]["short"])
            #     self.tabbed.AddPage(self.tabs[tab_id]["tab"].grid, self.tabs[tab_id]["title"])

            elif self.tabs[tab_id]["type"] == "e":
                self.tabs[tab_id]["tab"] = RowTable(self, tab_id, self.tabbed, self.tabs[tab_id]["short"])
                self.tabbed.AddPage(self.tabs[tab_id]["tab"].grid, self.tabs[tab_id]["title"])

            elif self.tabs[tab_id]["type"] == "t":
                self.tabs[tab_id]["tab"] = wx.Panel(self.tabbed, -1)
                self.tabs[tab_id]["text"] = wx.TextCtrl(self.tabs[tab_id]["tab"], size=(-1,-1), style=self.tabs[tab_id]["style"])
                self.tabbed.AddPage(self.tabs[tab_id]["tab"], self.tabs[tab_id]["title"])
                boxS = wx.BoxSizer(wx.VERTICAL)
                boxS.Add(self.tabs[tab_id]["text"], 1, wx.ALIGN_CENTER | wx.TOP | wx.EXPAND)
                self.tabs[tab_id]["tab"].SetSizer(boxS)

            elif self.tabs[tab_id]["type"] == "z":
                self.vizm = VizManager(self, tab_id, self.tabbed, self.tabs[tab_id]["title"]) 
                #self.tabs[tab_id]["tab"] = wx.Panel(self.tabbed, -1)
                self.tabs[tab_id]["tab"] = self.vizm.getSW()
                self.tabbed.AddPage(self.tabs[tab_id]["tab"], self.tabs[tab_id]["title"])

        self.tabbed.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChanged)
        ### splitter
        if self.hasSplit():
            self.splitter.Initialize(self.tabbed)
            self.splitter.SetSashGravity(1.)
            self.splitter.SetMinimumPaneSize(0)
            self.splitter.Bind(wx.EVT_SPLITTER_UNSPLIT, self.OnSplitchange)
        # # self.splitter.Initialize(self.tabbed)
        # # self.splitter.SplitHorizontally(self.tabbed, self.tabs["viz"]["tab"])
        self.tabbed.Bind(wx.EVT_SIZE, self.OnSize)

    def makeStatus(self, frame):
        ### status bar
        self.statusbar = frame.CreateStatusBar()
        self.statusbar.SetFieldsCount(5)
        self.statusbar.SetStatusWidths([25, 300, 150, -1, 200])
        self.updateDataInfo()
        rect = self.statusbar.GetFieldRect(0)
        self.buttViz = wx.StaticBitmap(self.statusbar, wx.NewId(), self.icons["split_frame"], size=(rect.height+4, rect.height+4))
        # self.buttViz.SetMargins(0, 0)
        # self.buttViz = wx.ToggleButton(self.statusbar, wx.NewId(), "s", style=wx.ALIGN_CENTER|wx.TE_RICH, size=(rect.height+4, rect.height+4))
        # self.buttViz.SetForegroundColour((0,0,0))
        # self.buttViz.SetBackgroundColour((255,255,255))
        self.buttViz.SetPosition((rect.x, rect.y))
        self.buttViz.Bind(wx.EVT_LEFT_UP, self.OnSplitchange)

        self.progress_bar = wx.Gauge(self.statusbar, -1, style=wx.GA_HORIZONTAL|wx.GA_SMOOTH)
        rect = self.statusbar.GetFieldRect(2)
        self.progress_bar.SetPosition((rect.x+2, rect.y+2))
        self.progress_bar.SetSize((rect.width-2, rect.height-2))
        self.progress_bar.Hide()

    def matchTabType(self, ffilter, tab=None):
        if tab is None:
            tab = self.selectedTab
        if "tab" in tab:
            return tab["type"].lower() in ffilter  
        return False
    def getTabsMatchType(self, ffilter):
        return [(ti,tab) for (ti, tab) in self.tabs.items() if self.matchTabType(ffilter, tab)]
    def getDefaultTabId(self, tabType):
        return self.main_tabs_ids[tabType]
    def getDefaultTab(self, tabType):
        if self.getDefaultTabId(tabType) in self.tabs and "tab" in self.tabs[self.getDefaultTabId(tabType)]:
            return self.tabs[self.getDefaultTabId(tabType)]["tab"]

        
######################################################################
###########     MENUS
######################################################################

    def updateProgressBar(self):
        if not self.plant.getWP().isActive():
            work_estimate, work_progress = (0, 0)
        else:
            work_estimate, work_progress = self.plant.getWP().getWorkEstimate()
        # print "PROGRESS", work_estimate, work_progress, type(work_estimate)
        if work_estimate > 0:
            self.progress_bar.SetRange(10**5)
            self.progress_bar.SetValue(math.floor(10**5*(work_progress/float(work_estimate))))
            self.progress_bar.Show()
        else:
            self.progress_bar.SetRange(1)
            self.progress_bar.SetValue(0)
            self.progress_bar.Hide()

    def makePopupMenu(self, frame):
        """
        Create and display a popup menu on right-click event
        """
        if self.matchTabType("t"):
            return
        # Popup the menu.  If an item is selected then its handler
        # will be called before PopupMenu returns.
        ct = 0
        menuCon = self.makeConSMenu(frame)
        if menuCon.GetMenuItemCount() > ct:
            ct = menuCon.GetMenuItemCount()
            menuCon.AppendSeparator()
        if self.matchTabType("ev"):
            self.makeRedMenuEV(frame, menuCon)
        elif self.matchTabType("r"):
            self.makeRedMenuR(frame, menuCon)
        if menuCon.GetMenuItemCount() > ct:
            ct = menuCon.GetMenuItemCount()
            menuCon.AppendSeparator()
        if self.matchTabType("evr"):
            self.makeVizMenu(frame, menuCon)

        frame.PopupMenu(menuCon)
        menuCon.Destroy()

    def makeConSMenu(self, frame, menuCon=None):
        if menuCon is None:
            menuCon = wx.Menu()
        if self.matchTabType("e"):            
            ID_FOC = wx.NewId()
            m_foc = menuCon.Append(ID_FOC, "Expand/Shrink column", "Expand/Shrink current column.")
            frame.Bind(wx.EVT_MENU, self.OnFlipExCol, m_foc)
        return menuCon

    def makeRedMenuEV(self, frame, menuRed=None):
        if menuRed is None:
            menuRed = wx.Menu()
        
        if self.matchTabType("ev"):
            if "tab" in self.selectedTab and self.selectedTab["tab"].GetNumberRows() > 0:

                if self.matchTabType("v") and self.selectedTab["tab"].hasFocusItemsL():
                    ID_DETAILS = wx.NewId()
                    m_details = menuRed.Append(ID_DETAILS, "View details", "View variable values.")
                    frame.Bind(wx.EVT_MENU, self.OnShowCol, m_details)

                if self.matchTabType("e"):
                    ID_HIGH = wx.NewId()
                    m_high = menuRed.Append(ID_HIGH, "Highlight in views", "Highlight the entity in all opened views.")
                    if self.viewsm.getNbViews() < 1:
                        menuRed.Enable(ID_HIGH, False)
                    frame.Bind(wx.EVT_MENU, self.OnHigh, m_high)

                if self.matchTabType("v") and self.selectedTab["tab"].nbItems() > 0:
                    ID_SVA = wx.NewId()
                    m_svla = menuRed.Append(ID_SVA, "Save Plots As...", "Save Plots as...")
                    frame.Bind(wx.EVT_MENU, self.OnExportListFigs, m_svla)


                if self.matchTabType("e") or ( self.matchTabType("v") and self.selectedTab["tab"].nbItems() > 0):
                    ID_FIND = wx.NewId()
                    m_find = menuRed.Append(ID_FIND, "Find\tCtrl+F", "Find by name.")
                    frame.Bind(wx.EVT_MENU, self.OnFind, m_find)

                if self.matchTabType("e") or ( self.matchTabType("v") and self.selectedTab["tab"].hasFocusItemsL()):
                    ID_ENABLED = wx.NewId()
                    m_enabled = menuRed.Append(ID_ENABLED, "En&able/Disable\tCtrl+D", "Enable/Disable current item.")
                    frame.Bind(wx.EVT_MENU, self.OnFlipEnabled, m_enabled)

                ID_ENABLEDALL = wx.NewId()
                m_enabledall = menuRed.Append(ID_ENABLEDALL, "&Enable All", "Enable all items.")
                frame.Bind(wx.EVT_MENU, self.OnEnabledAll, m_enabledall)

                ID_DISABLEDALL = wx.NewId()
                m_disabledall = menuRed.Append(ID_DISABLEDALL, "&Disable All", "Disable all items.")
                frame.Bind(wx.EVT_MENU, self.OnDisabledAll, m_disabledall)

        return menuRed

    def makeRedMenuR(self, frame, menuRed=None):
        if menuRed is None:
            menuRed = wx.Menu()
            
        if self.matchTabType("r"):
            if self.selectedTab["tab"].hasFocusContainersL():
                ID_NN = wx.NewId()
                m_newl = menuRed.Append(ID_NN, "New List\tCtrl+N", "New List.")
                frame.Bind(wx.EVT_MENU, self.OnNewList, m_newl)

                ID_ADP = wx.NewId()
                m_adp = menuRed.Append(ID_ADP, "Add to package", "Add/remove package list.")
                if self.getData() is not None and self.dw.isFromPackage and self.dw.getPackageSaveFilename() is not None and not self.selectedTab["tab"].hasFocusCLIsHist():
                    frame.Bind(wx.EVT_MENU, self.OnAddDelListToPack, m_adp)
                    if self.selectedTab["tab"].hasFocusCLInPack():
                        m_adp.SetText("Remove from package")
                else:
                    menuRed.Enable(ID_ADP, False)
                
                ID_SV = wx.NewId()
                m_svl = menuRed.Append(ID_SV, "Save List", "Save List.")
                frame.Bind(wx.EVT_MENU, self.OnSaveList, m_svl)
                if not self.selectedTab["tab"].hasFocusCLFile():
                    menuRed.Enable(ID_SV, False)

                ID_SVA = wx.NewId()
                m_svla = menuRed.Append(ID_SVA, "Save List As...", "Save List as...")
                frame.Bind(wx.EVT_MENU, self.OnSaveListAs, m_svla)
                
                ID_SVA = wx.NewId()
                m_svla = menuRed.Append(ID_SVA, "Save Plots As...", "Save Plots as...")
                frame.Bind(wx.EVT_MENU, self.OnExportListFigs, m_svla)
                

            if  self.selectedTab["tab"].hasFocusItemsL() and self.selectedTab["tab"].nbSelectedItems() == 1:
                ID_DETAILS = wx.NewId()
                m_details = menuRed.Append(ID_DETAILS, "View details", "View variable values.")
                frame.Bind(wx.EVT_MENU, self.OnShowCol, m_details)

            ID_DEL = wx.NewId()
            m_del = menuRed.Append(ID_DEL, "De&lete", "Delete current redescription.")
            frame.Bind(wx.EVT_MENU, self.OnDelete, m_del)                            
            m_cut = menuRed.Append(wx.ID_CUT, "Cu&t", "Cut current redescription.")
            frame.Bind(wx.EVT_MENU, self.OnCut, m_cut)
            m_copy = menuRed.Append(wx.ID_COPY, "&Copy", "Copy current redescription.")
            frame.Bind(wx.EVT_MENU, self.OnCopy, m_copy)
            m_paste = menuRed.Append(wx.ID_PASTE, "&Paste", "Paste current redescription.")
            frame.Bind(wx.EVT_MENU, self.OnPaste, m_paste)

            if self.selectedTab["tab"].nbItems() == 0:
                menuRed.Enable(wx.ID_CUT, False)
                menuRed.Enable(wx.ID_COPY, False)
            if self.selectedTab["tab"].isEmptyBuffer():
                menuRed.Enable(wx.ID_PASTE, False)
                
            if self.selectedTab["tab"].nbItems() > 0:
                ID_FIND = wx.NewId()
                m_find = menuRed.Append(ID_FIND, "Find\tCtrl+F", "Find by name.")
                frame.Bind(wx.EVT_MENU, self.OnFind, m_find)
                    
                ID_ENABLED = wx.NewId()
                m_enabled = menuRed.Append(ID_ENABLED, "En&able/Disable\tCtrl+D", "Enable/Disable current item.")
                frame.Bind(wx.EVT_MENU, self.OnFlipEnabled, m_enabled)

                ID_ENABLEDALL = wx.NewId()
                m_enabledall = menuRed.Append(ID_ENABLEDALL, "&Enable All", "Enable all items.")
                frame.Bind(wx.EVT_MENU, self.OnEnabledAll, m_enabledall)

                ID_DISABLEDALL = wx.NewId()
                m_disabledall = menuRed.Append(ID_DISABLEDALL, "&Disable All", "Disable all items.")
                frame.Bind(wx.EVT_MENU, self.OnDisabledAll, m_disabledall)

                ID_DELETEDALL = wx.NewId()
                m_deletedall = menuRed.Append(ID_DELETEDALL, "Delete Disabled", "Delete all disabled items.")
                frame.Bind(wx.EVT_MENU, self.OnDeleteDisAll, m_deletedall)


            if self.selectedTab["tab"].hasFocusItemsL() and self.selectedTab["tab"].nbSelectedItems() == 1:
                ID_NOR = wx.NewId()
                m_nor = menuRed.Append(ID_NOR, "&Normalize", "Normalize current redescription.")
                frame.Bind(wx.EVT_MENU, self.OnNormalize, m_nor)

                ID_EXPAND = wx.NewId()
                m_expand = menuRed.Append(ID_EXPAND, "E&xpand\tCtrl+E", "Expand redescription.")
                frame.Bind(wx.EVT_MENU, self.OnExpand, m_expand)

            if self.selectedTab["tab"].nbItems() > 0:
                ID_FILTER_ONE = wx.NewId()
                m_filter_one = menuRed.Append(ID_FILTER_ONE, "&Filter redundant to one\tCtrl+R", "Disable redescriptions redundant to current downwards.")
                frame.Bind(wx.EVT_MENU, self.OnFilterToOne, m_filter_one)

                ID_FILTER_ALL = wx.NewId()
                m_filter_all = menuRed.Append(ID_FILTER_ALL, "Filter red&undant all\tShift+Ctrl+R", "Disable redescriptions redundant to previous encountered.")
                frame.Bind(wx.EVT_MENU, self.OnFilterAll, m_filter_all)

                ID_FILTER_FLD = wx.NewId()
                m_filter_fld = menuRed.Append(ID_FILTER_FLD, "Filter on folds cover", "Disable redescriptions that do not adequately cover several folds.")
                frame.Bind(wx.EVT_MENU, self.OnFilterFolds, m_filter_fld)

                ID_PROCESS = wx.NewId()
                m_process = menuRed.Append(ID_PROCESS, "&Process redescriptions\tCtrl+P", "Sort and filter current redescription list.")
                frame.Bind(wx.EVT_MENU, self.OnProcessAll, m_process)

        if menuRed.GetMenuItemCount() == 0:
            self.appendEmptyMenuEntry(menuRed, "No items", "There are no items to edit.")

        return menuRed

    def makeVizMenu(self, frame, menuViz=None):
        if menuViz is None:
            countIts = 0
            menuViz = wx.Menu()
            
            #### not for popup menu
            if self.getVizm() is not None and self.getVizm().hasVizIntab():
                ID_CHECK = wx.NewId()
                m_check = menuViz.AppendCheckItem(ID_CHECK, "Plots in external windows", "Plot in external windows rather than the visualization tab.")
                if not self.getVizm().showVizIntab():
                    m_check.Check()
                frame.Bind(wx.EVT_MENU, self.OnVizCheck, m_check)
                menuViz.AppendSeparator()

        countIts = menuViz.GetMenuItemCount()
            
        queries = None
        if self.matchTabType("e") or ( self.matchTabType("r") and self.selectedTab["tab"].hasFocusItemsL() and self.selectedTab["tab"].nbSelectedItems() == 1 ) or ( self.matchTabType("v") and self.selectedTab["tab"].hasFocusItemsL() and self.selectedTab["tab"].nbSelectedItems() == 1 ):
            if self.matchTabType("r"):
                typ = "r"
                queries = self.selectedTab["tab"].getSelectedQueries()
            elif self.matchTabType("e"):
                typ = "e"
            elif self.matchTabType("v"):
                typ = "v"

            for item in self.viewsm.getViewsItems("R", typ, what=queries):
                ID_NEWV = wx.NewId()
                m_newv = menuViz.Append(ID_NEWV, "%s" % item["title"],
                                          "Plot %s." % item["title"])
                if not item["suitable"]:
                    m_newv.Enable(False)
                frame.Bind(wx.EVT_MENU, self.OnNewV, m_newv)
                self.ids_viewT[ID_NEWV] = item["viewT"]

        if self.matchTabType("r") and self.selectedTab["tab"].hasFocusContainersL()  and self.selectedTab["tab"].nbSelectedLists() > 0:
            for item in self.viewsm.getViewsItems("L", "r"):
                ID_NEWV = wx.NewId()
                m_newv = menuViz.Append(ID_NEWV, "%s" % item["title"],
                                          "Plot %s." % item["title"])
                if not item["suitable"]:
                    m_newv.Enable(False)
                frame.Bind(wx.EVT_MENU, self.OnNewVList, m_newv)
                self.ids_viewT[ID_NEWV] = item["viewT"]

        if ( self.matchTabType("r") and self.selectedTab["tab"].hasFocusItemsL() and self.selectedTab["tab"].nbSelectedItems() > 1 ):
            for item in self.viewsm.getViewsItems("L", "r"):
                ID_NEWV = wx.NewId()
                m_newv = menuViz.Append(ID_NEWV, "%s" % item["title"],
                                          "Plot %s." % item["title"])
                if not item["suitable"]:
                    m_newv.Enable(False)
                frame.Bind(wx.EVT_MENU, self.OnNewVItems, m_newv)
                self.ids_viewT[ID_NEWV] = item["viewT"]

                
        if menuViz.GetMenuItemCount() == countIts:
            self.appendEmptyMenuEntry(menuViz, "No visualization", "There are no visualizations.")

        return menuViz

    def makeProcessMenu(self, frame, menuPro=None):
        if menuPro is None:
            menuPro = wx.Menu()
        ID_MINE = wx.NewId()
        m_mine = menuPro.Append(ID_MINE, "&Mine redescriptions\tCtrl+M", "Mine redescriptions from the dataset according to current constraints.")
        if self.getData() is None:
            menuPro.Enable(ID_MINE, False)
        else:
            frame.Bind(wx.EVT_MENU, self.OnMineAll, m_mine)

        ct = menuPro.GetMenuItemCount()
        menuPro = self.makeStoppersMenu(frame, menuPro)
        if ct < menuPro.GetMenuItemCount():
            menuPro.InsertSeparator(ct)
        return menuPro

    def makeStoppersMenu(self, frame, menuStop=None):
        if menuStop is None:
            menuStop = wx.Menu()
        if self.plant.getWP().nbWorkers() == 0:
            ID_NOP = wx.NewId()
            m_nop = menuStop.Append(ID_NOP, "No process running", "There is no process currently running.")
            menuStop.Enable(ID_NOP, False)

        else:
            for wdt in self.plant.getWP().getWorkersDetails(): 
                ID_STOP = wx.NewId()
                self.ids_stoppers[ID_STOP] = wdt["wid"] 
                m_stop = menuStop.Append(ID_STOP, "Stop %s #&%s" % (wdt["wtyp"], wdt["wid"]), "Interrupt %s process #%s." % (wdt["wtyp"], wdt["wid"]))
                frame.Bind(wx.EVT_MENU, self.OnStop, m_stop)
        if self.plant.getWP().isActive():
            menuStop.AppendSeparator()
            ID_PLT = wx.NewId()
            m_plt = menuStop.Append(ID_PLT, self.plant.getWP().infoStr(), "Where processes are handled.")
            menuStop.Enable(ID_PLT, False)

        return menuStop

    def makeViewsMenu(self, frame, menuViews=None):
        if menuViews is None:
            menuViews = wx.Menu()

        menuViews.AppendMenu(wx.NewId(), "&Tabs", self.makeTabsMenu(frame))
        # self.makeTabsMenu(frame, menuViews)
        # if menuViews.GetMenuItemCount() > 0:
        #     menuViews.AppendSeparator()

        self.viewsm.makeViewsMenu(frame, menuViews) 
        return menuViews

    def makeTabsMenu(self, frame, menuTabs=None):
        if menuTabs is None:
            menuTabs = wx.Menu()

        for tab_id in self.tabs_keys:
            tab_prop = self.tabs[tab_id]
            ID_CHECK = wx.NewId()
            self.check_tab[ID_CHECK] = tab_id 
            m_check = menuTabs.Append(ID_CHECK, "%s" % tab_prop["title"], "Switch to %s tab." % tab_prop["title"])
            frame.Bind(wx.EVT_MENU, self.OnTabW, m_check)
        return menuTabs

    def makeFileMenu(self, frame, menuFile=None):
        if menuFile is None:
            menuFile = wx.Menu()
        # #### WARNING FOR DEBUG ONLY
        # ID_RELOAD = wx.NewId()
        # m_reload = menuFile.Append(ID_RELOAD, "&Reload\tShift+Ctrl+O", "Reload view code.")
        # frame.Bind(wx.EVT_MENU, self.OnReload, m_reload)

        m_open = menuFile.Append(wx.ID_OPEN, "&Open\tCtrl+O", "Open a project.")
        frame.Bind(wx.EVT_MENU, self.OnOpen, m_open)

        ## Save  
        m_save = menuFile.Append(wx.ID_SAVE, "&Save\tCtrl+S", "Save the current project.")
        if self.getData() is not None and self.dw.isFromPackage and self.dw.getPackageSaveFilename() is not None:
            frame.Bind(wx.EVT_MENU, self.OnSave, m_save)
        else:
            menuFile.Enable(wx.ID_SAVE, False)

        ## Save As...
        m_saveas = menuFile.Append(wx.ID_SAVEAS, "Save &As...\tShift+Ctrl+S", "Save the current project as...")
        if self.getData() is None:
            menuFile.Enable(wx.ID_SAVEAS, False)
        else:
            frame.Bind(wx.EVT_MENU, self.OnSaveAs, m_saveas)

        ## Import submenu
        submenuImport = wx.Menu()
        #submenuImportData = wx.Menu()
        ID_IMPORT_DATA_CSV = wx.NewId()
        m_impDataCSV = submenuImport.Append(ID_IMPORT_DATA_CSV, "Import Data", "Import data in CSV format.")
        frame.Bind(wx.EVT_MENU, self.OnImportDataCSV, m_impDataCSV)
        # ID_IMPORT_DATA_XML = wx.NewId()
        # m_impDataXML = submenuImport.Append(ID_IMPORT_DATA_XML, "Import Data from XML", "Import data in XML format.")
        # frame.Bind(wx.EVT_MENU, self.OnImportDataXML, m_impDataXML)
        # ID_IMPORT_DATA_TRIPLE = wx.NewId()
        # m_impDataTriple = submenuImport.Append(ID_IMPORT_DATA_TRIPLE, "Import Data from separate files", "Import data from separate files")
        # frame.Bind(wx.EVT_MENU, self.OnImportData, m_impDataTriple)
        
        # ID_IMPORT_DATA = wx.NewId()
        # m_impData = submenuImport.AppendMenu(ID_IMPORT_DATA, "Import &Data", submenuImportData)
        #m_impData = submenuImport.Append(ID_IMPORT_DATA, "Import &Data", "Import data into the project.")
        #frame.Bind(wx.EVT_MENU, self.OnImportData, m_impData)

        ID_IMPORT_PREFERENCES = wx.NewId()
        m_impPreferences = submenuImport.Append(ID_IMPORT_PREFERENCES, "Import &Preferences", "Import preferences into the project.")
        frame.Bind(wx.EVT_MENU, self.OnImportPreferences, m_impPreferences)
        
        ID_IMPORT_REDESCRIPTIONS = wx.NewId()
        m_impRedescriptions = submenuImport.Append(ID_IMPORT_REDESCRIPTIONS, "Import &Redescriptions", "Import redescriptions into the project.")
        if self.getData() is not None:
            frame.Bind(wx.EVT_MENU, self.OnImportRedescriptions, m_impRedescriptions)
        else:
            submenuImport.Enable(ID_IMPORT_REDESCRIPTIONS, False)

        ID_IMPORT = wx.NewId()
        m_import = menuFile.AppendMenu(ID_IMPORT, "&Import", submenuImport)

        
        ## Export submenu
        submenuExport = wx.Menu() # Submenu for exporting
        
        ID_EXPORT_REDESCRIPTIONS = wx.NewId()
        m_exportRedescriptions = submenuExport.Append(ID_EXPORT_REDESCRIPTIONS, "&Export Redescriptions\tShift+Ctrl+E", "Export redescriptions.")
        if self.getNbToExportReds() < 1:
            submenuExport.Enable(ID_EXPORT_REDESCRIPTIONS, False)
        else:
            frame.Bind(wx.EVT_MENU, self.OnExportRedescriptions, m_exportRedescriptions)

        ID_EXPORT_PREF = wx.NewId()
        m_exportPref = submenuExport.Append(ID_EXPORT_PREF, "&Export Preferences", "Export preferences.")
        frame.Bind(wx.EVT_MENU, self.OnPrintoutPreferences, m_exportPref)

        ID_TMPL_PREF = wx.NewId()
        m_tmplPref = submenuExport.Append(ID_TMPL_PREF, "Export Prefs Template", "Print out preferences template.")
        frame.Bind(wx.EVT_MENU, self.OnPrintoutPreferencesTmpl, m_tmplPref)
        ID_DEF_PREF = wx.NewId()
        m_defPref = submenuExport.Append(ID_DEF_PREF, "Export Prefs Default", "Print out default preferences.")
        frame.Bind(wx.EVT_MENU, self.OnPrintoutPreferencesDef, m_defPref)
        
        ID_EXPORT = wx.NewId()
        m_export = menuFile.AppendMenu(ID_EXPORT, "&Export", submenuExport)

        if True: #test_code:
            ID_RUN_TEST = wx.NewId()
            m_runTest = menuFile.Append(ID_RUN_TEST, "&Run test feature", "Trying a new feature.")
            frame.Bind(wx.EVT_MENU, self.OnRunTest, m_runTest)
            
        ## Preferences
        menuFile.AppendSeparator()
        m_preferencesdia = menuFile.Append(wx.ID_PREFERENCES, "P&references...\tCtrl+,", "Set preferences.")
        frame.Bind(wx.EVT_MENU, self.OnPreferencesDialog, m_preferencesdia)

        ## Worker setup
        if True:
                ID_CONN = wx.NewId()
                m_conndia = menuFile.Append(ID_CONN, "Wor&ker setup...\tCtrl+k", "Setup worker's connection.")
                frame.Bind(wx.EVT_MENU, self.OnConnectionDialog, m_conndia)

        ## Split setup
        if True:
                ID_SPLT = wx.NewId()
                m_spltdia = menuFile.Append(ID_SPLT, "Sp&lits setup...\tCtrl+l", "Setup learn/test data splits.")
                frame.Bind(wx.EVT_MENU, self.OnSplitDialog, m_spltdia)
                if not self.hasDataLoaded():
                        menuFile.Enable(ID_SPLT, False)

        ## Export submenu
        submenuFields = wx.Menu() # Submenu for exporting
        
        ID_FLD_GUI = wx.NewId()
        m_fldgui = submenuFields.Append(ID_FLD_GUI, "GUI", "Fields for GUI.")
        frame.Bind(wx.EVT_MENU, self.OnDefRedsFieldsGUI, m_fldgui)
        ID_FLD_TXT = wx.NewId()
        m_fldtxt = submenuFields.Append(ID_FLD_TXT, "Text export", "Fields for text export.")
        frame.Bind(wx.EVT_MENU, self.OnDefRedsFieldsOutTxt, m_fldtxt)
        ID_FLD_TEX = wx.NewId()
        m_fldtex = submenuFields.Append(ID_FLD_TEX, "LaTeX export", "Fields for LaTeX export.")
        frame.Bind(wx.EVT_MENU, self.OnDefRedsFieldsOutTex, m_fldtex)

        
        ID_FLD = wx.NewId()
        m_fld = menuFile.AppendMenu(ID_FLD, "Fields setup", submenuFields)

        menuFile.AppendSeparator()
        ## Quit
        m_quit = menuFile.Append(wx.ID_EXIT, "&Quit", "Close window and quit program.")
        frame.Bind(wx.EVT_MENU, self.OnQuit, m_quit)
        return menuFile

    def makeHelpMenu(self, frame, menuHelp=None):
        if menuHelp is None:
            menuHelp = wx.Menu()
        m_help = menuHelp.Append(wx.ID_HELP, "C&ontent", "Access the instructions.")
        frame.Bind(wx.EVT_MENU, self.OnHelp, m_help)
        
        m_about = menuHelp.Append(wx.ID_ABOUT, "&About", "About...")
        frame.Bind(wx.EVT_MENU, self.OnAbout, m_about)

        ID_LICENSE = wx.NewId()
        m_license = menuHelp.Append(ID_LICENSE, "&License", "View the license(s).")
        frame.Bind(wx.EVT_MENU, self.OnLicense, m_license)
        return menuHelp


    def appendEmptyMenuEntry(self, menu, entry_text, entry_leg=""):
        ID_NOR = wx.NewId()
        menu.Append(ID_NOR, entry_text, entry_leg)
        menu.Enable(ID_NOR, False)
    def makeMenuEmpty(self, menu_action, menuEmpty=None):
        if menuEmpty is None:
            menuEmpty = wx.Menu()
        self.appendEmptyMenuEntry(menuEmpty, "Nothing to %s" % menu_action, "Nothing to %s." % menu_action)
        return menuEmpty
        
    def makeMenu(self, frame):
        menuBar = wx.MenuBar()
        menuBar.Append(self.makeFileMenu(frame), "&File")
        me = None
        if self.matchTabType("ev"):
            me = self.makeRedMenuEV(frame)
        elif self.matchTabType("r"):
            me = self.makeRedMenuR(frame)
        if me is None or me.GetMenuItemCount() == 0:
            me = self.makeMenuEmpty("edit", menuEmpty=me)
        menuBar.Append(me, "&Edit")
        ## if self.matchTabType("evr"):
        menuBar.Append(self.makeVizMenu(frame), "&View")
        ## else:
        ##     menuBar.Append(self.makeMenuEmpty("visualize"), "&View")
        menuBar.Append(self.makeProcessMenu(frame), "&Process")
        menuBar.Append(self.makeViewsMenu(frame), "&Windows")
        menuBar.Append(self.makeHelpMenu(frame), "&Help")
        frame.SetMenuBar(menuBar)
        frame.Layout()
        frame.SendSizeEvent()

    def updateMenus(self):
        self.ids_viewT = {}
        self.ids_stoppers = {}
        self.check_tab = {}
        self.makeMenu(self.toolFrame)
        self.viewsm.makeMenusForViews()

######################################################################
###########     ACTIONS
######################################################################

    def hasSplit(self):
        return True #False

    def OnSplitchange(self, event):
        if self.hasSplit():
            ## if event.GetEventType() == 10022 or event.GetEventType() == 10029 or \
            ##        event.GetWindowBeingRemoved().GetParent().GetId() == self.splitter.GetId():
            if type(event) == wx.MouseEvent or \
               ( type(event) == wx.SplitterEvent and event.GetWindowBeingRemoved().GetParent().GetId() == self.splitter.GetId()):
                 self.getVizm().OnSplitchange()

    # def OnReload(self, event):
    #     ViewFactory.reloadCode()
        
    def OnOpen(self, event):
        if not self.checkAndProceedWithUnsavedChanges():
            return
        wcd = 'All files|*|Siren packages (*.siren)|*.siren'

        if self.dw.getPackageSaveFilename() is not None:
            dir_name = os.path.dirname(self.dw.getPackageSaveFilename())
        else:
            dir_name = os.path.expanduser('~/')
        path = dir_name            
        open_dlg = wx.FileDialog(self.toolFrame, message='Choose a file', defaultDir=dir_name,  
			wildcard=wcd, style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.LoadFile(path)
        open_dlg.Destroy()
        self.changePage(self.getDefaultTabId("r"))
        # DEBUGGING
        #wx.MessageDialog(self.toolFrame, 'Opened package from '+path).ShowModal()

    def LoadFile(self, path):
        try:
            self.dw.openPackage(path)
        except:
            raise
            return False
        else:
            self.resetConstraints()
            self.loadAll(refresh_rfields=True)
            return True
        
    def expand(self, params={}):
        if params is None:
            params = {}
        self.progress_bar.Show()
        if "red" in params and params["red"] is not None and params["red"].length(0) + params["red"].length(1) > 0:
            self.plant.getWP().addWorker("expander", self, params,
                                 {"results_track":0,
                                  "batch_type": "partial",
                                  "results_tab": "exp"})
        else:
            self.plant.getWP().addWorker("miner", self, params,
                                 {"results_track":0,
                                  "batch_type": "final",
                                  "results_tab": "exp"})
        self.checkResults(menu=True)

    def project(self, proj=None, vid=None):
        self.progress_bar.Show()
        if proj is not None and vid is not None:
            out = self.proj_cache.queryPC(proj, vid)
            if out == 0:
                # print "Found previous proj"
                self.readyProj(None, vid, proj)
            elif out < 0:
                # print "No previous proj"
                wid = self.plant.getWP().findWid([("wtyp", "projector"), ("vid", vid)])
                if wid is None:
                    self.plant.getWP().addWorker("projector", self, proj,
                                         {"vid": vid})
                    self.checkResults(menu=True)
            # else:
            #     print "Waiting previous proj"

    def checkResults(self, menu=False, once=False):
        # print "Check results\tnb working", self.plant.getWP().infoStr(), self.plant.getWP().nbWorking()
        updates = self.plant.getWP().checkResults(self)
        if menu:
            updates["menu"] = True
        if not once: 
            if self.plant.getWP().nbWorking() > 0:
                if self.call_check is None:
                    self.call_check = wx.CallLater(Siren.results_delay, self.checkResults)
                else:
                    self.call_check.Restart(Siren.results_delay)
            else:
                self.call_check = None

        if once or not self.plant.getWP().isActive():
            updates["progress"] = True
        self.doUpdates(updates) ## To update the worker stoppers

    def doUpdates(self, updates=None):
        if updates is None:
            updates={"menu":True }
        if "error" in updates:
            self.errorBox(updates["error"])
            self.appendLog("\n\n***" + updates["error"] + "***\n")
        if "menu" in updates:
            self.updateMenus()
        if "progress" in updates:
            self.updateProgressBar()
        if "status" in updates:
            self.statusbar.SetStatusText(updates["status"], 1)
        if "log" in updates:
            self.appendLog(updates["log"])

    def loggingDWError(self, output, message, type_message, source):
        self.errorBox(message)
        self.appendLog("\n\n***" + message + "***\n")

    def loggingLogTab(self, output, message, type_message, source):
        text = "%s" % message
        header = "@%s:\t" % source
        text = text.replace("\n", "\n"+header)
        self.appendLog(text+"\n")

    def appendLog(self, text):
        if "log" in self.tabs:
            self.tabs["log"]["text"].AppendText(text)


    def OnStop(self, event):
        if event.GetId() in self.ids_stoppers:
            self.plant.getWP().layOff(self.ids_stoppers[event.GetId()])
            self.checkResults(menu=True)
            
    def OnSave(self, event):
        if not (self.dw.isFromPackage and self.dw.getPackageSaveFilename() is not None):
            wx.MessageDialog(self.toolFrame, 'Cannot save data that is not from a package\nUse Save As... instead', style=wx.OK|wx.ICON_EXCLAMATION, caption='Error').ShowModal()
            return
        try:
            self.updatePackReds()
            self.dw.savePackage()
            self.markPackRedsWritten()
        except:
            pass
            
    def OnSaveAs(self, event):
        if self.dw.getPackageSaveFilename() is not None:
            dir_name = os.path.dirname(self.dw.getPackageSaveFilename())
        else:
            dir_name = os.path.expanduser('~/')

        save_dlg = wx.FileDialog(self.toolFrame, message="Save as", defaultDir=dir_name,
                                 style=wx.SAVE|wx.CHANGE_DIR)
        if save_dlg.ShowModal() == wx.ID_OK:
            path = save_dlg.GetPath()
            try:
                self.updatePackReds()
                self.dw.savePackageToFile(path)
                self.markPackRedsWritten()
                self.updateMenus()
            except:
                pass
        save_dlg.Destroy()

    def OnAddDelListToPack(self, event):
        if self.matchTabType("r"):
            lid = self.selectedTab["tab"].getFLLid()
            if lid is not None:
                self.selectedTab["tab"].OnAddDelListToPack(lid)
        
    def OnSaveList(self, event):
        info = self.getSaveListInfo()
        if info is not None:
            if info.get("nb_items", 0) < 1:
                wx.MessageDialog(self.toolFrame, 'Cannot save list: no redescriptions loaded',
                                 style=wx.OK|wx.ICON_EXCLAMATION, caption='Error').ShowModal()
                return

            elif info.get("path") is None:
                wx.MessageDialog(self.toolFrame, 'Cannot save list, no output\nUse Save List As... instead', style=wx.OK|wx.ICON_EXCLAMATION, caption='Error').ShowModal()
                return

            try:
                self.dw.exportRedescriptions(info["path"], info["items"])
                self.selectedTab["tab"].markSavedSrc(lid=info["lid"])
            except:
                pass


    def OnSaveListAs(self, event):
        info = self.getSaveListInfo()
        ## print "Save List As", info
        if info is not None and info.get("nb_items", 0) < 1:
            wx.MessageDialog(self.toolFrame, 'Cannot export redescriptions: no redescriptions loaded',
                             style=wx.OK|wx.ICON_EXCLAMATION, caption='Error').ShowModal()
            return

        if info.get("path") is not None:
            dir_name = os.path.dirname(info["path"])
        elif self.dw.getPackageSaveFilename() is not None:
            dir_name = os.path.dirname(self.dw.getPackageSaveFilename())
        else:
            dir_name = os.path.expanduser('~/')

        save_dlg = wx.FileDialog(self.toolFrame, message='Export redescriptions to:', defaultDir = dir_name, style = wx.SAVE|wx.CHANGE_DIR)
        if save_dlg.ShowModal() == wx.ID_OK:
            new_path = save_dlg.GetPath()
            try:
                self.dw.exportRedescriptions(new_path, info["items"])
                self.selectedTab["tab"].markSavedSrc(lid=info["lid"], path=new_path)
            except:
                pass
        save_dlg.Destroy()

    def OnExportListFigs(self, event):     
        """Shows a custom dialog to open the export files"""
        info = self.getSaveListInfo()
        if info.get("path") is not None:
            dir_name = os.path.dirname(info["path"])
        elif self.dw.getPackageSaveFilename() is not None:
            dir_name = os.path.dirname(self.dw.getPackageSaveFilename())
        else:
            dir_name = os.path.expanduser('~/')        
        if info is not None and info.get("nb_items", 0) < 1:
            wx.MessageDialog(self.toolFrame, 'Cannot export items: no items active',
                             style=wx.OK|wx.ICON_EXCLAMATION, caption='Error').ShowModal()
            return
        dlg = ExportFigsDialog(self, self.viewsm, info["items"], dir_name)
        dlg.showDialog()

    def OnDefRedsFieldsOutTex(self, event):
        self.OnDefRedsFields("tex")
    def OnDefRedsFieldsOutTxt(self, event):
        self.OnDefRedsFields("txt")        
    def OnDefRedsFieldsGUI(self, event):
        flk = "gui"
        tt = self.OnDefRedsFields(flk)
        if tt is not None:
            self.resetFields()            
        
    def OnDefRedsFields(self, flk):
        rp = Redescription.getRP()
        modifiers = rp.getModifiersForData(self.dw.getData())

        choice_list = []
        for k in rp.getAllFields(flk, modifiers):
            choice_list.append((k, ChoiceElement(k, "%s" % k)))
        selected_ids = rp.getCurrentListFields(flk, modifiers)
        
        dlg = MultiSelectorDialog(self, choice_list, selected_ids)
        fields = dlg.showDialog()
        if fields is not None:
            rp.setCurrentListFields(fields, flk, modifiers)
        return fields
    
    def quitFind(self):
        if self.findDlg is not None:
            self.selectedTab["tab"].quitFind()
            self.findDlg = None

    def OnFind(self, event):
        """Shows a custom dialog to open the three data files"""
        if self.findDlg is None:
            self.findDlg = FindDialog(self, self.selectedTab["tab"])
            self.findDlg.showDialog()
        else:
            self.findDlg.doNext()

    # def OnFindNext(self, event):
    #     """Shows a custom dialog to open the three data files"""
    #     self.selectedTab["tab"].updateFind()

    # def OnImportData(self, event):
    #     """Shows a custom dialog to open the three data files"""
    #     if self.dw.getData() is not None:
    #         if not self.checkAndProceedWithUnsavedChanges():
    #             return
    #     if len(self.dw.reds) > 0:
    #         sure_dlg = wx.MessageDialog(self.toolFrame, 'Importing new data erases old redescriptions.\nDo you want to continue?', caption="Warning!", style=wx.OK|wx.CANCEL)
    #         if sure_dlg.ShowModal() != wx.ID_OK:
    #             return
    #         sure_dlg.Destroy()

    #     dlg = ImportDataDialog(self)
    #     dlg.showDialog()

    def OnImportDataCSV(self, event):
        """Shows a custom dialog to open the two data files"""
        if not self.checkAndProceedWithUnsavedChanges():
            return 

        dlg = ImportDataCSVDialog(self)
        dlg.showDialog()
        self.changePage(self.getDefaultTabId("e"))
                            
    def OnImportPreferences(self, event):
        if not self.checkAndProceedWithUnsavedChanges(self.dw.preferences.isChanged):
            return
        dir_name = os.path.expanduser('~/')
        open_dlg = wx.FileDialog(self.toolFrame, message='Choose file', defaultDir = dir_name,
                                 style = wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            try:
                self.dw.importPreferencesFromFile(path)
            except:
                pass
        open_dlg.Destroy()
        self.resetConstraints()        
        if self.dw.needsReload():
            self.recomputeAll()
            # self.reloadReds()

        
    def OnImportRedescriptions(self, event):
        if not self.checkAndProceedWithUnsavedChanges():
            return
        reds, sortids, path  = (None, None, "")
        wcd = 'All files|*|Query files (*.queries)|*.queries|'
        dir_name = os.path.expanduser('~/')

        open_dlg = wx.FileDialog(self.toolFrame, message='Choose file', defaultDir = dir_name,
                                 style = wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            try:
                reds, sortids = self.dw.loadRedescriptionsFromFile(path)
            except:
                pass
        open_dlg.Destroy()        
        self.loadReds(reds, sortids, path=path)
        self.changePage(self.getDefaultTabId("r"))
        
    def OnExportRedescriptions(self, event):
        if self.getNbToExportReds() < 1:
            wx.MessageDialog(self.toolFrame, 'Cannot export redescriptions: no redescriptions loaded',
                             style=wx.OK|wx.ICON_EXCLAMATION, caption='Error').ShowModal()
            return
        
        if self.dw.getPackageSaveFilename() is not None:
            dir_name = os.path.dirname(self.dw.getPackageSaveFilename())
        else:
            dir_name = os.path.expanduser('~/')
        reds = self.getToExportReds()
        save_dlg = wx.FileDialog(self.toolFrame, message='Export redescriptions to:', defaultDir = dir_name, style = wx.SAVE|wx.CHANGE_DIR)
        if save_dlg.ShowModal() == wx.ID_OK:
            path = save_dlg.GetPath()
            try:
                self.dw.exportRedescriptions(path, reds)
            except:
                pass
        save_dlg.Destroy()

    def OnPrintoutPreferences(self, event):
        self.OnExportPreferences(mess="preferences")
    def OnPrintoutPreferencesTmpl(self, event):
        self.OnExportPreferences(mess="preferences template", conf_def=self.conf_defs+self.conf_defs_io)
    def OnPrintoutPreferencesDef(self, event):
        self.OnExportPreferences(mess="default preferences", inc_def=True)
        
    def OnExportPreferences(self, mess="preferences", inc_def=False, conf_def=None):
        if self.dw.getPackageSaveFilename() is not None:
            dir_name = os.path.dirname(self.dw.getPackageSaveFilename())
        else:
            dir_name = os.path.expanduser('~/')

        save_dlg = wx.FileDialog(self.toolFrame, message='Save %s template to:' % mess, defaultDir = dir_name, style = wx.SAVE|wx.CHANGE_DIR)
        if save_dlg.ShowModal() == wx.ID_OK:
            path = save_dlg.GetPath()
            try:
                self.dw.exportPreferences(path, inc_def, conf_def)
            except Exception: 
                pass
        save_dlg.Destroy()

        
    def getSelectionTab(self):
        return self.tabbed.GetSelection()
    
    def updateDataInfo(self):
        if self.hasDataLoaded():            
            # self.statusbar.SetStatusText("%s data rows" % self.dw.getData().rowsInfo(), 4)
            self.statusbar.SetStatusText("%s" % self.dw.getData(), 4)
        else:
            self.statusbar.SetStatusText("No data loaded", 4)

    
    def changePage(self, tabn):
        if tabn in self.tabs: # and not self.tabs[tabn]["hide"]:
            self.tabbed.ChangeSelection(self.tabs_keys.index(tabn))            
            self.OnPageChanged(-1)

    def OnListChanged(self, event):
        if self.findDlg is not None:
            self.findDlg.resetFind(self.selectedTab["tab"])
            
    def OnPageChanged(self, event):
        if not self.sysTLin() and type(event) == wx.BookCtrlEvent:
            self.tabbed.ChangeSelection(event.GetSelection())
            ## self.selectedTab["tab"].Show()
        
        tsel = self.tabbed.GetSelection()
        if tsel >= 0 and tsel < len(self.tabs_keys): 
            self.selectedTab = self.tabs[self.tabs_keys[tsel]]
            # print "Do updates", self.selectedTab, event
            self.doUpdates()
            # print "Done updates"
            
        self.OnListChanged(event)

    def OnSaveSuppAsVar(self, suppVect, name):
        self.dw.getData().addSuppCol(suppVect, name)
        self.reloadVars()
    def OnSaveSelAsVar(self, lids, name):
        self.dw.getData().addSelCol(lids, name)
        self.reloadVars()
        
    def OnNewVList(self, event):
        if self.matchTabType("r"):
            self.selectedTab["tab"].viewListData(viewT=self.ids_viewT[event.GetId()])
    def OnNewVItems(self, event):
        if self.matchTabType("r"):
            self.selectedTab["tab"].viewListData(lid=-1, viewT=self.ids_viewT[event.GetId()])

    def OnNewV(self, event):
        if self.matchTabType("evr"):
            self.selectedTab["tab"].viewData(viewT=self.ids_viewT[event.GetId()])

    def OnExpand(self, event):
        if self.matchTabType("r"):
            red = self.selectedTab["tab"].getSelectedItem()
            if red is not None:
                params = {"red": red.copy()}
                self.expand(params)

    def expandFromView(self, params=None):
        self.expand(params)

    def OnHigh(self, event):
        if self.matchTabType("e"):
            self.viewsm.setAllEmphasizedR([self.selectedTab["tab"].getSelectedPos()], show_info=False, no_off=True)

    def OnShowCol(self, event):
        shw = False
        if self.matchTabType("vr") and self.getDefaultTabId("e") in self.tabs:
            row = self.tabs[self.getDefaultTabId("e")]["tab"].showRidRed(self.tabs[self.getDefaultTabId("e")]["tab"].getSelectedRow(), self.selectedTab["tab"].getSelectedItem())
            shw = True

        if shw:
            self.showTab(self.getDefaultTabId("e"))

    def OnMineAll(self, event):
        self.expand()

    def OnFilterToOne(self, event):
        if self.matchTabType("r"):
            self.selectedTab["tab"].filterToOne(self.constraints.getFilterParams("redundant"))

    def OnFilterAll(self, event):
        if self.matchTabType("r"):
            self.selectedTab["tab"].filterAll(self.constraints.getFilterParams("redundant"))

    def OnFilterFolds(self, event):
        if self.matchTabType("r"):
            self.constraints.setFolds(self.dw.getData())
            self.selectedTab["tab"].processAll(self.constraints.getActions("folds"), True)
            ## self.selectedTab["tab"].filterAll(self.constraints.getFilterParams("redundant"))

    def OnProcessAll(self, event):
        if self.matchTabType("r"):
            self.selectedTab["tab"].processAll(self.constraints.getActions("final"), True)

    def OnFlipExCol(self, event):
        if self.matchTabType("e"):
            self.selectedTab["tab"].flipFocusCol(self.selectedTab["tab"].getSelectedCol())

    def OnFlipEnabled(self, event):
        if self.matchTabType("evr"):
            self.selectedTab["tab"].flipEnabled(self.selectedTab["tab"].getSelectedRow())

    def OnEnabledAll(self, event):
        if self.matchTabType("evr"):
            self.selectedTab["tab"].setAllEnabled()

    def OnDisabledAll(self, event):
        if self.matchTabType("evr"):
            self.selectedTab["tab"].setAllDisabled()

    def OnDelete(self, event):
        if self.matchTabType("r"):
            self.selectedTab["tab"].onDeleteAny()
    def OnDeleteDisAll(self, event):
        if self.matchTabType("r"):
            self.selectedTab["tab"].onDeleteDisAll()
            
    def OnNewList(self, event):
        if self.matchTabType("r"):
            self.selectedTab["tab"].onNewList()
            
    def OnCut(self, event):
        if self.matchTabType("r"):
            self.selectedTab["tab"].onCutAny()
            self.doUpdates({"menu":True})

    def OnCopy(self, event):
        if self.matchTabType("r"):
            self.selectedTab["tab"].onCopyAny()
            self.doUpdates({"menu":True})

    def OnPaste(self, event):
        if self.matchTabType("r"):
            self.selectedTab["tab"].onPasteAny()
            self.doUpdates({"menu":True})

    def OnNormalize(self, event):
        if self.matchTabType("r"):
            rid = self.selectedTab["tab"].getSelectedItemIid()
            red = self.selectedTab["tab"].getSelectedItem()
            if red is not None:
                redn, changed = red.getNormalized(self.dw.getData())
                if changed:
                    self.viewsm.dispatchEdit(redn, ikey=(self.selectedTab["tab"].tabId, rid))
                    # self.selectedTab["tab"].substituteItem(rid, redn)
                    # self.selectedTab["tab"].appendItemToHist(red)
                    self.doUpdates({"menu":True})

    def appendToHist(self, red):
        if self.getDefaultTabId("r") in self.tabs:
            self.tabs[self.getDefaultTabId("r")]["tab"].appendItemToHist(red)


    def flipRowsEnabled(self, rids):
        if self.getDefaultTabId("e") in self.tabs and len(rids)> 0:
            self.tabs[self.getDefaultTabId("e")]["tab"].flipAllEnabled(rids)

    def resetFields(self):
        for ti, tab in self.getTabsMatchType("r"):
            tab["tab"].resetFields()

    def recomputeAll(self):
        restrict = self.dw.getData().nonselectedRows()
        for ti, tab in self.getTabsMatchType("r"):
            tab["tab"].recomputeAll(restrict)
        self.viewsm.recomputeAll()

    def OnVizCheck(self, event):
        self.getVizm().setVizCheck(not event.IsChecked())
            
    def OnTabW(self, event):
        if event.GetId() in self.check_tab:
            tab_id = self.check_tab[event.GetId()]
            # if self.toolFrame.FindFocus() is not None and self.toolFrame.FindFocus().GetGrandParent() is not None \
            #        and self.toolFrame.FindFocus().GetGrandParent().GetParent() == self.toolFrame:
            #     if self.sysTLin():
            #         if not event.IsChecked():
            #             self.hideTab(tab_id)
            #             return
            # else: self.toTop()
            self.showTab(tab_id)

    def showTab(self, tab_id):
        # self.tabs[tab_id]["hide"] = False
        # self.tabs[tab_id]["tab"].Show()
        self.changePage(tab_id)

    # def hideTab(self, tab_id):
    #     self.tabs[tab_id]["hide"] = True
    #     self.tabs[tab_id]["tab"].Hide()

    def toTop(self):
        self.toolFrame.Raise()
        self.toolFrame.SetFocus()

    def OnPreferencesDialog(self, event):
        d = PreferencesDialog(self.toolFrame, self.dw)
        d.ShowModal()
        d.Destroy()
        self.resetConstraints()
        if self.dw.needsReload():
            self.recomputeAll()
            # self.reloadReds()

        
    def OnConnectionDialog(self, event):
        d = ConnectionDialog(self.toolFrame, self.dw, self.plant, self)
        tt = d.ShowModal()
        d.Destroy()

    def OnSplitDialog(self, event):
        d = SplitDialog(self.toolFrame, self.dw, self)
        d.ShowModal()
        d.Destroy()


    def OnHelp(self, event):
        self._onHelpOutside()
        # wxVer = map(int, wx.__version__.split('.'))
        # new_ver = wxVer[0] > 2 or (wxVer[0] == 2 and wxVer[1] > 9) or (wxVer[0] == 2 and wxVer[1] == 9 and wxVer[2] >= 3)
        # if new_ver:
        #     try:                                                      
        #         self._onHelpInside()
        #     except NotImplementedError:
        #         new_ver = False
        # if not new_ver:
        #    self._onHelpOutside()

    def _onHelpInside(self):
        import wx.html2
        import urllib
        import platform
        # DEBUG
        #self.toolFrame.Bind(wx.html2.EVT_WEBVIEW_ERROR, lambda evt: wx.MessageDialog(self.toolFrame, str(evt), style=wx.OK, caption='WebView Error').ShowModal())
        #self.toolFrame.Bind(wx.html2.EVT_WEBVIEW_LOADED, lambda evt: wx.MessageDialog(self.toolFrame, 'Help files loaded from '+evt.GetURL(), style=wx.OK, caption='Help files loaded!').ShowModal())
        if self.helpURL is None:
            self._onHelpOutside()
            return
        if self.helpFrame is None:
            self.helpFrame = wx.Frame(self.toolFrame, -1, self.titleHelp)
            self.helpFrame.Bind(wx.EVT_CLOSE, self._helpInsideClose)
            sizer = wx.BoxSizer(wx.VERTICAL)
            url = 'file://'+os.path.abspath(self.helpURL)
            if sys.platform == "darwin":
                # OS X returns UTF-8 encoded path names, decode to Unicode
                #url = url.decode('utf-8')
                # URLLIB doesn't like unicode strings, so keep w/ UTF-8 encoding
                # make the URL string URL-safe for OS X
                url = urllib.quote(url)
            browser = wx.html2.WebView.New(self.helpFrame, url=url)
            #browser.LoadURL('file://'+os.path.abspath(self.helpURL))
            sizer.Add(browser, 1, wx.EXPAND, 10)
            self.helpFrame.SetSizer(sizer)
            self.helpFrame.SetSize((900, 700))
        self.helpFrame.Show()
        self.helpFrame.Raise()

    def _helpInsideClose(self, event):
        self.helpFrame.Destroy()
        self.helpFrame = None

    def _onHelpOutside(self):
        import webbrowser
        try:
            ##webbrowser.open("file://"+ self.helpURL, new=1, autoraise=True)
            webbrowser.open(self.helpInternetURL, new=1, autoraise=True)
        except webbrowser.Error as e:
            self.logger.printL(1,'Cannot show help file: '+str(e)
                                   +'\nYou can find help at '+self.helpInternetURL+'\nor '+self.helpURL, "error", "help")        

    def OnAbout(self, event):
        wx.AboutBox(self.info)        

    def showCol(self, side, col):
        if self.getDefaultTabId("e") in self.tabs:
            self.tabs[self.getDefaultTabId("e")]["tab"].showCol(side, col)
            
    def showDetailsBox(self, rid, red):
        row = None
        if self.getDefaultTabId("e") in self.tabs:
            row = self.tabs[self.getDefaultTabId("e")]["tab"].showRidRed(rid, red)
        if row is not None:
            self.showTab(self.getDefaultTabId("e"))
        else:
            dlg = wx.MessageDialog(self.toolFrame,
                                   self.prepareDetails(rid, red),"Point Details", wx.OK|wx.ICON_INFORMATION)
            result = dlg.ShowModal()
            dlg.Destroy()

    def prepareDetails(self, rid, red):
        dets = "%d:\n" % rid 
        for side,pref in [(0,""), (1,"")]:
            dets += "\n"
            for lit in red.queries[side].listLiterals():
                dets += ("\t%s=\t%s\n" % (self.dw.getData().col(side,lit.colId()).getName(), self.dw.getData().getValue(side, lit.colId(), rid)))
        return dets


    def OnRunTest(self, event):
        pm = PolyMap()
        PointsMap, PointsIds = pm.getPointsFromData(self.dw.getData())
    
        final_polys, final_details, edges, nodes = pm.compute_polys(PointsMap)
        pm.plot_final_polys(final_polys, PointsMap)

        # if test_code is not None:
        #     test_code.test_run(self)    
    
    def OnSize(self, event):
        if self.getVizm() is not None:
            self.getVizm().resizeViz()
        event.Skip()
                
    def OnQuit(self, event):
        if self.plant.getWP().isActive():  
            if isinstance(self.plant.getWP(), WorkClient) and self.plant.getWP().nbWorking()>0:
                collectLater = self.returnLater(what=self.plant.getWP().hid)
                if collectLater == 0:
                    return
                self.plant.getWP().closeDown(self, collectLater > 0)   
            else :
                self.plant.getWP().closeDown(self)        
        if not self.checkAndProceedWithUnsavedChanges(what="quit"):
            return
        self.viewsm.deleteAllViews()
        if self.vizm is not None:
            self.vizm.OnQuit()
        self.toolFrame.Destroy()
        sys.exit()

    def OnLicense(self, event):
        import codecs # we want to be able to read UTF-8 license files
        license_text = None
        try:
            f = codecs.open(self.license_file, 'r', encoding='utf-8', errors='replace')
            license_text = f.read()
        except:
            wx.MessageDialog(self.toolFrame, 'No license found.', style=wx.OK, caption="No license").ShowModal()
            return

        external_license_texts = ''
        for ext in self.external_licenses:
            lic = 'LICENSE_'+ext
            try:
                lfile = Siren.searchData(lic, 'licenses')
                if lfile is not None:
                    f = codecs.open(lfile, 'r', encoding='utf-8', errors='replace')
                    external_license_texts += '\n\n***********************************\n\n' + f.read()
                    f.close()
            except:
                pass # We don't care about errors here

        if len(external_license_texts) > 0:
            license_text += "\n\nSiren comes bundled with other software for your convinience.\nThe licenses for this bundled software are below." + external_license_texts

        # Show dialog
        try:
            dlg = wx.lib.dialogs.ScrolledMessageDialog(self.toolFrame, license_text, "LICENSE")
        except Exception as e:
            wx.MessageDialog(self.toolFrame, 'Cannot show the license: '+str(e), style=wx.OK, caption="ERROR").ShowModal()
            sys.stderr.write(str(e))
        else:
            dlg.ShowModal()
            dlg.Destroy()
         
    def returnLater(self, test=None, what="continue"):
        reponse = -1
        dlg = wx.MessageDialog(self.toolFrame,'Some computations are underway (client id %s).\nDo you intend to collect the results later on?' % what, style=wx.YES_NO|wx.CANCEL|wx.NO_DEFAULT|wx.ICON_EXCLAMATION, caption='Computations in progress')
        tt = dlg.ShowModal()
        if tt == wx.ID_CANCEL: 
            reponse = 0
        if tt == wx.ID_YES: 
            reponse = 1
        return reponse

    def checkAndProceedWithUnsavedChanges(self, test=None, what="continue"):
        """Checks for unsaved changes and returns False if they exist and user doesn't want to continue
        and True if there are no unsaved changes or user wants to proceed in any case.
        If additional parameter 'test' is given, asks the question if it is true."""
        if test is None:
            test = self.hasRedsChanged({"srcTypOut": 'history'}) | self.dw.isChanged
        if test:
            dlg = wx.MessageDialog(self.toolFrame, 'Unsaved changes might be lost.\nAre you sure you want to %s?' % what, style=wx.YES_NO|wx.NO_DEFAULT|wx.ICON_EXCLAMATION, caption='Unsaved changes!')
            if dlg.ShowModal() == wx.ID_NO:
                return False
        return True

    def resetCoordinates(self):
        self.viewsm.deleteAllViews()

    def resetConstraints(self):
        self.constraints = Constraints(self.dw.getData(), self.dw.getPreferences())

    def resetLogger(self):
        verb = 1
        if self.dw.getPreferences() is not None and self.dw.getPreference('verbosity') is not None:
            verb = self.dw.getPreference('verbosity')

        self.logger.resetOut()
        if self.plant.getWP().isActive() and self.plant.getWP().getOutQueue() is not None:
            self.logger.addOut({"log": verb, "error": 0, "status":1, "time":0, "progress":2, "result":1}, self.plant.getWP().getOutQueue(), self.plant.getWP().sendMessage)
        else:
            self.logger.addOut({"*": verb,  "error":1,  "status":0, "result":0, "progress":0}, None, self.loggingLogTab)
        self.logger.addOut({"error":1, "dw_error":1}, "stderr")
        self.logger.addOut({"dw_error":1}, None, self.loggingDWError)

    def loadAll(self, refresh_rfields=False):
        if self.getVizm() is not None:
            self.getVizm().reloadVizTab()
        if self.plant is not None:
            self.plant.getWP().closeDown(self)
        self.reloadVars(review=False)
        self.reloadRows()
        self.clearReds()
        self.loadReds(refresh_fields=refresh_rfields)
        self.updateDataInfo()
        
    def reloadAll(self):
        if self.getVizm() is not None:
            self.getVizm().reloadVizTab()
        if self.plant is not None:
            self.plant.getWP().closeDown(self)
        self.reloadVars(review=False)
        self.reloadRows()
        self.reloadReds()
        self.updateDataInfo()
        
    def reloadVars(self, review=True):
        ## Initialize variable lists data
        if self.dw.getData() is not None:
            if self.getDefaultTabId("v") in self.tabs:
                self.tabs[self.getDefaultTabId("v")]["tab"].resetData(src='pack', data=self.dw.getData())
            details = {"names": self.dw.getData().getNames()}
            for ti, tab in self.getTabsMatchType("r"):
                tab["tab"].resetDetails(details, review)
                tab["tab"].resetFields()
            for ti, tab in self.getTabsMatchType("e"):
                tab["tab"].resetFields(self.dw, review)


    def reloadRows(self):
        ## Initialize variable lists data
        for ti, tab in self.getTabsMatchType("e"):
            if self.dw.getData() is not None:
                tab["tab"].resetData(self.dw.getDataRows())
            else:
                tab["tab"].resetData()

    def clearReds(self):
        for ti, tab in self.getTabsMatchType("r"):
            tab["tab"].resetData()
                        
    def reloadReds(self, all=True):
        ## Initialize red lists data
        if self.dw.getData() is not None:
            for ti, tab in self.getTabsMatchType("r"):
                tab["tab"].recomputeAll() #resetData()
        ## self.viewsm.deleteAllViews()
        self.doUpdates({"menu":True})
        self.dw.reloaded()
        
    def loadReds(self, reds=None, sortids=None, path=None, refresh_fields=False):
        ## Initialize red lists data
        tab = None
        if self.matchTabType("r"):
            tab = self.selectedTab
        elif self.getDefaultTabId("r") in self.tabs:
            tab = self.tabs[self.getDefaultTabId("r")]
        if tab is not None:
            if refresh_fields:
                tab["tab"].resetFields()

            if reds is None:
                for ri, rs in enumerate(self.dw.getReds()):
                    tab["tab"].addData(rs)
            else:
                tab["tab"].addData(src=('file', path, 0), data=reds)
                # tab["tab"].addData(src='pack', data=self.dw.getReds(), sord=self.dw.getShowIds())
                self.doUpdates({"menu":True})
        self.dw.reloaded()

    def hasRedsChanged(self, patts):
        return any([tab["tab"].hasRedsChanged(patts) for ti, tab in self.getTabsMatchType("r")])            

    def startFileActionMsg(self, msg, short_msg=''):
        """Shows a dialog that we're reading a file"""
        self.statusbar.SetStatusText(short_msg, 1)
        self.toolFrame.Enable(False)
        self.busyDlg = wx.BusyInfo(msg, self.toolFrame)
        #self.busyDlg = CBusyDialog.showBox(self.toolFrame, msg, short_msg, None)
        #self.busyDlg = PBI.PyBusyInfo(msg, parent=self.toolFrame, title=short_msg)
        # DEBUG
        # time.sleep(5)
        

    def stopFileActionMsg(self, msg=''):
        """Removes the BusyInfo dialog"""
        if self.busyDlg is not None:
            self.busyDlg.Destroy()
            # del self.busyDlg # Removes the dialog
            self.busyDlg = None
        self.toolFrame.Enable(True)
        self.statusbar.SetStatusText(msg, 1)
        
    def errorBox(self, message):
        if self.busyDlg is not None:
            del self.busyDlg
            self.busyDlg = None
        dlg = wx.MessageDialog(self.toolFrame, message, style=wx.OK|wx.ICON_EXCLAMATION|wx.STAY_ON_TOP, caption="Error")
        dlg.ShowModal()
        dlg.Destroy()

    # def warningBox(self, message):
    #     # if self.busyDlg is not None:
    #     #     del self.busyDlg
    #     #     self.busyDlg = None
    #     dlg = wx.MessageDialog(self.toolFrame, message, style=wx.OK|wx.ICON_INFORMATION|wx.STAY_ON_TOP, caption="Warning")
    #     dlg.ShowModal()
    #     dlg.Destroy()

    ##### receiving results
    def readyReds(self, wid, reds, tab):
        if len(reds) > 0 and self.getDefaultTabId("r") in self.tabs:
            for red in reds:
                red.recompute(self.getData())
                red.setRestrictedSupp(self.getData())
            src = ('run', wid)
            self.tabs[self.getDefaultTabId("r")]["tab"].insertItems(src, reds)

    def readyProj(self, wid, vid, proj):
        adjunct_vps = self.proj_cache.incomingPC(proj, vid)

        vv = self.viewsm.accessViewX(vid)
        if vv is not None:
            vv.readyProj(proj)
        for (avid, aproj) in adjunct_vps:
            vv = self.viewsm.accessViewX(avid)
            if vv is not None:
                vv.readyProj(aproj)
                
    # def readyProj(self, vid, proj):        
    #     adjunct_vps = self.proj_cache.incomingPC(proj, vid)
    #     # print "Ready proj", vid, proj, adjunct_vps
    #     # print "View IDS", self.view_ids.keys()
    #     if vid in self.view_ids:
    #         self.view_ids[vid].readyProj(proj)
    #     for (avid, aproj) in adjunct_vps:
    #         if avid in self.view_ids:
    #             self.view_ids[avid].readyProj(aproj)

    def getERCache(self):
        return self.er_cache
