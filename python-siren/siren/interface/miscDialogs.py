import sys
import wx
### from wx import ALIGN_RIGHT, ALL, CANCEL, CHANGE_DIR, EXPAND, HORIZONTAL, VERTICAL, OK, OPEN, TE_READONLY
### from wx import BoxSizer, Button, Choice, Dialog, FileDialog, FlexGridSizer, GridSizer, NewId, TextCtrl, StaticText
### from wx import EVT_BUTTON, EVT_CLOSE, EVT_KEY_UP
### from wx import ID_ANY, ID_APPLY, ID_FIND, ID_OK

import os.path, re
from ..reremi.classData import DataError, NA_str_def

import pdb

# class ImportDataDialog(object):
#     """Helper class to show the dialog for importing data file triplets"""
#     def __init__(self, parent):
#         self.parent = parent
#         self.dlg = wx.Dialog(self.parent.toolFrame, title="Import data")

#         LHStext = wx.StaticText(self.dlg, label='Left-hand side variables file:')
#         RHStext = wx.StaticText(self.dlg, label='Right-hand side variables file:')
#         Cootext = wx.StaticText(self.dlg, label='Coordinates file:')
#         RNamestext = wx.StaticText(self.dlg, label='Entities file:')

#         self.LHSfile = None
#         self.RHSfile = None
#         self.Coofile = None
#         self.RNamesfile = None

#         self.LHSfileTxt = wx.TextCtrl(self.dlg, value='', size=(500,10), style=wx.TE_READONLY)
#         self.RHSfileTxt = wx.TextCtrl(self.dlg, value='', style=wx.TE_READONLY)
#         self.CoofileTxt = wx.TextCtrl(self.dlg, value='', style=wx.TE_READONLY)
#         self.RNamesfileTxt = wx.TextCtrl(self.dlg, value='', style=wx.TE_READONLY)

#         LHSbtn = wx.Button(self.dlg, label='Choose', name='LHS')
#         RHSbtn = wx.Button(self.dlg, label='Choose', name='RHS')
#         Coobtn = wx.Button(self.dlg, label='Choose', name='Coordinates')
#         RNamesbtn = wx.Button(self.dlg, label='Choose', name='Entities')

#         LHSbtn.Bind(wx.EVT_BUTTON, self.onButton)
#         RHSbtn.Bind(wx.EVT_BUTTON, self.onButton)
#         Coobtn.Bind(wx.EVT_BUTTON, self.onButton)
#         RNamesbtn.Bind(wx.EVT_BUTTON, self.onButton)

#         gridSizer = wx.FlexGridSizer(rows = 4, cols = 3, hgap = 5, vgap = 5)
#         gridSizer.AddGrowableCol(1, proportion=1)
#         gridSizer.SetFlexibleDirection(wx.HORIZONTAL)
#         gridSizer.AddMany([(LHStext, 0, wx.ALIGN_RIGHT), (self.LHSfileTxt, 1, wx.EXPAND), (LHSbtn, 0),
#                            (RHStext, 0, wx.ALIGN_RIGHT), (self.RHSfileTxt, 1, wx.EXPAND), (RHSbtn, 0),
#                            (Cootext, 0, wx.ALIGN_RIGHT), (self.CoofileTxt, 1, wx.EXPAND), (Coobtn, 0),
#                            (RNamestext, 0, wx.ALIGN_RIGHT), (self.RNamesfileTxt, 1, wx.EXPAND), (RNamesbtn, 0)])

#         btnSizer = self.dlg.CreateButtonSizer(wx.OK|wx.CANCEL)
#         topSizer = wx.BoxSizer(wx.VERTICAL)
#         topSizer.Add(gridSizer, flag=wx.ALL, border=5)
#         topSizer.Add(btnSizer, flag=wx.ALL, border=5)

#         self.dlg.SetSizer(topSizer)
#         self.dlg.Fit()

#         self.open_dir = os.path.expanduser('~/')
#         self.wcd = 'All files|*|Numerical Variables (*.densenum / *.datnum)|*.densenum/*.datnum|Boolean Variables (*.sparsebool / *.datbool)|*.sparsebool/*.datbool'
#         self.names_wcd = 'All files|*|Information files (*.names)|*.names'


#     def showDialog(self):
#         if self.dlg.ShowModal() == wx.ID_OK:
#             try:
#                 if self.RHSfile is None:
#                     self.parent.dw.importDataFromMulFile(self.LHSfile)
#                 else:
#                     self.parent.dw.importDataFromMulFiles([self.LHSfile, self.RHSfile], None, self.Coofile, self.RNamesfile)
#             except:
#                 pass
#                 ##raise
#             else:
#                 self.parent.reloadAll()
#             finally:
#                 self.dlg.Destroy()
#             return True
#         else:
#             return False
                
            

#     def onButton(self, e):
#         button = e.GetEventObject()
#         btnName = button.GetName()
#         if btnName == 'Coordinates' or  btnName == 'Entities':
#             wcd = self.names_wcd
#         else:
#             wcd = self.wcd
#         open_dlg = wx.FileDialog(self.parent.toolFrame, message="Choose "+btnName+" file",
#                                  defaultDir=self.open_dir, wildcard=wcd,
#                                  style=wx.OPEN|wx.CHANGE_DIR)
#         if open_dlg.ShowModal() == wx.ID_OK:
#             path = open_dlg.GetPath()
#             self.open_dir = os.path.dirname(path)
#             if btnName == 'LHS':
#                 self.LHSfileTxt.ChangeValue(path)
#                 self.LHSfile = path
#                 # Both TextCtrl and variable hold the same info, but if the latter is empty is None,
#                 # making it compatible with dw.importDataFromMulFiles
#             elif btnName == 'RHS':
#                 self.RHSfileTxt.ChangeValue(path)
#                 self.RHSfile = path
#             elif btnName == 'Coordinates':
#                 self.CoofileTxt.ChangeValue(path)
#                 self.Coofile = path
#             elif btnName == 'Entities':
#                 self.RNamesfileTxt.ChangeValue(path)
#                 self.RNamesfile = path



class ImportDataCSVDialog(object):
    """Helper class to show the dialog for importing data file csv pairs"""
    def __init__(self, parent):
        self.parent = parent
        self.dlg = wx.Dialog(self.parent.toolFrame, title="Import data")

        LHStext = wx.StaticText(self.dlg, label='Left-hand side variables file:')
        RHStext = wx.StaticText(self.dlg, label='Right-hand side variables file:')

        self.dialect_options = {'delimiter': {"label": "Delimiter", "opts": [(None, "(auto)"), ('\t', 'TAB'), (';', ';'), (',', ','), (' ', 'SPC')]}}
        self.LHSfile = None
        self.RHSfile = None

        self.LHSfileTxt = wx.TextCtrl(self.dlg, value='', size=(500,10), style=wx.TE_READONLY)
        self.RHSfileTxt = wx.TextCtrl(self.dlg, value='', style=wx.TE_READONLY)

        so_sizer = wx.GridSizer(rows=1, cols=2*(1+len(self.dialect_options)), hgap=5, vgap=5)

        ctrl_id = wx.NewId()
        label = wx.StaticText(self.dlg, wx.ID_ANY, "Missing:")
        self.missing_ctrl = wx.TextCtrl(self.dlg, ctrl_id, NA_str_def)
        so_sizer.Add(label, 0, wx.ALIGN_RIGHT)
        so_sizer.Add(self.missing_ctrl, 0)

        self.dialect_ctrl = {}
        for item, details in self.dialect_options.items():
            ctrl_id = wx.NewId()
            label = wx.StaticText(self.dlg, wx.ID_ANY, details['label']+":")
            self.dialect_ctrl[item] = wx.Choice(self.dlg, ctrl_id)
            self.dialect_ctrl[item].AppendItems(strings=[v[1] for v in details['opts']])
            self.dialect_ctrl[item].SetSelection(0)
            so_sizer.Add(label, 0, wx.ALIGN_RIGHT)
            so_sizer.Add(self.dialect_ctrl[item], 0)
                        

        LHSbtn = wx.Button(self.dlg, label='Choose', name='LHS')
        RHSbtn = wx.Button(self.dlg, label='Choose', name='RHS')

        LHSbtn.Bind(wx.EVT_BUTTON, self.onButton)
        RHSbtn.Bind(wx.EVT_BUTTON, self.onButton)

        gridSizer = wx.FlexGridSizer(rows = 2, cols = 3, hgap = 5, vgap = 5)
        gridSizer.AddGrowableCol(1, proportion=1)
        gridSizer.SetFlexibleDirection(wx.HORIZONTAL)
        gridSizer.AddMany([(LHStext, 0, wx.ALIGN_RIGHT), (self.LHSfileTxt, 1, wx.EXPAND), (LHSbtn, 0),
                           (RHStext, 0, wx.ALIGN_RIGHT), (self.RHSfileTxt, 1, wx.EXPAND), (RHSbtn, 0)])

        btnSizer = self.dlg.CreateButtonSizer(wx.OK|wx.CANCEL)
        topSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer.Add(gridSizer, flag=wx.ALL, border=5)
        topSizer.Add(so_sizer, flag=wx.EXPAND|wx.ALL, border=5)
        topSizer.Add(btnSizer, flag=wx.ALL, border=5)

        self.dlg.SetSizer(topSizer)
        self.dlg.Fit()

        self.open_dir = os.path.expanduser('~/')
        self.wcd = 'All files|*|CSV files|*.csv'


    def showDialog(self):
        na = None
        dialect_dict = {}
        if self.dlg.ShowModal() == wx.ID_OK:
            tmp = self.missing_ctrl.GetValue()
            na = tmp
                
            for item, ctrl_single in self.dialect_ctrl.items():
                tmp = self.dialect_options[item]['opts'][ctrl_single.GetCurrentSelection()][0]
                if tmp is not None:
                    dialect_dict[item] = tmp
            try:
                self.parent.dw.importDataFromCSVFiles([self.LHSfile, self.RHSfile, dialect_dict, na])
            except:
                pass
                raise
            else:
                self.parent.loadAll()
            finally:
                self.dlg.Destroy()
            return True
        else:
            return False
                
    def onButton(self, e):
        button = e.GetEventObject()
        btnName = button.GetName()
        wcd = self.wcd
        open_dlg = wx.FileDialog(self.parent.toolFrame, message="Choose "+btnName+" file",
                                 defaultDir=self.open_dir, wildcard=wcd,
                                 style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.open_dir = os.path.dirname(path)
            if btnName == 'LHS':
                self.LHSfileTxt.ChangeValue(path)
                self.LHSfile = path
            elif btnName == 'RHS':
                self.RHSfileTxt.ChangeValue(path)
                self.RHSfile = path

class ExportFigsDialog(object):
    """Helper class to show the dialog for importing data file csv pairs"""
    def __init__(self, parent, vm, items, ddir=None):
        self.parent = parent
        self.vm = vm
        self.items = items
        self.dlg = wx.Dialog(self.parent.toolFrame, title="Export figures")

        Extext = wx.StaticText(self.dlg, label='Export file pattern:')
        self.exfile = None
        self.exfileTxt = wx.TextCtrl(self.dlg, value='', size=(500,10), style=wx.TE_READONLY)

        self.format_options = {'format': {"label": "Format", "order":1, "opts": [(None, ""), ('png', 'png'), ('eps', 'eps'), ('pdf', 'pdf')]},
                               'stamp': {"label": "Stamp", "order": 2, "opts": [(True, 'Yes'), (False, 'No')]},
                               'with_disabled': {"label": "Disabled", "order":3, "opts": [(False, 'Exclude'), (True, 'Include')]},
                               'viewT': {"label": "View type", "order":4, "opts": []}}

        viewTdef = vm.getDefaultViewT(typv="R")
        for v in vm.getViewsItems(typv="R"):
            self.format_options["viewT"]["opts"].append((v["viewT"], v["short_title"]))
            if v["viewT"] == viewTdef:
                self.format_options["viewT"]["opts"].insert(0, (v["viewT"], v["short_title"]))
        
        so_sizer = wx.FlexGridSizer(rows=2, cols=(2+len(self.format_options)), hgap=1, vgap=1)

        ctrl_id = wx.NewId()
        label = wx.StaticText(self.dlg, wx.ID_ANY, "Height:")
        self.height_ctrl = wx.TextCtrl(self.dlg, ctrl_id, "600")
        so_sizer.Add(label, 0, wx.ALIGN_RIGHT)
        so_sizer.Add(self.height_ctrl, 0)
        label = wx.StaticText(self.dlg, wx.ID_ANY, "Width:")
        self.width_ctrl = wx.TextCtrl(self.dlg, ctrl_id, "800")
        so_sizer.Add(label, 0, wx.ALIGN_RIGHT)
        so_sizer.Add(self.width_ctrl, 0)

        self.format_ctrl = {}
        iks = sorted(self.format_options.keys(), key= lambda x: self.format_options[x].get("order", 1))
        for item in iks:
            details = self.format_options[item]
            ctrl_id = wx.NewId()
            label = wx.StaticText(self.dlg, wx.ID_ANY, details['label']+":")
            self.format_ctrl[item] = wx.Choice(self.dlg, ctrl_id)
            self.format_ctrl[item].AppendItems(strings=[v[1] for v in details['opts']])
            self.format_ctrl[item].SetSelection(0)
            so_sizer.Add(label, 0, wx.ALIGN_RIGHT)
            so_sizer.Add(self.format_ctrl[item], 0)
                        

        Filebtn = wx.Button(self.dlg, label='Choose', name='ExFile')
        Filebtn.Bind(wx.EVT_BUTTON, self.onButton)

        gridSizer = wx.FlexGridSizer(rows = 1, cols = 3, hgap = 5, vgap = 5)
        gridSizer.AddGrowableCol(1, proportion=1)
        gridSizer.SetFlexibleDirection(wx.HORIZONTAL)
        gridSizer.AddMany([(Extext, 0, wx.ALIGN_RIGHT), (self.exfileTxt, 1, wx.EXPAND), (Filebtn, 0)])

        btnSizer = self.dlg.CreateButtonSizer(wx.OK|wx.CANCEL)
        topSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer.Add(gridSizer, flag=wx.ALL, border=5)
        topSizer.Add(so_sizer, flag=wx.EXPAND|wx.ALL, border=5)
        topSizer.Add(btnSizer, flag=wx.ALL, border=5)

        self.dlg.SetSizer(topSizer)
        self.dlg.Fit()

        if ddir is None:
            self.open_dir = os.path.expanduser('~/')
        else:
            self.open_dir = ddir
        self.wcd = 'All files|*'


    def showDialog(self):
        format_dict = {}
        if self.dlg.ShowModal() == wx.ID_OK:
            try:
                height = int(self.height_ctrl.GetValue())
            except:
                height = -1
            try:
                width = int(self.width_ctrl.GetValue())
            except:
                width = -1
                
            for item, ctrl_single in self.format_ctrl.items():
                tmp = self.format_options[item]['opts'][ctrl_single.GetCurrentSelection()][0]
                if tmp is not None:
                    format_dict[item] = tmp
            try:
                self.parent.dw.exportItemsFigs(self.vm, self.exfile, self.items, (width, height), format_dict)
            except:
                pass
                raise
            finally:
                self.dlg.Destroy()
            return True
        else:
            return False
                
    def onButton(self, e):
        button = e.GetEventObject()
        btnName = button.GetName()
        wcd = self.wcd
        open_dlg = wx.FileDialog(self.parent.toolFrame, message="Choose "+btnName+" file",
                                 defaultDir=self.open_dir, wildcard=wcd,
                                 style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.open_dir = os.path.dirname(path)
            if btnName == 'ExFile':
                self.exfileTxt.ChangeValue(path)
                self.exfile = path


class FindDialog(object):
    """Helper class to show the dialog for finding items"""
    def __init__(self, parent, page):        
        self.parent = parent
        self.dlg = wx.Dialog(self.parent.toolFrame, title="Find")
        self.resetFind(page)
        
        self.findTxt = wx.TextCtrl(self.dlg, value='', size=(500,30))
        self.findTxt.Bind(wx.EVT_KEY_UP, self.OnKey)
        self.findTxt.SetFocus()

        nextBtn = wx.Button(self.dlg, id=wx.ID_FIND, name='next')
        nextBtn.Bind(wx.EVT_BUTTON, self.onButton)
        self.dlg.Bind(wx.EVT_CLOSE, self.OnQuit)
        # btnSizer = self.dlg.CreateButtonSizer(wx.OK)
        # btnSizer = self.dlg.CreateStdDialogButtonSizer(WX.OK)
        # btnSizer.AddButton(wx.Button(self.dlg, id=wx.ID_APPLY, label="Next"))
        # #btnSizer.AddButton(wx.Button(self.dlg, id=wx.ID_OK, label="OK"))
        # btnSizer.Realize()

        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(self.findTxt, flag=wx.EXPAND|wx.ALL)
        topSizer.Add(nextBtn, 0)
        #topSizer.Add(btnSizer, flag=wx.ALL, border=5)

        self.dlg.SetSizer(topSizer)
        self.dlg.Fit()

    def resetValues(self, list_values):
        self.list_values = list_values
    def resetCallback(self, callback):
        self.callback = callback
    def resetFind(self, page):
        self.resetValues(page.getNamesList())
        self.resetCallback(page.updateFind)
        
    def showDialog(self):
        self.dlg.Show()
        # if self.dlg.ShowModal() == wx.ID_OK:
        #     #self.doFind(self.findTxt.GetValue())
        #     

    def onButton(self, e):
        button = e.GetEventObject()
        if button.GetId() == wx.ID_FIND:
            self.doNext()
            #self.parent.toolFrame.SetFocus()

    def OnKey(self, event):
        if len(self.findTxt.GetValue()) > 0:
            self.doFind(self.findTxt.GetValue())

    def doFind(self, patt):
        matching = []
        non_matching = []
        
        try:
            re.search(patt, "", re.IGNORECASE)
        except re.error:
            return
        
        for (x, value) in self.list_values:
            if re.search(patt, value, re.IGNORECASE) is not None:
                matching.append(x)
            else:
                non_matching.append(x)
        if self.callback is not None:
            self.callback(matching, non_matching, -1)
            ## print self.callback, matching, non_matching

    def doNext(self):
        if self.callback is not None:
            self.callback(cid=None)

    def OnQuit(self, event):
        self.parent.quitFind()
        self.dlg.Destroy()



########################################################################
class ChoiceElement:
    """"""

    #----------------------------------------------------------------------
    def __init__(self, id, lbl):
        """Constructor"""
        self.id = id
        self.lbl = lbl
    def getLabel(self):
        return self.lbl
    def getId(self):
        return self.id

    
########################################################################
class MultiSelectorDialog(object):

    #----------------------------------------------------------------------
    def __init__(self, parent, choice_list, selected_ids=[]):

        self.parent = parent
        self.dlg = wx.Dialog(self.parent.toolFrame, title="Fields selection", style=wx.RESIZE_BORDER|wx.DEFAULT_DIALOG_STYLE)
        self.default_choice_list = choice_list
        self.default_selected_ids = selected_ids
        
        self.lists = {}
        self.options = {}
        
        for ll in ["LHS", "RHS"]:
            self.lists[ll] = wx.ListBox(self.dlg, size=(200, 150), choices=[], style=wx.LB_MULTIPLE)
        self.buttons = {"up": wx.Button(self.dlg,-1, "^", pos=(110, 180)),
                        "down": wx.Button(self.dlg,-1, "v", pos=(110, 180)),
                        "add": wx.Button(self.dlg,-1, ">", pos=(110, 180)),
                        "rm": wx.Button(self.dlg,-1, "<", pos=(110, 180)),
                        "reset": wx.Button(self.dlg,-1, "Reset", pos=(110, 180))}

        self.buttons["up"].Bind(wx.EVT_BUTTON, self.onBtnUp)
        self.buttons["down"].Bind(wx.EVT_BUTTON, self.onBtnDown)
        self.buttons["add"].Bind(wx.EVT_BUTTON, self.onBtnAdd)
        self.buttons["rm"].Bind(wx.EVT_BUTTON, self.onBtnRm)
        self.buttons["reset"].Bind(wx.EVT_BUTTON, self.onBtnReset)        

        self.lists["LHS"].Bind(wx.EVT_LISTBOX_DCLICK, self.onLstAdd)
        self.lists["RHS"].Bind(wx.EVT_LISTBOX_DCLICK, self.onLstRm)

        self.populate(choice_list, selected_ids)

        sizer_btns = wx.BoxSizer(wx.VERTICAL)
        sizer_btns.AddStretchSpacer()
        for btn in ["up", "add", "rm", "down", "reset"]:
            sizer_btns.Add(self.buttons[btn], 0, wx.ALL, 5)
        sizer_btns.AddStretchSpacer()
            
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.lists["LHS"], 2, wx.EXPAND|wx.ALL, 5)
        sizer.Add(sizer_btns, 0, wx.EXPAND|wx.ALL, 5)
        sizer.Add(self.lists["RHS"], 2, wx.EXPAND|wx.ALL, 5)

        btnSizer = self.dlg.CreateButtonSizer(wx.OK|wx.CANCEL)
        topSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer.Add(sizer, 2, wx.EXPAND|wx.ALL, 5)
        topSizer.Add(btnSizer, 0, wx.EXPAND|wx.ALL, 5)
        
        self.dlg.SetSizer(topSizer)
        self.dlg.Fit()


    def populate(self, choice_list, selected_ids=[]):
        self.options = dict(choice_list)
        lhs_ids = [cid for cid, name in choice_list if cid not in selected_ids]
        self.setElemsList("LHS", lhs_ids)
        self.setElemsList("RHS", selected_ids)

    def resetLists(self):
        self.populate(self.default_choice_list, self.default_selected_ids)
        
    def getOrdSelect(self):
        which = "RHS"
        return [self.lists[which].GetClientData(p)["id"] for p in range(self.lists[which].GetCount())]        
        
    def setElemsList(self, which, eids=[]):
        self.lists[which].Clear()
        self.appendElemsList(which, eids)
    def appendElemsList(self, which, eids=[]):
        for eid in eids:
            self.lists[which].Append(self.options[eid].getLabel(), {"id": eid})
    def removeElemsList(self, which, poss=[]):
        for r in sorted(poss, reverse=True):
            self.lists[which].Delete(r)

    def onBtnReset(self, event):
        self.resetLists()
    def onBtnUp(self, event):
        self.doUp()
    def onBtnDown(self, event):
        self.doDown()
    def onBtnAdd(self, event):
        self.doAdd()
    def onBtnRm(self, event):
        self.doRm()
    def onLstAdd(self, event):
        self.doAdd()
    def onLstRm(self, event):
        self.doRm()
    def doAdd(self):
        self.doSwap("LHS", "RHS")
    def doRm(self):
        self.doSwap("RHS", "LHS")

    def insert(self, which, insert_pos, selids, select=False):
        reselect = []
        for ii, eid in enumerate(selids):
            self.lists[which].Insert(self.options[eid].getLabel(), insert_pos, {"id": eid})
            reselect.append(insert_pos+ii)
        if select:
            for ii in reselect:
                self.lists[which].SetSelection(ii)
        
    def doUp(self):
        which = "RHS"
        selposs = self.lists[which].GetSelections()        
        if len(selposs) > 0:
            dd = sorted(selposs, reverse=True)
            if dd[-1] > 0:
                insert_pos = dd[-1]-1
                selids = [self.lists[which].GetClientData(p)["id"] for p in dd]
                self.removeElemsList(which, dd)
                self.insert(which, insert_pos, selids, select=True)
                    
    def doDown(self):
        which = "RHS"
        selposs = self.lists[which].GetSelections()
        if len(selposs) > 0:
            dd = sorted(selposs, reverse=True)
            if dd[-1] + len(dd) < self.lists[which].GetCount():
                insert_pos = dd[-1]+1
                selids = [self.lists[which].GetClientData(p)["id"] for p in dd]
                self.removeElemsList(which, dd)
                self.insert(which, insert_pos, selids, select=True)
        
    def doSwap(self, LFrom, LTo):        
        selposs = self.lists[LFrom].GetSelections()
        if len(selposs) > 0:
            selids = [self.lists[LFrom].GetClientData(p)["id"] for p in selposs]
            self.appendElemsList(LTo, selids)
            self.removeElemsList(LFrom, selposs)

    def showDialog(self):
        fields = None
        if self.dlg.ShowModal() == wx.ID_OK:
            fields = self.getOrdSelect()
        self.dlg.Destroy()
        return fields
