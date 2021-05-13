import wx
import re
### from wx import ALIGN_BOTTOM, ALIGN_CENTER_HORIZONTAL, ALL, EXPAND, HORIZONTAL, VERTICAL
### from wx import EVT_BUTTON, EVT_CHECKBOX, EVT_CHOICE, EVT_CLOSE, EVT_TEXT, ID_ANY
### from wx import BoxSizer, Button, NewId, Panel, StaticText
import pdb
from classPreferencesDialog import PreferencesDialog, ApplyResetCancelDialog

### USAGE this class provides a wx Modal dialog to modify a dictionary of preferences managed with the PreferenceManager
## It is launched with the following command:
##    # def OnConnectionDialog(self, event):
##    #     d = ConnectionDialog(main_frame, pref_handle)
##    #     d.ShowModal()
##    #     d.Destroy()


class ConnectionDialog(PreferencesDialog):
    """
    Creates a preferences dialog to setup a worker connection
    """
    SUCCESS_FC = "DARKGREEN"
    WARNING_FC = "ORANGE"
    FAIL_FC = "RED"
    IP_NUM_MATCH = "[0-9][0-9]?[0-9]?\.[0-9][0-9]?[0-9]?\.[0-9][0-9]?[0-9]?\.[0-9][0-9]?[0-9]?"
    IP_LOC_MATCH = "[Ll]ocalhost"
    IPKN_FIELDS = ["workserver_ip", "workserver_port", "workserver_authkey", "workserver_client"]
    BTN_LABELS = {"start": "Start", "stop": "Stop", "ping": "Ping", "reconnect": "Relink"}

    button_types = [{"name":"action", "label":BTN_LABELS["ping"], "funct":"self.onAction"},
                {"name":"reset", "label":"Reset", "funct": "self.onReset"},
                {"name":"rtod", "label":"ResetToDefault", "funct": "self.onResetToDefault"},
                {"name":"quit", "label":"Done", "funct": "self.onClose"}]

    def __init__(self, parent, pref_handle, wp_handle, boss):
        """
        Initialize the config dialog
        """
        wx.Dialog.__init__(self, parent, wx.ID_ANY, 'Worker setup') #, size=(550, 300))

        self.parent = parent
        self.pref_handle = pref_handle
        self.boss = boss
        self.info_box = None
        self.controls_map = {}
        self.objects_map = {}
        self.tabs = []
        self.wp_handle = wp_handle
        self.sec_id = None
        self.no_problem = True
        self.pinged = None
        self.status = 0
        self.cancel_change = False # Tracks if we should cancel a page change

        section_name = "Network"
        ti, section = self.pref_handle.getPreferencesManager().getSectionByName(section_name)
        if ti is not None:
            sec_id = wx.NewId()
            self.tabs.append(sec_id)
            self.controls_map[sec_id] = {"button": {}, "range": {},
                             "open": {}, "boolean": {}, "single_options": {},
                             "multiple_options": {}, "color_pick": {}}
            
            conf = self
            # conf = wx.Panel(self.nb, -1)
            top_sizer = wx.BoxSizer(wx.VERTICAL)
            self.dispGUI(section, sec_id, conf, top_sizer)
            self.dispInfo(conf, top_sizer)
            self.makeButtons(sec_id, conf, top_sizer)
            conf.SetSizer(top_sizer)
            top_sizer.Fit(conf)
            self.setSecValuesFromDict(sec_id, self.pref_handle.getPreferences())

            status, msg, id_clients = self.wp_handle.getWP().getDetailedInfos()
            color = ConnectionDialog.FAIL_FC
            if status == "OK":
                color = ConnectionDialog.SUCCESS_FC
            type_wp = self.wp_handle.getWP().infoStr()
            self.updateInfo(type_wp+" --- "+msg, color)
            self.controls_map[sec_id]["button"]["reset"].Disable()
            
            for txtctrl in self.controls_map[sec_id]["open"].itervalues():
                self.Bind(wx.EVT_TEXT, self.changeHappened, txtctrl)
            for txtctrl in self.controls_map[sec_id]["range"].itervalues():
                self.Bind(wx.EVT_TEXT, self.changeHappened, txtctrl)
            for choix in self.controls_map[sec_id]["boolean"].itervalues():
                self.Bind(wx.EVT_CHOICE, self.changeHappened, choix)
            for choix in self.controls_map[sec_id]["single_options"].itervalues():
                self.Bind(wx.EVT_CHOICE, self.changeHappened, choix)
            for chkset in self.controls_map[sec_id]["multiple_options"].itervalues():
                for chkbox in chkset.itervalues():
                    self.Bind(wx.EVT_CHECKBOX, self.changeHappened, chkbox)
            self.sec_id = sec_id
            ipkn = self.getIPKN_WP()
            self.resetIPKN(ipkn)
            if type_wp == "Local":
                self.setactPing()
            else:
                self.setactStop()
            self.updateAction(ipkn)

        self.Centre()
        self.Layout()
        self.Fit()
        self.Bind(wx.EVT_CLOSE, self.onClose)

    def dispInfo(self, frame, top_sizer):

        sec_sizer= wx.BoxSizer(wx.VERTICAL)

        ########## ADD INFO BOX
        text_sizer = wx.BoxSizer(wx.HORIZONTAL)

        ctrl_id = wx.NewId()
        self.info_box = wx.StaticText(frame, wx.NewId(), "")
        self.box_color = self.info_box.GetForegroundColour()
        text_sizer.Add(self.info_box, 0, wx.EXPAND|wx.ALL, 5)
        sec_sizer.Add(text_sizer, 0, wx.EXPAND|wx.ALL, 5)
        top_sizer.Add(sec_sizer, 0,  wx.EXPAND|wx.ALL, 5)


    def makeButtons(self, sec_id, frame, top_sizer):
        btn_sizer = wx.BoxSizer(wx.HORIZONTAL)

        for button in self.button_types:
            btnId = wx.NewId()
            btn = wx.Button(frame, btnId, button["label"])
            frame.Bind(wx.EVT_BUTTON, eval(button["funct"]), btn)
            btn_sizer.Add(btn, 0)
            self.controls_map[sec_id]["button"][button["name"]] = btn
            self.objects_map[btnId] = (sec_id, "button", button["name"])

        top_sizer.Add(btn_sizer, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5)

    def updateInfo(self, text, color=None):
        if color is None:
            color = self.box_color
        self.info_box.SetForegroundColour(color)
        self.info_box.SetLabel(text)              

    def changeHappened(self, event):
        ipkn_txtctrl = self.getIPKN_TxtCtrl()
        if ipkn_txtctrl != self.getIPKN_params():
            self.controls_map[self.sec_id]["button"]["rtod"].Enable()
        else:
            self.controls_map[self.sec_id]["button"]["rtod"].Disable()
        if ipkn_txtctrl != self.getIPKN_WP():
            self.controls_map[self.sec_id]["button"]["reset"].Enable()
        else:
            self.controls_map[self.sec_id]["button"]["reset"].Disable()
        self.updateAction(ipkn_txtctrl)

    def updateAction(self, ipkn_txtctrl):
        if re.match(self.IP_LOC_MATCH, ipkn_txtctrl[0]) or re.match(self.IP_NUM_MATCH, ipkn_txtctrl[0]):
            self.controls_map[self.sec_id]["button"]["action"].Enable()
        else:
            self.controls_map[self.sec_id]["button"]["action"].Disable()
        action = self.controls_map[self.sec_id]["button"]["action"].GetLabel()
        if action in [self.BTN_LABELS["start"], self.BTN_LABELS["reconnect"]] and self.pinged is not None:
            match = True
            for ii in [0,1,2]:
                match &= (ipkn_txtctrl[ii] == self.pinged[ii])
            if match:
                if ipkn_txtctrl[3] in self.pinged[3]:
                    self.controls_map[self.sec_id]["button"]["action"].SetLabel(self.BTN_LABELS["reconnect"])
                else:
                    self.controls_map[self.sec_id]["button"]["action"].SetLabel(self.BTN_LABELS["start"])
            else:
                self.controls_map[self.sec_id]["button"]["action"].SetLabel(self.BTN_LABELS["ping"])
        
    def getIPKN_WP(self):
        return self.wp_handle.getWP().getParameters()

    def getIPKN_TxtCtrl(self):
        vdict = self.getSecValuesDict(self.sec_id)
        return tuple([vdict[f]["value"] for f in self.IPKN_FIELDS])
    
    def getIPKN_params(self):
        return tuple([self.pref_handle.getPreference(f) for f in self.IPKN_FIELDS])
    
    def setIPK_params(self):
        vdict = self.getSecValuesDict(self.sec_id)
        self.pref_handle.updatePreferencesDict(vdict)

    def resetIPKN(self, ipkn):
        for i,f in enumerate(self.IPKN_FIELDS):
            if i < len(ipkn) and ipkn[i] is not None:
                self.controls_map[self.sec_id]["open"][f].SetValue(str(ipkn[i]))
        
    def onReset(self, event):       
        self.resetIPKN(self.getIPKN_WP())
        self.controls_map[self.sec_id]["button"]["reset"].Disable()

    def onClose(self, event=None):
        if self.pinged is not None:
            self.pinged[-1].shutdown()
        #### return correct code
        self.EndModal(self.status)
        
    def onAction(self, event):
        action = self.controls_map[self.sec_id]["button"]["action"].GetLabel()
        # print "onAction", action
        if action == self.BTN_LABELS["ping"]:
            self.actPing()
        elif action == self.BTN_LABELS["start"]:
            self.actStart()
        elif action == self.BTN_LABELS["reconnect"]:
            self.actReconnect()
        elif action == self.BTN_LABELS["stop"]:
            self.actStop()

    def setactStop(self):
        self.controls_map[self.sec_id]["button"]["action"].SetLabel(self.BTN_LABELS["stop"])
        for i,f in enumerate(self.IPKN_FIELDS):
            self.controls_map[self.sec_id]["open"][f].Disable()

    def setactPing(self):
        self.controls_map[self.sec_id]["button"]["action"].SetLabel(self.BTN_LABELS["ping"])
        for i,f in enumerate(self.IPKN_FIELDS):
            self.controls_map[self.sec_id]["open"][f].Enable()            
            
    def actReconnect(self):
        ipkn_txtctrl = self.getIPKN_TxtCtrl()
        self.pinged[-1].resetNumClient(ipkn_txtctrl[-1])
        self.status = 2
        self.actStart()
        msg = self.wp_handle.getWP().infoStr() + " --- "
        self.updateInfo(msg + "Reconnecting, please wait...", ConnectionDialog.WARNING_FC)
        self.wp_handle.getWP().reconnection(self.boss)
        self.boss.checkResults(menu=True)
        self.updateInfo(msg + "Reconnected.", ConnectionDialog.SUCCESS_FC)
        
    def actStart(self):
        self.status += 1
        self.wp_handle.setWP(self.pinged[-1])
        self.pinged = None
        self.setactStop()
        msg = self.wp_handle.getWP().infoStr() + " --- "
        self.updateInfo(msg + "Started.", ConnectionDialog.SUCCESS_FC)

        
    def actPing(self):
        (ip, portnum, authkey, numclient) = self.getIPKN_TxtCtrl()
        self.no_problem = False
        tmpwp = self.wp_handle.setupWorkPlant(ip, portnum, authkey, numclient)
        color = ConnectionDialog.FAIL_FC
        try:
            status, info, id_clients = tmpwp.getDetailedInfos()            
            msg = tmpwp.infoStr() + " --- " + info
            if status == "OK":
                self.controls_map[self.sec_id]["button"]["action"].SetLabel(self.BTN_LABELS["start"])
                color = ConnectionDialog.SUCCESS_FC
                self.no_problem = True
                self.pinged = (ip, portnum, authkey, id_clients, tmpwp)
                
        except Exception as e:
            self.pinged = None
            self.no_problem = False
            msg = "Failed, check the parameters and try again (%s)" % e

        self.updateInfo(msg, color)
        self.Layout()
        self.Fit()
            
            
    def actStop(self):
        color = ConnectionDialog.FAIL_FC
        self.status = -1
        # self.wp_handle.getWP().shutdown()
        self.wp_handle.getWP().closeDown(self.boss)
        if self.no_problem:
            self.wp_handle.setWP(self.wp_handle.setupWorkPlant("local"))
            color = ConnectionDialog.SUCCESS_FC
            self.status = 0
        status, info, id_clients = self.wp_handle.getWP().getDetailedInfos()
        msg = self.wp_handle.getWP().infoStr() + " --- " + info
        self.updateInfo(msg, color)
        self.setactPing()
