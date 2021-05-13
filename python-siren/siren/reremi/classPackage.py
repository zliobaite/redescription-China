import tempfile
import os.path
import plistlib
import shutil
import zipfile
import codecs
import re

import pdb

from classRedescription import Redescription
from classData import Data, NA_str_def
from classQuery import Query
from classPreferencesManager import PreferencesReader, getPM
import toolRead as toolRead


class Package(object):
    """Class to handle the zip packages that contain data, preferences, results, etc. for redescription mining.
    """

    # CONSTANTS
    # Names of the files in the package
    DATA_FILENAMES = ['data_LHS.csv',
                     'data_RHS.csv']
    REDESCRIPTIONS_FILENAME = 'redescriptions.csv'
    PREFERENCES_FILENAME = 'preferences.xml'
    PLIST_FILE = 'info.plist'
    PACKAGE_NAME = 'siren_package'

    FILETYPE_VERSION = 5
    XML_FILETYPE_VERSION = 3

    NA_str = NA_str_def
    RED_FN_SEP = ";"

    CREATOR = "ReReMi/Siren Package"
    DEFAULT_EXT = ".siren"
    DEFAULT_TMP = "siren"

    def __init__(self, filename, callback_mess=None, mode="r"):
        if filename is not None:
            filename = os.path.abspath(filename)
            if mode !="w" and not os.path.isfile(filename):
                raise IOError('File does not exist')
            if mode !="w" and not zipfile.is_zipfile(filename):
                raise IOError('File is of wrong type')
        self.filename = filename
        self.callback_mess = callback_mess
        self.plist = dict(creator = self.CREATOR,
                          filetype_version = self.FILETYPE_VERSION)

    def __str__(self):
        return "PACKAGE: %s" % self.filename

    def raiseMess(self):
        if self.callback_mess is not None:
            self.callback_mess()

    def getFilename(self):
        return self.filename

    def getPackagename(self):
        return self.plist.get('package_name')

    def getFormatV(self):
        return self.plist.get('filetype_version', -1)
    def isOldXMLFormat(self):
        return self.getFormatV() <= self.XML_FILETYPE_VERSION
    def isLatestFormat(self):
        return self.getFormatV() == self.FILETYPE_VERSION

    def getSaveFilename(self):
        svfilename = self.filename
        if self.isOldXMLFormat():
            parts = self.filename.split(".")
            if len(parts) == 1:
                svfilename += "_new"
            elif len(parts) > 1:
                svfilename = ".".join(parts[:-1]) + "_new."+ parts[-1]
        return svfilename

    def getNamelist(self):
        return self.package.namelist()

    def closePack(self):
        if self.package is not None:
            self.package.close()
            self.package = None
            
    def openPack(self):
        try:
            self.package = zipfile.ZipFile(self.filename, 'r')
            plist_fd = self.package.open(self.PLIST_FILE, 'r')
            self.plist = plistlib.readPlist(plist_fd)
        except Exception:
            self.package = None
            self.plist = {}
            self.raiseMess()
            raise

######### READING ELEMENTS
##########################

    def read(self, pm, options_args=None):
        elements_read = {}
        self.openPack()

        try:
            preferences = self.readPreferences(pm, options_args)
            if preferences is not None:
                elements_read["preferences"] = preferences
            data = self.readData()
            if data is not None:
                elements_read["data"] = data
                reds = self.readRedescriptions(data)
                if reds is not None and len(reds) > 0:
                    elements_read["reds"] = reds
        finally:
            self.closePack()
        return elements_read

    def readPreferences(self, pm, options_args=None):
        # Load preferences
        preferences = None
        if 'preferences_filename' in self.plist:
            fd = None
            try:
                fd = self.package.open(self.plist['preferences_filename'], 'r')
                preferences = PreferencesReader(pm).getParameters(fd, options_args)
            except Exception:
                self.raiseMess()
                raise
            finally:
                if fd is not None:
                    fd.close()
        return preferences


    def readData(self):
        data = None
        # Load data
        if self.isOldXMLFormat():
            ################# START FOR BACKWARD COMPATIBILITY WITH XML
            if 'data_filename' in self.plist:
                try:
                    fd = self.package.open(self.plist['data_filename'], 'r')
                    data = Data.readDataFromXMLFile(fd)
                except Exception as e:
                    print e
                    self.raiseMess()
                    raise
                finally:
                    fd.close()

        elif 'data_LHS_filename' in self.plist:
            try:
                fdLHS = self.package.open(self.plist['data_LHS_filename'], 'r')
                if self.plist.get('data_RHS_filename', self.plist['data_LHS_filename']) != self.plist['data_LHS_filename']:
                    fdRHS = self.package.open(self.plist['data_RHS_filename'], 'r')
                else:
                    fdRHS = None
                NA_str = self.plist.get('NA_str', None)
                if NA_str is None and self.getFormatV() <= 4:
                    NA_str = Package.NA_str
                # pdb.set_trace()    
                data = Data([fdLHS, fdRHS, {}, NA_str], "csv")
            except Exception:
                data = None
                self.raiseMess()
                raise
            finally:
                fdLHS.close()
                if fdRHS is not None: 
                    fdRHS.close()
        return data

    def readRedescriptions(self, data):
        reds = []
        # Load redescriptions
        rp = Redescription.getRP()
        if 'redescriptions_filename' in self.plist:
            for file_red in self.plist['redescriptions_filename'].split(self.RED_FN_SEP):            
                try:
                    fd = self.package.open(file_red, 'r')
                    if self.isOldXMLFormat():
                        rs, rshowids = readRedescriptionsXML(fd, data)
                    else:
                        rs = []
                        rp.parseRedList(fd, data, rs)
                        rshowids = range(len(rs))
                except Exception:
                    self.raiseMess()
                    raise
                finally:
                    fd.close()
                reds.append({"items": rs, "src": ('file', file_red, 1), "rshowids": rshowids})
        return reds

######### WRITING ELEMENTS
##########################
    def getTmpDir(self):
        return tempfile.mkdtemp(prefix=self.DEFAULT_TMP)
            
    ## The saving function
    def writeToFile(self, filename, contents):
        # Store old package_filename
        old_package_filename = self.filename
        self.filename = os.path.abspath(filename)
        # Get a temp folder
        tmp_dir = self.getTmpDir()
        #package_dir = os.path.join(tmp_dir, filename)
        #os.mkdir(package_dir)

        # Write plist
        plist, filens = self.makePlistDict(contents)
        try:
            plistlib.writePlist(plist, os.path.join(tmp_dir, self.PLIST_FILE))
        except IOError:
            shutil.rmtree(tmp_dir)
            self.filename = old_package_filename
            self.raiseMess()
            raise

        # Write data files
        if "data" in contents:
            try:
                filenames = [os.path.join(tmp_dir, plist['data_LHS_filename']), None]
                if plist.get('data_RHS_filename', plist['data_LHS_filename']) != plist['data_LHS_filename']:
                    filenames[1] = os.path.join(tmp_dir, plist['data_RHS_filename'])
                writeData(contents["data"], filenames, toPackage = True)
            except IOError:
                shutil.rmtree(tmp_dir)
                self.filename = old_package_filename
                self.raiseMess()
                raise

        # Write redescriptions
        if "redescriptions" in contents:
            for rs in contents["redescriptions"]:            
                try:
                    writeRedescriptions(rs.get("items", []), os.path.join(tmp_dir, os.path.basename(rs["src"][1])),
                                        rs.get("rshowids"), names=False, with_disabled=True, toPackage = True)
                except IOError:
                    shutil.rmtree(tmp_dir)
                    self.filename = old_package_filename
                    self.raiseMess()
                    raise

        # Write preferences
        if "preferences" in contents:
            try:
                writePreferences(contents["preferences"], contents["pm"],
                                 os.path.join(tmp_dir, plist['preferences_filename']), toPackage = True)
            except IOError:
                shutil.rmtree(tmp_dir)
                self.filename = old_package_filename
                self.raiseMess()
                raise

        # All's there, so pack
        try:
            package = zipfile.ZipFile(self.filename, 'w')
            package.write(os.path.join(tmp_dir, self.PLIST_FILE),
                          arcname = os.path.join('.', self.PLIST_FILE))
            for eln, element in filens.items():                
                package.write(os.path.join(tmp_dir, element),
                              arcname = os.path.join('.', element),
                              compress_type = zipfile.ZIP_DEFLATED)
        except Exception:
            shutil.rmtree(tmp_dir)
            self.filename = old_package_filename
            self.raiseMess()
            raise
        finally:
            package.close()

        # All's done, delete temp file
        shutil.rmtree(tmp_dir)

    
    def makePlistDict(self, contents):
        """Makes a dict to write to plist."""
        d = dict(creator = self.CREATOR,
            filetype_version = self.FILETYPE_VERSION)
        
        if self.filename is None:
            d['package_name'] = self.PACKAGE_NAME
        else:
            (pn, suffix) = os.path.splitext(os.path.basename(self.filename))
            if len(pn) > 0:
                d['package_name'] = pn
            else:
                d['package_name'] = self.PACKAGE_NAME

        fns = {}              
        if "data" in contents:
            d['NA_str'] = contents["data"].NA_str
            fns['data_LHS_filename'] = self.DATA_FILENAMES[0]
            if not contents["data"].isSingleD():
                fns['data_RHS_filename'] = self.DATA_FILENAMES[1]
                                
        if "preferences" in contents:
            fns['preferences_filename'] = self.PREFERENCES_FILENAME
        d.update(fns)
        
        if "redescriptions" in contents and len(contents["redescriptions"]) > 0:
            base_names = [os.path.basename(c["src"][1]) for c in contents["redescriptions"]]
            d['redescriptions_filename'] = self.RED_FN_SEP.join(base_names)
            for ci, c in enumerate(base_names):
                fns['redescriptions_filename_%d' %ci] = c
        return d, fns


######### ADDITIONAL FUNCTIONS
###############################

def readRedescriptionsXML(filep, data):
    red = []
    show_ids = None
    
    doc = toolRead.parseXML(filep)
    if doc is not None:
        tmpreds = doc.getElementsByTagName("redescriptions")
        if len(tmpreds) == 1:
            reds_node = tmpreds[0]
            for redx in reds_node.getElementsByTagName("redescription"):
                tmp = Redescription()
                tmp.fromXML(redx)
                tmp.recompute(data)
                reds.append(tmp)
            tmpsi = reds_node.getElementsByTagName("showing_ids")
            if len(tmpsi) == 1:
                show_ids = toolRead.getValues(tmpsi[0], int)
                if len(show_ids) == 0 or min(show_ids) < 0 or max(show_ids) >= len(reds):
                    show_ids = None
    if show_ids is None:
        show_ids = range(len(reds))
    return reds, rshowids


def writeRedescriptions(reds, filename, rshowids=None, names = [None, None], with_disabled=False, toPackage = False, style="", full_supp=False, nblines=1, supp_names=None):
    if names is False:
        names = [None, None]
    if rshowids is None:
        rshowids = range(len(reds))
    red_list = [reds[i] for i in rshowids if reds[i].getEnabled() or with_disabled]
    if toPackage:
        fields_supp = [-1, "status_enabled"]
    else:
        fields_supp = None
        # with codecs.open(filename, encoding='utf-8', mode='w') as f:
    with open(filename, mode='w') as f:
        rp = Redescription.getRP()
        if style == "tex":
            f.write(codecs.encode(rp.printTexRedList(red_list, names, fields_supp, nblines=nblines), 'utf-8','replace'))
        else:
            f.write(codecs.encode(rp.printRedList(red_list, names, fields_supp, full_supp=full_supp, supp_names=supp_names, nblines=nblines), 'utf-8','replace'))
            
def writePreferences(preferences, pm, filename, toPackage=False, inc_def=False, core=False):
    sections = True
    helps = False
    if toPackage:
        sections = False
    if preferences is None or inc_def:
        helps = True
    with open(filename, 'w') as f:
        f.write(PreferencesReader(pm).dispParameters(preferences, sections, helps, inc_def, core))

def writeData(data, filenames, toPackage = False):
    data.writeCSV(filenames)

def saveAsPackage(filename, data, preferences=None, pm=None, reds=None):
    package = Package(None, None, mode="w")

    (filename, suffix) = os.path.splitext(filename)
    contents = {}
    if data is not None:
        contents['data'] = data                                
    if reds is not None and len(reds) > 0:
        contents['redescriptions'] = (self.REDESCRIPTIONS_FILENAME, reds, range(len(reds)))
    if preferences is not None:
        if pm is None:
            pm = getPM()
        contents['preferences'] = preferences
        contents['pm'] = pm

    package.writeToFile(filename+suffix, contents)

    
def getPrintParams(filename, data=None):
    params = {"with_disabled": False, "style": "", "full_supp":False, "nblines":1,
              "names": [None, None], "supp_names": None}

    named = re.search("[^a-zA-Z0-9]named[^a-zA-Z0-9]", filename) is not None
    supp_names = ( re.search("[^a-zA-Z0-9]suppnames[^a-zA-Z0-9]", filename) is not None ) or \
                 ( re.search("[^a-zA-Z0-9]suppids[^a-zA-Z0-9]", filename) is not None )

    params["with_disabled"] = re.search("[^a-zA-Z0-9]all[^a-zA-Z0-9]", filename) is not None
    params["full_supp"] = ( re.search("[^a-zA-Z0-9]support[^a-zA-Z0-9]", filename) is not None ) or supp_names
            
    if re.search(".tex$", filename):
        params["style"] = "tex"

    tmp = re.search("[^a-zA-Z0-9](?P<nbl>[1-3]).[a-z]*$", filename)
    if tmp is not None:
        params["nblines"] = int(tmp.group("nbl"))

    if named and data is not None:
        params["names"] = data.getNames()            
    if supp_names:
        params["supp_names"] = data.getRNames()            
    return params
