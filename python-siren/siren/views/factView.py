import classViewBasis
# import classViewBasis
# import classGView
# import classTDView
# import classMapView
# import classParaView
# import classViewEntitiesProj
# import classViewRedTree
#from classVProjView import VProjView

# import classViewList
# import classOverView
# import classInterView
# import classSimView

import pdb

def all_subclasses(cls):
    return cls.__subclasses__() + [g for s in cls.__subclasses__()
                                   for g in all_subclasses(s)]    

class ViewFactory(object):

    # typ_views = {"R": classGView.GView,
    #              "L": classViewList.ViewList}
    typ_views = {"R": classViewBasis.ViewRed,
                 "L": classViewBasis.ViewList}

        
    details_views_typs = {}
    viewsT_typs_map = {}
    for (typ, parent_class) in typ_views.items(): 
        details_views_typs[typ] = {}
        for cls in all_subclasses(parent_class):
            details_views_typs[typ].update(cls.getViewsDetails())
        for vt in details_views_typs[typ].keys():
            viewsT_typs_map[vt] = typ


    # @classmethod
    # def reloadCode(tcl):
    #     global classOverView
    #     classOverView = reload(classOverView)
    #     tcl.details_views_typs['L']['INT']['class'] = classInterView.InterView
    #     print "Reloaded InterView", tcl.details_views_typs['L']['INT']['class'].cversion
        
    @classmethod
    def getClasses(tcl, typv="R"):
        if typv in tcl.details_views_typs:
            return tcl.details_views_typs[typv]
        return {}

    @classmethod
    def getTypV(tcl, viewT):
        return tcl.viewsT_typs_map.get(viewT)

    @classmethod
    def getViewsInfo(tcl, typv="R", tabT=None, geo=False, what=None, excludeT=None):
        infos = [{"viewT": viewT, "title": details["title"], "short_title": details.get("short_title", details["title"]),
                  "ord": details["ord"], "suitable":details["class"].suitableView(geo, what, tabT)} \
                 for viewT, details in tcl.getClasses(typv).items() if (excludeT is None or viewT not in excludeT)]
        infos.sort(key=lambda x: (x["ord"], x["title"]))
        return infos

    @classmethod
    def getView(tcl, viewT, parent, vid):
        if tcl.getTypV(viewT) is not None:
            viewDet = tcl.details_views_typs[tcl.getTypV(viewT)][viewT]
            return viewDet["class"](parent, vid, viewDet["more"])

    @classmethod
    def getDefaultViewT(tcl, typv="R", geo=False, type_tab="r"):
        if typv == "L":
            return classViewBasis.ViewClustMap.TID
        elif typv == "R":
            if type_tab == "e":
                return classViewBasis.ViewEntitiesProj.defaultViewT
            elif geo:
                return classViewBasis.ViewRedMap.TID
            else:
                return classViewBasis.ViewRedPara.TID
