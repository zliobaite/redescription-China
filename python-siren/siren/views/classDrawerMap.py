import wx, numpy, re
# The recommended way to use wx with mpl is with the WXAgg backend. 
import matplotlib
matplotlib.use('WXAgg')

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import mpl_toolkits.basemap

from classDrawerBasis import DrawerEntitiesTD, DrawerBasis
from classDrawerClust import DrawerClustTD

import pdb


class MapBase:

    # circ_equ=2*numpy.pi*6378137.
    # circ_pol=2*numpy.pi*6356752.
    # circ_avg=2*numpy.pi*6371000.
    circ_def=2*numpy.pi*6370997.

    marg_f = 100.0
    proj_def = "mill"
    proj_names = {"None": None,
                  "Gnomonic": "gnom",
                  "Mollweide": "moll",
                  "Gall Stereographic Cylindrical": "gall",
                  "Miller Cylindrical": "mill",
                  "Mercator": "merc",
                  "Hammer": "hammer",
                  "Geostationary": "geos",
                  "Near-Sided Perspective": "nsper",
                  "van der Grinten": "vandg",
                  "McBryde-Thomas Flat-Polar Quartic": "mbtfpq",
                  "Sinusoidal": "sinu",
                  "Lambert Conformal": "lcc",
                  "Equidistant Conic": "eqdc",
                  "Cylindrical Equidistant": "cyl",
                  "Oblique Mercator": "omerc",
                  "Albers Equal Area": "aea",
                  "Orthographic": "ortho",
                  "Cassini-Soldner": "cass",
                  "Robinson": "robin",
                  ######
                  "Azimuthal Equidistant": "aeqd",
                  "Lambert Azimuthal Equal Area": "laea",
                  "Stereographic": "stere",
                  #############
                  "Cylindrical Equal Area": "cea", 
                  "Eckert IV": "eck4", 
                  "Kavrayskiy VII": "kav7",
                  "Polyconic": "poly",
                  "N/S-Polar Lambert Azimuthal": "nplaea",
                  "N/S-Polar Stereographic": "npstere",
                  "N/S-Polar Azimuthal Equidistant": "npaeqd",
                  # "Rotated Pole": "rotpole",
                  # "Transverse Mercator": "tmerc",
                      }

    proj_pk = {"aeqd": ["lat_0", "lon_0", "width", "height"],
               "laea": ["lat_0", "lon_0","width", "height"],
               "stere": ["lat_0", "lon_0","width", "height"],
               ###
               "npaeqd": ["lon_0", "boundinglat"],
               "nplaea": ["lon_0", "boundinglat"],
               "npstere": ["lon_0", "boundinglat"],
               ###
               "spaeqd": ["lon_0", "boundinglat"],
               "splaea": ["lon_0", "boundinglat"],
               "spstere": ["lon_0", "boundinglat"],              
               ##############
               "geos": ["lon_0"],
               "vandg": ["lon_0"],
               "moll": ["lon_0"],
               "hammer": ["lon_0"],
               "robin": ["lon_0"],
               "mbtfpq": ["lon_0"],
               "sinu": ["lon_0"],
               "eck4": ["lon_0"], 
               "kav7": ["lon_0"],
               "ortho": ["lat_0", "lon_0"],
               "nsper": ["lat_0", "lon_0", "satellite_height"],
               "poly": ["lat_0", "lon_0","width", "height"],
               "gnom": ["lat_0", "lon_0", "width", "height"],
               "cass": ["lat_0", "lon_0", "width", "height"],
               #############
               "eqdc": ["lat_0", "lon_0", "lat_1", "lat_2","width", "height"],
               "aea": ["lat_0", "lon_0", "lat_1", "lat_2","width", "height"],
               "omerc": ["lat_0", "lon_0","lat_1", "lon_1", "lat_2", "lon_2","width", "height"],               
               "lcc": ["lat_0", "lon_0","lat_1", "lon_1", "lat_2", "lon_2","width", "height"],
               #############
               "cea": ["llcrnrlat", "llcrnrlon", "urcrnrlat", "urcrnrlon"],
               "cyl": ["llcrnrlat", "llcrnrlon", "urcrnrlat", "urcrnrlon"],
               "merc": ["llcrnrlat", "llcrnrlon", "urcrnrlat", "urcrnrlon"],
               "mill": ["llcrnrlat", "llcrnrlon", "urcrnrlat", "urcrnrlon"],
               "gall": ["llcrnrlat", "llcrnrlon", "urcrnrlat", "urcrnrlon"],
                }        
        
    bounds_def = {"llon": -180., "ulon": 180., "llat": -90., "ulat": 90.}
    # bounds_try = {"llon": -180., "ulon": 180., "llat": -90., "ulat": 90.}

    @classmethod    
    def getBasemapProjSetts(tcl, prefs):
        proj = tcl.proj_def 
        if "map_proj" in prefs:
            tpro = re.sub(" *\(.*\)$", "", prefs["map_proj"]["data"])
            if tpro in tcl.proj_names:
                proj = tcl.proj_names[tpro]
        resolution = "c"
        if "map_resolution" in prefs:
            resolution = prefs["map_resolution"]["data"][0]

        return proj, resolution

    @classmethod
    def getBasemapBackSetts(tcl, prefs):
        draws = {"rivers": False, "coasts": False, "countries": False,
                 "states": False, "parallels": False, "meridians": False,
                 "continents":False, "lakes": False, "seas": False}
        ### DEBUG
        draws = {"rivers": False, "coasts": True, "countries": False,
                 "states": False, "parallels": False, "meridians": False,
                 "continents":False, "lakes": False, "seas": False}

        colors = {"line_color": "gray", "sea_color": "#F0F8FF", "land_color": "white", "none":"white"}
        more = {}
        
        for typ_elem in ["map_elem_area", "map_elem_natural", "map_elem_geop", "map_elem_circ"]:
            if typ_elem in prefs:
                for elem in prefs[typ_elem]["data"]:                    
                    draws[elem] = True
                    
        for k in ["map_back_alpha", "map_back_scale"]:
            if k in prefs:
                more[k] = prefs[k]["data"]/100.0
            else:
                more[k] = 1.
        for k in ["map_back"]:
            if k in prefs:
                more[k] = prefs[k]["value"]
            else:
                more[k] = 0
            
        for color_k in colors.keys():
            if color_k in prefs:
                colors[color_k] = "#"+"".join([ v.replace("x", "")[-2:] for v in map(hex, prefs[color_k]["data"])]) 
        return draws, colors, more

    @classmethod
    def getParallelsRange(tcl, bm_args):
        span = float(bm_args["urcrnrlat"] - bm_args["llcrnrlat"])
        # if bm_args["llcrnrlat"] < bm_args["urcrnrlat"]:
        #     span = float(bm_args["urcrnrlat"] - bm_args["llcrnrlat"])
        # else:
        #     span = (180. - bm_args["llcrnrlon"]) + (bm_args["urcrnrlon"] + 180.)
        opts = [60, 30, 10, 5, 1]
        p = numpy.argmin(numpy.array([((span/k)-5.)**2 for k in opts]))
        step = opts[p]
        # if bm_args["llcrnrlon"] < bm_args["urcrnrlon"]:
        return numpy.arange(int(bm_args["llcrnrlat"]/step)*step, (int(bm_args["urcrnrlat"]/step)+1)*step, step)
        # else:
        #     return numpy.concatenate([numpy.arange(int(bm_args["llcrnrlon"]/step)*step, (int(180./step)+1)*step, step),
        #                                   numpy.arange(int(-180./step)*step, (int(bm_args["urcrnrlon"]/step)+1)*step, step)])

    @classmethod
    def getMeridiansRange(tcl, bm_args):
        if bm_args["llcrnrlon"] < bm_args["urcrnrlon"]:
            span = float(bm_args["urcrnrlon"] - bm_args["llcrnrlon"])
        else:
            span = (180. - bm_args["llcrnrlon"]) + (bm_args["urcrnrlon"] + 180.)
        opts = [60, 30, 10, 5, 1]
        p = numpy.argmin(numpy.array([((span/k)-5.)**2 for k in opts]))
        step = opts[p]
        if bm_args["llcrnrlon"] < bm_args["urcrnrlon"]:
            return numpy.arange(int(bm_args["llcrnrlon"]/step)*step, (int(bm_args["urcrnrlon"]/step)+1)*step, step)
        else:
            return numpy.concatenate([numpy.arange(int(bm_args["llcrnrlon"]/step)*step, (int(180./step)+1)*step, step),
                                          numpy.arange(int(-180./step)*step, (int(bm_args["urcrnrlon"]/step)+1)*step, step)])
        
    @classmethod
    def getBasemapCorners(tcl, prefs, cextrema):
        coords = ["llon", "ulon", "llat", "ulat"]
        allundef = True
        ## try_bounds = {"llon": -30., "ulon": 30., "llat": 30., "ulat": 110.}

        mbounds = dict([("c_"+c, v) for (c,v) in tcl.bounds_def.items()])
        mbounds.update(dict([("margc_"+c, 1./tcl.marg_f) for (c,v) in tcl.bounds_def.items()]))
        for c in coords:
            ### get corners from settings
            if c in prefs:
                mbounds["c_"+c] = prefs[c]["data"]
            allundef &=  (mbounds["c_"+c] == -1)
        if allundef:
            ### if all equal -1, set corners to def, globe wide
            mbounds.update(tcl.bounds_def)
        else:
            ### get corners from data
            mbounds["llon"], mbounds["ulon"], mbounds["llat"], mbounds["ulat"] = cextrema #self.getParentCoordsExtrema()
            for coord in coords:
                ### if corners coords from settings lower than 180,
                ### replace that from data, and drop margin
                if numpy.abs(mbounds["c_"+coord]) <= 180: #numpy.abs(tcl.bounds_def[coord]):
                    mbounds[coord] = mbounds["c_"+coord]                    
                    mbounds["margc_"+coord] = 0.

        for coord in ["lon", "lat"]:
            mbounds["marg_l"+coord] = mbounds["margc_l"+coord] * (mbounds["u"+coord]-mbounds["l"+coord]) 
            mbounds["marg_u"+coord] = mbounds["margc_u"+coord] * (mbounds["u"+coord]-mbounds["l"+coord]) 
        return mbounds

    @classmethod
    def greatCircleAngle(tcl, x0, x1):
        ### compute great circle distance for degree coordinates (lat, long)
        rd0 = numpy.radians(x0)
        rd1 = numpy.radians(x1)
        return numpy.arccos(numpy.sin(rd0[0]) * numpy.sin(rd1[0]) \
            + numpy.cos(rd0[0]) * numpy.cos(rd1[0]) * numpy.cos(rd0[1] - rd1[1]))
    @classmethod
    def greatCircleAngleXY(tcl, x0, x1):
        ### compute great circle distance for degree coordinates (x,y)
        rd0 = numpy.radians(x0)
        rd1 = numpy.radians(x1)
        return numpy.arccos(numpy.sin(rd0[1]) * numpy.sin(rd1[1]) \
            + numpy.cos(rd0[1]) * numpy.cos(rd1[1]) * numpy.cos(rd0[0] - rd1[0]))

    @classmethod
    def greatCircleDist(tcl, x0, x1):
        ### compute great circle distance for degree coordinates (lat, long)
        return tcl.circ_def * tcl.greatCircleAngle(x0, x1)
    @classmethod
    def greatCircleDistXY(tcl, x0, x1):
        ### compute great circle distance for degree coordinates (x,y)
        return tcl.circ_def * tcl.greatCircleAngleXY(x0, x1)
    
    @classmethod
    def makeBasemapProj(tcl, prefs, cextrema):
        proj, resolution = tcl.getBasemapProjSetts(prefs)
        if proj is None:
            return None, None
        mbounds = tcl.getBasemapCorners(prefs, cextrema)
        ## print "MBOUNDS", "\n".join(["%s:%s" % (k,v) for (k,v) in mbounds.items()])

        llcrnrlon = numpy.max([-180., mbounds["llon"]-mbounds["marg_llon"]])
        urcrnrlon = numpy.min([180., mbounds["ulon"]+mbounds["marg_ulon"]])
        if urcrnrlon <= llcrnrlon:
            if "lon_0" in tcl.proj_pk[proj]:
                span_lon = (360+urcrnrlon-llcrnrlon)
            else:
                urcrnrlon = tcl.bounds_def["ulon"]
                llcrnrlon = tcl.bounds_def["llon"]
                span_lon = (urcrnrlon-llcrnrlon)
        else:
            span_lon = (urcrnrlon-llcrnrlon)

        lon_0 = llcrnrlon + span_lon/2.0
        if lon_0 > 180:
            lon_0 -= 360
            
        llcrnrlat = numpy.max([-90., mbounds["llat"]-mbounds["marg_llat"]])
        urcrnrlat = numpy.min([90., mbounds["ulat"]+mbounds["marg_ulat"]])
        if urcrnrlat <= llcrnrlat:
            urcrnrlat = tcl.bounds_def["ulat"]
            llcrnrlat = tcl.bounds_def["llat"]
        if "lat_0" in tcl.proj_pk[proj]:
            llcrnrlatT = numpy.max([-180., mbounds["llat"]-mbounds["marg_llat"]])
            urcrnrlatT = numpy.min([180., mbounds["ulat"]+mbounds["marg_ulat"]])
        else:
            llcrnrlatT = llcrnrlat
            urcrnrlatT = urcrnrlat 
        span_lat = (urcrnrlatT-llcrnrlatT)
        lat_0 = llcrnrlatT + span_lat/2.0
        if numpy.abs(lat_0) > 90:
            lat_0 = numpy.sign(lat_0)*(180 - numpy.abs(lat_0))

        boundinglat = 0
        height = span_lat/360.
        if numpy.sign(urcrnrlat) == numpy.sign(llcrnrlat):
            width = numpy.cos((numpy.pi/2.)*numpy.min([numpy.abs(urcrnrlat),numpy.abs(llcrnrlat)])/90.)*span_lon/360.
            if urcrnrlat > 0:
                boundinglat = llcrnrlat
            else:
                boundinglat = urcrnrlat
        else: ### contains equator, the largest, factor 1
            width = span_lon/360.
        height = numpy.min([height, 0.5])
        width = numpy.min([width, 1.])
        args_all = {"width": tcl.circ_def*width, "height": tcl.circ_def*height,
                    "lon_0": lon_0, "lat_0": lat_0,
                    "lon_1": lon_0, "lat_1": lat_0,
                    "lon_2": lon_0+5, "lat_2": lat_0-5, 
                    "llcrnrlon": llcrnrlon, "llcrnrlat": llcrnrlat,
                    "urcrnrlon": urcrnrlon, "urcrnrlat": urcrnrlat,
                    "boundinglat": boundinglat, "satellite_height": 30*10**6}
        if args_all["lat_1"] == -args_all["lat_2"]:
            args_all["lat_2"] = +4

        args_p = {"projection": proj, "resolution":resolution}
        if proj in ["npaeqd", "nplaea", "npstere"]:
            if boundinglat > 0:                
                args_p["projection"] = "np"+proj[2:]
            elif boundinglat < 0:
                args_p["projection"] = "sp"+proj[2:]
            else:
                args_p["projection"] = proj[2:]
        for param_k in tcl.proj_pk[args_p["projection"]]:
            args_p[param_k] = args_all[param_k]
        # print "Proj", args_p["projection"], "H", height, "W", width, "Corners", (llcrnrlon, llcrnrlat), (urcrnrlon, urcrnrlat) #, "args", args_all
        # print "--- ARGS ALL\n", "\n".join(["%s:%s" % (k,v) for (k,v) in args_all.items()])
        # print "--- ARGS P\n", "\n".join(["%s:%s" % (k,v) for (k,v) in args_p.items()])
        try:
            bm = mpl_toolkits.basemap.Basemap(**args_p)
            # print "<< Basemap init succeded!", args_p
        except ValueError:
            # print ">> Basemap init failed!", args_p
            # print "H", height, "W", width, "Corners", (llcrnrlon, llcrnrlat), (urcrnrlon, urcrnrlat), "args", args_all
            bm = None 
        ### print "BM Corners", (bm.llcrnrlon, bm.llcrnrlat), (bm.urcrnrlon, bm.urcrnrlat)
        return bm, args_all

    @classmethod
    def makeBasemapBack(tcl, prefs, bm_args, bm=None):
        if bm is None:
            return
        draws, colors, more = tcl.getBasemapBackSetts(prefs)
        bounds_color, sea_color, contin_color, lake_color = colors["none"], colors["none"], colors["none"], colors["none"]
        if draws["rivers"]:
            bm.drawrivers(color=colors["sea_color"])
        if draws["coasts"]:
            bounds_color = colors["line_color"]
            bm.drawcoastlines(color=colors["line_color"])
        if draws["countries"]:
            bounds_color = colors["line_color"]
            bm.drawcountries(color=colors["line_color"])
        if draws["states"]:
            bounds_color = colors["line_color"]
            bm.drawstates(color=colors["line_color"])
        if draws["continents"]:
            contin_color = colors["land_color"]
        if draws["seas"]:
            sea_color = colors["sea_color"]
        if draws["lakes"]:
            lake_color = colors["sea_color"]
            
        if draws["parallels"]:
            tt = tcl.getParallelsRange(bm_args)
            # print "parallels", tt
            bm.drawparallels(tt, linewidth=0.5, labels=[1,0,0,1])
        if draws["meridians"]:
            tt = tcl.getMeridiansRange(bm_args)
            # print "meridians", tt
            bm.drawmeridians(tt, linewidth=0.5, labels=[0,1,1,0])

        func_map = {1: bm.shadedrelief, 2: bm.etopo, 3: bm.bluemarble}
        bd = False
        if more.get("map_back") in func_map:
            ### HERE http://matplotlib.org/basemap/users/geography.html
            try:
                func_map[more.get("map_back")](alpha=more["map_back_alpha"], scale=more["map_back_scale"])
                bd = True
            except IndexError:
                bd = False
                print "Impossible to draw the image map background!"
        if not bd:
            if bounds_color != colors["none"] or sea_color != colors["none"]:
                bm.drawmapboundary(color=bounds_color, fill_color=sea_color)
            if contin_color != colors["none"] or lake_color != colors["none"] or sea_color != colors["none"]:
                bm.fillcontinents(color=contin_color, lake_color=lake_color)
            # bm.drawlsmask(land_color=contin_color,ocean_color=sea_color,lakes=draws["lakes"])

class DrawerMap(DrawerBasis):
    
    MAP_POLY = True
    def initPlot(self):
        self.bm, self.bm_args = MapBase.makeBasemapProj(self.view.getParentPreferences(), self.getPltDtH().getParentCoordsExtrema())
        
        if self.bm is not None:
            self.getPltDtH().setBM(self.bm)
            self.setAxe(self.getFigure().add_axes([0, 0, 1, 1]))
            self.bm.ax = self.getAxe()
        else:
            llon, ulon, llat, ulat = self.getPltDtH().getParentCoordsExtrema()
            midlon, midlat = (llon + ulon)/2, (llat + ulat)/2
            mside = max(abs(llon-midlon), abs(llat-midlat))
            self.setAxe(self.getFigure().add_subplot(111,
              xlim=[midlon-1.05*mside, midlon+1.05*mside],
              ylim=[midlat-1.05*mside, midlat+1.05*mside]))

    def getAxisLims(self):
        xx = self.axe.get_xlim()
        yy = self.axe.get_ylim()
        return (xx[0], xx[1], yy[0], yy[1])

    def makeFinish(self, xylims=(0,1,0,1), xybs=(.1,.1)):
        DrawerBasis.makeFinish(self, xylims, xybs)
        self.drawCondionArea()
        
    def drawCondionArea(self):
        if self.getParentData() is not None and self.getParentData().isGeoConditional():
            red = self.getPltDtH().getRed()
            if red is not None and red.hasCondition():
                qC = red.getQueryC()
                if len(qC) == 0:
                    return
                ax_lims = list(self.getAxisLims())
                cond_lims = list(self.getPltDtH().getParentCoordsExtrema())
                for term in qC.invTerms():
                    cid = term.colId()
                    if term.isLowbounded():
                        cond_lims[2*cid] = term.getLowb()
                    if term.isUpbounded():
                        cond_lims[2*cid+1] = term.getUpb()
                # corners = [self.getCoordXYtoP(cond_lims[xi], cond_lims[2+yi]) for xi, yi in [(0,0), (0,1), (1,1), (1,0)]]
                edges = []
                stp = 0.05
                edges.extend([self.getCoordXYtoP(cond_lims[0], cond_lims[2] + x*(cond_lims[3]-cond_lims[2])) for x in numpy.arange(0.,1., stp)])
                edges.extend([self.getCoordXYtoP(cond_lims[0] + x*(cond_lims[1]-cond_lims[0]), cond_lims[3]) for x in numpy.arange(0.,1., stp)])
                edges.extend([self.getCoordXYtoP(cond_lims[1], cond_lims[2] + x*(cond_lims[3]-cond_lims[2])) for x in numpy.arange(1., 0., -stp)])                
                edges.extend([self.getCoordXYtoP(cond_lims[0] + x*(cond_lims[1]-cond_lims[0]), cond_lims[2]) for x in numpy.arange(1., 0., -stp)])

                self.axe.add_patch(Polygon(edges, closed=True, fill=True, fc="yellow", ec="yellow", zorder=10, alpha=.3))
                # cxs, cys = zip(*corners)
                # self.axe.plot(cxs, cys, "o", color="red", zorder=10)
                # cxs, cys = zip(*edges)
                # self.axe.plot(cxs, cys, "s", color="blue", zorder=10)

                
    def makeBackground(self):
        MapBase.makeBasemapBack(self.view.getParentPreferences(), self.bm_args, self.bm)

    def drawPoly(self):
        return self.getPltDtH().hasPolyCoords() & self.getSettV("map_poly", self.MAP_POLY)

    def getPosInfo(self, x, y):
        if self.bm is None:
            return (x, y)
        else:
            return self.bm(x, y, inverse=True)
    def getCoordXYtoP(self, x, y):
        if self.bm is None:
            return (x, y)
        else:
            return self.bm(x, y)

        
    
class DrawerEntitiesMap(DrawerMap, DrawerEntitiesTD): pass
    
class DrawerClustMap(DrawerMap, DrawerClustTD): pass
