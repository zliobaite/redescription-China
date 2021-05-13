import pdb
import sys
import numpy as np
import voronoi_poly
import matplotlib.pyplot as plt
import scipy.spatial.distance
from shapely import geometry
from shapely.ops import polygonize
from shapely.geos import TopologicalError

def makePolys(pdp, boundaries):
    PointsMap=dict([(p, (c1,c2)) for (p, c1, c2) in pdp])
    sx, sy = zip(*PointsMap.values())
    D=scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(np.array([sx,sy]).T, metric="euclidean"))
    nnvs = []
    nds = []

    for i in range(D.shape[0]):
        nnvs.append(np.argsort(D[:,i])[1:3])
    for i in range(D.shape[0]):
        nds.append(min([D[i,j] for j in nnvs[i]]+[D[nnvs[j][0],j] for j in nnvs[i]]+[D[nnvs[j][1],j] for j in nnvs[i]]))

    vl=voronoi_poly.VoronoiPolygonsMod(PointsMap, BoundingBox=boundaries)
    orgs = {}
    ready = {}
    ready[None] = []
    for s, obj in vl.items():
        pos = obj["info"]
        dst = nds[pos]
        fact = 1.41
        ready[pos] = []
        orgs[pos] = []
        tmc = getContours(obj['obj_polygon'])
        if pos in [622, 1509]:
            plt.figure()
            plt.plot([boundaries[tp] for tp in [2,1,1,2]], [boundaries[tp] for tp in [0,0,3,3]], 'g-')
            for edge in obj['obj_polygon']:
                xx,yy=zip(*edge)
                plt.plot(xx,yy, 'r+:')
            for dd in tmc:
                xc,yc=zip(*dd)
                plt.plot(xc,yc, 'bo--')
            plt.text(xc[0],yc[0], '%d' % pos)
            plt.show(block=False)
        
        for data_ct in tmc:
            orgs[pos].append(data_ct)
            inter_p = None
            data_poly = geometry.Polygon(data_ct)
            lfact = fact*(1+obj["coordinate"][1]/90.)**2
            iterf = 1.
            while inter_p is None:
                rst = [(obj["coordinate"][0]+x*dst*lfact*iterf/2, obj["coordinate"][1]+y*dst*fact*iterf/2) for (x,y) in [(-1,-1), (-1,1), (1,1), (1,-1), (-1,-1)]]
                # if iterf == 1:
                #     print "----", iterf, rst, dst
                restrict_poly = geometry.Polygon(rst)
                if pos in [622, 1509]:
                    xx,yy=zip(*data_ct)
                    xc,yc=zip(*rst)
                    plt.figure()
                    plt.plot(xx,yy, 'r')
                    plt.plot(xc,yc, 'b')
                    plt.text(xc[0],yc[0], '%d' % pos)
                    plt.show(block=False)
                try:
                    inter_p = restrict_poly.intersection(data_poly)
                except TopologicalError as e:
                    print "TOP ERROR", pos, obj["coordinate"]
                    # inter_p = None
                    inter_p = -1
                    # iterf += .5
                    # if iterf > 6.:
                    #     inter_p = -1

            if inter_p == -1:
                mdi = np.max(scipy.spatial.distance.cdist(np.array([obj["coordinate"]]), data_ct, metric="euclidean"))
                if mdi > fact*dst:
                    rst = [(obj["coordinate"][0]+x*dst*lfact/2, obj["coordinate"][1]+y*dst*fact/2) for (x,y) in [(-1,-1), (-1,1), (1,1), (1,-1), (-1,-1)]]                
                    ready[pos].append(rst)
                else:
                    ready[pos].append(data_ct)

            else:
                ready[pos].append(list(inter_p.exterior.coords))

                out_p = data_poly.difference(restrict_poly)
                if type(out_p) in [geometry.multipolygon.MultiPolygon, geometry.collection.GeometryCollection]:
                    for op in out_p:
                        ready[None].append((pos, list(op.exterior.coords)))
                else:
                    ready[None].append((pos, list(out_p.exterior.coords)))
    return ready, orgs



def getContours(obj_polygon):
    ends_map = {}
    for edge in obj_polygon:
        for end in [0,1]:
            if not ends_map.has_key(edge[end]):
                ends_map[edge[end]] = set([edge[1-end]])
            else:
                ends_map[edge[end]].add(edge[1-end])

    vertices = ends_map.keys()
    ### remove hanging nodes
    for v in vertices:
        if len(ends_map[v]) == 0:
            del ends_map[v]
        elif len(ends_map[v]) == 1:
            ends_map[ends_map[v].pop()].remove(v)
            del ends_map[v]

    contours = []
    loop = False
    contour = []
    while len(ends_map) > 0:
        if len(contour) == 0:
            contour = [ends_map.keys()[0]]
        while not loop and len(ends_map) > 0:
            if len(ends_map[contour[-1]].intersection(contour)) == 1:
                nv = ends_map[contour[-1]].intersection(contour).pop()
                ends_map[contour[-1]].remove(nv)
                loop = True
            else:
                nv = ends_map[contour[-1]].pop()
            if len(ends_map[contour[-1]]) == 0:
                del ends_map[contour[-1]]
            ends_map[nv].remove(contour[-1])
            contour.append(nv)

        if len(ends_map[nv]) == 0:
            del ends_map[nv]
        prev_p = contour.index(nv)
        contours.append(contour[prev_p:]+[nv])
        del contour[prev_p:]
        loop = False
    return contours


def main(argv=[]):
    style = "plain"
    marg_f = 100
    if len(argv) > 2:
        marg_f = float(argv[2])
    if len(argv) > 3 and argv[3] == "xml":
        style = "xml"
    coords = np.loadtxt(argv[1], unpack=True, usecols=(1,0))
    llon, ulon, llat, ulat = [min(coords[0]), max(coords[0]), min(coords[1]), max(coords[1])]
    blon, blat = (ulon-llon)/marg_f, (ulat-llat)/marg_f
    polys = makePolys(zip(range(len(coords[0])), coords[0], coords[1]), [ulat+blat, llon-blon, ulon+blon, llat-blat])
    for i in range(10): # len(coords[0])):
        if len(polys.get(i, [])) > 0:
            if style == "xml":
                print " ".join([":".join(map(str,co)) for co in zip(*polys[i][0])])
            else:
                print " ".join([",".join(map(str,co)) for co in polys[i][0]])
        else:
            print ""

if __name__ == '__main__':
    main(sys.argv)
