"""Propagating 2D dynamics on the muller potential using OpenMM.

Currently, we just put a harmonic restraint on the z coordinate,
since OpenMM needs to work in 3D. This isn't really a big deal, except
that it affects the meaning of the temperature and kinetic energy. So
take the meaning of those numbers with a grain of salt.
"""
from simtk.unit import kelvin, picosecond, femtosecond, nanometer, dalton
import matplotlib as mpl
import matplotlib.pyplot as pp
import numpy as np
import simtk.openmm as mm
import voronoi


class MullerForce(mm.CustomExternalForce):
    """OpenMM custom force for propagation on the Muller Potential. Also
    includes pure python evaluation of the potential energy surface so that
    you can do some plotting"""
    aa = [-1, -1, -6.5, 0.7]
    bb = [0, 0, 11, 0.6]
    cc = [-10, -10, -6.5, 0.7]
    AA = [-200, -100, -170, 15]
    XX = [1, 0, -0.5, -1]
    YY = [0, 0.5, 1.5, 1]
    strength = 0.5

    def __init__(self):
        # start with a harmonic restraint on the Z coordinate
        expression = '1000.0 * z^2'
        for j in range(4):
            # add the muller terms for the X and Y
            fmt = dict(aa=self.aa[j], bb=self.bb[j], cc=self.cc[j], AA=self.AA[j], XX=self.XX[j], YY=self.YY[j])
            expression += '''+ {AA}*exp({aa} *(x - {XX})^2 + {bb} * (x - {XX}) 
                               * (y - {YY}) + {cc} * (y - {YY})^2)'''.format(**fmt)
        
        # Include scaling expression
        expression = "{strength}*(".format(strength=self.strength) + expression + ")"
        super(MullerForce, self).__init__(expression)

    @classmethod
    def potential(cls, x, y):
        "Compute the potential at a given point x,y"
        value = 0
        for j in range(4):
            value += cls.AA[j] * np.exp(cls.aa[j] * (x - cls.XX[j]) ** 2 + \
                cls.bb[j] * (x - cls.XX[j]) * (y - cls.YY[j]) + cls.cc[j] * (y - cls.YY[j]) ** 2)
        return value

    @classmethod
    def plot(cls, ax=None, minx=-1.5, maxx=1.2, miny=-0.2, maxy=2, **kwargs):
        "Plot the Muller potential"
        grid_width = max(maxx - minx, maxy - miny) / 200.0
        ax = kwargs.pop('ax', None)
        xx, yy = np.mgrid[minx : maxx : grid_width, miny : maxy : grid_width]
        V = cls.potential(xx, yy)
        # clip off any values greater than 200, since they mess up
        # the color scheme
        if ax is None:
            ax = pp
        ax.contourf(xx, yy, V.clip(max=200), 40, **kwargs)
        return minx, maxx, miny, maxy

    @classmethod
    def plot_voronoi(cls, points, limits=None, color='black'):
        PL = []
        for p in points:
            PL.append(voronoi.Site(x=p[0], y=p[1]))

        v, eqn, edges, wtf = voronoi.computeVoronoiDiagram(PL)

        xmin = wtf[0]
        ymin = wtf[1]
        xmax = wtf[2]
        ymax = wtf[3]

        if not limits is None:
            xmin = limits[0]
            xmax = limits[1]
            ymin = limits[2]
            ymax = limits[3]

        edge_points = []
        for (l, x1, x2) in edges:
            if x1 >= 0 and x2 >= 0:
                edge_points.append((v[x1], v[x2]))
            elif x1 == -1:
                a, b, c = eqn[l]
                x = xmin
                y = (c - x * a) / b
                if y < ymin:
                    y = ymin
                    x = (c - y * b) / a
                elif y > ymax:
                    y = ymax
                    x = (c - y * b) / a
                edge_points.append(((x, y), v[x2]))
            elif x2 == -1:
                a, b, c = eqn[l]
                x = xmax
                y = (c - x * a) / b
                if y < ymin:
                    y = ymin
                    x = (c - y * b) / a
                elif y > ymax:
                    y = ymax
                    x = (c - y * b) / a
                edge_points.append((v[x1], (x, y)))

        lines = mpl.collections.LineCollection(edge_points, linewidths=1, color=color)

        pp.gca().add_collection(lines)

        return wtf
