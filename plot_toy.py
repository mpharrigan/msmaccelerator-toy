import mdtraj
from matplotlib import pyplot as pp
import os
import mullerforce as mf

def plot_from_trajlist(trajlist):
    i = 0
    for t in trajlist:
        xs = t.xyz[:,0,0]
        ys = t.xyz[:,0,1]
        co = range(len(xs))
        pp.plot(xs,ys,'o')
        i+=1
    pp.show()

def get_trajlist(directory):
    trajlist = list()
    files = os.listdir(directory)
    for f in files:
        if f.endswith(".h5"):
            trajlist.append(mdtraj.load(directory + "/" + f))
    return trajlist

def load_and_plot(directory='trajs/'):
    mf.MullerForce.plot()
    plot_from_trajlist(get_trajlist(directory))

if __name__=="__main__":
    load_and_plot()