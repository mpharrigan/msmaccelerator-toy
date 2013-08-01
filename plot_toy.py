from matplotlib import pyplot as pp
from msmbuilder import io
import mdtraj
import mullerforce as mf
import os
import sqlite3 as sql

def plot_from_trajlist(trajlist):
    i = 0
    for t in trajlist:
        xs = t.xyz[:, 0, 0]
        ys = t.xyz[:, 0, 1]
        co = range(len(xs))
        pp.plot(xs, ys, 'o')
        i += 1

def get_trajlist_from_dir(directory):
    trajlist = list()
    files = os.listdir(directory)
    for f in files:
        if f.endswith(".h5"):
            trajlist.append(mdtraj.load(directory + "/" + f))
    return trajlist

def get_trajlist_from_sqlresult(sqlresult):
    trajlist = list()
    for row in sqlresult:
        filename = row[3]
        trajlist.append(mdtraj.load(filename))
    return trajlist

def get_model_from_sql(sqlresult):
    model = io.loadh(sqlresult[3])
    return model

def database_stuff(dbfile='db.sqlite'):
    # Connect
    conn = sql.connect(dbfile)
    c = conn.cursor()
    all_trajs = list()

    # Get a list of models
    c.execute("select id, time, protocol, path from models order by time asc")
    models = c.fetchall()
    for model in models:
        # For each model, find the trajectories that happened before it
        time = (model[1],)
        c.execute("""select id, time, protocol, path from trajectories
                where time < Datetime(?) order by time asc""", time)
        trajs = c.fetchall()
        all_trajs.append(trajs)
    return all_trajs, models

class ToyPlotter:

    def __init__(self, dbfile='db.sqlite'):
        self.all_trajs, self.all_models = database_stuff(dbfile)

    def plot_clustering_to(self, i):
        # Plot potential
        mf.MullerForce.plot()
        
        # Get trajectories
        trajs_sql = self.all_trajs[i]
        trajs = get_trajlist_from_sqlresult(trajs_sql)
        
        # Plot trajectories
        plot_from_trajlist(trajs)
        
        # Get generators
        model_sql = self.all_models[i]
        model = get_model_from_sql(model_sql)
        gen_i = model['generator_indices']
        points = list()
        for gen in gen_i:
            point = trajs[gen[0]][gen[1]].xyz[0,0]
            print(point[0:2])
            points.append(point[0:2]) # Just x,y
        
        mf.MullerForce.plot_voronoi(points)
        pp.show()

    def plot_trajs_to(self, i):
        # Plot potential
        mf.MullerForce.plot()
        
        # Get trajectories
        trajs_sql = self.all_trajs[i]
        trajs = get_trajlist_from_sqlresult(trajs_sql)
        
        # Plot trajectories and display
        plot_from_trajlist(trajs)
        pp.show()

    def plot_starting_structures(self, i):
        # Plot potential
        mf.MullerForce.plot()
        
        # Get generators
        model_sql = self.all_models[i]
        model = get_model_from_sql(model_sql)
        
        seed = io.loadh("seed_structures.h5")
        xs = seed['coordinates'][:,0,0]
        ys = seed['coordinates'][:,0,1]
        pp.plot(xs,ys,'o')
        pp.show()
        

def load_and_plot(directory='trajs/'):
    mf.MullerForce.plot()
    plot_from_trajlist(get_trajlist_from_dir(directory))
    
def test_cluster():
    p = ToyPlotter()
    p.plot_starting_structures(0)
    length = len(p.all_models)
    for i in range(length):
        p.plot_clustering_to(i)
    

if __name__ == "__main__":
    test_cluster()
