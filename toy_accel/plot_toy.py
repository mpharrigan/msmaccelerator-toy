from matplotlib import pyplot as pp
from msmbuilder import io
from msmbuilder import msm_analysis as msma
from msmaccelerator.core import markovstatemodel as msmm
import mdtraj
import mullerforce as mf
import os
import sqlite3 as sql
import numpy as np

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

def get_trajlist_from_fnlist(fnlist):
    trajlist = list()
    for fn in fnlist:
        trajlist.append(mdtraj.load(fn))
    return trajlist

def get_model_from_sql(sqlresult):
    # model = io.loadh(sqlresult[3])
    model = msmm.MarkovStateModel.load(sqlresult[3])
    return model

def get_starting_state_from_sql(sqlresult, top_fn):
    topology = mdtraj.load(top_fn)
    sslist = list()
    for row in sqlresult:
        filename = row[3]
        sslist.append(mdtraj.load(filename, top=topology))
    return sslist

def database_stuff(dbfile='db.sqlite'):
    # Connect
    conn = sql.connect(dbfile)
    c = conn.cursor()
    all_trajs = list()
    all_starting = list()

    # Get a list of models
    c.execute("select id, time, protocol, path from models order by time asc")
    models = c.fetchall()
    old_time = None
    for model in models:
        # Find each thing that is new to this particular model
        time = model[1]
        
        if old_time == None:
            condition = "where time < Datetime(?)"
            parameters = (time,)
        else:
            condition = "where time < Datetime(?) and time >= Datetime(?)"
            parameters = (time, old_time)
        
        # Get trajectories new for this model
        c.execute("""select id, time, protocol, path from trajectories
                %s order by time asc""" % condition, parameters)
        trajs = c.fetchall()
        all_trajs.append(trajs)
        
        # Get starting states new for this model
        c.execute("""select id, time, protocol, path from starting_states
                %s order by time asc""" % condition, parameters)
        sstates = c.fetchall()
        all_starting.append(sstates)
        
        old_time = time

                        
    return all_trajs, models, all_starting

class ToyPlotter:

    def __init__(self, dbfile, top_fn, fig_out_dir):
        self.all_trajs, self.all_models, self.all_starting = database_stuff(dbfile)
        self.top_fn = top_fn
        self.fig_out_dir = fig_out_dir + '/'
        self.fig_index = 0
        
        epsilon = 0.1
        limits = mf.MullerForce.plot()
        pp.clf()
        self.xlim = (limits[0], limits[1] - epsilon)
        self.ylim = (limits[2], limits[3] - epsilon)
        
    def save_fig(self):
        fn = self.fig_out_dir + 'frame%d.png' % self.fig_index
        
        if self.xlim is not None:
            pp.xlim(self.xlim)
        if self.ylim is not None:
            pp.ylim(self.ylim)
        
        pp.savefig(fn)
        pp.clf()
        self.fig_index += 1
        
    def show_fig(self):
        if self.xlim is not None:
            pp.xlim(self.xlim)
        if self.ylim is not None:
            pp.ylim(self.ylim)
            
        pp.show()
        

    def plot_clustering_to(self, i, show=False, save=True):
        # Plot potential
        limits = mf.MullerForce.plot()
        
        # Get the model
        model_sql = self.all_models[i]
        model = get_model_from_sql(model_sql)
        
        # Get and plot trajectories from the model
        traj_fns = model.traj_filenames
        trajs = get_trajlist_from_fnlist(traj_fns)
        for t in trajs:
            xs = t.xyz[:, 0, 0]
            ys = t.xyz[:, 0, 1]
            pp.plot(xs, ys, 'wo')
        
        # Get generators
        points = list()
        gen_i = model.generator_indices
        
        for gen in gen_i:
            point = trajs[gen[0]][gen[1]].xyz[0, 0]
            points.append(point[0:2])  # Just x,y
                
        
        mf.MullerForce.plot_voronoi(points)
        
        if show: self.show_fig()
        if save: self.save_fig()
        return limits
        
        
    def plot_trajs_movie_at(self, i, show=False, save=True):
        trajs_sql = self.all_trajs[i]
        trajs = get_trajlist_from_sqlresult(trajs_sql)
        
        n_frames = trajs[0].n_frames
        for frame in xrange(n_frames):
            # Plot the potential
            mf.MullerForce.plot()

            # Plot each trajectory
            for t in trajs:
                xs = t.xyz[0:frame, 0, 0]
                ys = t.xyz[0:frame, 0, 1]
                pp.plot(xs, ys, 'o-')
                
            # Put starting states on top
            self.plot_starting_states(i, show=False)
            
            if show: self.show_fig()
            if save: self.save_fig()

    def plot_starting_states(self, i, show=False, save=False):
        # Get states
        sstate_sql = self.all_starting[i]
        sstates = get_starting_state_from_sql(sstate_sql, self.top_fn)
        
        for sstate in sstates:
            x = sstate.xyz[0, 0, 0]
            y = sstate.xyz[0, 0, 1]
            pp.plot(x, y, 'o', markersize=12)
            
        
        if show: self.show_fig()
        if save: self.save_fig()
        
    def plot_starting_after(self, i, show=False, save=True):
        self.plot_clustering_to(i, show=False, save=False)
        self.plot_starting_states(i + 1, show=False, save=False)
        if show: self.show_fig()
        if save: self.save_fig()
        
    def get_implied_timescales(self, i, num_vals=2):
        # Get the model
        model_sql = self.all_models[i]
        model = get_model_from_sql(model_sql)
        
        transition_mat = model.transition_matrix
        lag_time = model.lag_time
        
        print('\n\n')
        print(i)
        
        if not msma.is_transition_matrix(transition_mat):
            print ("%d is not a transition matrix!" % i)
            
        eigenvals, _ = msma.get_eigenvectors(transition_mat, num_vals + 1, epsilon=1.0)
        eigenvals = eigenvals[1:]
        
        oo_lambda = [-1.0 / np.log(ev) for ev in eigenvals]
        
        print(oo_lambda)
        
        return oo_lambda
        
    def plot_implied_timescales(self, gold_vals):
        length = len(self.all_models)
        
        oo_lambdas = np.array([self.get_implied_timescales(i) for i in range(length)])
        num_ev = oo_lambdas.shape[1]
        
        for j in range(num_ev):
            ys = oo_lambdas[:,j]
            pp.plot(ys, 'o-')
    
        pp.hlines(gold_vals, 0, oo_lambdas.shape[0])
    
        pp.yscale('log')
        pp.show()

def view_starting_states(db_fn, top_fn, fig_out_dir):
    p = ToyPlotter(db_fn, top_fn, fig_out_dir)
    length = len(p.all_models)
    for i in range(length):
        p.plot_starting_states(i)  
    
def view_clustering(db_fn, top_fn, fig_out_dir):    
    p = ToyPlotter(db_fn, top_fn, fig_out_dir)
    length = len(p.all_models)
    for i in range(length):
        p.plot_clustering_to(i)
        
def view_movie(db_fn, top_fn, fig_out_dir):
    p = ToyPlotter(db_fn, top_fn, fig_out_dir)
    length = len(p.all_models)
    # for i in [0, 1, 2, length - 1]:
    for i in range(length):
        print("Plotting model %d" % i)
        p.plot_trajs_movie_at(i)
        p.plot_clustering_to(i)
        if i + 1 < length:
            p.plot_starting_after(i)
            
def quant_implied_timescales(db_fn, fig_out_dir):
    p = ToyPlotter(db_fn, None, fig_out_dir)
    p.plot_implied_timescales(gold_vals=[1551.418538, 23.88297])

