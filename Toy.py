
from msmbuilder import arglib
from toy_accel import app


def run_main(args):
    app.run_main(args)

def view_main(args):
    app.view_main(args)
    
def quant_main(args):
    app.quant_main(args)

if __name__ == "__main__":
    
    # Main parser
    parser = arglib.ArgumentParser(description="""Do everything with a 
            toy system.""")
    
    parser.add_argument('out_dir', help="""The directory to do everything in.""")
    parser.add_argument('top_fn', help="Topology file", default='single.pdb')
    
    subparsers = parser.add_subparsers(title="Available functions")
    
    # Main -> Run
    p_run = subparsers.add_parser('run')
    p_run.set_defaults(func=run_main)
    p_run.add_argument('-nstart', dest='num_starting_struct',
                       help="Number of random starting structures",
                       type=int,
                       default=10)
    p_run.add_argument('-s', dest='type_start',
                       help="What type of starting structures? random, fixed",
                       default='fixed')
    p_run.add_argument('-nstep', dest='n_steps',
                       help="Number of steps to simulate",
                       type=int,
                       default=10000)
    p_run.add_argument('-nrep', dest='rep_int',
                       help="Reporting interval",
                       type=int,
                       default=300)
    p_run.add_argument('-nround', dest='n_rounds',
                       help="Number of rounds",
                       type=int,
                       default=10)
    p_run.add_argument('-lt', dest='lag_time',
                       help="Lag time for MSM",
                       type=int,
                       default=1)
    p_run.add_argument('-b', dest='beta',
                       help="Beta for counts sampler. Need this for initial config as well",
                       type=float,
                       default=0.0)  

    run_subparsers = p_run.add_subparsers(title="Ways to run")
    
    # Main -> Run -> Constant
    p_run_const = run_subparsers.add_parser('constant')
    p_run_const.set_defaults(runtype='constant')
    p_run_const.add_argument('-neng', dest='n_engines',
                       help="Number of engines",
                       type=int,
                       default=5)  
    
    # Main -> Run -> Exponential
    p_run_exp = run_subparsers.add_parser('exponential')
    p_run_exp.set_defaults(runtype='exponential')
        
    # Main -> Run -> Varybeta
    p_run_varyb = run_subparsers.add_parser('varybeta')
    p_run_varyb.set_defaults(runtype='varybeta')
    p_run_varyb.add_argument('-neng', dest='n_engines',
                       help="Number of engines",
                       type=int,
                       default=5)
        
    
    # Main -> View
    p_view = subparsers.add_parser('view')

    p_view.add_argument('-db', dest='db_fn',
                        help="Database filename (relative)",
                        default='db.sqlite')
    p_view.add_argument('-fo', dest='fig_out_dir',
                        help="Figure output directory",
                        default='figs/')
    
    view_subparsers = p_view.add_subparsers(title="Do various things with data")
    
    # Main -> View -> Movie
    p_movie = view_subparsers.add_parser('movie')
    p_movie.set_defaults(func=view_main)
    p_movie.add_argument('-short', dest='is_short',
                         help="Make a short movie",
                         type=bool, default=False)
    p_movie.add_argument('-stride', dest='movie_stride',
                         help="Stride the movie frames by this",
                         type=int, default=10)

    
    # Main -> View -> Quant        
    p_quant = view_subparsers.add_parser('quant')
    p_quant.set_defaults(func=quant_main) 
    
    args = parser.parse_args()
    args.func(args)
    
