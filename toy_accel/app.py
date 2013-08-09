import toy_accel.plot_toy as plot_toy
import toy_accel.generate_system as generate_system
import toy_accel.run_accel as run_accel
import toy_accel.euclidean as euclidean
import os

def generate_main(args):
    # Parse arguments
    top_fn = args.top_fn
    out_dir = args.out_dir
    num_starting_struct = args.num_starting_struct
    type_start = args.type_start
    x = 0.5
    y = 0.0
    rep_int = args.rep_int
    n_steps = args.n_steps
    beta = args.beta
    lag_time = args.lag_time


    # Set up filenames
    os.mkdir(out_dir)
    seed_rel_fn = 'seed_structures.h5'
    seed_out_fn = out_dir + '/' + seed_rel_fn
    sys_fn = out_dir + '/' + 'system.xml'
    int_fn = out_dir + '/' + 'integrator.xml'
    config_fn = out_dir + '/' + 'msmaccelerator_config.py'
    rel_top_fn = "../" + top_fn
    metric_out_fn = out_dir + '/' + 'metric.pickl'

    # Generate structures
    if type_start == 'random':
        generate_system.generate_starting_structures_random(top_fn,
                num_starting_struct, seed_out_fn)
    elif type_start == 'fixed':
        generate_system.generate_starting_structures_fixed(
                top_fn, x, y, seed_out_fn)

    # Generate OpenMM stuff
    generate_system.generate_openmm(sys_fn, int_fn)

    # Generate config file
    generate_system.generate_config_file(rep_int, n_steps,
            rel_top_fn, seed_rel_fn, beta, config_fn, lag_time)

    # Generate custom metric
    euclidean.make(metric_out_fn)



def run_main(args):
    generate_main(args)
    run_accel.run_accel(args)

def view_main(args):
    db_fn = args.db_fn
    out_dir = args.out_dir
    top_fn = args.top_fn
    fig_out_dir = args.fig_out_dir
    is_short = args.is_short
    stride = args.movie_stride
    format_str = args.format

    if format_str == 'accel':
        suptitle = "MSMAccelerator"
        show_beta = True
    elif format_str == 'onetraj':
        suptitle = "One trajectory"
        show_beta = False
    else:
        print("Please choose a valid format option")
        return

    # Set up filenames
    db_out_fn = out_dir + '/' + db_fn
    fig_out_dir = out_dir + '/' + fig_out_dir
    os.mkdir(fig_out_dir)

    if is_short:
        plot_toy.view_clustering(db_out_fn, top_fn, fig_out_dir, suptitle)
    else:
        plot_toy.view_movie(db_out_fn, top_fn, fig_out_dir, stride, suptitle,
                            show_beta)
    # plot_toy.view_starting_states(db_out_fn, top_fn)
    # plot_toy.view_clustering(db_out_fn, top_fn)

def quant_main(args):
    """Perform quantitative analysis on the model(s).

    This will probably involve calculating two implied timescales
    and comparing them to analytical solution and/or really long run
    markov model gold standard.
    """
    db_fn = args.db_fn
    out_dir = args.out_dir
    fig_out_dir = args.fig_out_dir

    # Set up filenames
    db_out_fn = out_dir + '/' + db_fn
    fig_out_dir = out_dir + '/' + fig_out_dir
    # os.mkdir(fig_out_dir)
    plot_toy.quant_implied_timescales(db_out_fn, fig_out_dir)
