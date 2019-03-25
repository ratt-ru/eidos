import multiprocessing
from multiprocessing import Lock

def fun(f, q_in, q_out):
    """
    A helper function for the parmap function.

    Parameters
    ----------
    f : function
        Function to be evaluated using arguments in the queue q_in.
    q_in : multiprocessing Queue object
        A queue for the input arguments.
    q_out : multiprocessing Queue object
        A queue for the output of the function evaluation.
    """
    while True:
        i, x = q_in.get()
        if i is None:
            break
        q_out.put((i, f(x)))

def parmap(f, X, proc_power=1):
    """
    A parallelized implementation of the map function. See
    https://www.w3schools.com/python/ref_func_map.asp .

    Parameters
    ----------
    f : function
        Function onto which to map the input arguments in X.
    X : array_like
        The arguments to be fed to the function f. This can only handle a
        single argument for each evaluation of f.

    Returns
    -------
    out : list
        A list of the outputs for each function evaluation corresponding to the
        input arguments in X.
    """
    if proc_power<=1 and proc_power>0:
        nprocs=int(proc_power*multiprocessing.cpu_count())
    else:
        nprocs = multiprocessing.cpu_count()

    q_in = multiprocessing.Queue(1)
    q_out = multiprocessing.Queue()

    proc = [multiprocessing.Process(target=fun, args=(f, q_in, q_out)) for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()

    sent = [q_in.put((i, x)) for i, x in enumerate(X)]
    [q_in.put((None, None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(sent))]

    [p.join() for p in proc]

    return [x for i, x in sorted(res)]
