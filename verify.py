from fairsquare import fairProve
import argparse

def parseArgs():

    parser = argparse.ArgumentParser(
        description='''\u25B6\u25B6\u25B6 The Glorious Madison Algorithmic Fairness Prover.''')

    parser.add_argument("-f", "--file", required=True,
                        help="Input python bench")

    parser.add_argument("-e", "--epsilon", required=False, default=.1,
                        help="Epsilon fairness guarantee")

    parser.add_argument("-s", "--simulate", required=False, action='store_true', default=False,
                        help="Simulate program to check property")

    parser.add_argument("-sc", "--simTimes", required=False, default=1000,
                        help="Number of simulation times")

    parser.add_argument("-mf", "--finiteMaximize", required=False, action='store_true', default=False,
                        help="attempt to maximize the finite bounds of samples")

    parser.add_argument("-r", "--randomize", required=False, default=None,
                        help="If using finiteMaximize, randomly permutes the order of the variables to try maximizing using the given random seed")

    parser.add_argument("-mi", "--infiniteMaximize", required=False, action='store_true', default=False,
                        help="attempt to drop bounds for unbounded sides of samples")

    parser.add_argument("-p", "--plot", required=False, action='store_true', default=False,
                        help="Plots progress of sampling")

    parser.add_argument("-z", "--z3qe", required=False, action='store_true', default=False,
                        help="Use z3 quantifier elimination instead of redlog")

    parser.add_argument("-nh", "--numHists", required=False, default=5,
                        help="Number of buckets to use in histogram approximations (preferably odd)")

    parser.add_argument("-hb", "--histBound", required=False, default=3,
                        help="Sets the histogram bound to +/- variance times this number")

    parser.add_argument("-o", "--output", required=False, default=None,
                        help="File to write results")

    parser.add_argument("-t", "--timeout", required=False, default=None,
                        help="Use a timeout period in seconds after which to terminate (not including time for qelim) instead of an epsilon") 

    parser.add_argument("-a", "--adapt", required=False, action='store_true', default=False,
                        help="Dynamically increase the histogram bound for approximations when samples are small (in a rudimentary way")

    parser.add_argument("-ro", "--rotate", required=False, default='identity',
                        help="Enable rotations of the formula for improved convergence in some circumstances")

    parser.add_argument("-v", "--verbose", required=False, action='store_true', default=False,
                        help="Output a lot")

    return parser.parse_args()


def main():

    args = parseArgs()

    if args.simulate:
        simulate(args.file, args.simulate, args.simTimes)
        return

    #e = parse_fr(args.file)

    #print("\n\n== Population and program encoding == ")
    #print("Population model: ", e.model)
    #print("Program: ", e.program)

    #print("Probability dists: ", e.vdist)
    #print("fairness target", e.fairnessTarget)
    #print("sensitiveAttribute", e.sensitiveAttribute)
    #print("qualified", e.qualified)
    #print("\n\n")

    #print("ROTATE COMMAND LINE ARGUMENT:\n", args.rotate)

    timeoutarg = int(args.timeout) if args.timeout is not None else None
    randarg = int(args.randomize) if args.randomize is not None else None
    fairProve.proveFairness(args.file, args.output, float(args.epsilon), args.finiteMaximize, randarg, args.infiniteMaximize, 
            args.plot, args.z3qe, int(args.numHists), int(args.histBound), timeoutarg, args.adapt,
            args.rotate, args.verbose)

if __name__ == "__main__":
    main()
