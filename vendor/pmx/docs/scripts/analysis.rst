.. _script_analyse:

Analyse
-------

The ``pmx analyse`` script can be used to analyse the xvg output from Gromacs.
Effectively, it allows to perform all the steps described in the :ref:`API example <example_analysis>`
and more using a convenient command line tool. ::

    $ pmx analyse -h
    usage: pmx [-h] -fA dgdl [dgdl ...] -fB dgdl [dgdl ...]
       [-m method [method ...]] [-t temperature] [-o result file]
       [-b nboots] [-n nblocks] [--integ_only] [-iA work input]
       [-iB work input] [-oA work output] [-oB work output] [--reverseB]
       [--skip  | --slice   | --rand  | --index  [...]] [--prec] [--units]
       [--pickle] [--no_ks] [--cgi_plot] [--nbins] [--dpi] [--quiet]

    Calculates free energies from fast growth thermodynamic integration
    simulations. Available methods for free energy estimation: Crooks Gaussian
    Intersection (CGI); Benett Acceptance Ratio (BAR); Jarzinski equality (JARZ).

    optional arguments:
    -h, --help            show this help message and exit
    -fA dgdl [dgdl ...]   dgdl.xvg files for the A->B simulations. Use wildcard
                          to select multiple xvg files: e.g. "-fa
                          ./forward_results/dgdl*.xvg"
    -fB dgdl [dgdl ...]   dgdl.xvg files for the B->A simulations Use wildcard
                          to select multiple xvg files: e.g. "-fb
                          ./backward_results/dgdl*.xvg"
    -m method [method ...]
                          Choose one or more estimators to use from the
                          available ones: CGI, BAR, JARZ. Default is all.
    -t temperature        Temperature in Kelvin. Default is 298.15.
    -o result file        Filename of output result file. Default is
                          "results.txt."
    -b nboots             Number of bootstrap samples to use for the bootstrap
                          estimate of the standard errors. Default is 0 (no
                          bootstrap).
    -n nblocks            Number of blocks to divide the data into for an
                          estimate of the standard error. You can use this when
                          multiple independent equilibrium simulationshave been
                          run so to estimate the error from the repeats. Default
                          is 1 (i.e. no repeats). It assumes the dgdl files for
                          each repeat are read in order and are contiguous, e.g.
                          dgdl_0 to dgdl_9 is the first repeat, dgdl_10 to
                          dgdl_19 is the second one, etc.
    -w plot               Name of image file showing the distribution of work
                          values. Default is "wplot.png". If you want to avoid
                          saving this plot, pass "none" to this flag. If you
                          choose to calculate the free energy with multiple
                          estimators, the dG values shown on the plot will be
                          chosen following the hierarchy BAR > CGI > JARZ.
    --nbins               Number of histograms bins for the plot. Default is 20.
    --dpi                 Resolution of the plot. Default is 300.
    --reverseB            Whether to reverse the work values for the backward
                          (B->A) transformation. This is useful when in Gromacs
                          both forward and reverse simulations were run from
                          lambda zero to one. Default is False.
    --integ_only          Whether to do integration only; the integrated values
                          are computed and saved, and the program terminated.
                          Default is False.
    -iA work input        Two-column dat file containing the list of input files
                          and their respective integrated work values for the
                          forward (A->B) tranformation.
    -iB work input        Two-column dat file containing the list of input files
                          and their respective integrated work values for the
                          reverse (B->A) tranformation.
    -oA work output       File where to save the list of input dgdl files and
                          their respective integrated work values for the
                          forward (A->B) tranformation. Default is "integA.dat"
    -oB work output       File where to save the list of input dgdl files and
                          their respective integrated work values for the
                          reverse (B->A) tranformation. Default is "integB.dat"
    --skip                Skip files, i.e. pick every nth work value. Default is
                          1 (all); with 2, every other work value is discarded,
                          etc.
    --slice               Subset of trajectories to analyze. Provide list slice,
                          e.g. "10 50" will result in selecting
                          dgdl_files[10:50]. Default is all.
    --rand                Take a random subset of trajectories. Default is None
                          (do not take random subset)
    --index  [ ...]       Zero-based index of files to analyze (e.g. 0 10 20 50
                          60). It keeps the dgdl.xvg files according to their
                          position in the list, sorted according to the
                          filenames. Default is None (i.e. all dgdl are used).
    --prec                The decimal precision of the screen/file output.
                          Default is 2.
    --units               The units of the output. Choose from "kJ", "kcal",
                          "kT". Default is "kJ."
    --pickle              Whether to save the free energy results from the
                          estimators in pickled files. Default is False.
    --no_ks               Whether to do a Kolmogorov-Smirnov test to check
                          whether the Gaussian assumption for CGI holds. Default
                          is True; this flag turns it to False.
    --quiet               Minimal screen output.
