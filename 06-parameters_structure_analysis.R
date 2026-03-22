install.packages("remotes")
remotes::install_github("labroo2/deStructure")
library(deStructure)

output_file <- "structure_analysis/mainparams.txt"
mainparams(maxpops = NULL, burnin = 50000, numreps = 50000, infile = "structure_output.txt",
           outfile = "my_output_file", numinds = "449", numloci = "12595",
           ploidy = 4, missing = -9, onerowperind = 0, label = 1,
           popdata = 0, popflag = 0, locdata = 0, phenotype = 0, extracols = 0,
           markernames = 1, recessivealleles = 0, mapdistances = 0, phased = 0,
           phaseinfo = 0, markovphase = NULL, notambiguous = NULL, outpath = output_file)

output_file <- "structure_analysis/extraparams.txt"
extraparams(noadmix = 0, linkage = 0, usepopinfo = 0, locprior = 0, 
            freqscorr = 1, onefst = 0, inferalpha = 1, popalphas = 0,
            alpha = 1.0, inferlambda = 0, popspecificlambda = NULL,
            lambda = 1.0, fpriormean = 0.01, fpriorsd = 0.05,
            unifprioralpha = 1, alphamax = 10.0, log10rmin = NULL,
            log10rmax = NULL, log10propsd = NULL, log10rstart = NULL,
            gensback = NULL, migrprior = NULL, pfrompopflagonly = 0,
            locispop = NULL, locpriorinit = NULL, maxlocprior = NULL,
            printnet = NULL, printlambda = NULL, printqsum = NULL,
            sitebysite = NULL, printqhat = NULL, updatefreq = NULL,
            printlikes = NULL, intermedsave = NULL, echodata = NULL,
            ancestdist = 0, computeprob = 1, admburnin = NULL,
            alphapropsd = 0.025, startatpopinfo = 0, randomize = NULL,
            seed = NULL, metrofreq = 10, reporthitrate = NULL, outpath = output_file)

