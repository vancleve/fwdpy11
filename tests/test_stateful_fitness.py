#We use cppimport to build the module.
#For the sake of unit testing, we force
#a rebuild every time, but that clearly
#wouldn't be needed for normal use.
import cppimport
cppimport.force_rebuild()
cppimport.set_quiet(False)
snowdrift = cppimport.imp("snowdrift")
import unittest
import numpy as np
import fwdpy11 as fp11
import fwdpy11.temporal_samplers as fp11ts
import fwdpy11.wright_fisher
import matplotlib.pyplot as plt

def evolve_snowdrift(args):
    """
    We write the function taking a tuple 
    out of habit, simplifying later
    integration with multiprocessing or
    concurrent.futures.
    """
    N, ngens, pheno0, b1, b2, c1, c2, seed = args
    #Construct as single-deme object
    #with N diploids
    pop = fp11.Spop(N)
    #Initialize a random number generator
    rng = fp11.GSLrng(seed)
    sregions = [fp11.GaussianS(0, 1, 1, 0.01, 1.0)]
    recregions = [fp11.Region(0, 1, 1)]
    fitness = snowdrift.SpopSnowdrift(rng, b1, b2, c1, c2, 1, pheno0)

    nlist = np.array([N]*ngens,dtype=np.uint32)
    recorder = fp11ts.RecordNothing()
    fp11.wright_fisher.evolve_regions_sampler_fitness(
        rng, pop, nlist,
        0.0, 0.05,
        0.0,
        [], sregions, recregions, fitness, recorder, prune_all_fixations=False)
    #return our pop
    return pop, fitness

def plotPhenoHists(N, x0, b1, b2, c1, c2, seed = 42, times = [10, 100, 500, 1000, 2000]):
    # grayscale palette
    grays = [(x,x,x) for x in np.flip(np.linspace(0.1,1,len(times)),0)]

    # run sims
    ss = [evolve_snowdrift((N, i, x0, b1, b2, c1, c2, seed))[1] for i in times]

    # plot sims
    [plt.hist(np.array(ss[i].phenotypes), color=grays[i], alpha=0.5, histtype='step', fill=True)
         for i in range(len(ss))]
    plt.show()


# ### Test Fig 1A in Doebeli and Hauert (2004, Science)
# ### The phenotype distirbution starts at 0.1 and should evolve towards 0.6
# ### and then branch due to disruptive selection

# plotPhenoHists(500, 0.1, 6, -1.4, 4.56, -1.6, 42)

# ### Test Fig 1B in Doebeli and Hauert (2004, Science)
# ### The phenotype distirbution starts at 0.1 and should evolve towards 0.6

# plotPhenoHists(500, 0.1, 7, -1.5, 4.6, -1, 42)

# ### Test Fig 1D in Doebeli and Hauert (2004, Science)
# ### The phenotype distirbution starts at 0.9 and should evolve towards 0.0

# plotPhenoHists(500, 0.9, 7, -1.5, 8, -1, 42)

# ### Test Fig 1E in Doebeli and Hauert (2004, Science)
# ### The phenotype distirbution starts at 0.1 and should evolve towards 1.0

# plotPhenoHists(500, 0.1, 7, -1.5, 2, -1, 42)

def evolve_snowdrift_dist(ngens, N, x0, mu, sig, slope, b1, b2, c1, c2, grain = 1, seed = 42):

    pop = fp11.Spop(N)
    #Initialize a random number generator
    rng = fp11.GSLrng(seed)
    sregions = [fp11.GaussianS(0, 1, 1, sig, 1.0)]
    recregions = [fp11.Region(0, 1, 1)]
    fitness = snowdrift.SpopSnowdrift(rng, b1, b2, c1, c2, slope, x0)

    recorder = fp11ts.RecordNothing()
    phenos = list()
    
    nsteps = int(ngens / grain)
    nlist = np.array([N]*grain,dtype=np.uint32)
    for step in range(nsteps):
        fp11.wright_fisher.evolve_regions_sampler_fitness(
            rng, pop, nlist,
            0.0, mu,
            0.0,
            [], sregions, recregions, fitness, recorder, prune_all_fixations=False)

        # must append as 'list' otherwise weird buffer problem
        phenos.append(list(fitness.phenotypes))
        
    #return our pop
    return pop, fitness, phenos
    
def plotEvolveDist(vals, bins=100, interpolation='bicubic', vmin=None, vmax=None, tmult=1, phenorange=[0,1]):
    parr = np.array(vals)
    nruns = parr.shape[0]

    phist = np.zeros((nruns,bins))
    for i in range(nruns):
        ph, xe = np.histogram(parr[i].ravel(), bins=bins, range=phenorange)
        phist[i] = ph

    plt.imshow(phist, aspect='auto', interpolation=interpolation,
               origin='lower', cmap=plt.cm.binary, extent=(min(xe),max(xe),0,nruns),
               vmin=vmin, vmax=vmax)
    plt.yticks(np.around(np.linspace(0,nruns,11),-1), np.around(np.linspace(0,(nruns-1)*tmult+1,11),-1).astype('int'))
    plt.show()
