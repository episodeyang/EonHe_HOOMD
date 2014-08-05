import os, sys, subprocess
import numpy as np
sys.path.append("/mnt/s/_Lib/python/Projects/slab")
sys.path.append("/mnt/s/_Data/140312 - EonHe M007v5 Trident/analysis")
sys.path.append("/mnt/s/_Data/140312 - EonHe M007v5 Trident/experiment_M007v5_trident")
import util as util
# import analysis_util as analysis_util
from data_cache import dataCacheProxy

# Now the real hoomd stuff
# from hoomd_script import * # There are globals in hoomd_script that have to be imported as well.
from hoomd_script import *
#from hoomd_script import globals
import hoomd_script as hmd
import hoomd_plugins.evaluators_ext_conf as cylindrical_ext
import hoomd_plugins.evaluators_ext_ge as linear_ext

__author__ = 'yangg_000'

class step_builder():
    def __init__(self, T0):
        self.steps = [(0, T0)]
    def add(self, step, T):
        self.steps.append((self.last() + step, T))
    def last(self):
        return self.steps[-1][0]

class eheSimulation():
    def __init__(self, script, prefix='simulation', newFile=True):
        """script is the __file__ object os the simulation script that calls this method.
        """
        sim_path = os.path.dirname(os.path.realpath(script))
        print "=================================================="
        print sim_path
        print "=================================================="
        self.expt_path = sim_path
        if prefix != None: self.prefix = prefix
        self.cache = dataCacheProxy(self, newFile=newFile)
        self.filename = self.cache.filename
        self.dump_prefix = "system_0003"

    def bind_particles(self):
        self.all = group.all()

    def init_system(self, options, n=1000, L=50, dim=2, min_dist=0.001):
        globals.options = util.dict2obj(**options)
        #globals.options = lambda:None;
        #globals.options.mode = 'gpu'
        print "running in gpu mode by default! --Ge"
        #initialize systems, within a box with edge L
        self.system = hmd.init.create_random(N=n, name='A', min_dist=min_dist, box=hmd.data.boxdim(dimensions=dim, L=L))
        #hmd.init.create_random(N=n, name='A', min_dist=min_dist, box=hmd.data.boxdim(dimensions=dim, L=L))

    def pick_linear_potential(self, gamma=0.018, npower=2):
        self.confinement = linear_ext.external.ge()
        self.confinement.force_coeff.set('A', gamma=gamma, npower=npower)

    def pick_cylindrical_potential(self, gamma=0.018, npower=2):
        self.confinement = cylindrical_ext.external.confinement()
        self.confinement.force_coeff.set('A', gamma=gamma, npower=npower)

    def set_pair(self, k=0.0014, rmin=0.0005, rmax=19, table_width=1000000):
        # specify interactions between particle pairs

        def pwrlaw(r, rmin, rmax, k, mpower):
            """return the potential and the force when r apart, between rmin and rmax, using k as the 
            constant, w.r.t. power mpower"""
            V = k * (r ** mpower)
            F = -mpower * k * (r ** (mpower - 1))
            return V, F

        # now generate the table of interaction strength
        table = hmd.pair.table(width=table_width)
        table.pair_coeff.set('A', 'A', func=pwrlaw, rmin=rmin, rmax=rmax, coeff=dict(k=k, mpower=-1))

    def setup_integrator_analyzer_and_dump(self, dump_prefix=None, port=55000, period=300):
        hmd.dump.dcd(filename=self.dump_prefix + ".dcd", period=100)
        if dump_prefix!=None: self.dump_prefix= dump_prefix;
        # specify integration mode and timesteps dt in real time units
        self.integrator = hmd.integrate
        self.integrator.mode_standard(dt=0.005)
        xml = dump.xml(filename=self.dump_prefix + '.xml', vis=True)

        try:
            subprocess.call(['/bin/rm', self.expt_path + '/' + self.dump_prefix + '.dcd'])
        except:
            pass

        self.analyzer = hmd.analyze.imd(port=port, period=period)

    def run_brownian(self, T=0.01, runs=50000):
        # run brownian dynamics integration method with a step distance limit
        # for a bit to get rid of overlapping particles
        if not hasattr(self, 'bdnvt'):
            self.bdnvt = self.integrator.bdnvt(group=self.all, limit=2, T=T)
        hmd.run(runs)
        self.bdnvt.disable() # you have to stop the integration afterward.
        self.analyzer.disable()

    def set_Ts(self):
        Aneal = step_builder(0.1)
        Aneal.add(10000, 0.01)
        Aneal.add(10000, 0.001)
        Aneal.add(30000, 0.0001)
        Aneal.add(30000, 0.00001)
        Aneal.add(10000, 0.000001)
        Aneal.add(20000, 0.0000001)
        self.aneal = Aneal

    def run_aneal(self):
        #reenable normal integration mode
        self.integrator.mode_standard(dt=0.005)
        if not hasattr(self, 'nvt'):
            self.nvt = self.integrator.nvt(group=self.all, T=variant.linear_interp(self.aneal.steps), tau=0.5)
        self.analyzer.enable()
        hmd.run(self.aneal.last())
        self.analyzer.disable()
        self.nvt.disable()
        
    def save_xys(self):
         p = self.system.particles
         xys = np.array([pp.position for pp in p])

         self.cache.new_stack()
         self.cache.post('xys', xys)









