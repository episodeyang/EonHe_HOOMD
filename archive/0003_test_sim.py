# Now the real hoomd stuff
from hoomd_script import *
from hoomd_plugins import evaluators_ext_conf
from hoomd_plugins import evaluators_ext_ge
from pylab import *

# importing dataCache and other helpful stuff.
from ehe_simulation import eheSimulation
eheSim = eheSimulation(__file__, newFile=True)

#initialize systems, within a box with edge L
n = 1000
system = init.create_random(N=n, name='A', min_dist=0.001, box=data.boxdim(dimensions=2, L=50))
print globals.options

# specify interactions between particle pairs
def pwrlaw(r, rmin, rmax, k, mpower):
    """return the potential and the force when r apart, between rmin and rmax, using k as the constant, w.r.t. power mpower"""
    V = k * (r ** mpower)
    F = -mpower * k * (r ** (mpower - 1))
    return V, F

# now generate the table of interaction strength
table = pair.table(width=1000000)
table.pair_coeff.set('A', 'A', func=pwrlaw, rmin=.0005, rmax=19, coeff=dict(k=0.0014, mpower=-1))

#specify external potential
def eval_select(evalString):
    global confinement
    if evalString == 'linear':
        confinement = evaluators_ext_ge.external.ge()
    elif evalString == 'cylindrical':
        confinement = evaluators_ext_conf.external.confinement()
eval_select('linear')
confinement.force_coeff.set('A', gamma=0.018, npower=2)

all = group.all()
# specify integration mode and timesteps dt in real time units
integrate.mode_standard(dt=0.005)

file_prefix = "system_0003"
xml = dump.xml(filename=file_prefix + '.xml', vis=True)

import subprocess, os
# print '=============================================='
dirname = os.path.dirname(os.path.realpath(__file__))
print dirname
try:
    subprocess.call(['/bin/rm', dirname + '/' + file_prefix + '.dcd'])
except:
    pass
dump.dcd(filename=file_prefix + ".dcd", period=100)

# run brownian dynamics integration method with a step distance limit
# for a bit to get rid of overlapping particles
bdnvt = integrate.bdnvt(group=all, limit=2, T=.01)
analyzer = analyze.imd(port=55000, period=300)
run(5000)
bdnvt.disable() # you have to stop the integration afterward.
analyzer.disable()

#fire energy minimization
# analyze.log(quantities=['potential_energy'], period=1000, filename=prefix+'.log')
# fire = integrate.mode_minimize_fire(all, dt=0.01, ftol=1e-20, Etol=1e-20)
# analyzer.enable()
# run(5000)
# analyzer.disable()

class step_builder():
    def __init__(self, T0):
        self.steps = [(0, T0)]
    def add(self, step, T):
        self.steps.append((self.last() + step, T))
    def last(self):
        return self.steps[-1][0]

steps = step_builder(0.1)
steps.add(10000, 0.1)
steps.add(10000, 0.01)
steps.add(10000, 0.001)
steps.add(30000, 0.0001)
steps.add(30000, 0.00001)
steps.add(10000, 0.000001)
steps.add(10000, 0.0000001)

#reenable normal integration mode
integrate.mode_standard(dt=0.005)
nvt = integrate.nvt(group=all, T=variant.linear_interp(steps.steps), tau=0.5)
analyzer.enable()
run(steps.last())
analyzer.disable()
nvt.disable()

p = system.particles
xys = array([pp.position for pp in p])

print globals.options
eheSim.cache.new_stack()
eheSim.cache.post('xys', xys)
# scatter(xys[:, 0], xys[:, 1])
# gca().set_aspect('equal')
# show()
