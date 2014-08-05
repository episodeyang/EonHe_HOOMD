from hoomd_script import *
from hoomd_plugins import evaluators_ext_ge
from pylab import *

#initialize systems, within a box with edge L
n = 1000
system = init.create_random(N=n, name='A', min_dist=0.001, box=data.boxdim(dimensions=2, L=50))

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
confinement = evaluators_ext_ge.external.ge()
confinement.force_coeff.set('A', gamma=0.018, npower=2)

all = group.all()
#specify integration mode and timesteps dt in real time units
integrate.mode_standard(dt=0.005)
# run brownian dynamics integration method with a step distance limit
# for a bit to get rid of overlapping particles

file_prefix = "linear_parabolic"
xml = dump.xml(filename=file_prefix + '.xml', vis=True)
# dump.dcd(filename=file_prefix + ".dcd", period=100)

bdnvt = integrate.bdnvt(group=all, limit=2, T=.01)
analyzer = analyze.imd(port=54321, period=300)
run(5000)
bdnvt.disable() # you have to stop the integration afterward.
analyzer.disable()

#fire energy minimization
# analyze.log(quantities=['potential_energy'], period=1000, filename=prefix+'.log')
# fire = integrate.mode_minimize_fire(all, dt=0.01, ftol=1e-20, Etol=1e-20)
# analyzer.enable()
# run(5000)
# analyzer.disable()

#reenable normal integration mode
integrate.mode_standard(dt=0.005)
nvt = integrate.nvt(group=all, T=variant.linear_interp([(0, .01), (10000, .001), (30000, .0001), (40000, .00001), (50000, 0.000001), (60000, .0000001)]), tau=0.5)
analyzer.enable()
run(60000)
analyzer.disable()
nvt.disable()

p = system.particles
xys = array([pp.position for pp in p])
scatter(xys[:, 0], xys[:, 1])
gca().set_aspect('equal')
show()
