from hoomd_script import *
from pylab import *
from hoomd_plugins import evaluators_ext_conf

n = 100
system = init.create_random(N=n, name='A', box=data.boxdim(dimensions=2, L=200))

# specify interactions between particle pairs

def pwrlaw(r, rmin, rmax, k, mpower):
    V = k * (r ** mpower)
    F = -mpower * k * (r ** (mpower - 1))
    return V, F

table = pair.table(width=1000)
table.pair_coeff.set('A', 'A', func=pwrlaw, rmin=.05, rmax=3, coeff=dict(k=2, mpower=-1))

confinement = evaluators_ext_conf.external.confinement()
confinement.force_coeff.set('A', gamma=1, npower=2)

all = group.all()

integrate.mode_standard(dt=0.008)

dump.dcd(filename="test4", period=100)
bdnvt = integrate.bdnvt(group=all, limit=2, T=.01)
run(50000)
bdnvt.disable()

#analyze.log(quantities=['potential_energy'],period=1000,filename='logtest3.log')
#fire=integrate.mode_minimize_fire(all,dt=0.08,ftol=1e-20,Etol=1e-20)   
#run(5000000)

integrate.mode_standard(dt=0.008)
nvt = integrate.nvt(group=all, T=variant.linear_interp([(0, .004), (1000000, .0003)]), tau=0.5)
run(100000)
nvt.disable()

p = system.particles

#dump.dcd(filename="test_ge.dcd", period=100)
X = array([pp.position for pp in p])
scatter(X[:, 0], X[:, 1])
gca().set_aspect('equal')
show()
