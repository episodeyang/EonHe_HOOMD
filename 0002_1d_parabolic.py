from hoomd_script import *
from pylab import *
from hoomd_plugins import evaluators_ext_ge

#initialize system
n = 1000
system = init.create_random(N=n, name='A', box=data.boxdim(dimensions=2, L=40))

# specify interactions between particle pairs
def pwrlaw(r, rmin, rmax, k, mpower):
    V = k * (r ** mpower);
    F = -mpower * k * (r ** (mpower - 1));
    return (V, F)


table = pair.table(width=1000);
table.pair_coeff.set('A', 'A', func=pwrlaw, rmin=.05, rmax=3, coeff=dict(k=2, mpower=-1))

#specify external potential
confinement = evaluators_ext_ge.external.ge()
confinement.force_coeff.set('A', gamma=1, npower=2)

all = group.all()

#specify integration mode and timesteps dt in real time units
integrate.mode_standard(dt=0.005)

#run brownian dynamics integration method with a step distance limit
# for a bit to get rid of overlapping particles
dump.dcd(filename="test_ge5", period=100)
bdnvt = integrate.bdnvt(group=all, limit=2, T=.01)
run(50000)
bdnvt.disable()

#fire energy minimization
#analyze.log(quantities=['potential_energy'],period=1000,filename='logtest3.log')
#fire=integrate.mode_minimize_fire(all,dt=0.01,ftol=1e-20,Etol=1e-20)   
#run(5000)

#reenable normal integration mode
integrate.mode_standard(dt=0.005)
nvt = integrate.nvt(group=all, T=variant.linear_interp([(0, .01), (1000000, .001)]), tau=0.5)
run(100000)
nvt.disable()

p = system.particles

X = array([pp.position for pp in p])
scatter(X[:, 0], X[:, 1])
gca().set_aspect('equal')
show()
