from hoomd_script import *
from pylab import *
from hoomd_plugins import evaluators_ext_conf

n = 2000
system = init.create_random(N=n, name='A', box=data.boxdim(dimensions=2, L=50))

# specify interactions between particle pairs
def pwrlaw(r, rmin, rmax, k, mpower):
    V = k * (r ** mpower)
    F = -mpower * k * (r ** (mpower - 1))
    return V, F

table = pair.table(width=1000)
table.pair_coeff.set('A', 'A', func=pwrlaw, rmin=.05, rmax=10, coeff=dict(k=0.0014, mpower=-1))

confinement = evaluators_ext_conf.external.confinement()
confinement.force_coeff.set('A', gamma=0.018, npower=2)

all = group.all()
integrate.mode_standard(dt=0.005)
xml = dump.xml(filename='cylindrical_parabolic.xml', vis=True)

# dump.dcd(filename="cylindrical_parabolic.dcd", period=100)
bdnvt = integrate.bdnvt(group=all, limit=2, T=.1)
analyzer = analyze.imd(port=54321, period=300)
run(5000)
bdnvt.disable()
analyzer.disable()

#analyze.log(quantities=['potential_energy'],period=1000,filename='logtest3.log')
#fire=integrate.mode_minimize_fire(all,dt=0.08,ftol=1e-20,Etol=1e-20)   
#run(5000000)

integrate.mode_standard(dt=0.005)
nvt = integrate.nvt(group=all, T=variant.linear_interp([(0, .01), (10000, .001), (30000, .0001), (40000, .00001), (50000, 0.000001), (60000, .0000001)]), tau=0.5)
analyzer.enable()
run(60000)
analyzer.disable()
nvt.disable()

p = system.particles

dump.dcd(filename="test_ge.dcd", period=100)
xys = array([pp.position for pp in p])
scatter(xys[:, 0], xys[:, 1])
gca().set_aspect('equal')
show()
