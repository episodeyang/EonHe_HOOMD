from pylab import *
import argparse
from hoomd_script import *
print globals.options

parser = argparse.ArgumentParser(description='simulation script')
parser.add_argument('n', metavar='n', type=int, help='the number of particles simulated')
parser.add_argument('--resV', metavar='resV', type=int, help='the voltage of the traping potential')
parser.add_argument('--biasCorrection', dest='biasCorrection', type=float, default=0.30, help='set the ')

args = parser.parse_args()
n = args.n
resV = args.resV
biasCo = args.biasCorrection
C0 = resV * biasCo * 0.06 

print 'simulating n={} particles at resV={}V.'.format(n, resV)
print 'the bias correction factor is: {}. The gamma passed in is {} in eV/um'.format(biasCo, C0)

options = {'ny': None, 'linear': False, 'notice_level': 2, 'min_cpu': False, 'shared_msg_file': None, 'ignore_display': False, 'msg_file': None, 'nx': None, 'nrank': None, 'nz': None, 'user': [], 'onelevel': False, 'gpu': None, 'gpu_error_checking': False, 'mode': 'gpu', 'autotuner_enable': True, 'autotuner_period': 50000}

from ehe_simulation import eheSimulation

# You can't really run multiple times as of now...

eheSim = eheSimulation(__file__, newFile=True)
eheSim.init_system(options, n=1000, L=50, dim=2, min_dist=0.001)
eheSim.bind_particles();
eheSim.pick_linear_potential(gamma=C0, npower=2)
eheSim.set_pair(k=0.0014, rmin=0.0005, rmax=19, table_width=1000000)
eheSim.setup_integrator_analyzer_and_dump(port=55000, period=300)
# everthing below can be run multiple times.
# One can NOT change the confinement after the script has started.w

eheSim.run_brownian(T=0.01, runs=50000)
eheSim.set_Ts()
eheSim.run_aneal()
eheSim.run_brownian(T=0.01, runs=50000)
eheSim.save_xys()

