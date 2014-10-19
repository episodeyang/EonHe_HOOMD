from pylab import *
from hoomd_script import *
#print globals.options

import pickle, pprint
jobfilename = 'jobs.pkl'
with open(jobfilename, 'r') as f:
    jobs = pickle.load(f)
print 'Now reading simulation job file'
print "------------------------========list of jobs to simulate========-----------------------"
pprint.pprint(jobs)

import sys
sys.path.append("/mnt/s/_Lib/python/Projects/slab")
#sys.path.append("/mnt/s/_Data/140312 - EonHe M007v5 Trident/analysis")
sys.path.append("/mnt/s/_Data/140312 - EonHe M007v5 Trident/experiment_M007v5_trident")
import util as util

options = {'ny': None, 'linear': False, 'notice_level': 2, 'min_cpu': False, 'shared_msg_file': None, 'ignore_display': False, 'msg_file': None, 'nx': None, 'nrank': None, 'nz': None, 'user': [], 'onelevel': False, 'gpu': None, 'gpu_error_checking': False, 'mode': 'gpu', 'autotuner_enable': True, 'autotuner_period': 50000}

if jobs.__len__() == 0: 
    print '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
    print "job file is empty. simulation is completed!"
    print '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
    sys.exit('now exit')

current_job = jobs[0]
config = util.dict2obj(**current_job)

from ehe_simulation import eheSimulation
if hasattr(config, 'newFile'):
    eheSim = eheSimulation(__file__, newFile=config.newFile)
else:
    eheSim = eheSimulation(__file__)

k = 0.0014
current_job['k'] = k
C0 = config.resV * 0.06 * config.resV_correction
current_job['C0'] = C0

eheSim.cache.find_last_stack()
print eheSim.cache.current_stack
try:
    cached_job = eheSim.cache.get_dict('config')
    # not terribly robust way to compare but good enough for now
    if  cached_job != current_job:
        print "creating a new stack for new configuration"
        eheSim.cache.new_stack()
        eheSim.cache.set_dict('config', current_job)
except IOError as e:
    eheSim.cache.find_last_stack()
    print 'data file does not exist yet'
    print e
    eheSim.cache.set_dict('config', current_job)
    d = eheSim.cache.get_dict('config')
    print d
except KeyError as e:
    print 'current stack does not have configuration'
    print e
    eheSim.cache.set_dict('config', current_job)
    d = eheSim.cache.get_dict('config')
    print d

# You can't really run multiple simulation in a single script as of now...
eheSim.init_system(options, n=config.n, L=config.boxL, dim=2, min_dist=0.001)
eheSim.bind_particles()

eheSim.pick_linear_potential(gamma=C0, npower=2)
eheSim.set_pair(k=k, rmin=0.0005, rmax=19, table_width=1000000)
eheSim.setup_integrator_analyzer_and_dump(port=55000, period=300)
# everthing below can be run multiple times.
# One can NOT change the confinement after the script has started.w

eheSim.run_brownian(T=0.01, runs=50000)
eheSim.set_Ts()
eheSim.run_aneal()

# post is the highlevel append method with a stack pointer.
eheSim.cache.post('xys', eheSim.get_xys())
jobs = jobs[1:]
with open(jobfilename, 'wb') as f:
    pickle.dump(jobs, f, protocol=pickle.HIGHEST_PROTOCOL)

print 'simulated successfully', current_job
