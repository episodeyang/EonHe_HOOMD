import os, sys, subprocess, pickle
import numpy as np
sys.path.append("/mnt/s/_Lib/python/Projects/slab")
sys.path.append("/mnt/s/_Data/140312 - EonHe M007v5 Trident/analysis")
sys.path.append("/mnt/s/_Data/140312 - EonHe M007v5 Trident/experiment_M007v5_trident")
import util as util
# import analysis_util as analysis_util
from data_cache import dataCacheProxy

sim = lambda: None;
sim.expt_path = os.path.dirname(os.path.realpath(__file__))
sim.prefix = 'simulation'
# job_file = dataCacheProxy(sim, newFile=True, stack_prefix="job_")
# print 'the job file is located at: {}'.format(job_file.path)

resVs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]
ns = [6000, 8000, 10000, 14000, 18000, 25000, 40000]
# jobs = map(lambda n: {'n': n}, ns)
# job_file.note('simulating different numbers of particles')

jobs = []
for n in ns:
    jobs_n = map(lambda resV: {'n': n, 'resV': resV, 'resV_correction': 0.7, 'boxL': 50}, resVs)
    jobs.extend(jobs_n)
jobs[0]['newFile'] = True
import pprint, sys
pprint.pprint(jobs)
# sys.exit()

# jobs = [{'boxL': 50, 'n': 4000, 'resV': 1.4, 'resV_correction': 0.7}]
# jobs[0]['newFile'] = True


#this overwrites existing jobfiles.
jobfilename = os.path.join(sim.expt_path,  'jobs.pkl')
print 'now compiling the list of jobs to file: ' + jobfilename
with open(jobfilename, 'wb') as f:
    pickle.dump(jobs, f, protocol=pickle.HIGHEST_PROTOCOL)

for i in range(len(jobs)):
    worker_script = os.path.join(sim.expt_path, 'worker_script.py')
    print "worker_script", worker_script
    subprocess.call(['/home/ge/hoomd-install/bin/hoomd', worker_script])
