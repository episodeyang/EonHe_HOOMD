from pylab import *
from hoomd_script import *
print globals.options

globals.options = {'ny': None, 'linear': False, 'notice_level': 2, 'min_cpu': False, 'shared_msg_file': None, 'ignore_display': False, 'msg_file': None, 'nx': None, 'nrank': None, 'nz': None, 'user': [], 'onelevel': False, 'gpu': None, 'gpu_error_checking': False, 'mode': 'gpu', 'autotuner_enable': True, 'autotuner_period': 50000}

# importing dataCache and other helpful stuff.
from ehe_simulation import eheSimulation
eheSim = eheSimulation(__file__, newFile=True)

eheSim.init_system(globals.options, n=1000, L=50, dim=2, min_dist=0.001)
eheSim.pick_linear_potential(gamma=0.018, npower=2)
eheSim.set_pair(k=0.0014, rmin=0.0005, rmax=19, table_width=1000000)
eheSim.setup_analyzer_and_dump(port=55000, period=300)
eheSim.finish_setup();
eheSim.run_brownian(T=0.01, runs=50000)
eheSim.set_Ts()
eheSim.run_aneal()
eheSim.save_xys()

