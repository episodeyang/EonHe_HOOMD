To Run parametric sweeps of resonator bias voltage, and different particle numbers, we implemented the following architecture:

- jobCompiler: 0000_simulation_name.jobCompiler.py
- worker: sim_worker.py
- ehe_simulation class: ehe_simulation.py

The jobCompiler compiles the job parameter list into a list of dictionary objects and dump it into a pickle file:

```
dump file naming format: 0000_simulation_name.jobfile.pkl
```

To run a simulation:
```
hoomd xxx.jobCompiler.py
```

or to run a simuation that keeps going even if you log out of the secure shell:
```
nohup hoomd xxx.jobCompiler.py &
```

