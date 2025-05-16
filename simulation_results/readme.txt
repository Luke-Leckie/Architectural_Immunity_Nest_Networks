Simulation results read me

- Each file represents a simulation over 24-hours in 6-day nest network
- Filenames indicate the simulation parameter selection:
	- The number after "tun" indicates the denominator of the probability of transmission in tunnels. This has values at 1, 2, and 4.
	-The number after "ant_" indicates the number of agents in the simulation (1 or 20)
	-The final number, before ".csv", indicates the strength of self-isolation (0=no; 1=weak; 1.5=medium; 2=strong)
-Each file records more timesteps than are used for statistics. To avoid oversampling the data and to ensure we evenly sample to timesteps, we downsample this to be every 800 timesteps in our statistics scripts.
-In the files with self-isolation, self-isolation is performed in both the pathogen and control condition
