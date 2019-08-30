This file contains description of MonteCarlo simulation.

################
Simulated is equlibrium structure of NiAl intermetallic
for given chemical potential and temperature.
Publication:

################
The Nickiel Aluminum (NiAL) alloy with B2 phase modelled is
within simple Ising Hamiltonian approach.
Energy model assumes interraction within two coordination zones.
Simulation uses Semi Grand Cannonical SGC Monte Carlo alghoritm
to mimic a semi-open system.

################
Data to be collected are:
Internal energy of the system (dirE.dat)
Number of atoms (dirN.dat)
Snapshots of the structure (*pic.xyz)
Traces from the execution (control_file.dat)

################
To run the simulation start the script:
./run_sim

To plot no. of atoms, order parameter, energy run:
gnuplot plot_e.gp
gnuplot plot_n.gp
gnuplot plot_eta.gp

To view the structure start:
.\raswin.exe -script view.spt

################
Simulation files are with '*.in' extension.
conf.in - defines simulation alghortim and observability
structure.in - defines details of the used structure
energy.in - defines used energy model.
chem.in - defines chemical potentials of elements
