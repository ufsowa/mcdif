This file contains description of MonteCarlo simulation.

################
The Iron Platinum (FePt) alloy modelled is within simple Ising Hamiltonian approach.
Energy model assumes interraction within two coordination zones.
Simulation uses direct exchange alghoritm in canonical ensamble.

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
.\raswin.exe -script view.rsml

################
Simulation files are with '*.in' extension.
conf.in - defines simulation alghortim and observability
structure.in - defines details of the used structure
energy.in - defines used energy model.

