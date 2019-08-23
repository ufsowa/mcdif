This file contains description of MonteCarlo simulation.

################
Simulated is self-diffusion phenomena in NiAl intermetallic by
atomic jumps to neighbouring vacancies.
Equlibrium vacancy concentration is used.
Publication:

################
The Nickiel Aluminum (NiAL) alloy with B2 phase modelled is
within simple Ising Hamiltonian approach.
Energy model assumes interraction within two coordination zones.
Kinetic of atomic motion is realised by Resident-Time-Alghoritm on
cannonical ensamble.
Energy barriers are avaraged over atoms types per coordination zones
based on the results from Nudged-Elastic-Band calculation with
EAM potentials developed by Mishin.
Publication:

################
Data to be collected are:
Internal energy of the system (dirE.dat)
Number of atoms (dirN.dat)
Snapshots of the structure (*pic.xyz)
Traces from the execution (control_file.dat)
Squareroot of dispalcement for atoms (1dR2.dat)

################
Then run the simulations for all samples by the script:
./run_sim

View the structures during the simulations:
.\raswin.exe -script view.spt


To plot no. of atoms, order parameter, energy run:
gnuplot plot_r2.gp
gnuplot plot_d.gp
gnuplot plot_e.gp
gnuplot plot_n.gp
gnuplot plot_eta.gp


################
Simulation files are with '*.in' extension.
conf.in - defines simulation alghortim and observability
structure.in - defines details of the used structure
energy.in - defines used energy model.
barriers.in - defines energy barrier for the jumps
input.in - initial equlibrium configuration 
