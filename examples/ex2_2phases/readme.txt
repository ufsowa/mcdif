This file contains description of MonteCarlo simulation.

################
Simulated is thermodynamic equlibrium of vacancies in NiAl intermetallic.
Equlibrium between two phases occures when chemical potentials
of the constituents are equals.
Thus, to find equlibrium vacancy concentration, it is needed to find
such a point in chemical potentials space which produces two
different phases.
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
Before starting simulations prepare samples.
Copy sample directory for every point in samples.txt file
and replace dedicated chemical potentials in chem.in.

Then run the simulations for all samples by the script:
./run_sim

View the structures during the simulations:
.\raswin.exe -script view.spt


To plot no. of atoms, order parameter, energy run:
gnuplot plot_e.gp
gnuplot plot_n.gp
gnuplot plot_eta.gp


################
Simulation files are with '*.in' extension.
conf.in - defines simulation alghortim and observability
structure.in - defines details of the used structure
energy.in - defines used energy model.
chem.in - defines chemical potentials of elements
samples.txt - contains sampel points of chemical potentials space
