# mcdif
Monte Carlo based alghoritms used to simulate thermodynamic and kinetic properties in the alloys

author: Piotr Sowa
e-mail: piotrsowa87@gmail.com

NOTE:
In development.

---------------------------------

Dependancies:
GNUtoolkit: make, gcc

To compile:
make clean
make

To run:
./mcdif_{version} conf.in

  Definition of simulation is stored in files labeled *.in
  conf.in      - declaration of simulation parameters
  structure.in - definition of a crystall structure
  energy.in    - definition of energy interaction (Ising type)  
  barriers.in  - definition of energy barriers (Transition-State-Theory)
  chem.in      - input chemical potentials (Gran Cannonical Ensamble)
  stech_curve  - equlibrium vacancy concentration (for EQULIBRIUM option)

Output:
*.dat
*.xyz
control_output.dat - log file 

