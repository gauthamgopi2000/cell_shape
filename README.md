This is the source code for the deformable particle model (both adaptive and fixed $p_0$ versions).
The code is written in C++. To create the executable, run:

_g++ -o main main.cpp_

To run the executable in the same directory:

./main $\epsilon_a$ $\epsilon_b$ $\epsilon_l$ $\beta$ $Nt$ $f_0$ $\tau$ $\tau_p$

The inputs to the script are the following parameters:

$\epsilon_a$ = cell area conservation energy

$\epsilon_b$ = cell membrane bending energy

$\epsilon_l$ = cell perimeter elastic energy

$\beta$ = drag constant

$Nt$ = number of timesteps for simulation

$f_0$ = active brownian force

$\tau$ = persistence time for active brownian force

$\tau_p$ = perimeter adaptation timescale

The simulation is initalized with a packing of deformable particles at $\mathcal{A}_0=1.7$.

For some starting simulation parameters, you may use the following:
$\epsilon_a=10000$ 
$\epsilon_b=0.001$
$\epsilon_l=1$
$\beta=0.4$
$Nt=100000000$
$f_0=0.0009$
$\tau=500$
$\tau_p=1000$


The simulation will output  files called "coords.csv" and "log.csv". You may copy these files to the folder "matlab_analysis_files" and run "readcsv.m." This will create a file whose name starts with "dump...." that can be read by the freely available OVITO software for visualizing the simulation. 
