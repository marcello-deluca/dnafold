# origami_assembly_code

Software for simulating DNA origami self-assembly using overdamped Langevin (Brownian) dynamics and the dnafold model introduced in DeLuca et al., 2023 (currently deposited in biorxiv). 

Uses nlohmann's JSON parsing package from https://github.com/nlohmann/json (included as a subfolder for your convenience)

To build this software, modify the "parameters.hpp" file to the correct number of scaffold and staple beads, timestep and temperature, annealing, etc., then run "make" in the directory. The software will compile and produce an executable called "run_program" one folder back ("../run_program").

Arguments:

./run_program <json filename> <output filename>



