# DNAfold

Access our user guide [here](https://daniel-duke.github.io/DNAfold-docs/).

See our publication in Nature Communications [here](https://doi.org/10.1038/s41467-024-46998-y).

## Installation

Follow the steps below to install DNAfold:

- Ensure you are operating on either WSL, macOS, or a Linux OS (we tested Ubuntu and CentOS).

- Ensure you have the following installed:
	- Git
	- GNU Make
	- C++ Compiler

- Navigate to where you want to store DNAfold:

	`cd <PATH_TO_DNAFOLD>`

- Clone the repository:

	`git clone https://github.com/daniel-duke/DNAfold`

- Enter the clone, then the source directory:

	`cd DNAfold/src/`

- Build the software:

	`make`

This should create an executable in the root of the cloned repository, assuming compilation did not fail.

## Simulation

To get a feel for how the code works, feel free to run one of the examples. For instance, to simulate the folding of a smiley face, follow the steps below:

- Create and enter a simulation directory:

	`mkdir <MYDIR>; cd <MYDIR>`

- Copy the caDNAno file from the `examples` folder:

	`cp <PATH_TO_DNAFOLD>/examples/smile/smile.json .`

- Copy the input file from the `examples` folder:

	`cp <PATH_TO_DNAFOLD>/examples/smile/input.txt .`

- If desired, adjust the parameters in the input file.

- Run DNAfold:

	`<PATH_TO_DNAFOLD>/dnafold input.txt`

- To run the program in the background and save the standard output to a file, use the following command:

	`nohup <PATH_TO_DNAFOLD>/dnafold input.txt > report.out 2>&1 &`

The program will automatically parse the caDNAno file referenced in the input file and run a simulation of its self-assembly, in either a forced or unforced manner depending upon the setting in the input file.

## Visualization

To visualize the simulated folding, we recommend using [OVITO](https://ovito.org).

For rudimentary visualization (showing just bead positions):
- Load the output `.dat` file in OVITO (e.g., `smile.dat`)

For fancy visualization (showing backbone bonds):
- Open the topology file in OVITO (e.g. `smile.dat_TOPOLOGY.dat`), selecting "bond" or "molecular" style when prompted.
- One the right side of the window, under the "Visual elements" pane, select "Particles" and change the standard radius to 1.4 (also, ensure the radius scaling is set to 100%).
- Also under the "Visual elements" pane, select "Bonds" and increase the bond width to 2.8. If flat shading is enabled, turn it off.
- Under the "Add modification..." (right above "Visual elements"), select "Load trajectory" (under the "Modification" header).
- Navigate down to the panel titled "Trajectory Source: External file" and click the folder icon, then select the `.dat` file containing the trajectory (e.g. `smile.dat`).
- Voila! A beautiful visualization of DNA origami self-assembly.

## Examples

`donut`: a moderate 2D structure (336 beads / 2688 nucleotides) which self-assembles into a rectangle with a hole in the middle, unforced.

`smile`: a large 2D structure (1456 beads / 11,648 nucleotides) in the shape of a smiley face. Often gets topologically trapped with forced folding, although the resemblance to a smiley face sometimes still emerges.

`4HB`: a small 3D structure (64 bead, 512 nucleotides) based on the 4HB design studied in our manuscript. Assembles unforced within 5 minutes on most computers, with few, if any, defects.

`32HB`: a large 3D structure (1024 beads, 8192 nucleotides) identical to the 32HB design studied in our manuscript. Here, it assembles by force, leading to extensive defects and poor crystallinity of the product.

Note that in the case of the two smaller designs (`donut` and `4HB`), the number of simulation steps in the input file is much longer than the estimated folding time. We recommend checking the status of the folding using OVITO then quitting the simulation (`ctrl-c`, or `pkill dnafold` to kill all local DNAfold simulations) when it has reached a satisfactory state.
