dnafold: codes for simulating DNA origami folding

access our user guide at https://daniel-duke.github.io/dnafold-docs/


Quick start (bash instructions in parentheses, replace <PATH_TO_DNAFOLD> with wherever you put the files): 
1. Make sure you have git and gnu make installed, and a c++ compiler.
2. Clone this repository ("git clone https://github.com/marcello-deluca/dnafold")
3. in either WSL, a Linux OS (we have tested on Ubuntu and CentOS), or Mac OS, enter the "src" directory ("cd <PATH_TO_DNAFOLD>/src")
4. build the software by running make ("make")
5. There is now an executable in <PATH_TO_DNAFOLD>.
6. You can now run examples in the "examples" folder:
7. Navigate to the "examples" folder ("cd <PATH_TO_DNAFOLD>/examples/")
8. Navigate to one of the examples, e.g. 3D_4HB ("cd 3D_4HB")
9. Run the included input file: ("../../dnafold input.txt")
10. The program will automatically parse the design file and run a simulation of its self assembly, in either a forced or unforced manner depending upon the user setting.

Descriptions of each example:
- 3D_4HB: small 4HB structure, unguided. Assembles in under 5 minutes on most computers.
- laughing_face: a large structure (1456 beads / 11,648 nucleotides). Takes around 20 minutes. More diffusion will improve the final quality, but forced assembly can also make the structure topologically trapped.
- 2D_donut: a moderate structure which self-assembles into a rectangle with a hole in the middle by force. Assembles in a few minutes.
- 32HB: the 32HB structure that we folded in our manuscript, except assembled by force so that you can watch it fold quickly. Folds in 20-30 minutes.

To visualize simulation results using OVITO visualization software (https://ovito.org):
You have two options. 

The most rudimentary just shows particle locations:
- Just open the .dat file (e.g., in the 32HB folder, "32HB.dat")

The nicer version also shows the backbone of the ssDNA scaffold and staples:
- Open the topology file (e.g., in the 32HB folder, once the simulation is run: "32HB.dat_TOPOLOGY.dat")
- under the "add modification" selection pane on the top right, select "load file"
- navigate down to the panel titled "external file" and click the folder icon, then select the .dat file (e.g. "32HB.dat").
- You should now get a nice trajectory with all of the particle motion and with all of the backbone bonds filled in.
