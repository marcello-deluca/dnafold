## description
# this script runs dnafold for dsDNA with 384nb, then calculates the persistence length
# be sure to use "make" to create executable that stores two folders above the location of this script
# requires cadnano file, input file, and persistence_length.py in the same folder as this script

set -e
echo "Running simulation..."
rm -rf report.out
touch report.out
pwd >> report.out
echo "Job started: $(date)" >> report.out

printf "\n------ job file ------\n" >> report.out
cat job.sh >> report.out
printf "\n------ input file ------\n" >> report.out
cat input.txt >> report.out

printf "\n------ running ------\n" >> report.out
START_EPOCH=$(date +%s)
rm -rf output; mkdir output
../../dnafold input.txt >> report.out 2>&1
mv complements.dat simulation_metadata.dat output.dat* output

END_EPOCH=$(date +%s)
echo "Job completed: $(date) - $(($END_EPOCH - $START_EPOCH)) seconds elapsed" >> report.out

printf "\n------ calculate lp ------\n" >> report.out
python3 persistence_length_v5.py  --bases_skip 8 --steps_skip 10 >> report.out
