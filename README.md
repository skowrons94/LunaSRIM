# LunaSRIM

A collection of scripts used to run SRIM calculations with [costrouc/pysrim][https://github.com/costrouc/pysrim] docker image.

## Repository Structure

- **scripts/**: Contains the scripts used to make SRIM calculations both for implantation (```F_in_Fe.py```), energy loss (```H_in_C.py```) and angular straggling (```He_in_N.py```).
- **run.sh**: Bash script to execute the calculation.

## Usage

Put the new or modified script in the ```scripts/``` directory, modify the ```run.sh``` script to launch your ones. After the calculation has been endend, you should see the ```output/``` directory with the results from your SRIM calculation.
