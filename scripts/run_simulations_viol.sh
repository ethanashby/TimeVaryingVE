#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eashby@uw.edu
ml fhR
Rscript ./simulations_viol.R