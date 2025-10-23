#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eashby@fredhutch.org
ml fhR
Rscript ./simulations.R