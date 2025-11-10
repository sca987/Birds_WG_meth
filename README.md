
Chestnut-Crowned Babbler Methylation Analysis on the Lawrence Supercomputer

This repository contains the R scripts and documentation for filtering and modeling DNA methylation data in chestnut-crowned babbler (Pomatostomus ruficeps) fledglings.

The pipeline processes site-level methylation calls into a final, model-ready dataset, runs per-site beta regressions using GLMMs, and identifies significant CpG sites based on false discovery rate (FDR) and effect size (Δ% ≥ 25) thresholds.

All analyses were done on the Lawrence Supercomputer using a high-memory (himem) node.

Logging into the Lawrence Supercomputer

First, log in with your university credentials:

ssh your_username@login01.lawrence.usd.edu


Once logged in, you’ll be in your home directory.

Opening a High-Memory Interactive Node

To test or run your code interactively, open a himem node:

srun -p himem --mem=512G --cpus-per-task=48 --time=08:00:00 --pty bash

 What this means

srun — opens the node interactively

-p himem — selects the high-memory partition

--mem=512G — requests 512 GB of memory (up to 1.5 TB available)

--cpus-per-task=48 — number of CPUs (maximum = 48)

--time=08:00:00 — runtime limit (8 hours here; use 24 or 48 hours for long jobs)

--pty bash — keeps the session interactive

️ Things to keep in mind

The node will close if you close your laptop or lose connection.

It can also time out, so disable sleep mode on your computer.

Always save progress regularly and push results to shared storage or GitHub.

️ Basic Tools for the Lawrence Supercomputer
Action	Command / Tip
Copy	                  Highlight text (left-click)
Paste	                  Right-click (Ctrl + V will not work)
Cancel running job	    Ctrl + C
Cut a line in a script	Ctrl + K
View current directory	pwd
List files	            ls -lh
Change directories	    cd folder_name
Go up one directory	    cd ..

 Opening R and Fixing Module Errors

Once inside your working directory (for example: cd liebl_lab/shared/susan), load R:

module load R
R


If you get errors when loading R, try:

source /etc/profile.d/modules.sh
module load R


If you still get an error (often about gcc), load it manually first:

module load gcc/13.1.0-gcc-8.5.0
module load R
R


Now you are in the R environment and can run your scripts.
️ Note: You cannot edit lines interactively in R here — once you press Enter, that command runs and you cannot go back.

 Creating a SLURM Job File

Each job should include two files:

The R script you want to run.

A SLURM job script that tells the supercomputer how to run it.

Example SLURM job file:

#!/bin/bash
#SBATCH --job-name=glmm_beta
#SBATCH --partition=himem
#SBATCH --cpus-per-task=48
#SBATCH --mem=512G
#SBATCH --time=24:00:00
#SBATCH --output=glmm_beta_%j.log

module load R
Rscript 02_glmm_beta_parallel.R

Explanation

#!/bin/bash — must be the very first line

--job-name — name of your job (appears in the queue)

--partition — which node type to use (himem, short, etc.)

--cpus-per-task — number of CPUs requested

--mem — memory required (always over-estimate)

--time — time limit (also over-estimate to avoid timeouts)

--output — where the log file will be saved

You can view this log using:

cat glmm_beta_123456.log

 Submitting and Monitoring Jobs

Submit the job:

sbatch your_job_file.slurm


Check your job’s progress:

squeue -u your_username


View jobs by partition:

squeue -p himem

 Final Tips

Always confirm you’re in the correct directory before running scripts.

Use interactive nodes for testing, batch jobs for full analyses.

Save your results frequently.

When in doubt, over-estimate memory and time in SLURM.








