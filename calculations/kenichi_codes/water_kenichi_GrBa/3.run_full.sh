#!/bin/bash

cd h--o--h
sbatch 3.GrNrBaSq.slurm
cd ..

sleep 3

cd o-h--o
sbatch 3.GrNrBaSq.slurm
cd ..

sleep 3

cd h-o-h
sbatch 3.GrNrBaSq.slurm
cd ..


