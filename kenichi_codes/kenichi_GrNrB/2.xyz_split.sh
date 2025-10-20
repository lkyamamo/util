#!/bin/bash

num_atoms=$( head -1 ./calculation.xyz )
line_skip=$(( $num_atoms + 2))


mkdir pos_vel
cp calculation.xyz ./pos_vel/calculation_copy.xyz
cd pos_vel
split -l "$line_skip" -a 6 -d --additional-suffix=.xyz calculation_copy.xyz
rm calculation_copy.xyz
