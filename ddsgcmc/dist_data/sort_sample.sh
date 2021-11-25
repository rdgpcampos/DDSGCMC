#!/bin/bash

input=../energy_sampling/log.mc

awk '{if($1>0 && $11==1) {fname=sprintf("sample%dto%d.dat",$9,$10);print $6 >> fname; fname=sprintf("sample%dto%d.dat",$10,$9);print $7 >> fname;}}' $input
