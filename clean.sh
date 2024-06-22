#!/bin/bash
find "Results/plots" -type f -delete
find "Results/complete" -type f -delete
find "Results/completesemi" -type f -delete
find "Results/sound" -type f -delete
find "Results/systems" -type f -delete
rm result.txt
rm result3d.txt
echo "clean results"