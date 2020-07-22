#!/bin/bash

# use ". test.sh" to run, otherwise pts not recognized 

name="test1"
nwalkers=64
numPhotonsValue="1e7"
angleValue="0 deg"

# make directories if they don't already exist
mkdir -p $name/SKIRT_files
mkdir -p $name/Prospector_files

# move resources to current project
cp Resources/unassigned.ski $name/SKIRT_files/$name.ski
cp Resources/params.py $name/Prospector_files/params.py
cp Resources/compare.py $name/Prospector_files/compare.py
cp Resources/full_rf_wavelengths.npy $name/Prospector_files/full_rf_wavelengths.npy

cd $name
cd SKIRT_files

# change values in newly created .ski file 
sed -i '' "s/numPhotonsValue/$numPhotonsValue/g" $name.ski 
sed -i '' "s/angleValue/$angleValue/g" $name.ski

skirt $name.ski 
pts plot_seds .
cd ..
mv SKIRT_files/wave.npy Prospector_files/wave.npy
mv SKIRT_files/spec.npy Prospector_files/spec.npy
cd Prospector_files
python params.py --emcee --outfile=fit --nwalkers=$nwalkers
python compare.py --filename=fit_$nwalkers.h5

# go back to starting directory 
cd ..
cd ..