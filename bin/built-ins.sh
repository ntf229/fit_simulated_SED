#!/bin/bash

name="test2"

mainPath="$HOME/fit_simulated_SED"
projectPath="$HOME/fit_simulated_SED/projects/$name" # must contain config.txt (see ski_file_reader.py)

# make directories if they don't already exist
mkdir -p $projectPath/SKIRT_files
mkdir -p $projectPath/Prospector_files

# move ski file to current project
cp $mainPath/ski_files/built-ins.ski $projectPath/SKIRT_files/$name.ski

# change values in newly created .ski file 
python $mainPath/python/fit_simulated_sed/ski_file_reader.py --filename=$name.ski --projectPath=$projectPath


cd $projectPath/SKIRT_files # need to cd here for outputs to be in right place
skirt $name.ski 
python -m pts.do plot_seds . 
cd -

mv $projectPath/SKIRT_files/wave.npy $projectPath/Prospector_files/wave.npy
mv $projectPath/SKIRT_files/spec.npy $projectPath/Prospector_files/spec.npy

if [ -f $projectPath/Prospector_files/fit.h5 ] ; then
    rm $projectPath/Prospector_files/fit.h5
    echo removed fit file 
fi

python $mainPath/python/fit_simulated_sed/params.py --emcee --outfile=fit --path=$projectPath
python $mainPath/python/fit_simulated_sed/compare.py --filename=fit.h5 --path=$projectPath

