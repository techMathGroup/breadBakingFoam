########################################################################
# Instalation script for bread baking simulations

# Before run change the following variable to your path to solids4foam
# and in applications/solvers/breadBakingFoam/Make/options
pathToSolids4Foam="../solids4foam/"

########################################################################

# copy files to solids4foam 
rsync -r solids4foamAddOns/src $pathToSolids4Foam

# compile solids4foam
echo "Compiling solids4foam"
prevPath=$(pwd)
cd $pathToSolids4Foam
./Allwmake
cd $prevPath

# compile breadBakingFoam
echo "Compiling breadBakingFoam"
cd applications/solvers/breadBakingFoam
wmake
cd $prevPath

# compile boundary conditions
echo "Compiling boundary conditions"
cd src 
wmake -all 
cd $prevPath

########################################################################
