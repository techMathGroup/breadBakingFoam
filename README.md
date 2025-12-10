# OpenFOAM software for bread baking simulations

About breadBakingFoam
---------------------
breadBakingFoam is a free, open-source software based on OpenFOAM (https://openfoam.com) and solids4foam (https://www.solids4foam.com) capable of simulating bread baking process. The code is developed mostly by members of the techMathGroup of the Institute of Thermomechanics of the Czech Academy of Sciences (https://www.it.cas.cz/) as a part of the Natinal Agency for Agricultural research project on optimization of the energy consupmtion during bread baking.

<img alt="aaa_growth" width="75%" src="docs/aaa_growth.webp" />

[Installation](docs/installation.md)
--------------
1. install OpenFOAM-v2312 and solids4foam-v2.1 
2. change path to solids4foam in `install.sh` and `applications/solvers/breadBakingFoam/Make/options` files
3. source OpenFOAM
4. run `install.sh` script

Tutorials
---------
1. [2D Axisymmetrical geometry acording to https://doi.org/10.1002/aic.10518](docs/tutorial:-2D-axisymmetrical-bread.md)
2. [2D Axisymmetrical geometry acording to our custom experiments](docs/tutorial:-2D-our-custom-experiment.md)
3. [3D Geometry acording to our custom experiments](docs/tutorial:-3D-custom-experiment.md)

License
-------
breadBakingFoam is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software  Foundation, either version 3 of the License, or (at your option) any later version.  See http://www.gnu.org/licenses/, for a description of the GNU General Public License terms under which you can copy the files.
