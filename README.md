# fields-chimera
This repository hosts code for visualizing 1,2, or 3 dimensional vector fields in the Chimera biomolecular software suite

Upon installation with ```pip install .```, the following scripts will be available to use:

```
fields-chimera.py
```

The primary operating script is ```fields-chimera.py```, which will generate a ```.py``` chimera session file. The following options are available:
```
-i : Input vector field file
-o : Output file name (will end in .py by default, so don't include an extension)
-ch : Output file type (either 'chimera' or 'chimerax'); Default: chimerax
```

Input vector field files must be formatted in the following style:
```
#Sampling Density: M M M; Volume: Box: N N N
#Center: x y z
#Basis Matrix:
# B11 B12 B13
# B21 B22 B23
# B31 B32 B33
x y z vector_x vector_y vector_z
...
```

Here is an example:
```
#Sample Density: 10 10 10; Volume: Box: 1.500000 1.500000 1.500000
#Center: 34.806 37.859 18.465
#Basis Matrix:
#  0.233851   0.741795  -0.628532
#  0.952246 -0.0442104   0.302115
#   0.19632  -0.669167   -0.71671
-1.5 -1.5 -1.5 0.0205084 -0.167372 -0.191598
-1.5 -1.5 -1.35 0.011021 -0.170076 -0.184816
-1.5 -1.5 -1.2 0.00245585 -0.172498 -0.179772
-1.5 -1.5 -1.05 -0.00535891 -0.174655 -0.176227
....
```

The center, sample density, and basis lines are only required if you wish to rotate and translate the field to a different basis
