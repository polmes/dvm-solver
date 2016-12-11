# Discrete Vortex Method
Just another 4-digit NACA airfoil MATLAB solver

## Basic usage
```matlab
dvm % the script will ask for your input
```

## If you are feeling pro, however, or want to iterate over the script
```matlab
[Cl,Cm_LE,Cm_AC] = dvm('4412',5,500,'c','y',8,15) 
% The arguments are, in order of appearance:
% - 4-digit NACA (string)
% - Angle of attack (º)
% - Number of panels (integer)
% - Type of geometry discretization (a. Uniform, b. Full cosine, c. Optimal)
% - Include a flap? ['y'/'n']
% (-) Flap hinge position
% (-) Flap deflection angle
```
