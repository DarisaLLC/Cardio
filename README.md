# Cardio
**matlab** source files for calculation of contractile force of single iPSC-derived cardiomyocytes. 

Example of Use:
    `cardiomyocyte_3 (-1, dL, L, W, 0.0002, Cs)`
    
- dL = cell contraction (cm),
- L = original cell length (cm),
- W = cell width (cm),
- Cs = shear wave velocity of gel (cm/s) = 100 * sqrt(E/1.24) **E** is the Young's modulus in kPA. 

**The calculated contractile force is shown under "Total reaction force exerted by gel"**
