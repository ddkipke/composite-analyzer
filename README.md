# composite-analyzer

A matlab script to compute mechanical properties of a theoretical composite material

Run CompositeMaterialAnalyzer.m with MaterialData.xlsx in the same directory.
The program will offer a list of matrix material classes, and then individual matrix materials after you select a class.
Then it will do the same for fiber materials. Optionally, you can create your own material if you know
the modulus of elasticity, tensile strength, and density.

You then choose the volume fraction of the fiber, and whether the fibers are continuous.

The program will output the longitudinal and transvese moduli of elasticity, the longitudinal strength, 
and the specific strength of the composite.

It also displays graphs of stress-strain and modulus vs volume fraction of fiber.
