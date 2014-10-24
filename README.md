# SlowingDown

Simple Slowing Down Monte Carlo Code to calculate correct energy
distribution of diffusion coefficients. Currently works with
hydrogen-1 with simple elastic scattering physics.

## Building the code

To build the code, first navigate the source directory:
```bash
   cd src
```

Then run the build script that executes **cmake** and **make**:
```bash
   ./build.sh
```

This will create a **build** folder inside the source directory.
The executable can be found in **build/bin** called **slowdown**.
To run the code, navigate to the examples directory from the
source directory:
```bash
  cd ../examples
```

The input file is named **input.xml** and you can see some of the
different options. They should be self explanitory for now until
more instructions are written. To run the code, execute **slowdown**
from within the examples directory:
```bash
  ../src/build/bin/slowdown
```

After the code finishes, there will be numerous ASCII output files
in the examples directory. These are the individual tallies
as a function of energy group. To process these tallies to generate
diffusion coefficients and the transport-to-total ratio curve,
run the post process script from within the examples directory:
```bash
  python ../scripts/calculate_diffusion.py
```

This will open up a lot of plots.
