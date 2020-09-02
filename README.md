# thinAirfoilPy
Python code to compute lift and moment polars for airfoils based on the Thin Airfoil Theory.

## Usage
```bash
Usage: thinAirfoil.py MPXX [alphaDeg]

$ # Uncambered airfoil lift polars
$ thinAirfoil.py 0012
  CLa  (1/rad) =  6.283185
  CL0          =  0.000000
  alf0 (deg)   = -0.000000

$ # Cambered airfoil lift polars
$ thinAirfoil.py 5605
  CLa  (1/rad) =  6.283185
  CL0          =  0.710635
  alf0 (deg)   = -6.480218

$ # Uncambered airfoil lift at 10 deg
$ thinAirfoil.py 0012 10
  Alpha       CL         CM_le      CM_c/4     xcp
  10.000000   1.096623  -0.087266  -0.000000   0.250000

$ # Cambered airfoil lift at multiple angles
$ thinAirfoil.py 5605 0 2 4
  Alpha       CL         CM_le      CM_c/4     xcp/c
  0.000000   0.710635  -0.116148  -0.187232  -0.013471
  2.000000   0.929960  -0.133602  -0.187232   0.048667
  4.000000   1.149285  -0.151055  -0.187232   0.087088
```
