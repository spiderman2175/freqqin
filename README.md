# freqqin
Python2 script to process Gaussian frequency output: calculate quasi-RRHO contributions and/or displace molecule along a vibrational mode vector

Usage: freqqin.py gaussian_frequency_output_file [-c|a|o|d|i|h]

The script calculates entropy contributions for a molecular geometry, using an
output file of a Gaussian frequency calculation as input. If one imaginary
frequency exists in the output file, the option is also given to create Gaussian
input files for the reactant and product geometries corresponding to a
transition state characterized by this imaginary frequency.

Parameters that can be set by the user:<br />
-s: Select to [1] calculate entropy contributions, [2] displace the molecule
    along the Cartesian displacement vectors, or [3] do both.<br />
-c: Cutoff value (in cm**-1) for which frequencies are to be treated as
    low-lying modes (defaults to 100 cm**-1)<br />
-a: Value of the alpha parameter when the vibrational entropy is calculated
    according to DOI: 10.1002/chem.201200497 (defaults to 4)<br />
-o: Rotational symmetry number of the molecule (defaults to 1)<br />
-d: Displacement factor for the creation of Gaussian input files for reactant
    and product structures formed by vibrating along the imaginary frequencies
    (defaults to 0.1). The coordinates in those input files are calculated as
   [Original coordinates] +/- [Displacement factor] * [Vibrational displacement]<br />
-i: Include the imaginary frequencies v when calculating the entropy
    (set them to -i*v)<br />
-h: Display this message
