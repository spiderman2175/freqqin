#!/usr/bin/env python
import sys, re, math, numpy, getopt, os

usage = """
Usage: %s gaussian_frequency_output_file [-c|a|o|d|i|h]

The script calculates entropy contributions for a molecular geometry, using an
output file of a Gaussian frequency calculation as input. If one imaginary
frequency exists in the output file, the option is also given to create Gaussian
input files for the reactant and product geometries corresponding to a
transition state characterized by this imaginary frequency.

Parameters that can be set by the user:
-s: Select to [1] calculate entropy contributions, [2] displace the molecule
    along the Cartesian displacement vectors, or [3] do both.
-c: Cutoff value (in cm**-1) for which frequencies are to be treated as
    low-lying modes (defaults to 100 cm**-1)
-a: Value of the alpha parameter when the vibrational entropy is calculated
    according to DOI: 10.1002/chem.201200497 (defaults to 4)
-o: Rotational symmetry number of the molecule (defaults to 1)
-d: Displacement factor for the creation of Gaussian input files for reactant
    and product structures formed by vibrating along the imaginary frequencies
    (defaults to 0.1). The coordinates in those input files are calculated as
   [Original coordinates] +/- [Displacement factor] * [Vibrational displacement]
-i: Include the imaginary frequencies v when calculating the entropy
    (set them to -i*v)
-h: Display this message
""" % sys.argv[0]

# parameters set by the user

if len( sys.argv ) < 2:
    print usage
    sys.exit()

file_exists = os.path.isfile( sys.argv[1] )
if not file_exists:
    if sys.argv[1] == '-h':
        print usage
    else:
        print """
The file %s does not exist. Aborting
""" % sys.argv[1]
    sys.exit()

freqfile = sys.argv[1]
fo = open( freqfile, 'r' )
thefile = fo.read()

w0 = 100 #cm**-1
alpha = 4
o = 1
displm_factor = 0.1
include_imag = False
freqs_frozen = False

opts, args = getopt.getopt( sys.argv[2:], "s:c:a:o:d:if:h" )
for opt, arg in opts:
    if opt == '-h':
        print usage
        sys.exit()
    elif opt == '-s':
        usr_input = int(arg)
    elif opt == '-c':
        w0 = float(arg)
    elif opt == '-a':
        alpha = float(arg)
    elif opt == '-o':
        o = float(arg)
    elif opt == '-d':
        displm_factor = float(arg)
    elif opt == '-i':
        include_imag = True
    elif opt == '-f':
        freqs_frozen = True
        freqs_frozen_file = arg
        file_exists = os.path.isfile( freqs_frozen_file )
        if not file_exists:
            print """
The file %s does not exist. Aborting
""" % freqs_frozen_file
            sys.exit()

# define physical constants in atomic units
R = 8.3144621 #J mol**-1
molJ = 1 / (2625.49962 * 1E3) #Ha/(J/mol)
RHa = R * molJ #Ha
cal = 4.184 #J
kB = 3.16681520371153E-6 #1.3806488E-23 #m**2 kg s**-2 K**-1
h = 2 * math.pi #6.62606957E-34 #m**2 kg s**-1
atm = 3.398931578276133E-14 * 101325 #Pa
u = 1822.888479031408 #1.660538921E-27 #kg
A = 1.889726133921252 #1E-10 #m
cm1 = 2.418884324306202E-17 * 1E2 * 299792458 #s**-1

w0 = w0 * cm1

# catch all the vibrational frequencies
freqs_array = []
if freqs_frozen:
    ffo = open( freqs_frozen_file, 'r' )
    the_frozen_file = ffo.read()
    freqs = re.findall( "[ ]+[0-9]+\.[0-9]+[ ]+[0-9]+\.[0-9]+.*\n", the_frozen_file )
    for row in freqs:
        split_row = row.split()
        freq = float( split_row[1] )
        if len( split_row ) > 2:
            freq = -freq
        freqs_array.append( freq ) # unit: cm-1
    ffo.close()
else:    
    freqs = re.findall( "Frequencies --[ ]+-?[0-9]+\.[0-9]+.*\n", thefile )
    for row in freqs:
        split_row = row.split()
        for i in range( 2, len( split_row ) ):
            freqs_array.append( float( split_row[i] ) ) # unit: cm-1
freqs = numpy.array( freqs_array ) * cm1

# check for imaginary frequencies
no_imag = sum( freqs < 0 )
if no_imag > 0:
    print """
The file %s contains one or more imaginary frequencies. What do you want to do?
[1] Calculate entropic contributions
[2] Generate Gaussian input files for the reactants and products formed by
    vibrating along these vectors
[3] Do both""" % freqfile
    try:
        print "Your choice: " + str( usr_input )
    except:
        usr_input = 0
        while usr_input not in [ 1, 2, 3 ]:
            usr_input = input( 'Your choice: ' )
else:
    print """
The file %s contains no imaginary frequencies, proceeds to calculate entropic
contributions...""" % freqfile
    usr_input = 1

# catch the multiplicity
mult = re.search( "Multiplicity =[ ,0-9][0-9]*", thefile )
mult = int( mult.group(0).split("=")[1] )
print mult

# catch the cartesian coordinates for the molecule
coords = re.findall( "[ ]+[0-9]{1,3}[ ]+[0-9]{1,3}[ ]+-?[0-9][ ]+(-?[0-9]+\.[0-9]+)[ ]+(-?[0-9]+\.[0-9]+)[ ]+(-?[0-9]+\.[0-9]+)\n", thefile )
coordsA = numpy.array( coords, dtype=float )

# catch the atom number and mass for each atom in the molecule
masses = re.findall( "Atom[ ]+[0-9]{1,3} has atomic number[ ]+([0-9]{1,3}) and mass[ ]+([0-9]+\.[0-9]+)\n", thefile )
masses = numpy.array( masses, dtype = [ ( 'N', int ),( 'w', float ) ] )
masses['w'] = masses['w'] * u # unit: kg

# check so that there are no duplicates
no_atoms = len( masses[ 'w' ] )
coordsA_rows = coordsA.shape[ 0 ]
if coordsA_rows > no_atoms:
    coordsA = coordsA[ range( coordsA_rows - no_atoms, coordsA_rows ) ]

coords = coordsA * A

if usr_input in [ 1, 3 ]:

    # dispose of all imaginary frequencies
    freqs_cleaned = range( len( freqs ) )
    if not include_imag:
        freqs_cleaned = [ i for i, e in enumerate( freqs ) if e > 0 ]
    freqs_cleaned = freqs[ freqs_cleaned ]
    if include_imag:
        for i, e in enumerate( freqs ):
            if e < 0:
                freqs_cleaned[ i ] = - e
    
    # catch the temperature and the pressure
    TP = re.search( "Temperature[ ]*[0-9]*\.[0-9]* Kelvin.[ ]*Pressure[ ]*[0-9]*\.[0-9]* Atm.", thefile )
    T = float( TP.group(0).split()[1] )
    P = float( TP.group(0).split()[4] ) * atm

    # calculate molecular properties
    M = numpy.sum( masses['w'] )
    Rx = 1/M * sum( numpy.multiply( masses['w'], coords[:,0] ) ) #center
    Ry = 1/M * sum( numpy.multiply( masses['w'], coords[:,1] ) ) #of
    Rz = 1/M * sum( numpy.multiply( masses['w'], coords[:,2] ) ) #mass
    x = coords[:,0]
    y = coords[:,1]
    z = coords[:,2]
    m = masses['w']
    Ixx = sum( numpy.multiply( m, numpy.power(y,2) + numpy.power(z,2) ) ) #moment
    Ixy = -sum( numpy.multiply( m, numpy.multiply(x,y) ) )                #of
    Ixz = -sum( numpy.multiply( m, numpy.multiply(z,y) ) )                #inertia
    Iyy = sum( numpy.multiply( m, numpy.power(x,2) + numpy.power(z,2) ) )
    Iyz = -sum( numpy.multiply( m, numpy.multiply(y,z) ) )
    Izz = sum( numpy.multiply( m, numpy.power(x,2) + numpy.power(y,2) ) )
    I = numpy.array( [ [ Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz] ] )
    Ieig,Ivec = numpy.linalg.eig(I)
    OvK = h * freqs_cleaned / kB #(K)

    # calculate the partition function
    qt = (2 * math.pi * M * kB * T / h**2)**1.5 * kB * T / P
    qr = math.sqrt(math.pi) / o * (8 * math.pi**2 * kB * T / h**2)**1.5 * math.sqrt( numpy.prod(Ieig) )

    # calculate the contributions to entropy
    St = R * ( math.log(qt) + 1 + 1.5 )
    Se = R * math.log(mult)
    Sr = R * ( math.log(qr) + 1.5 )
    Svs = R * ( numpy.divide( OvK/T, numpy.exp(OvK/T) - 1 ) - numpy.log( 1 - numpy.exp(-OvK/T) ) )
    Sv = sum( Svs )
    Stot = St + Se + Sr + Sv

    # read the entropies calculated by Gaussian
    St_gauss = re.search( " Translational[ ]*[0-9]+\.[0-9]{3}[ ]*[0-9]+\.[0-9]{3}[ ]*([0-9]+\.[0-9]{3})\n", thefile )
    St_gauss = float( St_gauss.group(1) )
    St_gauss_str = ''
    if abs( St_gauss - St / cal ) > 0.0005:
        St_gauss_str = '(Gaussian: %7.3f)' % St_gauss
    Se_gauss = re.search( " Electronic[ ]*[0-9]+\.[0-9]{3}[ ]*[0-9]+\.[0-9]{3}[ ]*([0-9]+\.[0-9]{3})\n", thefile )
    Se_gauss = float( Se_gauss.group(1) )
    Se_gauss_str = ''
    if abs( Se_gauss - Se / cal ) > 0.0005:
        Se_gauss_str = '(Gaussian: %7.3f)' % Se_gauss
    Sr_gauss = re.search( " Rotational[ ]*[0-9]+\.[0-9]{3}[ ]*[0-9]+\.[0-9]{3}[ ]*([0-9]+\.[0-9]{3})\n", thefile )
    Sr_gauss = float( Sr_gauss.group(1) )
    Sr_gauss_str = ''
    if abs( Sr_gauss - Sr / cal ) > 0.0005:
        Sr_gauss_str = '(Gaussian: %7.3f)' % Sr_gauss
    Sv_gauss = re.search( " Vibrational[ ]*[0-9]+\.[0-9]{3}[ ]*[0-9]+\.[0-9]{3}[ ]*([0-9]+\.[0-9]{3})\n", thefile )
    Sv_gauss = float( Sv_gauss.group(1) )
    Sv_gauss_str = ''
    if abs( Sv_gauss - Sv / cal ) > 0.0005:
        Sv_gauss_str = '(Gaussian: %7.3f)' % Sv_gauss
    Stot_gauss = re.search( " Total[ ]*[0-9]+\.[0-9]{3}[ ]*[0-9]+\.[0-9]{3}[ ]*([0-9]+\.[0-9]{3})\n", thefile )
    Stot_gauss = float( Stot_gauss.group(1) )
    Stot_gauss_str = ''
    if abs( Stot_gauss - Stot / cal ) > 0.0005:
        Stot_gauss_str = '(Gaussian: %7.3f)' % Stot_gauss
    print """
Calculated in the same way as Gaussian does it:
St     = %8.3f cal/mol*K %s
Se     = %8.3f cal/mol*K %s
Sr     = %8.3f cal/mol*K %s
Sv     = %8.3f cal/mol*K %s
Stot   = %8.3f cal/mol*K %s""" % ( St / cal, St_gauss_str, Se / cal, Se_gauss_str, Sr / cal, Sr_gauss_str, Sv / cal, Sv_gauss_str, Stot /cal, Stot_gauss_str )

    # calculate the contributions to the internal energy
    Et = 1.5 * RHa * T
    Er = Et
    Ev = RHa * numpy.sum( numpy.multiply( OvK, 0.5 + 1 / ( numpy.exp(OvK / T) - 1 ) ) )
    zpCorr = RHa * numpy.sum( numpy.multiply( OvK, 0.5 ) )
    Etot = Et + Er + Ev
    Hcorr = Etot + kB * T
    Gcorr = Hcorr - T * Stot * molJ
    
    # read the corrections calculated by Gaussian
    zpCorr_gauss = re.search( " Zero-point correction=[ ]+([0-9]+\.[0-9]{6}).*\n", thefile )
    zpCorr_gauss = float( zpCorr_gauss.group(1) )
    zpCorr_gauss_str = ''
    if abs ( zpCorr_gauss - zpCorr ) > 0.0000005:
        zpCorr_gauss_str = '(Gaussian: %.6f)' % zpCorr_gauss
    Etot_gauss = re.search( " Thermal correction to Energy=[ ]+([0-9]+\.[0-9]{6}).*\n", thefile )
    Etot_gauss = float( Etot_gauss.group(1) )
    Etot_gauss_str = ''
    if abs ( Etot_gauss - Etot ) > 0.0000005:
        Etot_gauss_str = '(Gaussian: %.6f)' % Etot_gauss
    Hcorr_gauss = re.search( " Thermal correction to Enthalpy=[ ]+([0-9]+\.[0-9]{6}).*\n", thefile )
    Hcorr_gauss = float( Hcorr_gauss.group(1) )
    Hcorr_gauss_str = ''
    if abs ( Hcorr_gauss - Hcorr ) > 0.0000005:
        Hcorr_gauss_str = '(Gaussian: %.6f)' % Hcorr_gauss
    Gcorr_gauss = re.search( " Thermal correction to Gibbs Free Energy=[ ]+([0-9]+\.[0-9]{6}).*\n", thefile )
    Gcorr_gauss = float( Gcorr_gauss.group(1) )
    Gcorr_gauss_str = ''
    if abs ( Gcorr_gauss - Gcorr ) > 0.0000005:
        Gcorr_gauss_str = '(Gaussian: %.6f)' % Gcorr_gauss
    print """Thermal corrections at %.2f K and %.2f atm:
zpCorr = %11.6f Ha %s
Etot   = %11.6f Ha %s
Hcorr  = %11.6f Ha %s
Gcorr  = %11.6f Ha %s""" % ( T, P / atm, zpCorr, zpCorr_gauss_str, Etot, Etot_gauss_str, Hcorr, Hcorr_gauss_str, Gcorr, Gcorr_gauss_str )

    # calculate the vibrational entropy according to Grimme
    Bav = sum( Ieig ) / len( Ieig )
    mu = h / ( 8 * math.pi**2 * freqs_cleaned )
    mup = mu * Bav / ( mu + Bav )
    SRs = R * ( 0.5 + numpy.log( numpy.sqrt( 8 * math.pi**3 * mup * kB * T / h**2 ) ) )
    w = 1 / ( 1 + ( w0 / freqs_cleaned )**alpha )
    Ss = w * Svs + ( 1 - w ) * SRs
    S = sum( Ss )
    StotG = St + Se + Sr + S
    print """
Calculated with small vibrations (<~ %d cm**-1) treated like rotations:
(DOI: 10.1002/chem.201200497)
Sv     = %8.3f cal/mol*K
Stot   = %8.3f cal/mol*K""" % ( w0 / cm1, S / cal, StotG /cal )

    # calculate the contributions to the internal energy
    GcorrG = Hcorr - T * StotG * molJ
    print """Thermal corrections at %.2f K and %.2f atm:
Gcorr  = %11.6f Ha""" % ( T, P / atm, GcorrG )

    
    # set the entropic contributions of all vibrations smaller than the cutoff to a constant value
    freqs_w0 = [ i for i, e in enumerate( freqs_cleaned ) if e < w0 ]
    OvK_w0 = OvK
    OvK_w0[ freqs_w0 ] = h * w0 / kB * numpy.ones( len( freqs_w0 ) )
    Sv_w0 = R * sum( OvK_w0 / T / ( numpy.exp( OvK_w0 / T ) - 1 ) - numpy.log( 1 - numpy.exp( -OvK_w0 / T ) ) )
    Stot_w0 = St + Se + Sr + Sv_w0

    print """
Calculated with vibrations < %d cm**-1 set to %d cm**-1:
Sv     = %8.3f cal/mol*K
Stot   = %8.3f cal/mol*K""" % ( w0 / cm1, w0 / cm1, Sv_w0 / cal, Stot_w0 / cal )

    # calculate the contributions to the internal energy
    Ev_w0 = RHa * numpy.sum( numpy.multiply( OvK_w0, 0.5 + 1 / ( numpy.exp(OvK_w0 / T) - 1 ) ) )
    zpCorr_w0 = RHa * numpy.sum( numpy.multiply( OvK_w0, 0.5 ) )
    Etot_w0 = Et + Er + Ev_w0
    Hcorr_w0 = Etot_w0 + kB * T
    Gcorr_w0 = Hcorr_w0 - T * Stot_w0 * molJ
    print """Thermal corrections at %.2f K and %.2f atm:
zpCorr = %11.6f Ha
Etot   = %11.6f Ha
Hcorr  = %11.6f Ha
Gcorr  = %11.6f Ha""" % ( T, P / atm, zpCorr_w0, Etot_w0, Hcorr_w0, Gcorr_w0 )


if usr_input in [ 2, 3 ]:

    # using the imaginary frequency vibrational mode, create input files for reactant and product

    # catch the charge
    charge = re.search( "Charge =[ ]+(-?[0-9]+)", thefile )
    charge = int( charge.group(1) )

    # catch the displacement vectors for each vibrational mode
    displms = re.findall("[ ]*[0-9]{1,3}[ ]+[0-9]{1,3}[ ]+-?[0-9]+\.[0-9]{2}[ ]+-?[0-9]+\.[0-9]{2}[ ]+-?[0-9]+\.[0-9]{2}.*\n", thefile) #[ ]+-?[0-9]+\.[0-9]{2}[ ]+-?[0-9]+\.[0-9]{2}[ ]+-?[0-9]+\.[0-9]{2}[ ]+-?[0-9]+\.[0-9]{2}[ ]+-?[0-9]+\.[0-9]{2}[ ]+-?[0-9]+\.[0-9]{2}\n",thefile)
    displms_dict = {}
    i = -1
    for row in displms:
	split_row = row.split()
	reng = range(2,len(split_row),3)
	nreng = enumerate(reng)
	if split_row[0] == '1':
	    i += 1
	    vibr = [ x+i*3 for x in range( 0, len(reng) ) ]
	    for j,k in nreng:
		displms_dict[ vibr[j] ] = []
	    nreng = enumerate( reng )
	for j,k in nreng:
	    displms_dict[ vibr[j] ].append( [ float(x) for x in split_row[k:k+3] ] )

    fo.close()
    
    # select which imaginary frequencies to displace along
    imag_vibr = [ i for i, e in enumerate( freqs ) if e < 0 ]
    no_imag_vibr = len( imag_vibr )
    if no_imag_vibr > 1:
        vibr_input = raw_input( """
The file contains %d imaginary frequencies. Which ones do you want to displace along?
Give the answer as a space-separated list ('0' for all): """ % no_imag_vibr )
        if vibr_input is not '0':
            imag_vibr = vibr_input.split()
            imag_vibr = [ int(x)-1 for x in imag_vibr ]
    displm_dict = {}
    j = -1
    for i in imag_vibr:
        j += 1
        displm_dict[ j ] = numpy.array( displms_dict[ i ] )
    
    # function to convert to roman numerals (from http://code.activestate.com/recipes/81611-roman-numerals/)
    def int_to_roman(input):
        if type(input) != type(1):
            raise TypeError, "expected integer, got %s" % type(input)
        if not 0 < input < 4000:
            raise ValueError, "Argument must be between 1 and 3999"   
        ints = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
        nums = ('m',  'cm', 'd', 'cd','c', 'xc','l','xl','x','ix','v','iv','i')
        result = ""
        for i in range(len(ints)):
            count = int(input / ints[i])
            result += nums[i] * count
            input -= ints[i] * count
        return result

    # displace the geometries along the displacement vectors
    el = [ "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Uun", "Uuu", "Uub" ]
    displm_dict_len = len( displm_dict )
    displm_dict_keys = [ k for ( k, v ) in displm_dict.iteritems() ]
    # create an array of permutations of all possible displacement combinations
    no_of_permuts = 2**displm_dict_len
    permuts_signs = numpy.zeros( ( no_of_permuts, displm_dict_len ) )
    change_sign_every = no_of_permuts
    for i in range( displm_dict_len ):
        change_sign_every /= 2
        sign = 1
        sign_counter = 0
        for j in range( no_of_permuts ):
            if sign_counter == change_sign_every:
                sign = -sign
                sign_counter = 0
            sign_counter += 1
            permuts_signs[ j, i ] = sign
    permuts_signs *= displm_factor
    # for each displacement combination, create an input file
    print_vibr_no = False
    if range( len( permuts_signs[0] ) ) != imag_vibr:
        print_vibr_no = True
    for permut in permuts_signs:
        filename = freqfile[ :-4 ]
        for k, elem in enumerate( permut ):
            appendix = str( abs( int( elem*10 ) ) )
            if abs( elem ) < 1:
                appendix = "0" + appendix
            if elem < 0:
                appendix = "-" + appendix
            elif elem > 0:
                appendix = "+" + appendix
            if print_vibr_no:
                appendix = int_to_roman( imag_vibr[ k ] + 1 ) + appendix
            filename = filename + appendix
        filename = filename + ".com"
        while os.path.isfile( filename ):
            filename = input( 'The file %s does already exist. Enter a new filename (in quotes): ' % filename )
        fo = open( filename, 'w' )
        fo.write("""%%chk=%s.chk
%%mem=
%%nprocshared=
# route

%s

%d %d
""" % ( filename[ :-4 ], 'Title Card Required', charge, mult ) )
        displ_coords = coordsA
        for k, displm in enumerate( permut ):
            displ_coords = displ_coords + displm * displm_dict[ displm_dict_keys[ k ] ]
        for j, row in enumerate( displ_coords ):
            fo.write( '%2s              ' % el[ masses[ 'N' ][ j ] ] )
            for elem in row:
                fo.write( '%14.8f' % elem )
            fo.write( '\n' )
        fo.write( ' ' )
        fo.close()
        print """
The Gaussian input file %s has been created.""" % filename 
