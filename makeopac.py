#
# Import NumPy for array handling and system and os for file handling
#
import numpy as np
import math
import os, sys, glob, shutil
from subprocess import call
#
# Import natconst for natural constants in CGS units
#
import natconst as nc

# --------------------------------------------------------------------------------
# Store the opacity parameters into the file dust_param.inp for later use.
# --------------------------------------------------------------------------------
def makeopac(dustpar):
    """
    This subroutine stores the dust parameters into the file dust_param.inp.
    Whenever a new opacity calculation has to be done, this file contains
    the information. The store_dustpar() subroutine should not be called by
    the user. It is automatically called by the setup() routine. This means
    that once the dust opacity parameters are set in the setup() routine, they
    can (should) not be changed anymore, unless you make a new setup() call.

    The dustpar is a dictionary with parameters. They are all optional (with defaults):
    'material'       The name of the material out of which the dust is made.
                     This must be specified. A file containing the optical constants
                     named <material>.lnk must be given. Example: if dustpar['material']
                     equals 'pyrmg70', then a file pyrmg70.lnk must be present. This
                     file should contain 3 colums: col 1 = wavelength in micron, col 2 =
                     real part of the refractive index and col 3 = imaginary part of
                     the refractive index. You can find such files from the Jena
                     database: http://www.astro.uni-jena.de/index.php/laboratory/databases.html
                     Note that you also have to find on that database the material density
                     of the material. See 'matdens' below.
                     NOTE: The wavelength bins of that file do not have to be the same
                     as those that you use for the RADMC-3D model, BUT the smallest and
                     largest wavelength of this file should bracket the wavelengths you
                     use for the RADMC-3D model.
    'matdens'        Material density of the dust in gram/cm^3 (typically around 3.0 or so)
    'amin'           Minimum grain radius in cm
    'aminmic'        Minimum grain radius in micron
    'amax'           Maximum grain radius in cm
    'amaxmic'        Maximum grain radius in micron
    'pla'            Powerlaw of the grain size distribution (default -3.5)
    'na'             Number of grain size bins (default is 100, which is on the safe side;
                     smaller might also be ok)
    'nrang'          Nr of scattering angle bins to be considered (typically 181)
    'errtol'         The relative error tolerance for testing the following integral:
                        kappa_scat = 2 pi int_{-1}^{+1} Z_11(mu) dmu
    'refinefact'     If set, then use (instead of nrang), and produce, the scattering angular
                     grid in the file scatangle.inp. The refinefact is used in a simple
                     formula for refinement of the grid:
                        angle(i) = 90 * (exp(refinefact*(i-1)/nang)-1) / (exp(refinefact)-1)
                     where i=1...nang with nang the number of angle points between 0 and 90
                     degrees. Note that we only specify the angles between 0 and 90 degrees
                     in the file scatangle.inp because that is the way the BHMIE code of
                     Bohren & Hufmann / Draine is constructed. Recommended value e.g.
                     refinefact = 3.0.

    It returns: true if the dust parameters have changed compared to the previous
              values of dust_param.inp. If false, then the parameters given in
              dustpar are the same as found in dust_param.inp
    """
    #
    # The parameters
    #
    iformat    = 1               # For now the file format for dust_param.inp is always 1
    assert 'material' in dustpar, "Error: Without specifying 'material' in dustpar, I cannot proceed!"
    material   = dustpar['material'].strip()  # The material out of which the dust is made
    material   = material.strip()
    dustmodel  = 'powerlawdistr' # Default: Dust size distribution model
    amin       = 1e-5            # Default: (for the power law model) Minimum grain radius in cm
    amax       = 1e-2            # Default: (for the power law model) Maximum grain radius in cm
    pla        = -3.5            # Default: (for the power law model) index for size distribution
    na         = 100             # Default: Nr of grain size bins
    matdens    = 3.0             # Default: Material density of the dust in gram/cm^3
    nrang      = 181             # Default: Nr of scattering angles to be considered
    errtol     = 0.1             # Default: Error tolerance
    refinefact = 0.              # Default: Refinefact (0 means no refinement, just standard linear)
    a0         = 1e-4            # Default: (for the Gaussian model) peak size of the Gauss
    awidth     = 1e-5            # Default: (for the Gaussian model) width of the Gauss
    da         = 2e-5            # Default: (for the Gaussian model) width of the window
    wavelength = -1.0            # Needs to be specified. Wavelength in cm.
    singlewavelength = True      # If wavelength or wavelengthmic is specified, we calculate the opacity file for the single wavlength.
                                 # If not specified, we will calculate the opacity for the full wavelength range. But this makes the file size huge.
    if 'dustmodel' in dustpar:     dustmodel     = dustpar['dustmodel'].strip()
    if 'amin' in dustpar:          amin          = dustpar['amin']
    if 'amax' in dustpar:          amax          = dustpar['amax']
    if 'aminmic' in dustpar:       amin          = dustpar['aminmic']*1e-4
    if 'amaxmic' in dustpar:       amax          = dustpar['amaxmic']*1e-4
    if 'pla' in dustpar:           pla           = dustpar['pla']
    if 'na' in dustpar:            na            = dustpar['na']
    if 'matdens' in dustpar:       matdens       = dustpar['matdens']
    if 'nrang' in dustpar:         nrang         = dustpar['nrang']
    if 'errtol' in dustpar:        errtol        = dustpar['errtol']
    if 'refinefact' in dustpar:    refinefact    = dustpar['refinefact']
    if 'a0' in dustpar:            a0            = dustpar['a0']
    if 'awidth' in dustpar:        awidth        = dustpar['awidth']
    if 'da' in dustpar:            da            = dustpar['da']
    if 'wavelength' in dustpar:    wavelength    = dustpar['wavelength']
    if 'wavelengthmic' in dustpar: wavelength    = dustpar['wavelengthmic']*1e-4
    #assert wavelength > 0., "Error: wavelength is not specified."
    if wavelength < 0:             singlewavelength=False
    #
    # Check if the optical constants file exists
    #
    optconstfile = material+'.lnk'
    assert os.path.isfile(optconstfile), "Error: I could not find the optical constants file "+optconstfile+" in the local directory. Please read the manual of store_dustpar."
    #
    # If refinefact != 0 then make scatangle.inp
    #
    if(refinefact!=0):
        #
        # We make a refined angular grid
        #
        assert nrang%2!=0, "Error: nr of angle must be odd (nrang is nr of angles between 0 and 180 degrees)"
        nang  = int(((nrang)+1)/2)    # nang is nr of angles between and including 0 and 90
        angle = 90.*(np.exp(refinefact*np.linspace(0.,1.,nang))-1.0)/(math.exp(refinefact)-1.0)
        nrangout = 0  # Signal to the dust opacity code to use scatangle.inp
        #
        # We write the scattering angle file for the dust opacity
        # creation code. Here we write only the grid from 0 to 90
        # degrees, because the Mie code we use here mirrors the
        # grid from 0-90 degrees to 180-90 degrees.
        #
        with open('scatangle.inp','w+') as f:
            f.write('%d %d %13.6e\n' %( nang, 1, refinefact))
            np.savetxt(f,angle.T,fmt=['%13.6e'])
        #
        # The RADMC-3D code, on the other hand, can use any grid
        # of scattering angle. That file is formatted a bit different,
        # and includes the full 0-180 degrees.
        #
        with open('scattering_angular_grid.inp','w+') as f:
            f.write('%d\n' %(1))
            f.write('%d\n' %(nrang))
            anglemirror = 180.-angle[nang-2::-1]
            np.savetxt(f,angle.T,fmt=['%13.6e'])
            np.savetxt(f,anglemirror.T,fmt=['%13.6e'])
    else:
        nrangout = nrang
        if os.path.isfile("scatangle.inp"):
            os.remove("scatangle.inp")
        if os.path.isfile("scattering_angular_grid.inp"):
            os.remove("scattering_angular_grid.inp")
    #
    # Now make the dust_param.inp file
    #
    with open('dust_param.inp','w+') as f:
        f.write('%d\n' % 1)
        f.write(dustmodel+'\n')
        if dustmodel == 'powerlawdistr':
            f.write(material+'\n')
            f.write('%d\n' % na)
            f.write('%e\n' % amin)
            f.write('%e\n' % amax)
            f.write('%e\n' % pla)
            f.write('%e\n' % matdens)
            f.write('%d\n' % nrangout)
            f.write('%e\n' % errtol)
        elif dustmodel == 'gaussdistr':
            f.write(material+'\n')
            f.write('%d\n' % na)
            f.write('%e\n' % a0)
            f.write('%e\n' % awidth)
            f.write('%e\n' % da)
            f.write('%e\n' % matdens)
            f.write('%d\n' % nrangout)
            f.write('%e\n' % errtol)
        else:
            raise Exception('Unknown dust model')
    #
    # Make a string of the extension as a RADMC-3D opacity file
    #
    if singlewavelength:
        #
        # Make a string of the wavelength in units of micron
        #
        wlstring = '{:1.2e}'.format(wavelength*1e4)
        #
        #
        # Make the opacity file name
        #
        extension    = wlstring
    else:
        #
        # Make a string of the wavelength in units of micron
        #
        amaxmicstring = '{:3.2e}'.format(amax*1e4)
        #
        #
        # Make the opacity file name
        #
        extension    = material+'_amaxmic'+amaxmicstring
    opacfilename = 'dustkapscatmat_'+extension+'.inp'
    #
    # Check if the opacity file for this wavelength is already existing
    #
    if os.path.isfile(opacfilename): raise Exception('Error! The opacity file already exists.')
    #
    # Make the opacity file for this wavelength
    #
    #
    # Check if the executable is there
    #
    code = './makeopac_scatmat'
    assert os.path.isfile(code), "Error: The code makeopac_scatmat could not be found "
    #
    # Create the wavelength file for make_scatmat_distr.f90
    #
    if singlewavelength:
        with open('wavelength_micron_bh.inp','w+') as f:
            f.write('1\n')   # One wavelength
            f.write('%13.6e\n'%(wavelength*1e4))
    else:
        #
        # Write the wavelength_micron.inp file
        #
        lam1     = 0.1e0
        lam2     = 7.0e0
        lam3     = 25.e0
        lam4     = 1.0e4
        n12      = 20
        n23      = 100
        n34      = 30
        lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
        lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
        lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
        lam      = np.concatenate([lam12,lam23,lam34])
        nlam     = lam.size
        #
        # Write the wavelength file
        #
        # Note: If radmc3d mctherm is called, that will need and read the file wavelength_micron.inp.
        #
        with open('wavelength_micron_bh.inp','w+') as f:
            f.write('%d\n'%(nlam))
            np.savetxt(f,lam.T,fmt=['%13.6e'])
    #
    # Call the opacity code
    #
    call(code)
    #
    # Get the material name from the dust_param.inp
    #
    with open('dust_param.inp','r') as f:
        s        = f.readline()
        try: iformat  = int(s)
        except: print('Sorry, the first line of dust_param.inp must be the format number (an integer)')
        dustmodel = f.readline()
        dustmodel = dustmodel.strip()
        material = f.readline()
        material = material.strip()
    with open('dustopac.inp','w+') as f:
        f.write('2               Format number of this file\n')
        f.write('1               Nr of dust species\n')
        f.write('============================================================================\n')
        f.write('10              Way in which this dust species is read\n')  # 10 = read full scattering matrix
        f.write('0               0=Thermal grain\n')
        f.write(extension +  '   Extension of name of dustkappa_***.inp file\n')
        f.write('----------------------------------------------------------------------------\n')
##### end of the makeopac function ####


#
# Examples of execusion
#
#dustpar={'wavelengthmic':100.0,'material':'pyrmg70','nrang':181,'refinefact':3.0,'aminmic':0.01,'amaxmic':150.0,'na':200.}
#dustpar={'material':'dsharp','nrang':181,'refinefact':3.0,'aminmic':0.01,'amaxmic':150.0,'na':200.}
#makeopac(dustpar)
