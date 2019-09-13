# makeopac
Quick Start
- type 'make' (This creates the executable files, which will be called by the next line.)
- type 'python makeopac_runfile.py'

How to specify the parameters
dustpar={'wavelengthmic':870.0,'material':'dsharp','nrang':181,'refinefact':3.0,'aminmic':0.01,'amaxmic':150.0,'na':200.}
If you specify the wavelength, then this calculate the opacity file only at the single wavelength.
If the wavelength is not specified, this calculate the opacities at full wavelength range.
