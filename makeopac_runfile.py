import makeopac
#dustpar={'material':'dsharp','nrang':181,'refinefact':3.0,'aminmic':0.01,'amaxmic':150.0,'na':200.}
dustpar={'wavelengthmic':870.0,'material':'dsharp','nrang':181,'refinefact':3.0,'aminmic':0.01,'amaxmic':150.0,'na':200.,'matdens':1.675}
makeopac.makeopac(dustpar)
