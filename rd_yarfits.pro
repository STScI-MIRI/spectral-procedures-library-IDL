FUNCTION rd_yarfits,infile

; 23 Mar 15 created
;
; loads a YAAAR-style FITS table as a simple spectral data array

yar = mrdfits(infile,1,/silent)

len = n_elements(yar.wavelength)

out = fltarr(5,len)

out[0,*] = yar.wavelength
out[1,*] = yar.flux
out[2,*] = yar.error
out[3,*] = float(yar.order)
out[4,*] = float(yar.bit_flag)

RETURN,out
END
