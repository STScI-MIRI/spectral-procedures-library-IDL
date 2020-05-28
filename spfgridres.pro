PRO spfgridres,infile,lamfile,outfile,fwhm

; 23 Oct 15 modified to pass lamdata to spgridres, not just lam
;           this will ensure order information saved
; 25 Feb 14 created
;
; file wrapper for spgridres
;
; INPUT
; infile  - input spectral FITS file
; lamfile - spectral FITS file with new wavelength grid in col 0
; outfile - output spectral FITS file
; fwhm    - FWHM for gaussian when regridding
;           default = 0, which will force spgridres to compute it

; initializations

lcol=0

; check keywords

if (keyword_set(fwhm) eq 0) then fwhm=0

; load input file and header

input = readfits(infile,hdr,/silent)

; load wavelength grid from lamfile

lamdata = readfits(lamfile,/silent)
lam     = reform(lamdata[lcol,*])

; check wavelengths - if the same, then do NOT regrid

testlam = reform(input[lcol,*])
if (total(abs(testlam-lam)) gt 0) then gridflag=1 else gridflag=0

; regrid if the wavelength grids differ and write output

if (gridflag eq 1) then output = spgridres(input,lamdata,fwhm) else output=input

wr_spfits,outfile,output,-1,FITS_HDR=hdr

END
