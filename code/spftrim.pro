PRO spftrim,spin,spout,l0,l1,ORDARRAY=ordarray

; 19 May 06 created
;
; file wrapper for sptrim (cf)
;
; INPUT
;   spin     - file name of input spectral data array
;   spout    - file name for output
;   l0       - vector of starting wavelength for each order
;   l1       - vector of ending wavelength for each order
;   ordarray - vector of order numbers in order as given in l0, l1
;              this parameter helps the program in the event of missing orders
; OUTPUT - writes file spout
;

a=readfits(spin,hdr,/silent)
if (keyword_set(ordarray) eq 1) then b=sptrim(a,l0,l1,ORDARRAY=ordarray) $
  else b=sptrim(a,l0,l1)

wr_spfits,spout,b,-1,FITS_HDR=hdr

END
