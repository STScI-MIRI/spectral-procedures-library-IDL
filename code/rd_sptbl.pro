FUNCTION rd_sptbl,infile,header,fluxerror=fluxerror,order=order,$
  noerror=noerror,noorder=noorder

; 15 Mar 18 added noerror and noorder keywords to just set these to zero
;           also turned off passing header if it doesn't appear on command line
;  9 Aug 17 added fluxerror and segment keywords for optional column headers
; 28 Jul 14 added ability to handle header
; 20 Jul 14 created
;
; read a spectral data file in IPAC table format
; return result in the usual 4-column format:  lambda, F_nu, error, order

in=read_ipac_table(infile)
len=n_elements(in.flux)

out=fltarr(4,len)
out[0,*] = in.wavelength
out[1,*] = in.flux
if (keyword_set(noerror) eq 0) then begin
  if (keyword_set(fluxerror) eq 0) then out[2,*] = in.error $
    else out[2,*] = in.flux_error
endif
if (keyword_set(noorder) eq 0) then begin
  if (keyword_set(order) eq 0) then out[3,*] = in.segment $
    else out[3,*] = in.order
endif

if (keyword_set(header) ne 0) then $
  header=in.header_table_header ; to pass back header 

RETURN,out
END
