FUNCTION rd_cassis_tbl,infile,header,fluxerror=fluxerror,order=order,$
  noerror=noerror,noorder=noorder

; 23 Jul 18 created by modifying rd_sptbl.pro
;
; modifications are primarily for changes in column names
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
  if (keyword_set(order) eq 0) then out[3,*] = in.strip $
    else out[3,*] = in.order
endif

if (keyword_set(header) ne 0) then $
  header=in.header_table_header ; to pass back header 

RETURN,out
END
