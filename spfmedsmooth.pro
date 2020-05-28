PRO spfmedsmooth,spfin,spfout,smoothbox,PLOT=plot,_extra=e

;  1 Sep 10 created
;
; a wrapper for spmedsmooth
; default for smoothbox = 3

if (keyword_set(smoothbox) eq 0) then smoothbox=3

spin=readfits(spfin,hdr,/silent)

spout=spmedsmooth(spin,smoothbox,_extra=e)

wr_spfits,spfout,spout,-1,FITS_HDR=hdr

if (keyword_set(plot) ne 0) then begin
  spplot,spin,_extra=e
  spplot,spout,_extra=e
  zq=get_kbrd(1)
endif

END
