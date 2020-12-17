FUNCTION spfpos,spfile,FOV=fov,EXT=ext,POS=pos

; 27 Mar 12 rudimentary error handlin added
;           if /ext set and _EXT position not found, find _FOV position instead
; 20 Feb 12 added the EXT and POS keywords
;  5 Jul 11 added the FOV keyword 
; 25 Aug 09 created
;
; reads a spectral FITS file
; returns the extracted position, as defined by the EXTPOS keyword
;
; INPUT
;   spfile - spectral FITS file - but only the header matters here
;   fov    - keyword to extract RA_FOV and DEC_FOV
;   ext    - keyword to extract RA_EXT and DEC_EXT (default)
;   pos    - keyword to just extract EXTPOS from the header

data=readfits(spfile,hdr,/silent)

if (data[0] gt -1) then begin 
  if (keyword_set(pos) ne 0) then pos=sxpar(hdr,'EXTPOS') else begin
    if (keyword_set(fov) ne 0) then begin
      ra=sxpar(hdr,'RA_FOV', count=count)
      dec=sxpar(hdr,'DEC_FOV', count=count)
      pos=[ra,dec]
    endif else begin ; default
      ra=sxpar(hdr,'RA_EXT', count=count)
      if (count ne 1) then ra=sxpar(hdr,'RA_FOV', count=count)
      dec=sxpar(hdr,'DEC_EXT', count=count)
      if (count ne 1) then dec=sxpar(hdr,'DEC_FOV', count=count)
      pos=[ra,dec]
    endelse
  endelse
endif else pos='0'

RETURN,pos
END
