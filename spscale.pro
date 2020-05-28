FUNCTION spscale,spa,spb,l0,l1


; 26 Jun 04 spectra only need to have same number of wavelength elements
; 18 Jun 04 created
;
; spscale returns the scaling factor for spb needed to minimize
;   the summed square difference between spa and spb between
;   the wavelengths l0 and l1
; assumed that spa and spb have same wavelength scale

; check spa and spb to be sure they match

errcode=0

sza=size(spa) & szb=size(spb)

if (sza[2] ne szb[2]) then begin
  print,"Error in spscale.  Input spectra do not match."
  errcode=1
endif 

idx=where(spa[0,*] ge l0 and spa[0,*] le l1)
if (max(idx) lt 0) then begin
  print,"Error in spscale.  No data in range."
  errcode=1
endif

if (errcode eq 0) then begin
  num=total(reform(spa[1,idx])*reform(spb[1,idx]))
  den=total(reform(spb[1,idx])^2)
  scale=num/den
endif else scale=1.0

return,scale

END
