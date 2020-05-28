FUNCTION spzero,sp,l0,l1,INCLUSIVE=inclusive,EXCLUSIVE=exclusive

; 27 Jun 06 created
;
; spzero zeroes all data inside or outside the wavelength range l0-l1
; depending on whether /inclusive or /exclusive are set
; default is /inclusive
; "inclusive" deletes data up to and including the input wavelengths
; "exclusive" does not delete the data matching the input wavelengths

lcol=0 & fcol=1 & ecol=2 & ocol=3

if (keyword_set(exclusive) eq 1) then zeroflag=1 else zeroflag=0

spnew=sp
if (zeroflag eq 0) then idx=where(sp[lcol,*] ge l0 and sp[lcol,*] le l1) $
  else idx=where(sp[lcol,*] lt l0 or sp[lcol,*] gt l1) $

if (max(idx) gt -1) then begin
  spnew[fcol,idx]=0.0
  spnew[ecol,idx]=0.0
endif

RETURN,spnew
END
