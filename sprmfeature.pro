FUNCTION sprmfeature,sp,l0,l1,l2,l3,degree,sigratio

;  4 Aug 06 now calling spadd with /noerror set to not propagate errors
; 27 Jun 06 created
;
; sprmfeature fits and removes a feature by replacing it with a continuum
; fitted from l0-l1 to l2-l3 using a call to spcont

if (keyword_set(degree) eq 0) then begin
  spc=spcont(sp,l0,l1,l2,l3) 
endif else begin
  if (keyword_set(sigratio) eq 0) then begin
    spc=spcont(sp,l0,l1,l2,l3,degree) 
  endif else begin
    spc=spcont(sp,l0,l1,l2,l3,degree,sigratio)
  endelse
endelse

spdiff=spadd(sp,spc,/minus,/noerror)
spnew=spadd(sp,spdiff,/minus,/noerror)

RETURN,spnew
END
