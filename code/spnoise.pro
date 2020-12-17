FUNCTION spnoise,sparray,WIDTH=width

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 13 Apr 04 - Created by MJR
;
; Adds flux error to spfits images with no
; flux error
;
; sparray - an spfits array
; WIDTH - the width with which to calculate the median
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;default width is five pixels
if not( keyword_set(width) ) then width = 5

;find order range
if max(sparray[3,*]) gt 0 then begin
   omax = max(sparray[3,*])
endif else begin
   omax = 0
endelse

f = reform(sparray[1,*])
a = f

;calculate 'median', which is really the average of
;the two surrounding pixels. Proceed order by order
if omax eq 0 then begin
   for i=0,n_elements(f)-1 do begin
      if (i eq 0) then a[i] = (f[i] + f[i+2])/2
      if (i eq n_elements(f)-1) then a[i] = (f[i] + f[i-2])/2
      if ( i ne 0 and i ne n_elements(f)-1 ) then begin
         a[i] = (f[i-1] + f[i+1])/2
      endif
   endfor
endif else begin
   for j=0,omax-1 do begin
      idx = where(sparray[3,*] eq j+1)

      for k=0,n_elements(idx)-1 do begin
         if (k eq 0) then a[idx[k]] = (f[idx[k]] + f[idx[k]+2])/2
         if (k eq n_elements(idx)-1) then a[idx[k]] = (f[idx[k]] + f[idx[k]-2])/2 
         if (k ne 0 AND k ne n_elements(idx)-1) then begin
            a[idx[k]] = (f[idx[k]-1] + f[idx[k]+1])/2
         endif
      endfor
   endfor
endelse

sig = (f - a) / sqrt(2)
sig = median(sig,width)

sparray[2,*] = sig

return,sparray
END
