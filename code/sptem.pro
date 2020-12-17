FUNCTION sptem,sp_array,tem_array

;  8 Apr 14 addressing instance when tem_array has just one wavelength element
;           and also check that all where calls return valid indices
; 17 Mar 05 sort input template
; 14 Sep 04 replaced cspline with interpol (using /lsquadratic)
; 28 Dec 03 created
;
; sptem returns a spectral template array gridded to match the spec. data array
;
; INPUT
;   sp_array - n-column spectral data array
;              assumed sp_array[3,*] = order number
;   tem_array - n-column spectral template array
; OUTPUT
;   returns regridded template in an array with dimensions of sp_array

lcol=0 & fcol=1 & ocol=3
lam   = reform(sp_array[lcol,*])
order = reform(sp_array[ocol,*])
nlam  = n_elements(lam)
ncol  = n_elements(sp_array[*,0])
out_array=dblarr(ncol,nlam)

m0=min(order) & m1=max(order)
tem = tem_array
idx = sort(reform(tem[lcol,*]))                         ; sort input template
tem = tem[*,idx]

if (n_elements(tem[lcol]) gt 1) then begin              ; cut dupe wavelengths
  idx = where(tem[lcol,*] ne shift(tem[lcol,*],1)) 
  if (idx[0] gt -1) then tem = tem[*,idx]
endif

for m=m0,m1 do begin
  odx = where(order eq m)
  if (odx[0] gt -1) then begin
    if (n_elements(tem[lcol,*]) gt 1) then begin
      olam = lam[odx]
;     out_array[fcol,odx] = cspline(reform(tem[lcol,*]),reform(tem[fcol,*]),olam)
      out_array[fcol,odx] = interpol(reform(tem[fcol,*]),reform(tem[lcol,*]),$
        olam,/lsq)
    endif else begin ; just one element in template, so return that value
      out_array[fcol,odx] = tem[fcol,*]
    endelse
  endif
endfor

out_array[lcol,*]=lam
out_array[ocol,*]=order

RETURN,out_array
END 
