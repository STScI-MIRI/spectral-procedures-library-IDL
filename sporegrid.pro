FUNCTION sporegrid,sparray,lamarray

;  7 Dec 11 created
;
; sporegrid is a wrapper for spregrid
; regrids order by order, only regridding orders in both input arrays

; column definitions

lcol=0 & fcol=1 & ecol=2 & ocol=3
if (n_elements(lamarray[*,0]) ge 4) then olcol = 3 else olcol = 1

; determine min and max orders in sparray and lamarray

minorder=min([min(sparray[ocol,*]),min(lamarray[olcol,*])])
maxorder=max([max(sparray[ocol,*]),max(lamarray[olcol,*])])

count=0
for m=minorder,maxorder do begin
  idx=where(sparray[ocol,*] eq m)
  jdx=where(lamarray[olcol,*] eq m)

  if (idx[0] gt -1 and jdx[0] gt -1) then begin
    sporder = spregrid(sparray[*,idx],lamarray[*,jdx])
    if (count eq 0) then sp=sporder else sp=[[sp],[sporder]]
    count=count+1
  endif
endfor

RETURN,sp
END
