FUNCTION spdowngrid,sparray,lamarray

; 25 Jun 13 created
;
; spdowngrid regrids sparray to the wavelength grid in lamarray
;   using calls to spspot, which will use linear interpolation
; this routine is designed for cases where sparray is
; (1) on a very find grid w.r.t. lamarray
; (2) has high (or infinite) S/N
; (3) already has the appropriate resolution built into the data

; setup

lcol=0 & fcol=1 & ecol=2 & ocol=3 & mincol=4

; check dimensions of lamarray

lamcol=n_elements(lamarray[*,0])
lamrow=n_elements(lamarray[0,*])

; create outspec - minimum of four columns

outspec=fltarr(min([lamcol,mincol]),lamrow)
outspec[lcol,*] = lamarray[lcol,*]

; check dimensions of sparray

ncol=n_elements(sparray[*,0])
nrow=n_elements(sparray[0,*])

for i=0,lamrow-1 do begin

; flux - simple interpolation with spspot calls

  outspec[fcol,i] = spspot(sparray,lamarray[lcol,i],/silent)

; error - if it exists in sparray, interpolate it, too

  if (ncol gt 2) then $
    outspec[ecol,i] = spspot(sparray,lamarray[lcol,i],/error,/silent)

endfor

; order - copy directly from lamarray - assume col 1 if ncol=2, col 3 if ncol=4+

case lamcol of
  1 : lamocol=-1
  2 : lamocol=1
  3 : lamocol=2
  else : lamocol=3
endcase

if (lamocol gt 0) then outspec[ocol,*] = lamarray[lamocol,*]

RETURN,outspec
END
