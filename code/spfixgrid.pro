FUNCTION spfixgrid,sparray,lamarray

;  8 Dec 11 created
;
; checks to ensure that two spectra are on the same grid
;   sparray is the spectrum under consideration
;   lamarray is the comparison - if only 2 columns, they are assumed lam, order
; if the grids are the same, spfixgrid passes sparray back intact
; if not, it returns sparray after has been regridded to the lamarray grid 
;
; the code below was originally in /home/sloan/procs/irs/spccbat.pro
; but has been moved here as it's generic and useful

; column definitions

lcol=0 & ocol=3
if (n_elements(lamarray[*,0]) ge 4) then olcol = 3 else olcol = 1

lam=reform(lamarray[lcol,*])
ord=reform(lamarray[olcol,*])
nlam=n_elements(lam)

goodflag=1

lamtest=reform(sparray[lcol,*])
ordtest=reform(sparray[ocol,*])
nlamtest=n_elements(lamtest)
dellam=lamtest-shift(lamtest,1)
dellam=max(abs(dellam[1:nlamtest-2]))

if (max(abs(lamtest-lam)) gt 0.01*dellam) then goodflag=0
if (max(abs(ordtest-ord)) gt 0.0) then goodflag=0
if (nlamtest ne nlam) then goodflag=0

; if goodflag=0, then regrid

if (goodflag eq 0) then sp=sporegrid(sparray,lamarray) else sp=sparray

RETURN,sp

END
