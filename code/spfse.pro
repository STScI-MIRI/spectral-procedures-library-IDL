FUNCTION spfse,spfile,starfile,PLOT=plot,_extra=e

; 23 Dec 06, now passing extra keywords to spse and on to spplot
;  5 Nov 04, modified slightly
; 23 Aug 04, created, a wrapper for spse.pro

a = readfits(spfile,/silent)
b = readfits(starfile,/silent)

if (not keyword_set(plot)) then return,spse(a,b,_extra=e) $
  else return,spse(a,b,/plot,_extra=e)

end
