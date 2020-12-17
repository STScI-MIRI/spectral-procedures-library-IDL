FUNCTION spfsum,file,l0,l1

; 20 Dec 05 created 
;
; spfsum reads a spectral FITS file, calls spsum, returns result

sp=readfits(file,/silent)
result=spsum(sp,l0,l1)

RETURN,result
END
