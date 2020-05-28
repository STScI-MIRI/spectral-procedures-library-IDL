FUNCTION spcolor,sp,l0,l1,l2,l3,ERROR=error

;  6 Sep 05 corrected error calculation
; 29 Apr 05 created
;
; given a spectrum and four wavelengths, compute a color
; uses a Rayleigh-Jeans tail to calibrate the color

fa   = spavg(sp,l0,l1,err=erra)
fb   = spavg(sp,l2,l3,err=errb)
;err  = sqrt( (fa/fb)^2 * ( (erra/fa)^2 + (errb/fb)^2 ) )

color0 = -2.5*alog10((mean([l0,l1])/mean([l2,l3]))^2)
color  = -2.5*alog10(fa/fb)+color0
errcol = sqrt( (2.5/alog(10)) * ( (erra/fa)^2 + (errb/fb)^2 ) )

if (n_elements(error) eq 0) then return,color else return,[color,errcol]
END
