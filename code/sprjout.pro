FUNCTION sprjout,inspec,RJ=rj

; 17 Feb 11 added RJ flag - if 0 or not set, just return input spectrum
; 11 Feb 11 created from sprj
;
; sprjout returns a spectrum with the flux and errors divided by lambda^2
; assumes col 0 = wavelength, col 1 = flux, col 2 = error

sz=size(inspec)

lcol=0 & fcol=1 & ecol=2

outspec=inspec

if (keyword_set(rj) ne 0) then begin
  outspec[fcol,*]=inspec[fcol,*]/inspec[lcol,*]^2
  if (sz[1] gt 2) then outspec[ecol,*]=inspec[ecol,*]/inspec[lcol,*]^2
endif

return,outspec
END
