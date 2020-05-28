PRO spfplot,infile,outfile,landscape=landscape,_extra=e

; 27 Feb 07 modified so that files ending with "eps" are saved encapsulated
; 24 May 04 created
;
; reads a spectral fits file, calls spplot to plot it as a PS file
; see spplot.pro for an explanation of legal keywords

if (keyword_set(landscape) eq 1) then landflag=1 else landflag=0

a=readfits(infile,/silent)

set_plot,'ps'
device,file=outfile
if (landflag eq 1) then device,/landscape
len=strlen(outfile)
if (strmid(outfile,len-3,3) eq 'eps') then device,/encapsulated

spplot,a,_extra=e

device,/close
set_plot,'x'

END
