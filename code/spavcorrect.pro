FUNCTION spavcorrect,inspec,av,inex,FILE=file,NORM=norm

; 25 Nov 12 flipped order of inputs, added FILE and NORM keywords
;           corrected normalization at K (from 0.122 to 0.112)
;  2 Sep 10 modified by G. Sloan
; 28 Jul 10 created by Dominic A. Ludovici
;
; Corrects an input spectrum for extinction
; INPUT
;   inspec - input spectrum (spectral data array)
;   av     - A_v - extinction in visual magnitudes of the source
;   inex   - extinction correction - must cover wavelength range of inspec
;            expected that inex is in mag/mag and normalized at K
;   file   - keyword to specify extinction file (instead of passing the data)
;   NORM   - keyword, if not set or 0, normalization is at V, otherwise K
; OUTPUT   - returns the extinction-corrected spectrum

; setup

lcol=0 & fcol=1 & ecol=2 & ocol=3

; read extinction correction if file keyword set
;   otherwise it must be passed

if (keyword_set(file) ne 0) then inex=readfits(file,/silent)

; set normalization value to adust K-normalized extinction to V-normalization
; note that Chiar & Tielens suggest 0.09, not 0.112

if (keyword_set(norm) eq 0) then normalization=0.112 else normalization=1.0

; add cols to inex until it has four (to be a legal spectral data array) 

while(n_elements(inex[*,0]) lt 4) do inex=[inex,fltarr(1,n_elements(inex[0,*]))]

; outex = extinction scaled to A_k (A_k = 0.122 A_v -- Rieke & Lebofsky 1983)
; regrid inex to wavelengths of inspec
; convert from magnitudes to a throughput

outex = sptimes(inex,normalization*av)
outex = spregrid(outex,inspec)
outex[1,*] = 2.512^(-outex[1,*])

; correct input spectrum for extinction

outspec=spdivide(inspec,outex)

RETURN,outspec
END
