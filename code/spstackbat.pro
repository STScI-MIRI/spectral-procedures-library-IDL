PRO spstackbat,infilelist,outfile,NORMALIZATION=normalization,SLBONUS=slbonus,LLBONUS=llbonus,DIAG=diag,MULTIPLE=multiple,FIXORDER=fixorder

;  3 Sep 14 repaired bug in Aug 14 improvements
;  8 Aug 14 norm=3 created - normalize entire stacked spectrum to max=1
;  7 Aug 14 now ignoring any fields to the right of a semicolon
; 29 Jul 14 modified to check input file for IPAC table or FITS
; 25 Jul 14 removed a bug - comment lines can now appear anywhere
;  2 Jan 14 now handling multiple lengths of spectral data files
;           omitting orders as requested in the input file
;           and no longer using header from first input spectrum
;  5 Mar 08 fixed minor bug with non-multiple mode
; 10 Apr 06 added multiple mode
; 30 Mar 06 added diagnostic flag, check for comment line (begins with ;)
; 30 Oct 04 removed some diagnostic print commands
; 29 Oct 04 added normalization parameter (default is 0)
;           and keywords slbonus and llbonus
; 30 Sep 04 created
;
; takes an input file list, reads all data in
; then finds normalization for each order of each file
; applies these, coadds, computes uncertainty in mean at each wavelength
; NOTE:  all files must be on same wavelength grid
;
; INPUT
;   infilelist    - a list of spectral FITS files to be coadded (no header)
;   outfile       - name of output FITS file to write
;   normalization - 0 (default) for no normalization (i.e. simple avg)
;                 - 1 to normalize each order to its mean separately
;                 - 2 to normalize each order to its max separately
;                 - 3 to normalize entire spectrum to its max
;   slbonus       - keyword, if set, slbonus will use SL2 normalizations 
;   llbonus       - keyword, if set, llbonus will use LL2 normalizations 
;   diag          - diagnostic flag to turn on some output
;   multiple      - multiple coadd mode
;                   if set, line 1 = number of coadds
;                   each group preceded by a count of files and output file
;   fixorder      - fix normalization to the given order
; OUTPUT - writes coadded file to outfile as a spectral FITS file
;
; formatting of input lines:
; if first character is ";" then treated as a comment and skipped
; if more than one string, following strings list orders to remove
;
; USER NOTE:
;   if setting norm=1 or norm=2, it is strongly advisable to set fixorder, too
;
; code currently optimized for lo-res data, which should be normalized 
;   aperture by aperture

; open input file 

openr,fi,infilelist,/get_lun
line=' '

; parse keywords, set flags, and set for multiple output files if needed

if (keyword_set(normalization) gt 0) then normflag=normalization else normflag=0

if (keyword_set(multiple) ne 0) then begin
  multiflag=1     ; multiple output files
  readf,fi,nfout  ; read number of output spectra from input file
  nfout=fix(nfout)
  if (n_elements(outfile) eq 0) then outfile=' '
endif else begin
  multiflag=0    ; single output file
  nfout=1        ; set number of output spectra=1
  if (n_elements(outfile) eq 0) then $ 
    print,'Warning in spstackbat.  Writing output to generic file out.fits.'
endelse

; more setup

lcol=0 & fcol=1 & ecol=2 & ocol=3
fitsfile=' '
fmt='(10(f8.3))'

; set up stops for wavelength averages
; this should be part of input file eventually

for i=1,nfout do begin ; loop through output spectra

  print,i,nfout

; if multiflag set, then read in number of input files and output filename 
; else output filename already set, and set no. of input files to -1

  if (multiflag eq 1) then begin

;   check for comments - keep reading until we reach a non-comment

    goflag=0
    while (goflag eq 0) do begin
      readf,fi,line
      checkchar = strmid(line,0,1)
      if (checkchar ne ';') then begin
        split = strsplit (line,' ',/extract)
        nfin = fix(split[0])
        outfile = split[1]
        goflag=1 
      endif
    endwhile
;   readf,fi,nfin,outfile                    ; outfile includes a leading space
;   nfin=fix(nfin)
;   outfile=strcompress(outfile,/remove_all) ; remove all spaces from filename
  endif else nfin=-1

  fcount=0

; read list of input files
; for each, read file in, check that grids have same number of elements,
;   and load into a common flux array

; load spectral data into arrays l, o, f
; loop to end of file or fcount=nfin (latter can't happen if multiflag=0)

  while (not eof(fi) and fcount ne nfin) do begin 
    readf,fi,line
    split1   = strsplit(line,' ',/extract)
    filename = split1[0]
    split2   = strsplit(filename,'.',/extract)
    suffix   = split2[n_elements(split2)-1]
    print,fcount,nfin,' ',filename
    checkchar = strmid(filename,0,1) 
    if (checkchar ne ';') then begin

;     check file suffix and read file

      if (suffix eq 'tbl') then sp=rd_sptbl(filename,hdr) else $
        sp=readfits(filename,hdr,/silent)

;     if input list contains orders to remove, then remove them

      if (n_elements(split1) gt 1) then for n=1,n_elements(split1)-1 do begin 
        if (strmid(split1[n],0,1) eq ';') then break else begin
          order=float(split1[n])
          sp=spcut(sp,order=order)
        endelse
      endfor

      if (fcount eq 0) then begin

;       first file 
;       load wavelength, flux, and order vectors
;       load header and savesp array from first spectrum

        l = reform(sp[lcol,*])
        f = reform(sp[fcol,*])
        o = reform(sp[ocol,*])
        h = hdr                ; NOTE - not used anymore
        savesp = sp

;       also build list of allowed orders

        omin=min(o) & omax=max(o) & norders=0
        for m=omin,omax do begin
          idx=where(o eq m, ocnt)
          if (ocnt gt 0) then begin
            if (norders eq 0) then orders=m else orders=[orders,m]
            norders ++
          endif
        endfor

      endif else begin

;     if second or later file in group, ensure that new data aligns with old
;     pad with "-99" where data are missing

        new_l = reform(sp[lcol,*]) ; not yet used - should check for alignment
        new_f = reform(sp[fcol,*])
        new_o = reform(sp[ocol,*])

        load_f = l-l - 99.0          ; initialize flux to load with -99 values

;       iterate through allowed orders, loading load_f when matches are made

        for m=0,norders-1 do begin
          idx=where(o eq orders[m])
          odx=where(new_o eq orders[m], ocnt)
          if (ocnt gt 0) then begin
            load_f[idx] = new_f[odx]
          endif 
        endfor

;       update the f array

        f=[[f],[load_f]]
      endelse

      fcount ++                ; fcount=no. of input files at END of loop
    endif
  endwhile

; determine normalization factors for each order and file

  normfact=fltarr(omax+1,fcount)

  if (normflag gt 0) then begin
    for j=0,fcount-1 do begin
      for m=omin,omax do begin
        idx=where(o eq m)
        if (max(idx) gt -1) then begin
          olam=l[idx]
          meanlam=0.5*(min(olam)+max(olam))
          dellam=0.333*(max(olam)-min(olam))
          l0 = meanlam-dellam
          l1 = meanlam+dellam
          
          savesp[fcol,*] = f[*,j]
          normfact[m,j] = spavg(savesp,l0,l1)
        endif 
      endfor 
    endfor
  endif else normfact[*,*]=1.0

if (keyword_set(diag) ne 0) then begin
  for j=0,fcount-1 do print,normfact[*,j],format=fmt
  print,'---------'
endif

; normalize normalization factors to max or average for each order

  for m=omin,omax do begin
    idx=where(normfact[m,*] gt -90)
    maxnorm=max(normfact[m,idx])
    avgnorm=mean(normfact[m,idx])
    if (normflag eq 1 and avgnorm gt 0) then $
      normfact[m,idx]=normfact[m,idx]/avgnorm
    if (normflag eq 2 and maxnorm gt 0) then $
      normfact[m,idx]=normfact[m,idx]/maxnorm
  endfor

; if slbonus set, then copy normalization factors for order=2 into order=3
;   because bonus order *should* behave as order 2 (but doesn't...)
; if llbonus set, then repeat for orders 5 and 6

  if (n_elements(slbonus) eq 1) then normfact[3,*]=normfact[2,*]
  if (n_elements(llbonus) eq 1) then normfact[6,*]=normfact[5,*]

  for j=0,fcount-1 do begin

;   if fixorder set, then use normalizations for that order for all orders

    if (keyword_set(fixorder) ne 0) then begin
      idx=where(orders eq fixorder, ocnt)
      if (ocnt gt 0) then begin 
        for m=0,omax do if (m ne fixorder and normfact[m,j] gt -90) then $
          normfact[m,j] = normfact[fixorder,j]
      endif else print,'Warning, fixorder not an input order, ignoring.'
    endif

    if (keyword_set(diag) ne 0) then print,normfact[*,j],format=fmt

    for m=0,omax do begin
      idx = where(o eq m)
      normval = normfact[m,j]
      if (max(idx) gt 0 and normval gt -90) then f[idx,j] = f[idx,j]/normval
    endfor
  endfor

; build new spectral data array
; if only one file begin "coadded" then old already copied into new

  if (fcount gt 1) then for j=0,n_elements(l)-1 do begin
    idx = where(f[j,*] gt -90, lcnt)
    savesp[fcol,j] = mean(f[j,idx])
    if (lcnt gt 1) then savesp[ecol,j] = stddev(f[j,idx])/sqrt(float(lcnt)) $
      else savesp[ecol,j] = 0
  endfor

; implement normalization if norm=3

  if (normflag eq 3) then savesp = spdivide(savesp,spmax(savesp))

; write result with header of first file

  wr_spfits,outfile,savesp,-1
  print,'Wrote ',outfile,' from ',fcount,' files.'

endfor

close,/all
END
