FUNCTION spcombine,a,b,weight,TEMPLATE=template,ORDER=order,RJ=rj,BKSPACE=bkspace,NEW1=new1,NEW2=new2,PLOT=plot,_extra=e

;  2 Apr 12 fixing bug with new1 and new2 if order keyword set
; 23 Feb 12 added new1 and new2 keywords to return interim smoothed spectra
;  3 Jan 11 added bkspace keyword
; 31 Dec 10 added template keyword
; 30 Jul 10 adjusting error calculation (had a bug when in RJ units)
; 27 Jul 10 added order keyword to limit operations to a single order
; 26 Jul 10 created
;
; uses bspline_iterfit to combine two spectra
;   the weight vector determines which of the spectra are used for each element
; based on /home/don/Reduction/nods/tgrid.pro
; currently assumed that they are on the same wavelength grid
; input:
;   a        - first input spectrum
;   b        - second input spectrum
;   weight   - code for each pixel - 0 = use both, 1 = use a, 2 = use b
;   template - optional keyword to force spline to input a (1) or b (2)
;   order    - optional parameter to limit order
;   rj       - optional keyword to make calculations in Rayleigh-Jeans space
;   new1     - optional keyword to return the 1st input spectrum renormalized
;   new2     - optional keyword to return the 2nd input spectrum renormalized
;   plot     - optional keyword to force plotting (expects colors to be defined)
;   bkspace  - optional keyword to pass to bspline_fit (and bspline_bkpts)
;              sets spacing of break points for spline fit in x
;              can be one element or have one element per order
; output:
;   returns combined spectrum

; set up col definitions, check keywords

lcol=0 & fcol=1 & ecol=2 & ocol=3
if (keyword_set(plot) eq 0) then plotflag=0 else begin
  if (plot lt 0) then plotflag=0 else plotflag=1
endelse
if (keyword_set(rj) eq 0) then rjflag=0 else rjflag=1
if (keyword_set(template) eq 0) then template=0
if (keyword_set(bkspace) eq 0) then bkspace=0.2

; check weight array

len=n_elements(a[lcol,*])
if (n_elements(weight) ne len) then weight=fltarr(len) ; set to zero

; check orders and order keyword

omin=min(a[ocol,*]) & omax=max(a[ocol,*])

if (keyword_set(order) ne 0) then begin
  if (order ge omin and order le omax) then begin 
    omin=order & omax=order 
  end
endif

; expand bkspace keyword if it has one element and there are multiple orders

if (omin ne omax and n_elements(bkspace) eq 1) then begin
  dummy=bkspace
  bkspace=fltarr(omax-omin+1)
  bkspace[*] = dummy
endif

out=a
new1=out
new2=out

; plot set-up

if (plotflag eq 1) then begin
  !p.font=1 & pth=3 & lth=1.5 & psz=0.85 ; psz=1.00 & !p.charsize=1.5
  xtt='!9l!4 (!9m!4m)'
  if (rjflag eq 0) then ytt='!4F!d!9n!n!4 (Jy)' $
    else ytt='!9l!4!u2!n F!d!9n!n!4 (Jy !9m!4m!u2!n)'
endif

for o=omin,omax do begin ; build output iterating through the orders

  idx=where(a[ocol,*] eq o)
  olen=n_elements(idx)
  out_f=fltarr(olen) & out_e=fltarr(olen)

  if (idx[0] ne -1) then begin

;   determine the difference spectrum and spline-fit it with bspline_iterfit

    diff=spadd(a[*,idx],b[*,idx],/minus)
    l=reform(diff[lcol,*]) & f=reform(diff[fcol,*])

    af=reform(a[fcol,idx]) & bf=reform(b[fcol,idx])
    ae=reform(a[ecol,idx]) & be=reform(b[ecol,idx])
    w=reform(weight[idx])

    if (rjflag eq 1) then begin ; shift to RJ units if requested
      f=f*l*l & af=af*l*l & bf=bf*l*l
    endif

    cc=bspline_iterfit(l,f,bkspace=bkspace[o-omin],yfit=diff_f,$
      maxiter=10,lower=1,upper=1)

;   set anew_f and bnew_f - on which the splines will be laid
;   if template set to 1 or 2, then af or bf serve as the template
;   otherwise the new template is the average

    case template of
      1 : begin
        anew_f = af
        bnew_f = bf + diff_f
      end
      2 : begin
        anew_f = af - diff_f
        bnew_f = bf
      end
      else : begin ; default, template = 0 or never set
        anew_f = af - diff_f/2.0
        bnew_f = bf + diff_f/2.0
      end
    endcase

;   generate the spline-shifted spectra for a and b - anew and bnew

    for j=0,olen-1 do begin
      case w[j] of
        1    : begin & out_f[j]=anew_f[j] & out_e[j]=ae[j] & end
        2    : begin & out_f[j]=bnew_f[j] & out_e[j]=be[j] & end
        else : begin ; default, wt = 0 or wt > 2 (which defaults here to 0)
                 out_f[j] = mean([anew_f[j],bnew_f[j]])
                 err1 = sqrt((ae[j]^2+be[j]^2))/2.0
                 err2 = stddev([anew_f[j],bnew_f[j]])/2.0
                 if (rjflag eq 1) then err2=err2/l[j]^2
                 out_e[j] = max([err1,err2])
;                 if (j eq 100) then print,l[j],anew_f[j],bnew_f[j],$
;                   stddev([anew_f[j],bnew_f[j]])/2.0,out_e[j]
               end
      endcase
    endfor

;   assemble the smoothed spectra anew_f and bnew_f

;    if (o eq omin) then begin
;      newflux_a = anew_f             & newflux_b = bnew_f 
;    endif else begin
;      newflux_a = [newflux_a,anew_f] & newflux_b = [newflux_b,bnew_f]
;    endelse

;   plotting

    if (plotflag eq 1) then begin
      plot,l,af,_extra=e,xtit=xtt,xth=pth,ytit=ytt,yth=pth,/nodata
      oplot,l,af,psym=-4,symsize=psz,th=lth,col=2
      oplot,l,bf,psym=-1,symsize=psz,th=lth,col=8
      oplot,l,anew_f,th=1.5*lth,col=2
      oplot,l,bnew_f,th=1.5*lth,col=8
    endif

    out[fcol,idx]=out_f
    out[ecol,idx]=out_e

;   set new1 and new2 in case they were defined in the call to spcombine

    new1[fcol,idx] = anew_f
    new2[fcol,idx] = bnew_f

;   if rjflag set, then return fluxes to F_nu units

    if (rjflag eq 1) then begin
      out[fcol,idx]  = out[fcol,idx]  / out[lcol,idx]^2
      new1[fcol,idx] = new1[fcol,idx] / new1[lcol,idx]^2
      new2[fcol,idx] = new2[fcol,idx] / new2[lcol,idx]^2
    endif

  endif
endfor

RETURN,out
END
