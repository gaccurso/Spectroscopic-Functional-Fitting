FUNCTION raw_data_iram_my_reduction, ny, nx, SN, file, obsid, graphname

READCOL, 'LOWMASS_catalog.dat', id, ra, dec, redshift, logM, SFR, r50, r90, D25, colour, format='L,F,F,F,F,F,F,F,F,F', comment='#'
AD = WHERE(id EQ obsid)
z = redshift[AD]
ra= ra[AD]
dec = dec[AD]
lambda_em = 157.74  
PI = 3.1416
lambda_obs= (1+z)*lambda_em
data = read_ascii(file, data_start=4,  DELIMITER=',')  ; start on the 4th line  ;up to data.field01[49,*]

!P.MULTI = [0, 5, 6]
DEVICE, RETAIN=3, DECOMPOSED=0
entry_device = !D.NAME & help, entry_device 
set_plot, 'ps' 
device, file=graphname
DEVICE, /COLOR

vel = fltarr(n_elements(data.field01[0,*]))   
for i = 0, n_elements(data.field01[0,*])-1 do begin
    vel[i] = (((data.field01[0,i]-lambda_obs)/lambda_obs)*3e+05)   ;centering the data on the CII line
endfor
  
nrows = n_elements(vel)
rmsnoise = fltarr(25)
widths = fltarr(25)
B = fltarr(25)
xsize = 0.12
ysize = 0.18
ymin = -0.2
x0 = [0.05, 0.8]
fwhm1=22.0
fwhm2=11.0
spaxelwidth = 9.4
loadct, 13
for i=0, 24 do begin
    data = read_ascii(file, data_start=4,  DELIMITER=',')
   
    nz = nrows - 1.0 - nx    
    first1 = max(where(vel[ny:nz] LT 0.0))
    second1 = min(where(vel[ny:nz] GT 0.0))
    lambda_estimate = max(data.field01[2*i+1,(first1-4):(second1+4)])
    afit = gaussfit(vel[ny:nz], data.field01[2*i+1,ny:nz], A1, NTERMS=6, estimates=[lambda_estimate, 0.00, 50, 0.0, 0.0, 0.0]) ;;;;this is causing a bug
    new_k = data.field01[2*i+1, *] - A1[3] - (A1[4]*vel[*]) - (A1[5]*vel[*]*vel[*])
    g = reform(new_k[ny:((nrows/2)-4)])
    j = reform(new_k[((nrows/2)+4):nz])   	
    noise = stddev([g, j])
   
    y0 = [x0[0]+xsize, x0[1]+ysize]
    ymax = max(data.field01[25,ny:nz])
    yrange = [ymin, ymax]
    xrange = [-2000.0, 2000.0]
    remainder = i mod 5.0
   
    if (remainder EQ 0.0) AND (i LT 19) then begin
    	divisor = i/5.0
    	position = [x0[0], (x0[1]-(divisor*ysize)), y0[0], (y0[1]-(divisor*ysize))]
    	plot, vel[ny:nz], new_k[ny:nz], ytitle='Line Flux Emission (Jy)', position=position, XTICKFORMAT="(A1)" , xrange=xrange, yrange=yrange, XTICKS=2, XMINOR=4  	;xtitle=TeXtoIDL("Velocity (km s^{-1})")
    endif else begin
    	if (i EQ 20) then begin
    		divisor = i/5.0
    		position = [x0[0], (x0[1]-(divisor*ysize)), y0[0], (y0[1]-(divisor*ysize))]
       		plot, vel[ny:nz], new_k[ny:nz], xtitle=TeXtoIDL("Velocity (km s^{-1})"), ytitle='Line Flux Emission (Jy)' , position = position, xrange=xrange, yrange=yrange, XTICKV=[-1000, 0, 1000], XTICKS=2, XMINOR=4  
    	endif else begin
  			if (i GT 20) then begin
  			    divisor = i/5.0
  				position = [(x0[0]+(remainder*xsize)), (x0[1]-(4.0*ysize)), (y0[0]+(remainder*xsize)), (y0[1]-(4.0*ysize))]
  				plot, vel[ny:nz], new_k[ny:nz], xtitle=TeXtoIDL("Velocity (km s^{-1})"), position = position, xrange=xrange, yrange=yrange, YTICKFORMAT="(A1)",  XTICKV=[-1000, 0, 1000], XTICKS=2, XMINOR=4   	;, ytitle='Line Flux Emission (Jy)' ,
  			endif else begin
 				divisor = i/5
 				position = [(x0[0]+(remainder*xsize)), (x0[1]-(divisor*ysize)), (y0[0]+(remainder*xsize)), (y0[1]-(divisor*ysize))] 			  						
  				plot, vel[ny:nz], new_k[ny:nz], position=position, XTICKFORMAT="(A1)", xrange=xrange, yrange=yrange, YTICKFORMAT="(A1)" , XTICKS=2, XMINOR=4      ;, xtitle=TeXtoIDL("Velocity (km s^{-1})"), ytitle='Line Flux Emission (Jy)'
  			endelse	
  		endelse
  	endelse
  	
    signal= max(new_k[ny:nz])
    rmsnoise[i] = sqrt((noise^2.0 + (mean([g, j]))^2.0)) 
    ;print, noise, signal
    
IF signal/noise gt SN THEN BEGIN
    lambda_max = max(new_k[ny:nz])
    afit = gaussfit(vel[ny:nz], new_k[ny:nz], A2, NTERMS=3, estimates=[lambda_max, 0.00, 50])
    ;PRINT, A2[0], A2[1], A2[2], "   FWHM=", 2.35482*A2[2], "   TOTAL FLUX=", A2[0]*A2[2]*(SQRT(2*PI)), "   S/N ratio =", signal/noise
    oplot, vel[ny:nz], afit, color=80
    widths[i] = 2.35482*A2[2]
    B[i] = A2[0]*A2[2]*(SQRT(2*PI))

ENDIF ELSE BEGIN
    ;PRINT, "No fitting due to no detection", "S/N ratio =", signal/noise 
ENDELSE
endfor 


A2= 21.4/2.3548				;the IRAM beam
A0= 1.0/(SQRT(2.0*pi)*A2)
n=100.0
x = (findgen(n) * (47.0/(n-1)) - 23.5) # replicate(1.0, n)
y = transpose(x)
nk =6.0

f1 = A0*exp((-(x)^2.0)/(2.0*A2*A2))*A0*exp((-(y)^2.0)/(2.0*A2*A2))
print, "mean for the whole gaussian=", mean(f1)
print, "total volume,", int_tabulated_2d(x,y,f1)

m = fltarr(4,6)

m[0,*]=[0, 9.4, 9.4, 18.8, 18.8, 18.8]
m[1,*]=[0, 0, 9.4, 0, 9.4, 18.8]
m[2,*]=[-4.7, 4.7, 4.7, 14.1, 14.1, 14.1]
m[3,*]=[-4.7, -4.7, 4.7, -4.7, 4.7, 14.1]

c = fltarr(nk)
U = fltarr(nk)

for i=0, nk-1 do begin 

	q = (findgen(n) * (9.4/(n-1)) + m[2,i] ) # replicate(1.0, n)
	p = transpose((findgen(n) * (9.4/(n-1)) + m[3,i] ) # replicate(1.0, n))
	f2 = A0*exp((-(q)^2)/(2.0*A2*A2))*A0*exp((-(p)^2.0)/(2.0*A2*A2))
	print, "w_i", mean(f2)/(A0*A0)
	c[i]= mean(f2)/(A0*A0)

endfor

U[0]= exp((-(0)^2)/(2*A2*A2))*exp((-(0)^2)/(2*A2*A2))
U[1]= exp((-(9.4)^2)/(2*A2*A2))*exp((-(0)^2)/(2*A2*A2))
U[2]= exp((-(9.4)^2)/(2*A2*A2))*exp((-(9.4)^2)/(2*A2*A2))
U[3]= exp((-(18.8)^2)/(2*A2*A2))*exp((-(0)^2)/(2*A2*A2))
U[4]= exp((-(18.8)^2)/(2*A2*A2))*exp((-(9.4)^2)/(2*A2*A2))
U[5]= exp((-(18.8)^2)/(2*A2*A2))*exp((-(18.8)^2)/(2*A2*A2))
 

P = fltarr(nx,1)
P = [c[5], c[4], c[3], c[4], c[5], c[4], c[2], c[1], c[2], c[4], c[3], c[1], c[0], c[1], c[3], c[4], c[2], c[1], c[2], c[4], c[5], c[4], c[3], c[4], c[5]]

L = fltarr(nx,1)
L = [U[5], U[4], U[3], U[4], U[5], U[4], U[2], U[1], U[2], U[4], U[3], U[1], U[0], U[1], U[3], U[4], U[2], U[1], U[2], U[4], U[5], U[4], U[3], U[4], U[5]]

;units = 1.90186e+03        ; a unit conversion factor to go from mJ mu_m to J km/s

print, "TOTAL FLUX SEEN IN IRAM=", (total(B*L))    ; No need to divide by total(P) as this gaussian is already normalised



SCII= total(B*L)
SCIICENTRAL = B[12]


zero1 = where(widths GT 0.0)
newwidths = widths[zero1]
newrmsnoise = rmsnoise[zero1]
newerrors = transpose((newrmsnoise[*]*newwidths[*])/(sqrt(newwidths[*]/119.7)))
newweightings = transpose(L[zero1])
totalerrors = total(newerrors*L)


rmserror = mean(newrmsnoise)
dispersion = mean(newwidths)
SCIIERROR = (rmserror*dispersion)/(sqrt(dispersion/119.7))
SCIIERRORCENTRAL = (rmsnoise[12]*widths[12])/(sqrt(widths[12]/119.7))

print, rmserror
print, " "
print, SCIIERROR
PRINT, " "

SCIIERROR = totalerrors

LCII= (1.04e-03)*(SCII)*vrestcii*red*D*D
LCIIERROR= ((1.04e-03)*(SCIIERROR)*vrestcii*red*D*D) + (0.3*LCII)     ;this 30% cii error comes from pogloghisch paper, due to cross talk of the spaxels.

LCIICENTRAL = (1.04e-03)*(SCIICENTRAL)*vrestcii*red*D*D
LCIICENTRALERROR = (1.04e-03)*(SCIIERRORCENTRAL)*vrestcii*red*D*D + (0.3*LCIICENTRAL)



A = fltarr(9)
A= [obsid, z, total(B), total(B*E), total(B*L), SCII,  LCII, LCIICENTRAL, LCIICENTRALERROR]
return, A

end




















