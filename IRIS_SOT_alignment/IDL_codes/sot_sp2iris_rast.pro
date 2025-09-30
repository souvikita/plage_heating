;Aim: to co-align SOT/SP rasters to IRIS 2832 continuum raster
;date: Sep. 5, 2025

;######################## IRIS time to seconds conversion #################3
function time_iris_secs,time_str
  time_secs=fltarr(n_elements(time_str))
  for i=0,n_elements(time_str)-1 do begin

   time = strmid(time_str[i],11,8)
   str1=strsplit(time,':',/extract)
   time_secs[i]=float(str1[0])*60*60.+float(str1[1])*60+float(str1[2])
  endfor
return, time_secs
end
;##################################################################

dpath_sot = '/archive1/hinode/sot/level1hao/2017/07/29/SP3D/20170729_235421/'
dpath_iris_l2 = '/irisa/data/level2/2017/07/29/20170729_234952_3620010077/'

;----- Reading IRIS raster file and their meeta data below for each step----
iris_rast_files = file_search(dpath_iris_l2+'*raster*')
d = iris_obj(rast_files[0])
win_cont = d->getvar(5, /load); for the 2832 continuum
ypos_rast = d->getypos()
xpos_rast = d->getxpos()
dxpos_iris = mean(xpos_rast[1:*]-xpos_rast[0:3])
dypos_iris = mean(ypos_rast[1:*]-ypos_rast[0:3])
hdr = d->gethdr(iext, struct=struct)
t_UTC=d->ti2utc()
time_secs_iris=time_iris_secs(t_UTC)

;----- Reading SOT raster files and their meta data for each slit position----.
sot_raster_files = file_search(dpath_sot+'*fits*')
t_obs = strarr(n_elements(sot_raster_files))
y_sampling_sot = fltarr(n_elements(sot_raster_files))
x_fov = y_sampling_sot*0
y_fov = y_sampling_sot*0
x_cen = y_sampling_sot*0
y_cen = y_sampling_sot*0
x_sampling_sot = y_sampling_sot*0
dat_fe_core = fltarr(n_elements(sot_raster_files),512)

for rast_idx=0, n_elements(sot_raster_files)-1 do begin
	read_sot, sot_raster_files[rast_idx], idx_sot, dat_sot,/silent
	exp_time_sot = idx_sot.exptime
	t_obs[rast_idx] = idx_sot.date_obs ;reading the start time. You could also consider the end time. "String"
	y_sampling_sot[rast_idx]= idx_sot.CDELT2; pixel scale in the y-direction
	x_sampling_sot[rast_idx] = idx_sot.xscale ;arcsec/step in the slit scanning direction
	y_fov[rast_idx] = idx_sot.FOVY
	y_cen[rast_idx] = idx_sot.YCEN
	x_cen[rast_idx] = idx_sot.XCEN
	x_fov[rast_idx] = idx_sot.FOVX
	dat_fe_core[rast_idx,*] = reform(dat_sot[round(idx_sot.crpix1),*,3]/float(exp_time_sot)); the last index is for Stokes I,Q,U,V. Selecting only Stokes I
	print,string(13b)+'% finished: ',float(rast_idx)*100./(n_elements(sot_raster_files)),format='(a,f4.0,$)'

endfor
expanded_sot = fltarr(n_elements(t_UTC),round(512*y_sampling_sot[0]/dypos_iris) )

time_sot_secs = time_iris_secs(t_obs) ; converting the SOT time to seconds

print,'Now entering the IRIS time loop'
for time_rt =0, n_elements(t_UTC)-1 do begin
;	print,'Now entering the IRIS time loop'
	target_time=time_secs_iris[time_rt]
	minimizer = abs(time_sot_secs-target_time)
	if n_elements(where(minimizer eq min(minimizer))) gt 1 then begin
		index = min(where(minimizer eq min(minimizer)))
		temp_sot = reform(dat_fe_core[index,*])
		crval1 = x_cen[index]
	endif else begin
		temp_sot = reform(dat_fe_core[where(minimizer eq min(minimizer)),*])
		crval1 = x_cen[where(minimizer eq min(minimizer))]
	endelse
	expanded_sot[time_rt,*] = congrid(temp_sot,round(512*y_sampling_sot[time_rt]/dypos_iris),cubic=-0.5)
	delta_solar_x = crval1 - xpos_rast[time_rt]
	delta_pixel_x = x_sampling_sot[time_rt] ; the sampling does not change with time. So take any!
	print,string(13b)+'% finished: ',float(time_rt)*100./(n_elements(t_UTC)),format='(a,f4.0,$)'


endfor
		
end



