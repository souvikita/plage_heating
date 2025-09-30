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

dpath_sot = '/archive1/hinode/sot/level1hao/2017/03/04/SP3D/20170304_005005/' ;'/archive1/hinode/sot/level1hao/2017/07/29/SP3D/20170729_235421/'
dpath_iris_l2 = '/irisa/data/level2/2017/03/04/20170304_004113_3620106076/' ;'/irisa/data/level2/2017/07/29/20170729_234952_3620010077/'

;------------- READING IRIS SJI AND THEIR META DATA BELOW ----------;
read_iris_l2, dpath_iris_l2+'iris_l2_20170304_004113_3620106076_SJI_2832_t000.fits', idx_cont, dat_cont, /silent
date_time_cont=idx_cont.date_obs
time_secs_iris_cont=time_iris_secs(date_time_cont)
xcen_cont = idx_cont.XCEN
ycen_cont = idx_cont.YCEN
exptime_cont = idx_cont.exptime
dxpos_cont = idx_cont.CDELT1
dypos_cont = idx_cont.CDELT2

;----- Reading IRIS raster file and their meeta data below for each step----
iris_rast_files = file_search(dpath_iris_l2+'*raster*')
d = iris_obj(iris_rast_files[0])
win_cont = d->getvar(5, /load); for the 2832 continuum
ypos_rast = d->getypos()
xpos_rast = d->getxpos()
dxpos_iris = mean(xpos_rast[1:*]-xpos_rast[0:3])
dypos_iris = mean(ypos_rast[1:*]-ypos_rast[0:3])
hdr = d->gethdr(iext, struct=struct)
t_UTC=d->ti2utc()
time_secs_iris=time_iris_secs(t_UTC)
wav_cont = d->getlam(5)
wav_idx_cont = where(abs(wav_cont - 2832.79) eq min(abs(wav_cont - 2832.79))) ;will correspond to the first index of the IRIS raster cube

;----- Reading SOT raster files and their meta data for each slit position----.
sot_raster_files = file_search(dpath_sot+'*fits*')
t_obs = strarr(n_elements(sot_raster_files))
y_sampling_sot = fltarr(n_elements(sot_raster_files))
x_fov = y_sampling_sot*0
y_fov = y_sampling_sot*0
x_cen = y_sampling_sot*0
y_cen = y_sampling_sot*0
x_sampling_sot = y_sampling_sot*0
dat_fe_cont = fltarr(n_elements(sot_raster_files),512)

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
	cont_index = 1 ; first wavelength position out of 112 wavelength indices. 
;	dat_fe_core[rast_idx,*] = reform(dat_sot[round(idx_sot.crpix1),*,0]/float(exp_time_sot)); the last index is for Stokes I,Q,U,V. Selecting only Stokes I
	dat_fe_cont[rast_idx,*] = reform(dat_sot[1,*,0]/float(exp_time_sot)) 
	print,string(13b)+'\n% finished: ',float(rast_idx)*100./(n_elements(sot_raster_files)),format='(a,f4.0,$)'

endfor
window, 0
plot_image, reform(dat_fe_cont),  min=6500,max=1e4
;stop

time_sot_secs = time_iris_secs(t_obs) ; converting the SOT time to seconds
expanded_sot = fltarr(n_elements(time_secs_iris_cont),round(512*y_sampling_sot[0]/dypos_cont[0]) ); 512 is the pixels along the SOT slit N-S direction
expanded_iris_rast = fltarr(n_elements(time_secs_iris_cont),round(1095*dypos_iris/dypos_cont[0])); 1095 is the pixels along the SOT slit N-S direction
crval_sot = fltarr(n_elements(time_secs_iris_cont))
crval_iris_rast = crval_sot*0
index_sot_saved = []
index_iris_rast_saved = []

for time_rt=0, n_elements(time_secs_iris_cont)-1 do begin
	target_time=time_secs_iris_cont[time_rt]
	minimizer_sot = abs(time_sot_secs-target_time)
	minimizer_iris_rast = abs(time_secs_iris-target_time)
	
	if n_elements(where(minimizer_sot eq min(minimizer_sot))) gt 1 or n_elements(where(minimizer_iris_rast eq min(minimizer_iris_rast))) then begin
		index = min(where(minimizer_sot eq min(minimizer_sot)))
		index_iris = min(where(minimizer_iris_rast eq min(minimizer_iris_rast)))
		temp_sot = reform(dat_fe_cont[index,*])
		temp_iris_rast = reform(win_cont[wav_idx_cont,*,index_iris])
		crval_sot[time_rt] = x_cen[index] ;SOT
		crval_iris_rast[time_rt] = xpos_rast[index_iris]; IRIS rast
		index_sot_saved = [index_sot_saved,index]
		index_iris_rast_saved = [index_iris_rast_saved,index_iris]
	endif else begin
		temp_sot = reform(dat_fe_cont[where(minimizer_sot eq min(minimizer_sot)),*])
		temp_iris_rast = reform(win_cont[wav_idx_cont,*,where(minimizer_iris_rast eq min(minimizer_iris_rast))])
		crval_sot[time_rt] = x_cen[where(minimizer_sot eq min(minimizer_sot))]; SOT
		crval_iris_rast[time_rt] = xpos_rast[where(minimizer_iris_rast eq min(minimizer_iris_rast))]; IRIS rast
		index_sot_saved = [index_sot_saved,index]
		index_iris_rast_saved = [index_iris_rast_saved,index_iris]

	endelse
	size_temp_sot = size(temp_sot)
	size_temp_iris_rast = size(temp_iris_rast)
	expanded_sot[time_rt,*] = congrid(temp_sot,round(size_temp_sot[1]*y_sampling_sot[time_rt]/dypos_cont[time_rt]),cubic=-0.5)
	expanded_iris_rast[time_rt,*] = congrid(temp_iris_rast,round(size_temp_iris_rast[1]*dypos_iris/dypos_cont[time_rt]),cubic=-0.5)

	print,string(13b)+'% finished: ',float(time_rt)*100./(n_elements(time_secs_iris_cont)),format='(a,f4.0,$)'

endfor
window,1

;!p.multi=[0,2,1]
;plot_image, reform(expanded_sot), min=6500,max=1e4
cgplot, time_secs_iris_cont[index_sot_saved], crval_sot, color = 'dodger blue', label ='SOT', thick=2
cgoplot, time_secs_iris_cont[index_iris_rast_saved], crval_iris_rast, color='red', label ='IRIS rast', thick=2
;
;!p.multi=0

end
