;Aim: to co-align SOT/SP rasters to IRIS 2832 continuum raster. However, what if each hinode and IRIS has multiple rasters?
;date: Sep. 30, 2025

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

device, GET_DECOMPOSED=originalIDCstate
device, DECOMPOSED=0
dpath_sot =  '/archive1/hinode/sot/level1hao/2017/03/05/SP3D/20170305_164005/';'/archive1/hinode/sot/level1hao/2017/03/04/SP3D/20170304_005005/'
dpath_iris_l2 = '/irisa/data/level2/2017/03/05/20170305_164021_3620106076/' ;/irisa/data/level2/2017/03/04/20170304_004113_3620106076/'

;------------- READING IRIS SJI AND THEIR META DATA BELOW ----------;
read_iris_l2, dpath_iris_l2+'iris_l2_20170305_164021_3620106076_SJI_2832_t000.fits', idx_cont, dat_cont, /silent
date_time_cont=idx_cont.date_obs
time_secs_iris_cont=time_iris_secs(date_time_cont)
xcen_cont = idx_cont.XCEN
ycen_cont = idx_cont.YCEN
exptime_cont = idx_cont.exptime
dxpos_cont = idx_cont.CDELT1
dypos_cont = idx_cont.CDELT2

;----- Reading IRIS raster file and their meta data below for each step----
iris_rast_files = file_search(dpath_iris_l2+'*raster*')
d = iris_obj(iris_rast_files[0])
win_cont = d->getvar(5, /load); for the 2832 continuum
ypos_rast = d->getypos()
xpos_rast = d->getxpos()
dxpos_iris = mean(xpos_rast[1:*]-xpos_rast[0:3])
dypos_iris = mean(ypos_rast[1:*]-ypos_rast[0:3])
;hdr = d->gethdr(iext, struct=struct)
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

loadct,0
plot_image, reform(dat_fe_cont),  min=4000,max=9000
!p.multi=0
dim_orig_fe_cont = size(dat_fe_cont)
time_sot_secs = time_iris_secs(t_obs) ; converting the SOT time to seconds
expanded_sot = fltarr(n_elements(time_secs_iris_cont),round(dim_orig_fe_cont[2]*y_sampling_sot[0]/dypos_cont[0]) ); 512 is the pixels along the SOT slit N-S direction
;expanded_iris_rast = fltarr(n_elements(time_secs_iris_cont),round(1095*dypos_iris/dypos_cont[0])); 1095 is the pixels along the SOT slit N-S direction
crval_sot = fltarr(n_elements(time_secs_iris))
crval_iris_rast = crval_sot*0
index_sot_saved = []
index_iris_rast_saved = []

	for time_rt =0, n_elements(time_secs_iris)-1 do begin
		target_time=time_secs_iris[time_rt]
		minimizer_sot = abs(time_sot_secs-target_time)
		;if time_secs_iris[0] gt time_sot_secs[time_rt] then begin
	;		print, ' Exiting the loop'
	;		break
	;	endif
		minimizer_iris_rast = abs(time_secs_iris-target_time)
		if n_elements(where(minimizer_sot eq min(minimizer_sot))) gt 1 or n_elements(where(minimizer_iris_rast eq min(minimizer_iris_rast))) gt 1 then begin
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
			index_sot_saved = [index_sot_saved,where(minimizer_sot eq min(minimizer_sot))]
			index_iris_rast_saved = [index_iris_rast_saved,where(minimizer_iris_rast eq min(minimizer_iris_rast))]

		endelse
		size_temp_sot = size(temp_sot)
		size_temp_iris_rast = size(temp_iris_rast)
		print,string(13b)+'% finished: ',float(time_rt)*100./(n_elements(time_secs_iris)),format='(a,f4.0,$)'


	endfor

cgplot, time_secs_iris[index_sot_saved], crval_sot, color = 'dodger blue', thick=2, xtitle='time [secs]',ytitle='solar X (arcsec)',/window
cgoplot, time_secs_iris[index_iris_rast_saved], crval_iris_rast, color='red', thick=2, /addcmd

;stop
; judge by the above plot. It looks like we can consider all the steps of the IRIS data for the first raster.

iris_t_aligned = transpose(reform(win_cont[wav_idx_cont,*,*]))
sot_t_aligned = reform(dat_fe_cont)
size_iris_t = size(iris_t_aligned)
size_sot_t = size(sot_t_aligned)

iris_exp = congrid(iris_t_aligned, round(size_iris_t[1]*dxpos_iris/x_sampling_sot[0]),round(size_iris_t[2]*dypos_iris/y_sampling_sot[0]),cubic=-0.5)

window, 1
!p.multi=0
plot_image,iris_exp, min=0, max=800 

size_iris_exp = size(iris_exp)

sot_x_zoom = [70,170]
sot_y_zoom = [230,330]

iris_x_zoom = [20,120]
iris_y_zoom = [150,250]

shifts = trk(sot_t_aligned[70:170,230:330], iris_exp[20:120,150:250])
sot_shifted = frac_shift(sot_t_aligned, shifts[0], shifts[1])

window, 2
loadct,3
result1 = correl_images(sot_t_aligned[70:170,230:330],iris_exp[20:120,150:250],xshift=8,yshift=8)
result2 = correl_images(sot_shifted[70:170,230:300],iris_exp[20:120,150:250],xshift=8,yshift=8)

!p.multi=[0,2,1]

plot_image, result1, title='before cross correlation'

plot_image, result2, title='after cross correlation'
loadct,0
!p.multi=0

final_sot_x_crop = [50,50+376]
final_sot_y_crop = [80,80+406]


wset,0

!p.multi=0
plot_image, sot_shifted[final_sot_x_crop[0]:final_sot_x_crop[1],final_sot_y_crop[0]:final_sot_y_crop[1]]

sot_final_shifted_for_paper = sot_shifted[final_sot_x_crop[0]:final_sot_x_crop[1],final_sot_y_crop[0]:final_sot_y_crop[1]]
blink, [0,1]


print, ' Obtaining the disambiguated Bx, By and Bz from the L2.1 files....'

scanid = '20170305_164005' ; this is from the folder that contains the data
sotsp_getdata,scanid,level=2, outdir ='/sanhome/bose/HOP_307/'  ; downloads both level 2 and 2.1 data files in the outdir folder

flist=file_search(scanid+'/'+scanid+'.fits');  gets name of level 2 fits file
read_sotsp,flist[0],index,data,scan_info=scan_info,/xycoord

;data will contain (Nx,Ny,42) dimensions, where each of the 42 planes in the array corresponds to one of the physical quantities of interest.
;These physical quantities are identified by name in the ftype tag in the index structure array.
ftypes = index.ftype ; will be a 42 element array

Bx = reform(data[*,*,-4]) ;Magnetic Field Zonal Component (positive west)
By = reform(data[*,*,-3]) ;Magnetic Field Meridional Component (positive north)
Bz = reform(data[*,*,-2]) ;Magnetic Field Radial Component (positive up)

print, ' saving the necessary parameters and arrays for post processing.'

x_sampling_sot = x_sampling_sot[0]
y_sampling_sot = y_sampling_sot[0] 
save, Bx, By, Bz, filename='/sanhome/bose/HOP_307/'+scanid+'/bxbybz_20170305_164005.sav' 
save, final_sot_x_crop,  final_sot_y_crop, shifts, $
dypos_iris,dxpos_iris, x_sampling_sot, y_sampling_sot, filename='/sanhome/bose/HOP_307/'+scanid+'/alignment_params.sav'
save, result1, result2, sot_final_shifted_for_paper, iris_exp,  filename='/sanhome/bose/HOP_307/'+scanid+'/cross_correl_sot_iris_aligned_maps.sav'
save, t_UTC, t_obs, crval_sot, crval_iris_rast, index_sot_saved, filename = '/sanhome/bose/HOP_307/'+scanid+'/times_and_xcen.sav'

end
