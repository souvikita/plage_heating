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

device, GET_DECOMPOSED=originalIDCstate
device, DECOMPOSED=0
dpath_sot = '/archive1/hinode/sot/level1hao/2017/07/29/SP3D/20170729_235421/'; '/archive1/hinode/sot/level1hao/2017/03/04/SP3D/20170304_005005/'
dpath_iris_l2 ='/irisa/data/level2/2017/07/29/20170729_234952_3620010077/';  '/irisa/data/level2/2017/03/04/20170304_004113_3620106076/' 

;------------- READING IRIS SJI AND THEIR META DATA BELOW ----------;
read_iris_l2, dpath_iris_l2+'iris_l2_20170729_234952_3620010077_SJI_2832_t000.fits', idx_cont, dat_cont, /silent
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
	expanded_sot[time_rt,*] = congrid(temp_sot,round(size_temp_sot[1]*y_sampling_sot[time_rt]/dypos_cont[time_rt]),cubic=-0.5)
	expanded_iris_rast[time_rt,*] = congrid(temp_iris_rast,round(size_temp_iris_rast[1]*dypos_iris/dypos_cont[time_rt]),cubic=-0.5)

	print,string(13b)+'% finished: ',float(time_rt)*100./(n_elements(time_secs_iris_cont)),format='(a,f4.0,$)'

endfor
window,1

;!p.multi=[0,2,1]
;plot_image, reform(expanded_sot), min=6500,max=1e4
;cgplot, time_secs_iris_cont[index_sot_saved], crval_sot, color = 'dodger blue', label ='SOT', thick=2
;cgoplot, time_secs_iris_cont[index_sot_saved], crval_iris_rast, color='red', label ='IRIS rast', thick=2
cgplot, crval_sot, color = 'dodger blue', thick=2, xtitle='time indices', ytitle='solar X (arcsec)', psym=-1, /window
cgoplot,crval_iris_rast, color='red', thick=2, psym=-2, /addcmd
cgLegend, Title=['SOT rast', 'IRIS rast'], PSym=[-1,-2], linestyle=[1,1], Color=['dodger blue', 'red'], location=[65,-150], /data, /addcmd, length=0.0, /box
;
;!p.multi=0

end


; t_utc[index_iris_rast_saved[27-11]:index_iris_rast_saved[27+11]] ;IRIS raster 27th index is where the intersection takes place. We are considering 10 mins before and rough;ly after
;IDL>  t_utc[index_iris_rast_saved[27-11]],t_UTC[index_iris_rast_saved[27+11]]
;----> 2017-07-30T00:08:10.770
;----> 2017-07-30T00:32:34.460
;the IRIS rast xpos pointing changes from -44" to -10" during the above time interval

;t_obs[index_sot_saved[27-11]:index_sot_saved[27+11]] ## SOT times

;We will consider the SOT slit positions at times below. This is motivated by not changing unbderlying magnetic field analysis of HMI delta_B
;IDL> t_obs[index_sot_saved[27-11]],t_obs[index_sot_saved[27+11]]
;-----> 2017-07-30T00:08:09.526
;-----> 2017-07-30T00:32:33.982
;
;******* Beyond the times above the magnetic fields would surely change(?).
; Anyway, the SOT x_cen pointing changes from -86" to 31" during the above time interval
;y_sampling_sot = 0.317
;x_sampling_sot = 0.297140
;
;dypos_iris = 0.1663
;dxpos_iris = 0.388177

;sot_t_aligned = reform(dat_fe_cont[index_sot_saved[27-11]:index_sot_saved[27+11],*])
;iris_t_aligned = transpose(reform(win_cont[wav_idx_cont,*, index_iris_rast_saved[27-11]:index_iris_rast_saved[27+11]]))

;IDL> size_iris_t = size(iris_t_aligned)
;IDL> size_iris_t
;           2          89        1095           4       97455
;89 slits with 15s cadence

;IDL> size_sot_t = size(sot_t_aligned)
;IDL> size_sot_t
;'           2         384         512           4      196608
; 384 slits with ~3.2s cadence

; I resampled IRIS raster to SOT raster pixel scales below using cubic interpolation as shown below:

; iris_exp = congrid(iris_t_aligned, round(size_iris_t[1]*dxpos_iris/x_sampling_sot[0]),round(size_iris_t[2]*dypos_iris/y_sampling_sot[0]),cubic=-0.5)
; 
; IDL> help, iris_exp
; IRIS_EXP        FLOAT     = Array[116, 575]
; The above is of course re-sampled from IRIS rast --> SOT pixel scales.

; ***** Important:: Even though we have somewhat aligned in time because now the resampled IRIS and sot_t_aligned cover 24 mins duration, 
; it is important to remember that during these 24 minutes the IRIS raster pointing changes from -44" to -10" in solar X while SOT changes (massively) 
; from -86" to 31". I have to make sure that both the slits point at roughly the same spatial location during this time interval. For that:

; x_cen_iris_delta_t = xpos_rast[index_iris_rast_saved[27-11]:index_iris_rast_saved[27+11]] ; assume this is the master as this moves less
; x_cen_sot_delta_t = x_cen[index_sot_saved[27-11]:index_sot_saved[27+11]]

; cgplot, x_cen[index_sot_saved[27-11]:index_sot_saved[27+11]], color='red', thick=2
; cgoplot, xpos_rast[index_iris_rast_saved[27-11]:index_iris_rast_saved[27+11]], color='blue

; *** Check the above overplot and see what I was talking about *****
; Therefore, below is my attempt to spatially align it as well.
;

; target_x1 = x_cen_iris_delta_t[0]
; minimizing_x1 = abs(x_cen_sot_delta_t- target_x1)
; sot_x_idx_start = where(minimizing_x eq min(minimizing_x))
; now since the re-sampled IRIS raster array has the same pixel scale as the SOT, you can add simply add 115 to this index.
; sot_x_idx_end = sot_x_idx_start + 115

; *****So the final SOT roughly aligned array would be reform(sot_t_aligned[sot_x_idx_start:sot_x_idx_end,*])

; hard to find a nice match between the two continuum images but still trying to blink
; save_loc = '/sanhome/bose/HOP_307/'
; restore, save_loc+'bxbybz_29Jul2017_plage.sav', /v; restores: bx2, by2, bz2 all disambiguated
; blos_aligned = reform(bz2[index_sot_saved[27-11]:index_sot_saved[27+11],*])
; plot_image, reform(blos_aligned[141:141+115,*]), max=600, min=-600

;;; WHAT AN IRONY!!! AFTER SUCH DETAILED PROCEDURE, THE PART THAT'S MOST WELL ALIGNED IS THE LEAST INTERESTING REGION DEVOID OF MOST STRONG PLAGE. SOME NEGATIVE
;; POLARITY IF YOU CHECK BLOS_ALIGNED
;; FUCK!!!
;;
; NEXT I'll have perform a cross-correlation between the reform(sot_t_aligned[sot_x_idx_start:sot_x_idx_end,*]) IRIS raster. 
;; Not sure which two should I use. continuum ususally works best but still. Maybe I will use Mg II K core image and check the SOT

;; Anyway, starting again with Mg II k 2796.2 

; win_mgk =  d->getvar(7, /load)
; wav_mgk_index =where(abs(d->getlam(7) - 2796.20) eq min(abs(d->getlam(7) - 2796.20)))
; iris_mg_t_aligned = transpose(reform(win_mgk[wav_mgk_index,*,index_iris_rast_saved[27-11]:index_iris_rast_saved[27+11]]))
; iris_mg_exp = congrid(iris_mg_t_aligned, round(size_iris_t[1]*dxpos_iris/x_sampling_sot[0]),round(size_iris_t[2]*dypos_iris/y_sampling_sot[0]),cubic=-0.5)
;
; the crop is iris_mg_exp[*,40:40+511] for IRIS raster
; and the somewhat rough crop of SOT would be reform(blos_aligned[70:70+115,*])
; 
; iris_cont_to_align = iris_exp[*,40:40+511] ; the y offset is based on chi-by-eye of blos and mg ii K core
; sot_cont_to_align = reform(sot_t_aligned[60:70+125,*]); ------------""------------------------'
; blos_to_align = reform(blos_aligned[60:70+125,*])
; plot_image, shift_blos[0:115,*], min=-600, max=600

; now it is the time to go back to the two continuum rasters and perform a cross-correlation to align them perfectly. 

; shift = trk(sot_cont_to_align[0:115,100:250], iris_cont_to_align[*,100:250]) ; selecting a sub-field for the detailed cross-correlation. Ref is iris raster image
; shift_sot = frac_shift(sot_cont_to_align, shift[0], shift[1])
; plot_image, iris_mg_exp[*, 40:40+511], min=0,max=1.4e3
; 
