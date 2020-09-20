
;plot AMPERE current density

pro ampere_plot_j, $
  amp_d, $ AMPERE data from read_ampere_ncdf
  tstart, $ ;start time to plot in string or themis time double
  jmin=jmin, $
  jmax=jmax
  
  
  if keyword_set(jmin) then jmin eq jmin else jmin = 0.2
  if keyword_set(jmax) then jmax eq jmax else jmax = 1.0
  
  ;find plotting position
  yr = long(time_string(tstart, tformat='YYYY'))
  mo = long(time_string(tstart, tformat='MM'))
  dy = long(time_string(tstart, tformat='DD'))
  hr = long(time_string(tstart, tformat='hh'))
  mt = long(time_string(tstart, tformat='mm'))
  sc = long(time_string(tstart, tformat='ss'))


  gd = where(amp_d[*].start_yr eq yr and amp_d[*].start_mo eq mo and amp_d[*].start_dy eq dy $
    and amp_d[*].start_hr eq hr and amp_d[*].start_mt eq mt and amp_d[*].start_sc eq sc, pc)


end

;plot AMPERE fB vectors
; these are the fitted vectors
; from the AMPERE grd.ncdf files
; 
; need to add _extra key word

pro ampere_plot_fb, $
  amp_d, $ AMPERE data from read_ampere_ncdf
  tstart ;start time to plot in string or themis time double

  ;find plotting position
  yr = long(time_string(tstart, tformat='YYYY'))
  mo = long(time_string(tstart, tformat='MM'))
  dy = long(time_string(tstart, tformat='DD'))
  hr = long(time_string(tstart, tformat='hh'))
  mt = long(time_string(tstart, tformat='mm'))
  sc = long(time_string(tstart, tformat='ss'))

  
  gd = where(amp_d[*].start_yr eq yr and amp_d[*].start_mo eq mo and amp_d[*].start_dy eq dy $
    and amp_d[*].start_hr eq hr and amp_d[*].start_mt eq mt and amp_d[*].start_sc eq sc, pc)

  if pc gt 0 then begin
    colat_p = amp_d[gd].colat
    mlt_p   = amp_d[gd].mlt
    east    = amp_d[gd].dbeast1
    north   = amp_d[gd].dbnorth1
  endif else begin
    print, 'no data to plot'
    return
  endelse

  ;setup plotting region
  ;want option for extra key words here
  ampere_plots
  ;plot vectors
  ampere_vec_plot,colat_p, mlt_p, north, east
  ;plot mlt legend
  mlt_legend
end


;plot AMPERE dB vectors
; from the AMPERE invert.ncdf
; 
; need to add _extra key word
pro ampere_plot_db, $
  lat, $ ; latitude position
  mlt, $ ; mlt position
  e_vec, $ ; east db vectors
  n_vec, $ ; north db vectors
  t_vec, $ ; time of vectors
  tstart, $ ;start time to plot in string or themis time double
  colat = colat, $ ;if lat is colat set this word
  tdur = tdur ; duration in minutes

  ;get colatitude position
  if keyword_set(colat) then begin
    colat = lat
  endif else begin
    colat = 90-lat
  endelse
  ;get duration of for plotting vectors
  if keyword_set(tdur) then tdur=tdur else tdur=10
 
  ;find good data to plot
  gd = where(t_vec ge time_double(tstart) and t_vec le time_double(tstart)+tdur*60., pc)
  if pc gt 0 then begin
    colat_p = colat[gd]
    mlt_p   = mlt[gd]
    east    = e_vec[gd]
    north   = n_vec[gd]
  endif else begin
    print, 'no data to plot'
    return
  endelse
 
  ;setup plotting region
  ;want option for extra key words here
  ampere_plots
  ;plot vectors
  ampere_vec_plot,colat_p, mlt_p, north, east
  ;plot mlt legend
  mlt_legend
end

; plot contours for latitude axis
pro lat_axis, $
  lat_val, $ ; array of values to draw latitude contours
  lat_leg, $ ; position to draw latitude legend
  axis_ct = axis_ct, $
  axis_col = axis_col

  if keyword_set(axis_ct) then axis_ct = axis_ct else axis_ct = 0
  if keyword_set(axis_col) then axis_col = axis_col else axis_col = !p.color

  theta = findgen(360)*!dtor
  r = fltarr(360)
  loadct,axis_ct, /silent
  for j = 0L, n_elements(lat_val) - 1 do begin
    ;if lat_arr[j] gt xmax then continue
    r[*]  = long(reform(lat_val[j]))
    oplot, /polar, r, theta,linestyle = 1, color=axis_col

    xpos = r[0]*cos(lat_leg*!dtor)
    ypos = r[0]*sin(lat_leg*!dtor)
    xyouts, xpos, ypos, strtrim(long(90-r[0]),2)
    ;might put this in later
    ;    if hemisphere eq 'north' then begin
    ;      xyouts, xpos, ypos, strtrim(long(90-r[0]),2)
    ;    endif else begin
    ;      xyouts, xpos, ypos, strtrim(-1*long(90-r[0]),2)
    ;    endelse
  endfor
end


; plot contours for mlt axis
pro mlt_axis, $
  mlt_pos, $ ; array of values for MLT contours
  xval, $ ; max latitude for length of contours
  axis_ct = axis_ct, $ 
  axis_col = axis_col
  
  if keyword_set(axis_ct) then axis_ct = axis_ct else axis_ct = 0
  if keyword_set(axis_col) then axis_col = axis_col else axis_col = !p.color

  loadct,axis_ct, /silent
  for i=0L, n_elements(mlt_pos)-1 do begin
    pos_start = 0
    pos_end   = xval
    mlt_rad = (mlt_pos[i]*360./24.)*!dtor - !pi/2.
    xp = [pos_start*cos(mlt_rad),pos_end*cos(mlt_rad)]
    yp =[pos_start*sin(mlt_rad),pos_end*sin(mlt_rad)]
    oplot, xp, yp, linestyle=1, color=axis_col
  endfor
end

; add an mlt legend
pro mlt_legend 
  ;get max and min points from plot
  tr = convert_coord(!p.clip[2],!p.clip[3], /device,/data)
  bl = convert_coord(!p.clip[0],!p.clip[1], /device,/data)

  shift_y = (tr[1]-bl[1])/30.
  xyouts, 0, bl[1]-shift_y, '0 MLT', /data, alignment = 0.5
  xyouts, 0, tr[1]+shift_y/2., '12 MLT', /data, alignment = 0.5
  xyouts, bl[0], 0, '18 MLT ', /data, alignment = 1.0
  xyouts, tr[0], 0, ' 6 MLT', /data, alignment = 0.0
end


; plot ampere dB and fB vectors on a 
; at somepoint add color here
pro ampere_vec_plot, $
  v_colat, $
  v_mlt, $
  v_n, $
  v_e, $
  vsize = vsize
  
  ;vector size for normalization
  if keyword_set(vsize) then vsize=vsize else vsize=75.
  ;get max and min points from plot
  tr = convert_coord(!p.clip[2],!p.clip[3], /device,/data)
  bl = convert_coord(!p.clip[0],!p.clip[1], /device,/data)
  xmin = round(bl[0])
  xmax = round(tr[0])
  ymin = round(bl[1])
  ymax = round(tr[1])
  
  ;calculate vector start and end points
  ;start on the negative y axis and rotate
  ; everything based on MLT location
  ; north vectors
  nstart = -1*v_colat
  nend   = (nstart+v_n/vsize)
  ; east vectors
  estart = 0
  eend   = (estart+v_e/vsize)
  ; rotate vectors based on MLT
  ystart = sin(!dtor*v_mlt*360./24 )*estart + cos(!dtor*v_mlt*360./24 )*nstart
  yend   = sin(!dtor*v_mlt*360./24 )*eend + cos(!dtor*v_mlt*360./24 )*nend

  xstart = cos(!dtor*v_mlt*360./24 )*estart - sin(!dtor*v_mlt*360./24 )*nstart
  xend   = cos(!dtor*v_mlt*360./24 )*eend - sin(!dtor*v_mlt*360./24 )*nend

  ;plot vectors
  for i =0 , n_elements(v_n)-1 do begin
    if ystart[i] gt ymax or ystart[i] lt ymin then continue
    if xstart[i] gt xmax or xstart[i] lt xmin then continue
    if yend[i] gt ymax or yend[i] lt ymin then continue
    if xend[i] gt xmax or xend[i] lt xmin then continue

    plots, xstart[i],ystart[i], psym = sym(1), symsize = 0.5
    plots, xstart[i],ystart[i]
    plots, xend[i],yend[i], /continue
  endfor
  ;plot vector legend
  ;create a legend
  nstart = ymin+3
  nend   = nstart+0
  estart = 0+30
  eend   = estart+250/vsize

  loadct, 39, /silent
  plots, estart,nstart, psym = sym(1), symsize = 0.5
  plots, estart,nstart
  plots, eend,nend, /continue
  leg_pos = convert_coord(eend,nend,/data,/to_device)
  xyouts, leg_pos[0]+50./vsize,leg_pos[1],' 250 nT', /device
end

;setup AMPERE plotting area
pro ampere_plots
  
  lat_pos_leg = -35
  
  mlt_sec  = [19,20,21,22,23,0,1,2,3,4]
  mlt_arr = [0,3,6,9,12,15,18,21]

  ;radii for ploting constant lats
  lat_arr = [10,20,30,40,50]

  ;max plotting co-latitude
  ymax = 50
  ymin = -50
  xmax = 50
  xmin = -50
  
  loadct,0,/silent
  plot,[0,0],[0,0], /nodata, xrange=[xmin,xmax],yrange=[ymin,ymax], $
    /isotropic, xstyle = 1, ystyle = 1, xticks = 1, xminor =1, yticks = 1, yminor = 1, $
    xtickf='no_ticks', ytickf='no_ticks', thick = !p.thick

  theta = findgen(360)*!dtor
  r = dblarr(360)

  lat_axis, lat_arr, lat_pos_leg
  mlt_axis, mlt_arr, xmax
end


;main
fixplot

window, 0, xsize = 1450, ysize = 450
!p.multi = [0,3,1]
!x.margin = [15,15]

; get db data
db = AMPERE_db_rotate(ldate='20170831')
; get fittend and current data
read_ampere_ncdf,'D:\data\AMPERE\20170831.0000.86400.120.north.grd.ncdf', amp




ampere_plot_db,db.mlat_geo,db.mlt, db.mv_east1, db.mv_north1,db.t_th, pt
ampere_plot_fb, amp, pt

end

