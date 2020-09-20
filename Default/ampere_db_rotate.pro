

function AMPERE_db_rotate, $
  ldate = ldate, $
  lfile = lfile, $
  d_dir = d_dir, $
  aacgm_path = aacgm_path, $
  coef_year = coef_year

  if keyword_set(d_dir) then d_dir = d_dir else d_dir = 'D:\data\AMPERE\rawdB\'
  if keyword_set(ldate) then begin
    ncdfname = d_dir+time_string(ldate, tformat='YYYY')+'\'+time_string(ldate, tformat='YYYYMMDD')+'Amp_invert.ncdf'      
  endif else if keyword_set(lfile) then begin
    ncdfname = lfile
  endif else begin
    print, 'Either a date or a file name must be specified for loading'
  endelse


  radius_earth = 6378.1
  ;string date for geopack recalc
  d_str = strmid(ncdfname,22,8,/reverse) 
  yy = double(time_string(d_str,tformat='YYYY'))
  mm = double(time_string(d_str,tformat='MM'))
  dd = double(time_string(d_str,tformat='DD'))

  read_ampere_db, ncdfname, t_dec_hr, sv_num, p_num, p_eci, b_eci, sv_q, ds
  ;radius of earth in km
  r_earth = 6378.1
  
  ;aacgm_set_path, aacgm_path
  ;aacgm_load_coef, coef_year

  geopack_recalc,yy,mm,dd,12,00,00,/date
  ;convert eci/gei pos to geographic cooridinates
  geopack_conv_coord, reform(p_eci[0,*])/1000.,reform(p_eci[1,*])/1000.,reform(p_eci[2,*])/1000., $
    x_geo,y_geo,z_geo, /from_gei, /to_geo
  geopack_sphcar, x_geo,y_geo,z_geo, r_geo, the_geo, phi_geo, /to_sphere
  lat_geo = 90-(the_geo*180/!pi)
  lon_geo = phi_geo*180/!pi
    
  ;convert eci/gei b to geo
  geopack_conv_coord,reform(b_eci[0,*]),reform(b_eci[1,*]),reform(b_eci[2,*]), $
    bx_geo,by_geo,bz_geo, /from_gei, /to_geo
  
  ;convert geo b from cartersian to spherical
  geopack_bcarsp, x_geo, y_geo, z_geo, bx_geo, by_geo, bz_geo, $
    br_geo, bth_geo, bphi_geo
  
  ;calculate aacgm b field
  aacgm_conv_vec,lat_geo,lon_geo,r_geo-radius_earth,bth_geo,bphi_geo, $
    mlat,mlon,mvth1,mvph1,mvth2,mvph2,err, /to_aacgm

  ;calculate cdf epoch and mlt for
  ; position
  hh = long(t_dec_hr)
  dm = (t_dec_hr-hh)*60.
  mn = long(dm)
  ds = (dm-mn)*60
  ss = long(ds)
  ms = long((ds-ss)*1000)

  cdf_epoch, ce, yy, mm, dd, hh, mn, ss, ms, /compute_epoch
  mlt = aacgm_mlt(ce, mlon)
  
  t_th = time_double(d_str) + t_dec_hr*3600.
  
  p_geo=transpose([[x_geo],[y_geo],[z_geo]])
  b_geo=transpose([[bx_geo],[by_geo],[bz_geo]])
  
  r_dat = {d_str:d_str,t_dhr:t_dec_hr, t_th:t_th, $
           sv_num:sv_num, p_num:p_num, sv_q:sv_q, ds:ds, $
           p_eci:p_eci, b_eci:b_eci, $
           p_geo:p_geo, b_geo:b_geo, $
           br_geo:br_geo, bth_geo:bth_geo, bz_geo:bz_geo, $
           lat_geo:lat_geo, lon_geo:lon_geo, $
           mlat_geo:mlat, mlon_geo:mlon, mlt:mlt, $
           mv_east1:mvph1, mv_north1:(-1*mvth1), $
           mv_east2:mvph2, mv_north2:(-1*mvth2)}
            
           
  return, r_dat
end


;MAIN
;
;read in db data
ncdfname = 'D:\data\AMPERE\rawdB\2017\20170831Amp_invert.ncdf'



a = AMPERE_db_rotate(ldate='20170831')


end
