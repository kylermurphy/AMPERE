

function AMPERE_db_rotate, $


  aacgm_path = aacgm_path
  coef_year = coef_year

  ;radius of earth in km
  r_earth = 6378.1
  
  aacgm_set_path, aacgm_path
  aacgm_load_coef, coef_year




  return, 0
end
