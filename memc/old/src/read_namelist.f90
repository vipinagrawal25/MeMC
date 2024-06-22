module readnamelist
    use, intrinsic :: iso_c_binding
    implicit none
    integer, parameter :: char_len = 128
    contains

subroutine convert_cstr_fstr(c_str, f_str)
    character(kind=c_char, len=1), dimension(char_len), intent(in) :: c_str
    character(len=char_len), intent(out) :: f_str
    integer :: i

   f_str = " "
   loop_string: do i=1, char_len
      if ( c_str(i) == c_null_char ) then
         exit loop_string
      else
         f_str(i:i) = c_str (i)
      end if
   end do loop_string
end subroutine

subroutine convert_fstr_cstr(f_str, c_str)
    character(kind=c_char, len=1), dimension(char_len), intent(out) :: c_str
    character(len=char_len), intent(in) :: f_str
    integer :: i

   loop_string: do i=1, char_len
      if (i<len(trim(f_str))+1)then
         c_str (i) = f_str(i:i)
         else
         c_str (i) = c_null_char
      endif

   end do loop_string

end subroutine


subroutine Membrane_listread(N, coef_bend, YY, radius, bdry_type, parafile) bind(c, name='Membrane_listread')
    real(kind=c_double) :: coef_bend, YY
    real(kind=c_double) :: radius
    integer(kind=c_int) :: N, bdry_type
    character(kind=c_char, len=1), dimension(char_len), intent(in) :: parafile
    character(len=char_len) :: f_fname

    namelist /Membrane/ N, coef_bend, YY, radius, bdry_type
    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100,file=f_fname,status='old')
    read(unit=100,nml=Membrane)
    close(unit=100)
end subroutine

subroutine Stick_listread(do_stick, pos_bot_wall, sigma, epsilon, &
        theta,  parafile) bind(c, name='Stick_listread')
    real(kind=c_double) :: pos_bot_wall, sigma, theta, epsilon 
    character(kind=c_char, len=1), dimension(char_len), intent(in) :: parafile
    character(len=char_len) :: f_fname
    logical(kind=c_bool) :: do_stick;


    namelist /Stick/ do_stick, pos_bot_wall, sigma, epsilon, theta

    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100,file=f_fname,status='old')
    read(unit=100,nml=Stick)
    close(unit=100)
end subroutine

subroutine MC_listread(algo, dfac, kbt, is_restart,&
 tot_mc_iter, dump_skip, parafile) bind(c, name='MC_listread')
 integer(kind=c_int) :: tot_mc_iter, dump_skip
 logical (kind=c_bool) :: is_restart
    real(kind=c_double) :: dfac, kbt
    character(kind=c_char, len=1), dimension(char_len) :: algo;
    character(kind=c_char, len=1), dimension(char_len), intent(in) :: parafile
    character(len=char_len) :: f_fname
    character(len=char_len) :: Mcalgo

    namelist /mcpara/ Mcalgo, dfac, kbt, is_restart, tot_mc_iter, dump_skip

    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100,file=f_fname,status='old')
    read(unit=100,nml=mcpara)
    close(unit=100)
    call convert_fstr_cstr(Mcalgo, algo)

end subroutine

subroutine Spcurv_listread(which_spcurv, minC, maxC, theta, &
    parafile) bind(c, name='Spcurv_listread')
    real(kind=c_double) :: minC, maxC, theta
    character(kind=c_char, len=1), dimension(char_len) :: which_spcurv
    character(kind=c_char, len=1), dimension(char_len), intent(in) ::  parafile
    character(len=char_len) :: f_fname
    character(len=char_len) :: spcurv_which

    namelist /spcurvpara/ spcurv_which, minC, maxC, theta
    
    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100,file=f_fname,status='old')
    read(unit=100,nml=spcurvpara)
    close(unit=100)
    call convert_fstr_cstr(spcurv_which, which_spcurv)
end subroutine


subroutine Activity_listread(which_act, minA, maxA, &
    parafile) bind(c, name='Activity_listread')
    real(kind=c_double) :: minA, maxA 
    character(kind=c_char, len=1), dimension(char_len) :: which_act
    character(kind=c_char, len=1), dimension(char_len), intent(in) ::  parafile
    character(len=char_len) :: f_fname
    character(len=char_len) :: act_which

    namelist /Actpara/ act_which, maxA, minA
    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100,file=f_fname,status='old')
    read(unit=100,nml=Actpara)
    close(unit=100)
    call convert_fstr_cstr(act_which, which_act)

end subroutine

subroutine Afm_listread(do_afm, tip_rad, tip_pos_z, sigma, epsilon, &
             parafile) bind(c, name='Afm_listread')
         real(c_double) :: tip_rad, tip_pos_z, sigma, epsilon
    character(kind=c_char, len=1), dimension(char_len), intent(in) :: parafile
    character(len=char_len) :: f_fname
    logical(kind=c_bool) :: do_afm

    namelist /afmpara/ do_afm, tip_rad, tip_pos_z, sigma, epsilon 
    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100,file=f_fname,status='old')
    read(unit=100,nml=afmpara)
    close(unit=100)

end subroutine


subroutine Volume_listread(do_volume, is_pressurized, coef_vol_exp, pressure, &
         parafile) bind(c, name='Volume_listread')

     real(c_double) :: coef_vol_exp, pressure
     logical(kind=c_bool) :: do_volume, is_pressurized 
     character(kind=c_char, len=1), dimension(char_len), intent(in) ::  parafile
     character(len=char_len) :: f_fname

     namelist /Volpara/ do_volume, is_pressurized, coef_vol_exp, pressure
    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100,file=f_fname,status='old')
    read(unit=100,nml=Volpara)
    close(unit=100)
end subroutine


subroutine Fluid_listread(is_fluid,  min_allowed_nbr, fluidize_every, fac_len_vert, &
    parafile) bind(c, name='Fluid_listread')

     logical(kind=c_bool) :: is_fluid
     integer(kind=c_int) :: min_allowed_nbr, fluidize_every
     real(kind=c_double) :: fac_len_vert
     character(kind=c_char, len=1), dimension(char_len), intent(in) ::  parafile
     character(len=char_len) :: f_fname

    namelist /fluidpara/ is_fluid, min_allowed_nbr, fluidize_every, fac_len_vert
    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100,file=f_fname,status='old')
    read(unit=100,nml=fluidpara)
    close(unit=100)
end subroutine

subroutine Area_listread(do_area, coef_area_exp, &
         parafile) bind(c, name='Area_listread')

     real(c_double) :: coef_area_exp
     logical(kind=c_bool) :: do_area
     character(kind=c_char, len=1), dimension(char_len), intent(in) ::  parafile
     character(len=char_len) :: f_fname

     namelist /Areapara/ do_area, coef_area_exp
    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100,file=f_fname,status='old')
    read(unit=100,nml=Areapara)
    close(unit=100)
end subroutine

subroutine Spring_listread(do_spring, icompute, nPole_eq_z, sPole_eq_z, &
         parafile) bind(c, name='Spring_listread')

       real(c_double) :: nPole_eq_z, sPole_eq_z
       integer(c_int) :: icompute
       logical(kind=c_bool) :: do_spring
       character(kind=c_char, len=1), dimension(char_len), intent(in) :: parafile
    character(len=char_len) :: f_fname

    namelist /springpara/ do_spring, icompute, nPole_eq_z, sPole_eq_z 
    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100,file=f_fname,status='old')
    read(unit=100,nml=springpara)
    close(unit=100)
    end subroutine
end module