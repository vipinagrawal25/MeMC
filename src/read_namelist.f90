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


subroutine Membrane_listread(N, coef_bend, YY, vol_exp, &
        sp_curv, pressure, radius, isfluid, parafile) bind(c, name='Membrane_listread')
    real(kind=c_double) :: coef_bend, YY, vol_exp
    real(kind=c_double) :: sp_curv, pressure, radius
    integer(kind=c_int) :: N
    character(kind=c_char, len=1), dimension(char_len) :: parafile
    character(len=char_len) :: f_fname
    logical :: isfluid

    namelist /Membrane/ N, coef_bend, YY, vol_exp, sp_curv, pressure, radius, isfluid

    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100,file=f_fname,status='old')
    read(unit=100,nml=Membrane)
    close(unit=100)
end subroutine

subroutine Ljbot_listread(pos_bot_wall, sigma, epsilon, &
        theta,  parafile) bind(c, name='Ljbot_listread')
    real(kind=c_double) :: pos_bot_wall, sigma, theta, epsilon 
    character(kind=c_char, len=1), dimension(char_len) :: parafile
    character(len=char_len) :: f_fname


    namelist /Ljbot/ pos_bot_wall, sigma, epsilon, theta

    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100,file=f_fname,status='old')
    read(unit=100,nml=Ljbot)
    close(unit=100)
end subroutine

subroutine MC_listread(algo, dfac, kbt, is_restart,&
 tot_mc_iter, dump_skip, parafile) bind(c, name='MC_listread')
 integer(kind=c_int) :: tot_mc_iter, dump_skip
 logical (kind=c_bool) :: is_restart
    real(kind=c_double) :: dfac, kbt 
    character(kind=c_char, len=1), dimension(char_len) :: algo, parafile
    character(len=char_len) :: f_fname
    character(len=char_len) :: Mcalgo


    namelist /mcpara/ Mcalgo, dfac, kbt, is_restart, tot_mc_iter, dump_skip

    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100,file=f_fname,status='old')
    read(unit=100,nml=mcpara)
    close(unit=100)
    call convert_fstr_cstr(Mcalgo, algo)

end subroutine


subroutine Activity_listread(which_act, minA, maxA, &
    parafile) bind(c, name='Activity_listread')
    real(kind=c_double) :: minA, maxA 
    character(kind=c_char, len=1), dimension(char_len) :: which_act
    character(kind=c_char, len=1), dimension(char_len) ::  parafile
    character(len=char_len) :: f_fname
    character(len=char_len) :: act_which


    namelist /Actpara/ act_which, maxA, minA 

    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100,file=f_fname,status='old')
    read(unit=100,nml=Actpara)
    close(unit=100)
    call convert_fstr_cstr(act_which, which_act)

end subroutine

subroutine afm_listread(tip_rad, tip_pos_z, sigma, epsilon, &
             parafile) bind(c, name='afm_listread')
         real(c_double) :: tip_rad, tip_pos_z, sigma, epsilon
    character(kind=c_char, len=1), dimension(char_len) :: parafile
    character(len=char_len) :: f_fname

    namelist /afmpara/ tip_rad, tip_pos_z, sigma, epsilon 
    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100,file=f_fname,status='old')
    read(unit=100,nml=afmpara)
    close(unit=100)

end subroutine


 subroutine spring_listread(icompute, nPole_eq_z, sPole_eq_z, &
         parafile) bind(c, name='spring_listread')

       real(c_double) :: nPole_eq_z, sPole_eq_z
       integer(c_int) :: icompute
       character(kind=c_char, len=1), dimension(char_len) :: which_act, parafile
    character(len=char_len) :: f_fname

    namelist /springpara/ icompute, nPole_eq_z, sPole_eq_z 
    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100,file=f_fname,status='old')
    read(unit=100,nml=springpara)
    close(unit=100)

     end subroutine
end module
