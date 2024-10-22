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

subroutine MC_listread(algo, dfac, kbt, is_restart, &
    tot_mc_iter, dump_skip, is_fluid, min_allowed_nbr, &
    fluidize_every, fac_len_vertices, parafile) bind(c, name='MC_listread')

    use, intrinsic :: iso_c_binding
    implicit none

    integer(kind=c_int) :: tot_mc_iter
    integer(kind=c_int) :: dump_skip, min_allowed_nbr
    integer(kind=c_int) :: fluidize_every
    logical(kind=c_bool) :: is_restart, is_fluid
    real(kind=c_double) :: dfac, kbt, fac_len_vertices
    character(kind=c_char, len=1), dimension(char_len), intent(out) :: algo
    character(kind=c_char, len=1), dimension(char_len), intent(in) :: parafile
    character(len=char_len) :: f_fname
    character(len=char_len) :: Mcalgo

    namelist /mcpara/ Mcalgo, dfac, kbt, is_restart, tot_mc_iter, dump_skip, &
    is_fluid, min_allowed_nbr, fluidize_every, fac_len_vertices

    call convert_cstr_fstr(parafile, f_fname)
    open(unit=200, file=f_fname, status='old')
    read(unit=200, nml=mcpara)
    close(unit=200)

    ! Convert Mcalgo to algo
    call convert_fstr_cstr(Mcalgo, algo)

end subroutine

subroutine ElectroRead(charge1, charge2, conc, parafile) bind(c, name="ElectroRead")
  real (kind=c_double) :: charge1, charge2, conc
  character(kind=c_char, len=1), dimension(char_len), intent(in) ::  parafile
  character(len=char_len) :: f_fname

  namelist /electrostatpara/ charge1, charge2, conc
    call convert_cstr_fstr(parafile, f_fname)
    open(unit=200,file=f_fname,status='old')
    read(unit=200,nml=electrostatpara)
    close(unit=200)
end subroutine

subroutine SelfAvoidRead(isselfrepulsive, boxsize, cutoff, minl, sig, eps, parafile) bind(c, name='SelfAvoidRead')
    real(kind=c_double) :: boxsize, cutoff, sig, eps, minl
    logical(c_bool) :: isselfrepulsive
    character(kind=c_char, len=1), dimension(char_len), intent(in) :: parafile
    character(len=char_len) :: f_fname

    namelist /selfavoidpara/ isselfrepulsive, boxsize, cutoff, minl, sig, eps
        call convert_cstr_fstr(parafile, f_fname)
        open(unit=200,file=f_fname,status='old')
        read(unit=200,nml=selfavoidpara)
        close(unit=200)
end subroutine

subroutine BendRead(bend1, bend2, spC1, spC2, parafile) bind(c, name="BendRead")
    real (kind=c_double) :: bend1, bend2, spC1, spC2
    character(kind=c_char, len=1), dimension(char_len), intent(in) ::  parafile
    character(len=char_len) :: f_fname

    namelist /bendpara/ bend1, bend2, spC1, spC2 
        call convert_cstr_fstr(parafile, f_fname)
        open(unit=200,file=f_fname,status='old')
        read(unit=200,nml=bendpara)
        close(unit=200)
end subroutine

subroutine MeshRead(bdry_cdt, nghst, radius, ncomp, compfrac, parafile) bind(c, name="MeshRead")
    integer(kind=c_int) :: bdry_cdt, nghst, ncomp
    real(kind = c_double) :: radius, compfrac
    character(kind=c_char, len=1), dimension(char_len), intent(in) :: parafile
    character(len=char_len) :: f_fname

    namelist /meshpara/ bdry_cdt, radius, nghst, ncomp, compfrac
        call convert_cstr_fstr(parafile, f_fname)
        open(unit=200,file=f_fname,status='old')
        read(unit=200,nml=meshpara)
        close(unit=200)
end subroutine

subroutine LipidRead(iregsoln, kai, epssqby2, parafile) bind(c, name="LipidRead")
    real(kind=c_double) :: iregsoln, kai, epssqby2
    character(kind=c_char, len=1), dimension(char_len), intent(in) ::  parafile
    character(len=char_len) :: f_fname

    namelist /lipidpara/ iregsoln, kai, epssqby2
    call convert_cstr_fstr(parafile, f_fname)
    open(unit=200, file=f_fname, status='old')
    read(unit=200, nml=lipidpara)
    close(unit=200)
end subroutine


subroutine StretchRead(YY1, YY2, do_volume, is_pressurized, coef_vol_expansion, &
               pressure, coef_area_expansion, do_area, parafile) bind(c, name="StretchRead")
    logical (kind=c_bool) :: do_volume, is_pressurized, do_area
    real (kind=c_double) :: YY1, YY2, pressure, coef_area_expansion, coef_vol_expansion
    character(kind=c_char, len=1), dimension(char_len), intent(in) ::  parafile
    character(len=char_len) :: f_fname

    namelist /stretchpara/ YY1, YY2, do_volume, is_pressurized, coef_vol_expansion, &
                         pressure, coef_area_expansion, do_area 
        call convert_cstr_fstr(parafile, f_fname)
        open(unit=200,file=f_fname, status='old')
        read(unit=200,nml=stretchpara)
        close(unit=200)
end subroutine

end module