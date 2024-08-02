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

subroutine BendRead(coef_bend, minC, maxC, theta, spcurv, parafile) bind(c, name="BendRead")
  real (kind=c_double) :: coef_bend, minC, maxC, theta, spcurv
  character(kind=c_char, len=1), dimension(char_len), intent(in) ::  parafile
  character(len=char_len) :: f_fname

  namelist /Bendpara/ coef_bend, minC, maxC, theta, spcurv 
    call convert_cstr_fstr(parafile, f_fname)
    open(unit=200,file=f_fname,status='old')
    read(unit=200,nml=Bendpara)
    close(unit=200)
end subroutine

subroutine MeshRead(bdry_cdt, nghst, radius, parafile) bind(c, name="MeshRead")
  integer (kind=c_int) :: bdry_cdt, nghst
  real(kind = c_double) :: radius
  character(kind=c_char, len=1), dimension(char_len), intent(in) :: parafile
  character(len=char_len) :: f_fname

  namelist /Meshpara/ bdry_cdt, radius, nghst
    call convert_cstr_fstr(parafile, f_fname)
    open(unit=200,file=f_fname,status='old')
    read(unit=200,nml=Meshpara)
    close(unit=200)
end subroutine

subroutine StretchRead(YY, do_volume, is_pressurized, coef_vol_expansion, &
               pressure, coef_area_expansion, do_area, parafile) bind(c, name="StretchRead")
  logical (kind=c_bool) :: do_volume, is_pressurized, do_area
  real (kind=c_double) :: YY, pressure, coef_area_expansion, coef_vol_expansion
  character(kind=c_char, len=1), dimension(char_len), intent(in) ::  parafile
  character(len=char_len) :: f_fname

  namelist /Stretchpara/ YY, do_volume, is_pressurized, coef_vol_expansion, &
                         pressure, coef_area_expansion, do_area 
    call convert_cstr_fstr(parafile, f_fname)
    open(unit=200,file=f_fname,status='old')
    read(unit=200,nml=Stretchpara)
    close(unit=200)
end subroutine

subroutine MC_listread(algo, dfac, kbt, is_restart,&
 tot_mc_iter, dump_skip, is_fluid, min_allowed_nbr, &
                fluidize_every, fac_len_vertices, parafile) bind(c, name='MC_listread')
 integer(kind=c_int) :: tot_mc_iter, dump_skip, min_allowed_nbr, fluidize_every
 logical (kind=c_bool) :: is_restart, is_fluid
    real(kind=c_double) :: dfac, kbt, fac_len_vertices
    character(kind=c_char, len=1), dimension(char_len) :: algo;
    character(kind=c_char, len=1), dimension(char_len), intent(in) :: parafile
    character(len=char_len) :: f_fname
    character(len=char_len) :: Mcalgo

    namelist /mcpara/ Mcalgo, dfac, kbt, is_restart, tot_mc_iter, dump_skip, &
         is_fluid, min_allowed_nbr, fluidize_every, &
         fac_len_vertices

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

subroutine LipidRead(ncomp, &
        parafile) bind(c, name='Lipid_listread')
    integer(kind=c_int) :: ncomp
    character(kind=c_char, len=1), dimension(char_len), intent(in) ::  parafile
    character(len=char_len) :: f_fname

    namelist /lipidpara/ ncomp
    call convert_cstr_fstr(parafile, f_fname)
    open(unit=100, file=f_fname, status='old')
    read(unit=100, nml=lipidpara)
    close(unit=100)
end subroutine

end module