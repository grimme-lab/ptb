subroutine readline(line,floats,strings,logicals, cs,cf,cl)
   use iso_fortran_env, only: wp => real64,output_unit
   implicit none
   real(wp) :: floats(10)
   character(len=*),intent(in) :: line
   character(len=80)  :: strings(10)
   logical, dimension(10) :: logicals

   real(wp) :: num
   character(len=80) :: stmp,str
   character(len=1)  :: digit
   integer  :: i, ty, cs, cf, cl

   stmp=''
   cl = 1
   cs = 1
   cf = 1
   strings(:)=''
   do i=1,len(trim(line))
      digit=line(i:i)
      if(digit.ne.' '.and.digit.ne.char(9)) then  !should exclude tabstops and blanks, 9 is ascii code for tab
         stmp=trim(stmp)//trim(digit)
      elseif(stmp.ne.'')then
         call checktype(stmp,num,str,ty)      !get type of string, 0=number, 1=character
         if(ty.eq.0) then
            floats(cf)=num
            cf=cf+1
         elseif(ty.eq.1) then
            strings(cs)=str
            cs=cs+1
         elseif(ty.eq.2) then
            logicals(cl) = merge(.true., .false., int(num) .eq. 1) 
            cl = cl + 1
         else
            write(output_unit, &
               & '(''readline: problem in checktype, must abort'')')
            exit
         endif
         stmp=''
      endif
      if(i.eq.len(trim(line))) then  !special case: end of line
         call checktype(stmp,num,str,ty)
         if(ty.eq.0) then
            floats(cf)=num
            cf=cf+1
         elseif(ty.eq.1) then
            strings(cs)=str
            cs=cs+1
         elseif(ty.eq.2) then
            logicals(cl) = merge(.true., .false., int(num) .eq. 1) 
            cl = cl + 1
         else
            write(output_unit, &
               & '(''readline: problem in checktype, must abort'')')
            exit
         endif
         stmp=''
      endif
   enddo
   
   cs = cs - 1
   cf = cf - 1
   cl = cl - 1
end subroutine readline
