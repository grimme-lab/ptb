
!cuts the at blanks and tabstops and returns all floats and strings in order of occurence
subroutine cutline(line,floats,strings)
implicit none
real*8 floats(*),num
character*128 line,str,stmp
character*80 strings(3)
character*1 digit
integer i,ty,cs,cf

stmp=''
cs=1
cf=1
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
   else
      write(*,*)'Problem in checktype, must abort'
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
   else
      write(*,*)'Problem in checktype, must abort'
      exit
   endif
   stmp=''
   endif
enddo
end


!this checks the type of the string and returns it cast to real or as string.
subroutine checktype(field,num,str,ty)
implicit none
character*(*) field,str
real*8 num
integer e,ty
logical is_num, is_logical

ty=99
str=''
is_num=.false.
read(field,'(F10.5)',IOSTAT=e)num !cast string on real and get error code; 0 means success.
if(e.eq.0)is_num=.true.
if(is_num)then
   if(index(field,'.').ne.0) then  !check for integer/real
   read(field,'(F30.16)')num
   ty=0
   else                       !if integer, add .0 to string; otherwise cast to real does not work
   str=trim(field)//'.0'
   read(str,'(F30.16)')num
   str=''
   ty=0
   endif
else
   if (trim(adjustl(field)) == 'true' .or. trim(adjustl(field)) == 'false') then
      is_logical = .true.
      if (field == '.TRUE.') then
            num = 1.0_8
      else
            num = 0.0_8
      end if
      ty = 2  ! Assign type 2 for logical
   else
   
      str=field
      ty=1
   endif
endif

end


