program parascript
    implicit none
    integer :: myunit,i,j,elmtnum,ioerror
    integer :: carg,period,group,despar
    integer :: elemrange(6),arraysize
    logical :: per,gr,par,parname,help,fullparn
    character(len=128) :: username,pathname,atmp,parn
    character(len=8) :: x1
    real*8 :: efield
    real*8 :: globpar(10,2),allpar(150,86)
    real*8 :: expscal1(10,86),ener_par1(10,86),ener_par2(10,86), &
        floats(1:10), & 
        shell_xi(10,86),  &  
        shell_cnf1(10,86),  &
        shell_cnf2(10,86),  &
        shell_cnf3(10,86),  &  
        expscal3(10,86),  &  
        shell_cnf4(10,86),  &  
        shell_resp1(10,86),  & 
        shell_resp2(10,86),  &
        paramset(10,86)
   
    expscal1 = 0.0d0
    ener_par1 = 0.0d0
    ener_par2 = 0.0d0
    floats = 0.0d0
    shell_xi = 0.0d0
    shell_cnf1 = 0.0d0
    shell_cnf2 = 0.0d0
    shell_cnf3 = 0.0d0
    expscal3 = 0.0d0
    shell_cnf4 = 0.0d0
    shell_resp1 = 0.0d0
    shell_resp2 = 0.0d0
    paramset = 0.0d0

    per=.false.
    gr=.false.
    par=.false.
    parname=.false.
    fullparn=.false.
    help=.false.
    carg = command_argument_count()
    do i=1,carg
        call getarg(i,atmp)
        if(index(atmp,'-per').ne.0) then
            per=.true.
            call getarg(i+1,atmp)
            read(atmp,*) period
        endif
        if(index(atmp,'-g').ne.0) then
            gr=.true.
            call getarg(i+1,atmp)
            read(atmp,*) group
        endif
        if(index(atmp,'-par').ne.0) then
            par=.true.
            call getarg(i+1,atmp)
            read(atmp,*) despar
        endif
        if(index(atmp,'-name').ne.0) then
            parname=.true.
            call getarg(i+1,atmp)
            read(atmp,*) parn
        endif
        if(index(atmp,'-nfull').ne.0) then
            fullparn=.true.
            call getarg(i+1,atmp)
            read(atmp,*) parn
        endif
        if(index(atmp,'-help').ne.0) then
            help=.true.
        endif
    enddo
    
    if (parname.and.fullparn) then
        write(*,*) "Choose either -name OR -nfull"
        stop
    endif

    if (help) then
        write(*,'(a)') "Possible parameter set names after '-name' or '-nfull' are:"
        write(*,'(3x,a)') "expscal1"
        write(*,'(3x,a)') "ener_par1"
        write(*,'(3x,a)') "ener_par2"
        write(*,'(3x,a)') "shell_xi"
        write(*,'(3x,a)') "shell_cnf1"
        write(*,'(3x,a)') "shell_cnf2"
        write(*,'(3x,a)') "shell_cnf3"
        write(*,'(3x,a)') "expscal3"
        write(*,'(3x,a)') "shell_cnf4"
        write(*,'(3x,a)') "shell_resp1"
        write(*,'(3x,a)') "shell_resp2"
        stop
    endif
    
    if ((.not.par).and.(.not.parname).and.(.not.fullparn)) write(*,*) "DESIRED PARAMETER NOT GIVEN!"
    if ((.not.per).and.(.not.gr)) write(*,*) "DESIRED ELEMENT RANGE NOT GIVEN!"

    elemrange=0
    if (per) then
       SELECT CASE (period)
          CASE (1)
              WRITE(*,'(a,i1,a)')  "Period ",period," selected."
              elemrange(1)=1
              elemrange(2)=2
          CASE (2)
              WRITE(*,'(a,i1,a)')  "Period ",period," selected."
              elemrange(1)=3
              elemrange(2)=10
          CASE (3)
              WRITE(*,'(a,i1,a)')  "Period ",period," selected."
              elemrange(1)=11
              elemrange(2)=18
          CASE (4)
              WRITE(*,'(a,i1,a)')  "Period ",period," selected."
!             elemrange(1)=19
!             elemrange(2)=36
              elemrange(1)=21
              elemrange(2)=30
          CASE (5)
              WRITE(*,'(a,i1,a)')  "Period ",period," selected."
!             elemrange(1)=37
!             elemrange(2)=54
              elemrange(1)=39
              elemrange(2)=48
          CASE (6)
              WRITE(*,'(a,i1,a)')  "Period ",period," selected."
!             elemrange(1)=55
!             elemrange(2)=57
!             elemrange(3)=71
!             elemrange(4)=86
              elemrange(1)=57
              elemrange(2)=80
       END SELECT
   elseif (gr) then
       SELECT CASE (group)
          CASE (1)
              elemrange(1)=1
              elemrange(2)=3
              elemrange(3)=11
              elemrange(4)=19
              elemrange(5)=37
              elemrange(6)=55
              arraysize=6
          CASE (2)
              elemrange(1)=4
              elemrange(2)=12
              elemrange(3)=20
              elemrange(4)=38
              elemrange(5)=56
              arraysize=5
          CASE (3)
              elemrange(1)=21
              elemrange(2)=39
              elemrange(3)=57
              arraysize=3
          CASE (4)
              elemrange(1)=22
              elemrange(2)=40
              elemrange(3)=72
              arraysize=3
          CASE (5)
              elemrange(1)=23
              elemrange(2)=41
              elemrange(3)=73
              arraysize=3
          CASE (6)
              elemrange(1)=24
              elemrange(2)=42
              elemrange(3)=74
              arraysize=3
          CASE (7)
              elemrange(1)=25
              elemrange(2)=43
              elemrange(3)=75
              arraysize=3
          CASE (8)
              elemrange(1)=26
              elemrange(2)=44
              elemrange(3)=76
              arraysize=3
          CASE (9)
              elemrange(1)=27
              elemrange(2)=45
              elemrange(3)=77
              arraysize=3
          CASE (10)
              elemrange(1)=28
              elemrange(2)=46
              elemrange(3)=78
              arraysize=3
          CASE (11)
              elemrange(1)=29
              elemrange(2)=47
              elemrange(3)=79
              arraysize=3
          CASE (12)
              elemrange(1)=30
              elemrange(2)=48
              elemrange(3)=80
              arraysize=3
          CASE (13)
              elemrange(1)=5
              elemrange(2)=13
              elemrange(3)=31
              elemrange(4)=49
              elemrange(5)=81
              arraysize=5
          CASE (14)
              elemrange(1)=6
              elemrange(2)=14
              elemrange(3)=32
              elemrange(4)=50
              elemrange(5)=82
              arraysize=5
          CASE (15)
              elemrange(1)=7
              elemrange(2)=15
              elemrange(3)=33
              elemrange(4)=51
              elemrange(5)=83
              arraysize=5
          CASE (16)
              elemrange(1)=8
              elemrange(2)=16
              elemrange(3)=34
              elemrange(4)=52
              elemrange(5)=84
              arraysize=5
          CASE (17)
              elemrange(1)=9
              elemrange(2)=17
              elemrange(3)=35
              elemrange(4)=53
              elemrange(5)=85
              arraysize=5
          CASE (18)
              elemrange(1)=2
              elemrange(2)=10
              elemrange(3)=18
              elemrange(4)=36
              elemrange(5)=54
              elemrange(6)=86
              arraysize=6
       END SELECT
!      WRITE(*,'(a,i2,a)')  "Group ",group," selected."
   endif
   if (par) then
!      write(*,'(a,i3,a)') "Parameter ",despar," selected."
   endif


!   write(*,'(a)') "Reading in parameters..."
    !call GET_ENVIRONMENT_VARIABLE("USER",username)
    !pathname= "/home/" // trim(username) // "/.atompara"
    pathname="~/.atompara"
    open(newunit=myunit,file=pathname,status='old')
    do i=1,2
        read(myunit,*) globpar(1:10,i)
    enddo
    i=0
    do ! number element in PSE
        i=i+1
        read(myunit,*,iostat=ioerror) elmtnum
        if (ioerror.ne.0) exit
        if (elmtnum.ne.i) then
!           write(*,'(a,i2,a,i2,a)') "Element ",i," to ",elmtnum-1," is not defined."
            i=elmtnum
        endif
        !write(*,'(i2)') elmtnum
        read(myunit,*) expscal1(1:10,elmtnum)   ! 1 -10  E
        allpar(1:10,elmtnum)=expscal1(1:10,elmtnum)
        read(myunit,*) ener_par1(1:10,elmtnum)    ! 11-20  E
        allpar(11:20,elmtnum)=ener_par1(1:10,elmtnum)
        read(myunit,*) ener_par2(1:10,elmtnum)    ! 21-30  E
        allpar(21:30,elmtnum)=ener_par2(1:10,elmtnum)
        read(myunit,*) floats(1:10)      ! 31-40
        allpar(31:40,elmtnum)=floats(1:10)
        read(myunit,*) floats(1:10)      ! 41-50
        allpar(41:50,elmtnum)=floats(1:10)
        read(myunit,*) floats(1:10)      ! 51-60
        allpar(51:60,elmtnum)=floats(1:10)
        read(myunit,*) floats(1:10)      ! 61-70
        allpar(61:70,elmtnum)=floats(1:10)
        read(myunit,*) shell_xi(1:10,elmtnum)    ! 71-80  shellQ
        allpar(71:80,elmtnum)=shell_xi(1:10,elmtnum)
        read(myunit,*) shell_cnf1(1:10,elmtnum)    ! 81-90    "
        allpar(81:90,elmtnum)=shell_cnf1(1:10,elmtnum)
        read(myunit,*) shell_cnf2(1:10,elmtnum)    ! 91-100   "
        allpar(91:100,elmtnum)=shell_cnf2(1:10,elmtnum)
        read(myunit,*) shell_cnf3(1:10,elmtnum)    ! 101-110  "
        allpar(101:110,elmtnum)=shell_cnf3(1:10,elmtnum)
        read(myunit,*) expscal3(1:10,elmtnum)   ! 111-120  "
        allpar(111:120,elmtnum)=expscal3(1:10,elmtnum)
        read(myunit,*) shell_cnf4(1:10,elmtnum)    ! 121-130  "
        allpar(121:130,elmtnum)=shell_cnf4(1:10,elmtnum)
        read(myunit,*) shell_resp1(1:10,elmtnum)  ! 131-140  "
        allpar(131:140,elmtnum)=shell_resp1(1:10,elmtnum)
        read(myunit,*) shell_resp2(1:10,elmtnum)  ! 141-150  "
        allpar(141:150,elmtnum)=shell_resp2(1:10,elmtnum)
    enddo

   
    paramset=0.0d0
    parn=trim(parn)
    if (parname.or.fullparn) then
    SELECT CASE (parn)
       CASE ("expscal1")
           paramset(1:10,1:86)=expscal1(1:10,1:86)
       CASE ("ener_par1")
           paramset(1:10,1:86)=ener_par1(1:10,1:86)
       CASE ("ener_par2")
           paramset(1:10,1:86)=ener_par2(1:10,1:86)
       CASE ("shell_xi")
           paramset(1:10,1:86)=shell_xi(1:10,1:86)
       CASE ("shell_cnf1")
           paramset(1:10,1:86)=shell_cnf1(1:10,1:86)
       CASE ("shell_cnf2")
           paramset(1:10,1:86)=shell_cnf2(1:10,1:86)
       CASE ("shell_cnf3")
           paramset(1:10,1:86)=shell_cnf3(1:10,1:86)
       CASE ("expscal3")
           paramset(1:10,1:86)=expscal3(1:10,1:86)
       CASE ("shell_cnf4")
           paramset(1:10,1:86)=shell_cnf4(1:10,1:86)
       CASE ("shell_resp1")
           paramset(1:10,1:86)=shell_resp1(1:10,1:86)
       CASE ("shell_resp2")
           paramset(1:10,1:86)=shell_resp2(1:10,1:86)
    END SELECT
    WRITE(*,'(a,1x,a,1x,a)') "Parameter set",trim(parn),"selected."
    endif

    open(newunit=myunit,file='param_ana.dat',status='replace')
    if (par) then
      if (per) then
          do i=elemrange(1),elemrange(2)
              if (allpar(despar, i) .ne. 0.0d0) then
                  write(myunit,*) i, allpar(despar,i)
              endif
          enddo
          if (elemrange(3).ne.0) then
              do i=elemrange(3),elemrange(4)
                  if (allpar(despar, i) .ne. 0.0d0) then
                      write(myunit,*) i, allpar(despar,i)
                  endif
              enddo
          endif
      endif
      if (gr) then
          do i=1,arraysize
              if (allpar(despar, elemrange(i)) .ne. 0.0d0) then
                  write(myunit,*) elemrange(i), allpar(despar, elemrange(i))
              endif
          enddo
      endif
    else if (parname) then
        if (per) then
            do i=elemrange(1),elemrange(2)
                if (paramset(1,i) .ne. 0.0d0) then
                    write(myunit,'(i2,a)',advance='no') i," "
                    do j=1,8
                        if (paramset(j,i) .ne. 0.0d0) then
                            if (.not.j.eq.8) then
                                write(myunit,'(999(f12.8,1x))',advance='no') paramset(j,i)
                            else
                                write(myunit,'(f12.8,1x)') paramset(j,i)
                            endif
                        else
!                           write(myunit,'(/)',advance='no')
                            exit
                        endif
                    enddo
                endif
            enddo
            if (elemrange(3).ne.0) then
                do i=elemrange(3),elemrange(4)
                    if (paramset(1,i) .ne. 0.0d0) then
                        write(myunit,'(i2,a)',advance='no') i," "
                        do j=1,8
                            if (paramset(j,i) .ne. 0.0d0) then
                                if (.not.j.eq.8) then
                                    write(myunit,'(999(f12.8,1x))',advance='no') paramset(j,i)
                                else
                                    write(myunit,'(f12.8,1x)') paramset(j,i)
                                endif
                            else
!                               write(myunit,'(/)',advance='no')
                                exit
                            endif
                        enddo
                    endif
                enddo
            endif
        endif
        if (gr) then
            do i=1,arraysize
                if (paramset(1,elemrange(i)) .ne. 0.0d0) then
                    write(myunit,'(i2,a)',advance='no') elemrange(i)," "
                    do j=1,8
                        if (paramset(j,elemrange(i)) .ne. 0.0d0) then
                            if (.not.j.eq.8) then
                                write(myunit,'(999(f12.8,1x))',advance='no') paramset(j,elemrange(i))
                            else
                                write(myunit,'(f12.8,1x)') paramset(j,elemrange(i))
                            endif
                        else
!                           write(myunit,'(/)',advance='no')
                            exit
                        endif
                    enddo
                endif
            enddo
        endif
    else if (fullparn) then
        if (per) then
            do i=elemrange(1),elemrange(2)
                if (paramset(1,i) .ne. 0.0d0) then
                    write(myunit,'(i2,a)',advance='no') i," "
                    do j=1,9
                        write(myunit,'(9(f12.8,1x))',advance='no') paramset(j,i)
                    enddo
                    write(myunit,'(f12.8,1x,/)',advance='no') paramset(10,i)
                endif
            enddo
            if (elemrange(3).ne.0) then
                do i=elemrange(3),elemrange(4)
                    if (paramset(1,i) .ne. 0.0d0) then
                        write(myunit,'(i2,a)',advance='no') i," "
                        do j=1,9
                            write(myunit,'(9(f12.8,1x))',advance='no') paramset(j,i)
                        enddo
                        write(myunit,'(f12.8,1x,/)',advance='no') paramset(10,i)
                    endif
                enddo
            endif
        endif
        if (gr) then
            do i=1,arraysize
                if (paramset(1,elemrange(i)) .ne. 0.0d0) then
                    write(myunit,'(i2,a)',advance='no') elemrange(i)," "
                    do j=1,9
                        write(myunit,'(999(f12.8,1x))',advance='no') paramset(j,elemrange(i))
                    enddo
                    write(myunit,'(f12.8,1x,/)',advance='no') paramset(10,elemrange(i))
                endif
            enddo
        endif
    endif
    close(myunit)
!   write(*,*) "Printout of 'param_ana.dat':"
    call system("cat param_ana.dat")

end program parascript
