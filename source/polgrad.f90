      open(unit=21,file='polgrad')
      write(*,*) 'reading reference (TM) polgrad and calculating ...'
      tenmap(1)=1
      tenmap(2)=3
      tenmap(3)=6
      tenmap(4)=2
      tenmap(5)=4
      tenmap(6)=5
      x=0.015_wp
      test2=.true.
      read(21,'(a)') atmp
      if(index(atmp,'ORCA').ne.0) test2=.false.
      if(test2) read(21,'(a)') atmp ! TM read
      do m=1,6
      if(test2) read(21,'(a)') atmp ! TM read
      do i=1,n
         read(21,*) fdgrad_ref(1:3,i,tenmap(m))
      enddo
      enddo
      close(21)

      call modbas(n,at,4) 
      do i=1,n
         do j=1,3
            xyz(j,i)=xyz(j,i)+x    
            call calcrab(n,at,xyz,rab)
            call sint(n,ndim,at,xyz,rab,S,xnorm)       ! exact S
            call dipint(n,ndim,at,xyz,rab,xnorm,pnt,D3)! dipole integrals
            call pgtb(.false.,-2,n,ndim,nel,nopen,ihomo,at,chrg,xyz,z,rab,pnt,xnorm,S,T,D3,&
     &               efield,ML1,ML2,psh,q,P,F,eps,eel,ecoul,wbo,dip,alpr)
            xyz(j,i)=xyz(j,i)-2_wp*x
            call calcrab(n,at,xyz,rab)
            call sint(n,ndim,at,xyz,rab,S,xnorm)       ! exact S
            call dipint(n,ndim,at,xyz,rab,xnorm,pnt,D3)! dipole integrals
            call pgtb(.false.,-2,n,ndim,nel,nopen,ihomo,at,chrg,xyz,z,rab,pnt,xnorm,S,T,D3,&
     &               efield,ML1,ML2,psh,q,P,F,eps,eel,ecoul,wbo,dip,alpl)
            fdgrad(j,i,1:6)=(alpr(1:6)-alpl(1:6))/(2_wp*x)
            xyz(j,i)=xyz(j,i)+x   
         enddo
      enddo
      
