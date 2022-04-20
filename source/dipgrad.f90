
      open(unit=21,file='dipgrad')
      read(21,'(a)') atmp
      write(*,*) 'reading reference (TM) dipgrad and calculating ...'
      x=0.015_wp
      k=0
      do i=1,n
         do j=1,3
            k=k+1
            read(21,*) dipgrad_ref(1:3,k)
            xyz(j,i)=xyz(j,i)+x    
            call calcrab(n,at,xyz,rab)
            call sint(n,ndim,at,xyz,rab,S,xnorm)       ! exact S
            call dipint(n,ndim,at,xyz,rab,xnorm,pnt,D3)! dipole integrals
            call pgtb(.false.,-1,n,ndim,nel,nopen,ihomo,at,chrg,xyz,z,rab,pnt,xnorm,S,T,D3,&
     &                efield,ML1,ML2,psh,q,P,F,eps,eel,ecoul,wbo,dipr,alp)
            xyz(j,i)=xyz(j,i)-2_wp*x
            call calcrab(n,at,xyz,rab)
            call sint(n,ndim,at,xyz,rab,S,xnorm)       ! exact S
            call dipint(n,ndim,at,xyz,rab,xnorm,pnt,D3)! dipole integrals
            call pgtb(.false.,-1,n,ndim,nel,nopen,ihomo,at,chrg,xyz,z,rab,pnt,xnorm,S,T,D3,&
     &                efield,ML1,ML2,psh,q,P,F,eps,eel,ecoul,wbo,dipl,alp)
            dipgrad(1:3,k)=(dipr(1:3)-dipl(1:3))/(2_wp*x)
!           write(*,'(2i3,6F12.6)') i,j,dipr,dipl
            xyz(j,i)=xyz(j,i)+x   
         enddo
      enddo
      close(21)
