      open(unit=21,file='fdgrad')
      write(*,*) 'reading reference (TM) fdgrad and calculating ...',y
      y=0.0001_wp ! field strength
      eftmp=0
      eftmp(1,1)= y    
      eftmp(1,2)=-y    
      eftmp(2,3)= y    
      eftmp(2,4)=-y    
      eftmp(3,5)= y    
      eftmp(3,6)=-y    
      x=0.005_wp
      do m=1,6
      efield(1:3)=eftmp(1:3,m)
      read(21,'(a)') atmp
      k=0
      do i=1,n
         do j=1,3
            k=k+1
            read(21,*) fdgrad_ref(1:3,k,m)
            xyz(j,i)=xyz(j,i)+x    
            call calcrab(n,at,xyz,rab)
            call modbas(n,at,4) 
            call sint(n,ndim,at,xyz,rab,S,xnorm)       ! exact S
            call shellq(.false.,prop,n,ndim,nel,nopen,ihomo,at,chrg,xyz,z,rab,pnt,xnorm,S,&
     &                  efield,psh,q,P,wbo,dip,alp)
            fdgrad(1:3,k,m)=dip(1:3)
            xyz(j,i)=xyz(j,i)-2_wp*x
            call calcrab(n,at,xyz,rab)
            call modbas(n,at,4) 
            call sint(n,ndim,at,xyz,rab,S,xnorm)       ! exact S
            call shellq(.false.,prop,n,ndim,nel,nopen,ihomo,at,chrg,xyz,z,rab,pnt,xnorm,S,&
     &                  efield,psh,q,P,wbo,dip,alp)
            fdgrad(1:3,k,m)=(fdgrad(1:3,k,m)-dip(1:3))/(2_wp*x)
            xyz(j,i)=xyz(j,i)+x   
         enddo
      enddo
      read(21,'(a)') atmp
      enddo
      close(21)
      
      k=0
      do i=1,n
         do j=1,3
            k=k+1
            y=0.0001_wp  ! field strength
            fdgrad(1:3,k,7)=(fdgrad(1:3,k,1)-fdgrad(1:3,k,2))/(2_wp*y)
            fdgrad(1:3,k,8)=(fdgrad(1:3,k,3)-fdgrad(1:3,k,4))/(2_wp*y)
            fdgrad(1:3,k,9)=(fdgrad(1:3,k,5)-fdgrad(1:3,k,6))/(2_wp*y)
            y=0.0020_wp ! field strength TM script
            fdgrad_ref(1:3,k,7)=(fdgrad_ref(1:3,k,1)-fdgrad_ref(1:3,k,2))/(2_wp*y)
            fdgrad_ref(1:3,k,8)=(fdgrad_ref(1:3,k,3)-fdgrad_ref(1:3,k,4))/(2_wp*y)
            fdgrad_ref(1:3,k,9)=(fdgrad_ref(1:3,k,5)-fdgrad_ref(1:3,k,6))/(2_wp*y)
         enddo
      enddo
