
      open(unit=21,file='gradient')
      read(21,'(a)') atmp
      read(21,'(a)') atmp
      do i=1,n
         read(21,*) floats(1:3)
         if(abs(floats(1)-xyz(1,i)).gt.1.d-8 .or. &
     &      abs(floats(2)-xyz(2,i)).gt.1.d-8 .or. &
     &      abs(floats(3)-xyz(3,i)).gt.1.d-8) &
     &      stop 'inconsistent coordinates in gradient and coord'
      enddo
      write(*,*) 'reading reference (TM) grad ...'
      do i=1,n
         read(21,*) grad_ref(1:3,i)
      enddo
      close(21)

      ngrad = 0

      step=0.020_wp
      grad = 0
      if(n.eq.2.and.abs(xyz(1,1)).lt.1d-6.and.abs(xyz(2,1)).lt.1d-6.and.abs(xyz(1,2)).lt.1d-6.and.abs(xyz(2,2)).lt.1d-6)then
         xyz(3,1)=xyz(3,1)+step
!---------------------------    right step
         ngrad = ngrad + 1
         call calcrab(n,at,xyz,rab)
         if (calc_ptb_grad) then
          call sint(n,ndim,at,xyz,rab,S,xnorm)       ! exact S
          call pgtb(.false.,0,n,ndim,nel,nopen,ihomo,at,chrg,filter,xyz,z,rab,pnt,xnorm,S,D3,efield,ML1,ML2,psh,q,P,F,eps,wbo,dipr,alp)
          write(atmp,'(''mv ptb_dump ptb_dump_'',i0)') ngrad
          call system(atmp)
         endif
         write(atmp,'(''cp ptb_dump_'',i0,'' ptb_dump'')') ngrad
         call system(atmp)
         call ptb_energy(.false.,n,ndim,nopen,at,z,xyz,rab,q,psh,wbo,P,eref,er)
!---------------------------      
         xyz(3,1)=xyz(3,1)-step*2d0
!---------------------------    left one
         ngrad = ngrad + 1
         call calcrab(n,at,xyz,rab)
         if (calc_ptb_grad) then
          call sint(n,ndim,at,xyz,rab,S,xnorm)       ! exact S
          call pgtb(.false.,0,n,ndim,nel,nopen,ihomo,at,chrg,filter,xyz,z,rab,pnt,xnorm,S,D3,efield,ML1,ML2,psh,q,P,F,eps,wbo,dipr,alp)
          write(atmp,'(''mv ptb_dump ptb_dump_'',i0)') ngrad
          call system(atmp)
         endif
         write(atmp,'(''cp ptb_dump_'',i0,'' ptb_dump'')') ngrad
         call system(atmp)
         call ptb_energy(.false.,n,ndim,nopen,at,z,xyz,rab,q,psh,wbo,P,eref,el)
!---------------------------    back
         xyz(3,1)=xyz(3,1)+step
         grad(3,1)= (er-el)/(2d0*step)
         grad(3,2)=-(er-el)/(2d0*step)
      else
         call getsymmetry2(.true.,n,at,xyz,0.0001d0,ntrans,ict,trans) 
         if(ntrans.gt.1) then ! symmetric
            dgen = 0
            do i=1,n
               do j=1,ntrans
                  if(ict(i,j).eq.i) dgen(i)=dgen(i)+1   
               enddo
               dgen(i)=ntrans / dgen(i)
            enddo
         else
            dgen = 1
         endif
         do i=2,n ! first atom added later
            if(ntrans.gt.1.and.sum(abs(grad(:,i))).gt.1d-12) cycle  !  symmetry related atom already done
            do j=1,3
               if(ntrans.gt.1.and.abs(xyz(j,i)).lt.1d-12) cycle  ! =0 by symmetry
               if(ntrans.gt.1) write(*,*) 'atom ',i,' dir',j,' degen. by symmetry ',dgen(i)
               xyz(j,i)=xyz(j,i)+step
!---------------------------    right step
         ngrad = ngrad + 1
         call calcrab(n,at,xyz,rab)
         if (calc_ptb_grad) then
          call sint(n,ndim,at,xyz,rab,S,xnorm)       ! exact S
          call pgtb(.false.,0,n,ndim,nel,nopen,ihomo,at,chrg,filter,xyz,z,rab,pnt,xnorm,S,D3,efield,ML1,ML2,psh,q,P,F,eps,wbo,dipr,alp)
          write(atmp,'(''mv ptb_dump ptb_dump_'',i0)') ngrad
          call system(atmp)
         endif
         write(atmp,'(''cp ptb_dump_'',i0,'' ptb_dump'')') ngrad
         call system(atmp)
         call ptb_energy(.false.,n,ndim,nopen,at,z,xyz,rab,q,psh,wbo,P,eref,er)
               xyz(j,i)=xyz(j,i)-step*2d0
!---------------------------    left one
         ngrad = ngrad + 1
         call calcrab(n,at,xyz,rab)
         if (calc_ptb_grad) then
          call sint(n,ndim,at,xyz,rab,S,xnorm)       ! exact S
          call pgtb(.false.,0,n,ndim,nel,nopen,ihomo,at,chrg,filter,xyz,z,rab,pnt,xnorm,S,D3,efield,ML1,ML2,psh,q,P,F,eps,wbo,dipr,alp)
          write(atmp,'(''mv ptb_dump ptb_dump_'',i0)') ngrad
          call system(atmp)
         endif
         write(atmp,'(''cp ptb_dump_'',i0,'' ptb_dump'')') ngrad
         call system(atmp)
         call ptb_energy(.false.,n,ndim,nopen,at,z,xyz,rab,q,psh,wbo,P,eref,el)
!---------------------------    back
               xyz(j,i)=xyz(j,i)+step
               grad(j,i)= (er-el)/(2d0*step)
            enddo
            if(ntrans.gt.1) then
               dum = 0
               dum(1:3,i) = grad(1:3,i)
               grad(1:3,i) = 0
               call grdsym(dum,n,ntrans,ict,trans) ! symmetrize
               grad = grad + dum * dble(dgen(i))
            endif
         enddo
         grad(1,1)=-sum(grad(1,2:n)) ! use other forces to compute first atom
         grad(2,1)=-sum(grad(2,2:n))
         grad(3,1)=-sum(grad(3,2:n))
      endif
      write(*,*) 'norm of gradient                                ',sqrt(sum((grad         )**2))
      write(*,*) 'norm of reference gradient                      ',sqrt(sum((grad_ref     )**2))
      write(*,*) 'norm of gradient difference to reference vector ',sqrt(sum((grad-grad_ref)**2))
      call system('rm -rf ptb_dump')
