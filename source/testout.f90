         open(unit=112,file='.datap')  
         open(unit=113,file='.dataq')  
         open(unit=114,file='.datab')  
         open(unit=115,file='.datad')  
         do i=1,n
         write(113,*) z(i)-q_ref(i),z(i)-q(i) 
         do j=1,bas_nsh(at(i))
            write(112,*) psh_ref(j,i),psh(j,i) 
         enddo
         enddo
         close(112)
         close(113)
         do i=2,n
         do j=1,i-1 
            if(abs(wbo_ref(j,i)).gt.0.1) write(114,'(2F16.8,4i4)') wbo_ref(j,i), wbo(j,i), j, i, at(j),at(i)
         enddo
         enddo
         close(114)
         do i=1,3
         if(abs(dip_ref(i)).gt.0.01) write(115,*) dip_ref(i),dip(i)
         enddo
         close(115)
         if(raman)then
         open(unit=126,file='.datara')  
         do m=1,6
         do i=1,n
         do j=1,3
         if(abs(fdgrad_ref(j,i,m)).gt.1d-5) write(126,'(2F22.14)') fdgrad_ref(j,i,m), fdgrad(j,i,m)
         enddo
         enddo
         enddo
         close(126)
         endif
         if(dgrad)then
         open(unit=116,file='.datadg')  
         do i=1,3*n
         do j=1,3
         if(abs(dipgrad_ref(j,i)).gt.1d-5) write(116,'(2F22.14)') dipgrad_ref(j,i), dipgrad(j,i)
         enddo
         enddo
         close(116)
         if(alp_ref(1).gt.0)then
         open(unit=117,file='.datapol')  
         do i=1,6  
         if(abs(alp_ref(i)).gt.1d-2) write(117,'(2F22.14)') alp_ref(i), alp(i)
         enddo
         close(117)
         endif
         endif
         if(ekinref.gt.0.001)then
         open(unit=117,file='.datat')  
         write(117,*) ekinref,ekin
         close(117)
         endif
