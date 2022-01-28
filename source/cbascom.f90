module cbascom
      use iso_fortran_env, only : wp => real64
      implicit none
      private :: wp
      public

      integer  cncao, cnsao, cnpr
      integer  ncorelist
      integer  cbas_nsh    (86) 
      integer  cbas_lsh (13,86) 
      integer  cbas_npr (13,86)
      integer  cbas_npq (13,86)
      real(wp) cslexpo  (13,86)  
      real(wp) clev     (13,86)  

      integer, allocatable :: corelist(:)  
      integer, allocatable :: cprim_npr(:)  
      integer, allocatable :: cprim_count(:)  
      integer, allocatable ::  aoshellc(:,:)  
      integer, allocatable :: caoshellc(:,:)  
      real(wp),allocatable :: cprim_exp(:)  
      real(wp),allocatable :: cprim_cnt(:)  

end module cbascom       
