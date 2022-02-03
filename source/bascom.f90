module bascom
      use iso_fortran_env, only : wp => real64
      implicit none
      private :: wp
      public

      integer  ncao, nsao, npr, nsh
      integer  bas_nsh(86) 
      real(wp) bas_ec(2,5,10,86)  ! read data: exp/contr,max contr,max shells,elements
      integer  bas_lsh (10,86) 
      integer  bas_npr (10,86)

      integer, allocatable :: prim_npr(:)  
      integer, allocatable :: prim_count(:)  
      integer, allocatable ::  aoshell(:,:)  
      integer, allocatable :: caoshell(:,:)  
      integer, allocatable :: shell2ao(:)    
      integer, allocatable :: shmap(:,:)    
      integer, allocatable :: aoat(:)    
      real(wp),allocatable :: prim_exp(:)  
      real(wp),allocatable :: prim_cnt(:)  

end module bascom       
