module aescom
      use iso_fortran_env, only : wp => real64
      implicit none
      private :: wp
      public

      real(wp),allocatable :: qm  (  :)  
      real(wp),allocatable :: dipm(:,:)  
      real(wp),allocatable :: qp  (:,:)  
      real(wp),allocatable :: pint(:,:)  

end module aescom       
