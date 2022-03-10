!> Proxy module for exporting all BLAS interfaces and wrappers
module gtb_la
  use gtb_la_level1, only : &
    & la_axpy, la_copy, la_dot, la_dotc, la_dotu, la_dsdot, la_sdsdot, &
    & la_rotg, la_rotm, la_rotmg, la_scal, la_rot, la_swap, la_abs1, &
    & la_asum, la_nrm2, la_iamax
  use gtb_la_level2, only : &
    & la_gbmv, la_gemv, la_ger, la_gerc, la_geru, la_sbmv, la_spmv, &
    & la_hbmv, la_hemv, la_spr2, la_spr, la_syr2, la_syr, la_her2, la_her, &
    & la_symv, la_hpmv, la_hpr2, la_hpr, la_tbmv, la_tbsv, la_tpmv, la_tpsv, &
    & la_trmv, la_trsv
  use gtb_la_level3, only : &
    & la_gemm, la_hemm, la_her2k, la_herk, la_symm, la_syr2k, la_syrk, &
    & la_trsm, la_trmm
  use gtb_lapack_eig, only : &
    & la_syev, la_syevx, la_sygvx, la_sygvd
  implicit none
  public
end module gtb_la
