      function pgrid(i,j,k,ngrx,ngry)

! --- converts 3d array to 1D array
      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'
      integer(KIND=normal_int)::pgrid,i,j,k,ngrx,ngry
      pgrid = i + j * ngrx + k * ngry*ngrx
      return
      end
