      program getfrac
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Ordering the .res file to input coord for fort.77
!     Readin files: fort.33,source_filename(.res)
!     Format of fort.33, which is from source_filename.inp:
!     source_filename
!     nbead
!     C1 H2 C5 O4 ...
!     Output file: source_filename.frac
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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

      character(LEN=default_path_length)::infile,outfile
      character(LEN=default_string_length)::sname(500),celltype
      integer(KIND=normal_int)::nbead, i, j
      real::cell(6),bead(3,500),dum
      
      read(4,*) infile
      read(4,*) nbead

      outfile=trim(infile) // '.frac'

      OPEN(22,FILE=infile,STATUS='old')
      read(22,*)
      read(22,*) celltype,dum,cell
      read(22,*)
      read(22,*)
      do i=1,nbead
         read(22,*) sname(i),dum,bead(1:3,i)
      end do
      close(22)
      
      print*, nbead
      OPEN(23,FILE=outfile,STATUS='unknown')
      write(23,*) cell
      write(23,*) nbead
      do i=1,nbead
         write(23,*) sname(i), bead(1:3,i)
      end do
      close(23)
      print*, 'coordinates are in : ', outfile
      end program
