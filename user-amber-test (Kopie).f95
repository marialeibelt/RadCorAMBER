                 !!!!!!!!!!!!!!!!!!!!!
                     MODULE  USER
                 !!!!!!!!!!!!!!!!!!!!!

  use mcmule

  implicit none

!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!

  integer, parameter :: nrq = 2
  integer, parameter :: nrbins = 500
  real(kind=prec), parameter :: &
       min_val(nrq) = (/ .3e-3,  95.e3/) !rad, MeV
  real(kind=prec), parameter :: &
       max_val(nrq) = (/ 2.e-3, 101.e3/) !rad, MeV
  integer :: userdim = 0

!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!

!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!

    !! ============================================== !!
    !! DO NOT EVEN THINK ABOUT CHANGING ANYTHING HERE !!
    !! ============================================== !!

  integer :: namesLen=6
  integer :: filenamesuffixLen=10
  integer :: nq=nrq
  integer :: nbins=nrbins



!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!

            !! ----------------------------------------- !!
            !!     There are two versions of binning     !!
            !!     One for computing   d \sigma/ d Q     !!
            !!     One for computing  Q d \sigma/ d Q    !!
            !!  choose by setting the variable bin_kind  !!
            !! ----------------------------------------- !!
  integer :: bin_kind = 0       !!  0 for d \sig/dQ; +1 for Q d \sig/dQ;


  contains


  SUBROUTINE FIX_MU

  musq = mM**2

  END SUBROUTINE FIX_MU



  SUBROUTINE INITUSER
  print*, "This is Mary testing the McMule userfile <3"
  print*, " * E_mu = 150 GeV"
  print*, " * 0.3 < th_mu < 2 mrad"
  print*, " * E_mu > 70 GeV"
  
  call initflavour("mu-p", Mmu**2+Mproton**2+2*Mproton*100.e3)
  END SUBROUTINE


  FUNCTION QUANT(q1,q2,q3,q4,q5,q6,q7)

  real (kind=prec), intent(in) :: q1(4),q2(4),q3(4),q4(4), q5(4),q6(4),q7(4)
  real (kind=prec) :: ql1(4),ql2(4),ql3(4),ql4(4), ql5(4),ql6(4),ql7(4)  ! in lab frame
  real (kind=prec) :: th3,q3perp,Emu
  real (kind=prec) :: quant(nr_q)
  
  !! ==== keep the line below in any case ==== !!
  call fix_mu

  ql1 = boost_rf(q2,q1)
  ql2 = boost_rf(q2,q2)
  ql3 = boost_rf(q2,q3)
  ql4 = boost_rf(q2,q4)
  ql5 = boost_rf(q2,q5)
  ql6 = boost_rf(q2,q6)
  ql7 = boost_rf(q2,q7)

  q3perp=sqrt(ql3(1)**2+ql3(2)**2)
  th3 = atan2(q3perp,ql3(3))   ! scattering angle in rad
  Emu = ql3(4)
  pass_cut = .true.

  ! scattering angle cut
  if(th3 < 0.3e-3) pass_cut = .false.
  if(th3 > 2.e-3) pass_cut = .false.
  !write(*,*) 'Emu ', ql3(4)
  if(Emu < 70.e3) then 
  	pass_cut =.false.
  	!write(*,*) 'Emu cut ', Emu
  endif

  names(1) = 'th3'
  quant(1) = th3
  names(2) = 'Emu'
  quant(2) = Emu

  END FUNCTION QUANT


  SUBROUTINE USEREVENT(X, NDIM)
  integer :: ndim
  real(kind=prec) :: x(ndim)
  userweight = 1.
  END SUBROUTINE USEREVENT


                 !!!!!!!!!!!!!!!!!!!!!!!
                     END MODULE  USER
                 !!!!!!!!!!!!!!!!!!!!!!!
