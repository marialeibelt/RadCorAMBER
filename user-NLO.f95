                 !!!!!!!!!!!!!!!!!!!!!
                     MODULE  USER
                 !!!!!!!!!!!!!!!!!!!!!

  use mcmule

  implicit none

!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!

  integer, parameter :: nrq = 12			!th3,Emu,th5,Eph, th5x, th5y, +same in cms
  integer, parameter :: nrbins = 500
  real(kind=prec), parameter :: &
       min_val(nrq) = (/ .3e-3,  95.e3, -12.e-3, 1.e3, -12.e-3, -12.e-3, &
       			.3e-3,  95.e3, -12.e-3, 1.e3, -12.e-3, -12.e-3/)      	!rad, MeV, rad, MeV
  real(kind=prec), parameter :: &
       max_val(nrq) = (/ 2.e-3, 101.e3, 12.e-3, 100.e3, 12.e-3, 12.e-3, &
       			 2.e-3, 101.e3, 12.e-3, 100.e3, 12.e-3, 12.e-3/) 	!rad, MeV, rad, MeV
  integer :: userdim = 0

!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!

!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!

    !! ============================================== !!
    !! DO NOT EVEN THINK ABOUT CHANGING ANYTHING HERE !!
    !! ============================================== !!

  integer :: namesLen=12
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
  print*, " * 1.35 < th_mu < 1.65 mrad"
  print*, " * E_mu > 70 GeV"
  print*, " * E_ph > 200MeV"
  
  call initflavour("mu-p", Mmu**2+Mproton**2+2*Mproton*100.e3)
  END SUBROUTINE


  FUNCTION QUANT(q1,q2,q3,q4,q5,q6,q7)

  real (kind=prec), intent(in) :: q1(4),q2(4),q3(4),q4(4), q5(4),q6(4),q7(4)
  real (kind=prec) :: ql1(4),ql2(4),ql3(4),ql4(4), ql5(4),ql6(4),ql7(4)  ! in lab frame
  real (kind=prec) :: th3,q3perp,q5perp,th5,Emu,Eph
  real (kind=prec) :: phi5_x, phi5_y  ! in lab frame
  real (kind=prec) :: q3perp_cms, th3_cms, Emu_cms, Eph_cms, th5_cms, q5perp_cms, phi5_x_cms, phi5_y_cms
  real (kind=prec) :: quant(nrq)
  
  !! ==== keep the line below in any case ==== !!
  call fix_mu
  ! proton frame
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
  Eph_cut = 1.e3
  
  ! cms frame
  q3perp_cms = sqrt(q3(1)**2+q3(2)**2)
  th3_cms = atan2(q3perp_cms,q3(3))   ! scattering angle in rad
  Emu_cms = q3(4)
  
  pass_cut = .true.

  ! Muon cuts
  if(th3.lt.1.35e-3) pass_cut = .false.
  if(th3.gt.1.65e-3) pass_cut = .false.
  !write(*,*) 'Emu ', ql3(4)
  if(Emu.lt.70.e3) pass_cut =.false.
  
  ! initialize photon variables
  Eph = 0._prec
  th5 = 0._prec
  Eph_cms = 0._prec
  th5_cms = 0._prec

  ! apply photon cuts only if real photon exists
  if (ql5(4) > 0._prec) then
    Eph = ql5(4)
    q5perp = sqrt(ql5(1)**2 + ql5(2)**2)
    th5 = atan2(q5perp, ql5(3))
    phi5_x = atan2(ql5(2),ql5(1))
    phi5_y = atan2(ql5(1),ql5(2))
    Eph_cms = q5(4)
    q5perp_cms = sqrt(q5(1)**2 + q5(2)**2)
    th5_cms = atan2(q5perp_cms, q5(3))
    phi5_x_cms = atan2(q5(2),q5(1))
    phi5_y_cms = atan2(q5(1),q5(2))
    if ((Eph < Eph_cut) .or. (abs(th5).gt.12.e-3))pass_cut = .false.
  endif


  names(1) = 'th3'
  quant(1) = th3
  names(2) = 'Emu'
  quant(2) = Emu
  names(3) = 'th5'
  quant(3) = th5
  names(4) = 'Eph'
  quant(4) = Eph
  names(5) = 'phi5_x'
  quant(5) = phi5_x
  names(6) = 'phi5_y'
  quant(6) = phi5_y
  
  names(7) = 'th3_cms'
  quant(7) = th3_cms
  names(8) = 'Emu_cms'
  quant(8) = Emu_cms
  names(9) = 'th5_cms'
  quant(9) = th5_cms
  names(10) = 'Eph_cms'
  quant(10) = Eph_cms
  names(11) = 'phi5_x_cmEph_cms = q5(4)
    q5perp_cms = sqrt(q5(1)**2 + q5(2)**2)
    th5_cms = atan2(q5perp_cms, q5(3))
    phi5_x_cms = atan2(q5(2),q5(1))
    phi5_y_cms = atan2(q5(1),q5(2))s'
  quant(11) = phi5_x_cms
  names(12) = 'phi5_y_cms'
  quant(12) = phi5_y_cms
 

  END FUNCTION QUANT


  SUBROUTINE USEREVENT(X, NDIM)
  integer :: ndim
  real(kind=prec) :: x(ndim)
  userweight = 1.
  END SUBROUTINE USEREVENT


                 !!!!!!!!!!!!!!!!!!!!!!!
                     END MODULE  USER
                 !!!!!!!!!!!!!!!!!!!!!!!
