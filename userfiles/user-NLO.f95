                 !!!!!!!!!!!!!!!!!!!!!
                     MODULE  USER
                 !!!!!!!!!!!!!!!!!!!!!

  use mcmule

  implicit none

!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!

  integer, parameter :: nrq = 34			!th3,Emu,th5,Eph, phi5, +same in cms, 
							!x5,y5,ySlices(11x),xSlices(11x)
  integer, parameter :: nrbins = 500
  real(kind=prec), parameter :: &
       min_val(nrq) = (/1.3e-3,  95.e3, -12.e-3, 1.e3, -12.e-3,&		!rad,MeV,rad,MeV,rad
       			1.3e-3,  95.e3, -12.e-3, 1.e3, -12.e-3,&  		!rad,MeV,rad,MeV,rad
       			-0.5, -0.5,&						!m,m				
       			-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,&!m
       			-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5/)!m
  real(kind=prec), parameter :: &
       max_val(nrq) = (/ 1.7e-3, 101.e3, 12.e-3, 100.e3, 12.e-3,&		!rad,MeV,rad,MeV,rad
       			 1.7e-3, 101.e3, 12.e-3, 100.e3, 12.e-3,&		!rad,MeV,rad,MeV,rad
       			  0.5, 0.5,&						!m,m
       			  0.5, 0.5,0.5, 0.5,0.5, 0.5,0.5, 0.5,0.5, 0.5,0.5,&	!m
       			  0.5, 0.5,0.5, 0.5,0.5, 0.5,0.5, 0.5,0.5, 0.5,0.5/) 	!m
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
  real (kind=prec) :: th3,q3perp,q5perp,th5,Emu,Eph,Eph_cut
  real (kind=prec) :: phi5 ! in lab frame
  real (kind=prec) :: q3perp_cms,th3_cms,Emu_cms,Eph_cms,th5_cms,q5perp_cms,phi5_cms
  real (kind=prec) :: d_detec,x5,y5
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
  d_detec = 10 !m, dummy value!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! dummy value
  
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
    phi5 = atan2(ql5(2),ql5(1))
    Eph_cms = q5(4)
    q5perp_cms = sqrt(q5(1)**2 + q5(2)**2)
    th5_cms = atan2(q5perp_cms, q5(3))
    phi5_cms = atan2(q5(2),q5(1))
    
    x5 = d_detec*tan(th5)*cos(phi5)
    y5 = d_detec*tan(th5)*sin(phi5)
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
  names(5) = 'phi5'
  quant(5) = phi5
  
  names(6) = 'th3_cms'
  quant(6) = th3_cms
  names(7) = 'Emu_cms'
  quant(7) = Emu_cms
  names(8) = 'th5_cms'
  quant(8) = th5_cms
  names(9) = 'Eph_cms'
  quant(9) = Eph_cms
  names(10) = 'phi5_cms'
  quant(10) = phi5_cms
  
  names(11) = 'x5'
  quant(11) = x5
  names(12) = 'y5'
  quant(12) = y5
  
  
  !y Slices
  names(13) = 'x5_B1'
  pass_cut(13) = (-0.2 < y5).and.(y5 < -0.16)
  quant(13) = x5
  names(14) = 'x5_B2'
  pass_cut(14) = (-0.16 < y5).and.(y5 < -0.12)
  quant(14) = x5
  names(15) = 'x5_B3'
  pass_cut(15) = (-0.12 < y5).and.(y5 < -0.08)
  quant(15) = x5
  names(16) = 'x5_B4'
  pass_cut(16) = (-0.08 < y5).and.(y5 < -0.04)
  quant(16) = x5
  names(17) = 'x5_B5'
  pass_cut(17) = (-0.04 < y5).and.(y5 < 0.)
  quant(17) = x5
  names(18) = 'x5_B6'
  pass_cut(18) = (0. < y5).and.(y5 < 0.04)
  quant(18) = x5
  names(19) = 'x5_B7'
  pass_cut(19) = (0.04 < y5).and.(y5 < 0.08)
  quant(19) = x5
  names(20) = 'x5_B8'
  pass_cut(20) = (0.08 < y5).and.(y5 < 0.12)
  quant(20) = x5
  names(21) = 'x5_B9'
  pass_cut(21) = (0.12 < y5).and.(y5 < 0.14)
  quant(21) = x5
  names(22) = 'x5_B10'
  pass_cut(22) = (0.14 < y5).and.(y5 < 0.16)
  quant(22) = x5  
  names(23) = 'x5_B11'
  pass_cut(23) = (0.16 < y5).and.(y5 < 0.2)
  quant(23) = x5
  
  !x Slices
  names(24) = 'y5_B1'
  pass_cut(24) = (-0.2 < x5).and.(x5 < -0.16)
  quant(24) = y5
  names(25) = 'y5_B2'
  pass_cut(25) = (-0.16 < x5).and.(x5 < -0.12)
  quant(25) = y5
  names(26) = 'y5_B3'
  pass_cut(26) = (-0.12 < x5).and.(x5 < -0.08)
  quant(26) = y5
  names(27) = 'y5_B4'
  pass_cut(27) = (-0.08 < x5).and.(x5 < -0.04)
  quant(27) = y5
  names(28) = 'y5_B5'
  pass_cut(28) = (-0.04 < x5).and.(x5 < 0.)
  quant(28) = y5
  names(29) = 'y5_B6'
  pass_cut(29) = (0. < x5).and.(x5 < 0.04)
  quant(29) = y5
  names(30) = 'y5_B7'
  pass_cut(30) = (0.04 < x5).and.(x5 < 0.08)
  quant(30) = y5
  names(31) = 'y5_B8'
  pass_cut(31) = (0.08 < x5).and.(x5 < 0.12)
  quant(31) = y5
  names(32) = 'y5_B9'
  pass_cut(32) = (0.12 < x5).and.(x5 < 0.14)
  quant(32) = y5
  names(33) = 'y5_B10'
  pass_cut(33) = (0.14 < x5).and.(x5 < 0.16)
  quant(33) = y5  
  names(34) = 'y5_B11'
  pass_cut(34) = (0.16 < x5).and.(x5 < 0.2)
  quant(34) = y5
  
  

  

  END FUNCTION QUANT


  SUBROUTINE USEREVENT(X, NDIM)
  integer :: ndim
  real(kind=prec) :: x(ndim)
  userweight = 1.
  END SUBROUTINE USEREVENT


                 !!!!!!!!!!!!!!!!!!!!!!!
                     END MODULE  USER
                 !!!!!!!!!!!!!!!!!!!!!!!
