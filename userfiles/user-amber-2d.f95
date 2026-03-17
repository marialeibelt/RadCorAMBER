                 !!!!!!!!!!!!!!!!!!!!!
                     MODULE  USER
                 !!!!!!!!!!!!!!!!!!!!!

  use mcmule
  implicit none

!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!

  integer, parameter :: nrq = 21
  integer, parameter :: nrbins = 200
  real(kind=prec), parameter :: &
       min_val(nrq) = (/ 0., -175.,  0.,    0.,     0.,     0.,  0., &
                         0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,     /) 
  real(kind=prec), parameter :: &
       max_val(nrq) = (/ 2.,   25., 25., 2000., 50000., 50000.,  pi, &
                         pi, pi, pi, pi, pi, pi, pi, pi, pi, pi,     /)
  integer :: userdim = 0
  integer :: scenario

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
  read*, scenario
  write(filenamesuffix,'(I1)') scenario
  
  call initflavour("mu-p", Mmu**2+Mproton**2+2*Mproton*100.e3)
  
  print*, "This is a McMule userfile for AMBER"
  print*, " * E_mu = 100 GeV"
  print*, " * 0.3 < th_mu < 2 mrad"
  if(scenario==0) print*, " * S0 :: --"
  if(scenario==1) print*, " * S1 :: E_mu > 90 GeV, 50 MeV < E_y < 5 GeV"
  if(scenario==2) print*, " * S2 :: if E_y > 0.3 GeV then th_y > 16 mrad"
  
  END SUBROUTINE


  FUNCTION QUANT(q1,q2,q3,q4,q5,q6,q7)

  real (kind=prec), intent(in) :: q1(4),q2(4),q3(4),q4(4), q5(4),q6(4),q7(4)
  real (kind=prec) :: ql1(4),ql2(4),ql3(4),ql4(4), ql5(4),ql6(4),ql7(4)  ! in lab frame
  real (kind=prec) :: gah(4), gas(4), gahcm(4), gascm(4), ga2nd(4)
  real (kind=prec) :: th3,ene1,q3perp,th4,q4perp
  real (kind=prec) :: egah,cthyh,cthyyh,thyh,thyyh,phyh,phyyh
  real (kind=prec) :: quant(nr_q)
  !! ==== keep the line below in any case ==== !!
  call fix_mu

  pol1 = (/ 0._prec, 0._prec, 0._prec, 0._prec /)

  ql1 = boost_rf(q2,q1)
  ql2 = boost_rf(q2,q2)
  ql3 = boost_rf(q2,q3)
  ql4 = boost_rf(q2,q4)
  ql5 = boost_rf(q2,q5)
  ql6 = boost_rf(q2,q6)
  ql7 = boost_rf(q2,q7)

  ! define hard (gah) and soft (gas) photons
  if(q5(4) > q6(4)) then
     gahcm = q5;  gascm = q6
  else
     gahcm = q6;  gascm = q5
  endif
  if(ql5(4) > ql6(4)) then
     gah = ql5;  gas = ql6
  else
     gah = ql6;  gas = ql5
  endif

  ga2nd = 0._prec
  if(gas(4) > 10._prec) ga2nd = gas

  ene1 = ql1(4)
  egah = gah(4)

  q3perp=sqrt(ql3(1)**2+ql3(2)**2)
  th3 = 1000.*atan2(q3perp,ql3(3))   ! scattering angle in mrad

  q4perp=sqrt(ql4(1)**2+ql4(2)**2)
  th4 = 1000.*(atan2(q4perp,ql4(3))-0.5*pi)   ! scattering angle-pi/2 in mrad

  cthyh = cos_th(ql1,gah)    
  thyh = acos(cthyh)  ! polar angle of y in proton RF plane in mrad
  phyh = 180./pi*phi(gah)

  pass_cut = .true.

  ! scattering angle cut
  if(th3 < 0.3_prec) pass_cut = .false.
  if(th3 > 2._prec) pass_cut = .false.

  ! scenario S1
  if(scenario==1) then
     if(ql3(4) < 90.e3)  pass_cut =.false.
  endif

  if(scenario==2) then
    if(ql5(4) > 300._prec) then
      if(1000.*atan2(sqrt(ql5(1)**2+ql5(2)**2),ql5(3)) < 16.) pass_cut=.false.
    endif
    if(ql6(4) > 300._prec) then
      if(1000.*atan2(sqrt(ql6(1)**2+ql6(2)**2),ql6(3)) < 16.) pass_cut=.false.
    endif
  endif

  names(1) = 'th3'
  quant(1) = th3
  names(2) = 'th4'
  quant(2) = th4
  names(3) = 'Tp'         !  quant(3) = ene1 - Mm*ene1/(ene1*(1.-cth3)+Mm) 
  quant(3) = ql4(4) - Mm
  names(4) = 'dEl'
  quant(4) = ql1(4)-ql3(4)
  names(5) = 'Q2e'
  quant(5) = -sq(ql1-ql3)
  names(6) = 'Q2p'
  quant(6) = -sq(ql2-ql4)
  
  names(7) = 'thyh'
  quant(7) = thyh
  names(8) = 'thyhB1'
  pass_cut(8) = (50.+(5000.-50.)/10.*0. < egah).and.(egah < 50.+(5000.-50.)/10.*1.)
  quant(8) = thyh
  names(9) = 'thyhB2'
  pass_cut(9) = (50.+(5000.-50.)/10.*1. < egah).and.(egah < 50.+(5000.-50.)/10.*2.)
  quant(9) = thyh
  names(10) = 'thyhB3'
  pass_cut(10) = (50.+(5000.-50.)/10.*2. < egah).and.(egah < 50.+(5000.-50.)/10.*3.)
  quant(10) = thyh
  names(11) = 'thyhB4'
  pass_cut(11) = (50.+(5000.-50.)/10.*3. < egah).and.(egah < 50.+(5000.-50.)/10.*4.)
  quant(11) = thyh
  names(12) = 'thyhB5'
  pass_cut(12) = (50.+(5000.-50.)/10.*4. < egah).and.(egah < 50.+(5000.-50.)/10.*5.)
  quant(12) = thyh
  names(13) = 'thyhB6'
  pass_cut(13) = (50.+(5000.-50.)/10.*5. < egah).and.(egah < 50.+(5000.-50.)/10.*6.)
  quant(13) = thyh
  names(14) = 'thyhB7'
  pass_cut(14) = (50.+(5000.-50.)/10.*6. < egah).and.(egah < 50.+(5000.-50.)/10.*7.)
  quant(14) = thyh
  names(15) = 'thyhB8'
  pass_cut(15) = (50.+(5000.-50.)/10.*7. < egah).and.(egah < 50.+(5000.-50.)/10.*8.)
  quant(15) = thyh
  names(16) = 'thyhB9'
  pass_cut(16) = (50.+(5000.-50.)/10.*8. < egah).and.(egah < 50.+(5000.-50.)/10.*9.)
  quant(16) = thyh
  names(17) = 'thyhB10'
  pass_cut(17) = (50.+(5000.-50.)/10.*9. < egah).and.(egah < 50.+(5000.-50.)/10.*10.)
  quant(17) = thyh

  END FUNCTION QUANT


  SUBROUTINE USEREVENT(X, NDIM)
  integer :: ndim
  real(kind=prec) :: x(ndim)
  userweight = 1.
  END SUBROUTINE USEREVENT


                 !!!!!!!!!!!!!!!!!!!!!!!
                     END MODULE  USER
                 !!!!!!!!!!!!!!!!!!!!!!!
