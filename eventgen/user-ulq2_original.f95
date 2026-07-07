                 !!!!!!!!!!!!!!!!!!!!!
                     MODULE  USER
                 !!!!!!!!!!!!!!!!!!!!!

  use mcmule
  implicit none

!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!==!

  integer, parameter :: nrq = 5
  integer, parameter :: nrbins = 500
  real(kind=prec), parameter :: &
       min_val(nrq) = (/ 47.5, 46., 15.e2, 12.e2, 46./)
  real(kind=prec), parameter :: &
       max_val(nrq) = (/ 52.5, 51., 20.e2, 22.e2, 51./)

  integer :: userdim
  integer :: bspread = 0
  real(kind=prec) :: ebeam

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

  !! ==== Specify the scale mu AND musq==mu**2 ==== !!

  musq = me**2

  END SUBROUTINE FIX_MU



  SUBROUTINE INITUSER
  print*, "ULQ2 with electrons"
  print*, " * beam energy 50 MeV (H, C) or 60 MeV (D)"
  print*, " * momentum cut at 46 MeV (H, C) or 55 MeV (D)"
  print*, " * angular acceptance 50+-2 deg"
  print*, " * photon-inclusive"

  print*, " * e-p scattering"
  call initflavour("e-p-",Mel**2+Mproton**2+2*Mproton*50._prec)
  END SUBROUTINE




  FUNCTION QUANT(q1,q2,q3,q4,q5,q6,q7)

  real (kind=prec), intent(in) :: q1(4),q2(4),q3(4),q4(4), q5(4),q6(4),q7(4)
  real (kind=prec) :: q1Rest(4),q2Rest(4),q3Rest(4),q4Rest(4),q5Rest(4),q6Rest(4)
  real (kind=prec) :: quant(nr_q)
  real (kind=prec) :: thetal, cthetal, pfin, efin, qsql, qsqp, pmin

  !! ==== keep the line below in any case ==== !!
  call fix_mu

  pol1 = (/ 0._prec, 0._prec, 0._prec, 0._prec /)

  pass_cut = .true.

  q1Rest = boost_rf(q2,q1)  ! incomping lepton
  q2Rest = boost_rf(q2,q2)  ! proton at rest
  q3Rest = boost_rf(q2,q3)  ! outgoing lepton
  q4Rest = boost_rf(q2,q4)  ! recoiling proton
  q5Rest = boost_rf(q2,q5)  ! outgoing photon (if present)
  q6Rest = boost_rf(q2,q6)  ! outgoing photon (if present)

  cthetal = cos_th(q1Rest,q3Rest)
  thetal = acos(cos_th(q1Rest,q3Rest))
  pfin = sqrt(q3Rest(1)**2+q3Rest(2)**2+q3Rest(3)**2)
  efin = q3Rest(4)
  qsql = -sq(q1-q3)
  qsqp = -sq(q2-q4)

  if(which_piece(1:5) == "mp2mp" .OR. which_piece(1:5) == "ms2ms") then
    if (thetal.lt. 48*pi/180.) pass_cut=.false.  ! theta > 48 deg
    if (thetal.gt. 52*pi/180.) pass_cut=.false.  ! theta < 52 deg
    pmin = 46._prec                              ! minimum el momentum (MeV)
    if (sqrt(q3Rest(1)**2+q3Rest(2)**2+q3Rest(3)**2) < pmin) pass_cut = .false.
  elseif(which_piece(1:5) == "md2md") then
    if (thetal.lt. 48*pi/180.) pass_cut=.false.  ! theta > 48 deg
    if (thetal.gt. 52*pi/180.) pass_cut=.false.  ! theta < 52 deg
    pmin = 55._prec                              ! minimum el momentum (MeV)
    if (sqrt(q3Rest(1)**2+q3Rest(2)**2+q3Rest(3)**2) < pmin) pass_cut = .false.
  endif

  names(1) = "thetal"
  quant(1) = 180*thetal/pi
  names(2) = "pfin"
  quant(2) = pfin
  names(3) = "qsql"
  quant(3) = qsql
  names(4) = "qsqp"
  quant(4) = qsqp
  names(5) = "efin"
  quant(5) = efin

  END FUNCTION QUANT


  SUBROUTINE USEREVENT(X, NDIM)
  integer :: ndim
  real(kind=prec) :: x(ndim)
  userweight = 1.
  END SUBROUTINE USEREVENT

  !SUBROUTINE USEREVENT(X, NDIM)
  !use integrands, only: fluxfac, xinormcut, xinormcut1, xinormcut2
  !integer :: ndim
  !real(kind=prec) :: x(ndim), z1, z2
  !real(kind=prec), parameter :: e0 = 50.25_prec  ! was 50._prec
  !real(kind=prec), parameter :: sig = e0 * 0.1e-2_prec
  !
  ! using Box-Mueller algorithm to generate 2 normal numbers
  ! with mean=0 and sigma=1
  !z1 = sqrt(-2*log(x(1))) * cos(2*pi*x(2))
  !z2 = sqrt(-2*log(x(1))) * sin(2*pi*x(2))
  !
  !if(bspread==0) then
  !   userdim = 0
  !   userweight = 1.
  !elseif(bspread==1) then
  !   userdim = 2
  !   ebeam = e0 + sig * z1
  !   call initflavour("e-p", Mel**2+Mproton**2+2*Mproton*ebeam)
  !   fluxfac = 0.5/sq_lambda(scms,me,mm)
  !   xieik1 = xinormcut* (1.-(me+mm)**2/scms)
  !   xieik2 = xinormcut* (1.-(me+mm)**2/scms)
  !   xicut1 = xinormcut1*(1.-(me+mm)**2/scms)
  !   xicut2 = xinormcut2*(1.-(me+mm)**2/scms)
  !   userweight = 1.
  !endif  
  !END SUBROUTINE USEREVENT

  FUNCTION EVENT_DISTANCE(Q1,Q2,Q3,Q4,Q5,Q6,Q7, P1,P2,P3,P4,P5,P6,P7)
  real (kind=prec), intent(in) :: q1(4),q2(4),q3(4),q4(4), q5(4),q6(4),q7(4)
  real (kind=prec), intent(in) :: p1(4),p2(4),p3(4),p4(4), p5(4),p6(4),p7(4)
  real (kind=prec) :: event_distance
  real (kind=prec) :: Pfin, Qfin
  real (kind=prec) :: Peta, Qeta
  real (kind=prec) :: Q3lab(4), P3lab(4)

  q3lab = boost_rf(q2,q3)  ! outgoing lepton
  p3lab = boost_rf(p2,p3)  ! outgoing lepton

  Qfin = absvec(q3lab) ; Pfin = absvec(p3lab)
  Qeta = acos(Q3lab(3) / Qfin) * 180 / pi
  Peta = acos(P3lab(3) / Pfin) * 180 / pi

  event_distance = sqrt((Qfin-Pfin)**2 + (Qeta-Peta)**2)

  END FUNCTION

                 !!!!!!!!!!!!!!!!!!!!!!!
                     END MODULE  USER
                 !!!!!!!!!!!!!!!!!!!!!!!
