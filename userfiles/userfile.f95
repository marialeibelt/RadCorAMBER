module user
  use mcmule
  implicit none

  integer, parameter :: nrq = 9
  integer, parameter :: nrbins = 500
  
  real(kind=prec), parameter :: min_val(nrq) = (/ &
    0.07e-3_prec, 0._prec, -13e-3_prec, 0._prec, -pi, &    	      	! th3[rad], Emu[MeV], th5[rad], Eph[MeV], phi5[rad]
    -1._prec, &                                                	 	  ! costh3[]
    1._prec, &                                                  	  ! Qsq in MeV^2
    -0.2_prec,-0.2_prec /)! , &    !                             		! x5[m], y5[m]
  real(kind=prec), parameter :: max_val(nrq) = (/ &
    8.005e-3_prec, 101.e3_prec,  13.e-3_prec, 101.e3_prec,  pi, &   ! th3[rad], Emu[MeV], th5[rad], Eph[MeV], phi5[rad]
    1._prec, &                                                    	! costh3[]
    10.e4_prec, &	                                                 	! Qsq in MeV^2
    0.2_prec,0.2_prec /)                                  	      	! x5[m], y5[m]
    
  integer :: userdim = 0
  integer :: namesLen = 12
  integer :: filenamesuffixLen = 10
  integer :: nq = nrq
  integer :: nbins = nrbins
  integer :: bin_kind = 0       !! 0 for d\sigma/dQ; +1 for Q d\sigma/dQ;

  contains

  subroutine FIX_MU
    musq = mM**2
  end subroutine FIX_MU

  subroutine INITUSER
    print*, "Welcome to Mary's McMule userfile <3"
    print*, " * 0.316 < th_mu < 4.mrad"
    print*, " * no photon cut"
    print*, " * d_detec = 30m"
    print*, " * Ethresh = 5._prec"
    print*, " * Ththresh = 15._prec"
  
    call initflavour("mu-p", Mmu**2+Mproton**2+2*Mproton*100.e3)
  end subroutine INITUSER

  function QUANT(q1,q2,q3,q4,q5,q6,q7)
    real(kind=prec), intent(in) :: q1(4),q2(4),q3(4),q4(4),q5(4),q6(4),q7(4)
    real(kind=prec) :: ql1(4),ql2(4),ql3(4),ql4(4), ql5(4)
    real(kind=prec) :: thmu_low,thmu_up
    real(kind=prec) :: th3,costh3,q3perp,Emu
    real(kind=prec) :: q5perp,th5,phi5,Eph
    real(kind=prec) :: Qsq
    real(kind=prec) :: d_detec,x5,y5
    real(kind=prec) :: quant(nr_q)
    
    integer :: n_bands,i,offset_x,offset_y,last_hist_nr
    real(kind=prec) :: band_min, bin_width
    
    character(len=3) :: str_i    
    character(len=3) :: Qsq_window,thmu_window
  
    call fix_mu

    ! Proton rest frame / lab frame
    ql1 = boost_rf(q2,q1) ! muon in
    ql2 = boost_rf(q2,q2) ! proton in (ungenutzt)
    ql3 = boost_rf(q2,q3) ! muon out
    ql4 = boost_rf(q2,q4) ! proton out (ungenutzt)
    ql5 = boost_rf(q2,q5) ! photon
 
    d_detec = 30._prec
  
    ! Kinematische Variablen berechnen
    q3perp = sqrt(ql3(1)**2 + ql3(2)**2)
    th3    = atan2(q3perp, ql3(3))
    Emu    = ql3(4)
    costh3 = cos(th3)
    Qsq    = -sq(ql1-ql3)
  
    Eph    = ql5(4)
    q5perp = sqrt(ql5(1)**2 + ql5(2)**2)
    th5    = atan2(q5perp, ql5(3)) 
    phi5   = atan2(ql5(2), ql5(1))
    
    x5 = d_detec*tan(th5)*cos(phi5)
    y5 = d_detec*tan(th5)*sin(phi5) !>0 if phi5>0 and <0 if phi5<0
  
    ! --- CUTS AUSWERTEN ---
    pass_cut = .true.

    thmu_low = 0.316e-4_prec
    thmu_up  = 4.000e-3_prec

    if ((th3.lt.thmu_low) .or. (th3.gt.thmu_up)) pass_cut = .false. 
   

    ! --- OBSERVABLE SPEICHERN ---
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
    names(6) = 'costh3'
    quant(6) = costh3
    names(7) = 'Qsq'
    quant(7) = Qsq
    
    names(8) = 'x5'
    quant(8) = x5
    names(9) = 'y5'
    quant(9) = y5
    

  end function QUANT

  subroutine USEREVENT(X, NDIM)
    integer :: ndim
    real(kind=prec) :: x(ndim)
    userweight = 1._prec
  end subroutine USEREVENT

  FUNCTION EVENT_DISTANCE(Q1,Q2,Q3,Q4,Q5,Q6,Q7, P1,P2,P3,P4,P5,P6,P7)
  real (kind=prec), intent(in) :: q1(4),q2(4),q3(4),q4(4), q5(4),q6(4),q7(4)
  real (kind=prec), intent(in) :: p1(4),p2(4),p3(4),p4(4), p5(4),p6(4),p7(4)
  real (kind=prec) :: event_distance, Ediff, Thdiff, th3_q, th3_p, Ethresh, Ththresh
  real (kind=prec) :: q3perp, p3perp
  real (kind=prec) :: Q3lab(4), P3lab(4)

  q3lab = boost_rf(q2,q3)  ! outgoing lepton
  p3lab = boost_rf(p2,p3)  ! outgoing lepton
  q3perp = sqrt(q3lab(1)**2 + q3lab(2)**2)
  p3perp = sqrt(p3lab(1)**2 + p3lab(2)**2)

  th3_q = atan2(q3perp, q3lab(3))
  th3_p = atan2(p3perp, p3lab(3))

  Ediff = abs(q3lab(4)-p3lab(4))
  Thdiff = abs(th3_q - th3_p) *1e6! in μrad
  Ethresh = 5._prec
  Ththresh = 15._prec
  
  event_distance = sqrt((Ediff/Ethresh)**2 + (Thdiff/Ththresh)**2)

  END FUNCTION

end module user
