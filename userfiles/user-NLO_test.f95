module user
  use mcmule
  implicit none

  integer, parameter :: nrq = 7
  integer, parameter :: nrbins = 500
  
  real(kind=prec), parameter :: min_val(nrq) = (/ &
    0.295e-3_prec, 95.e3_prec, -13e-3_prec, 50._prec, -pi, &    ! th3[rad], Emu[MeV], th5[rad], Eph[MeV], phi5[rad]
    -1._prec, &                                                 ! costh3[]
    1._prec /)                                                  ! Qsq in MeV^2

  real(kind=prec), parameter :: max_val(nrq) = (/ &
    2.005e-3_prec, 101.e3_prec,  13.e-3_prec, 101.e3_prec,  pi, &   ! th3[rad], Emu[MeV], th5[rad], Eph[MeV], phi5[rad]
    1._prec, &                                                  ! costh3[]
    10.e4_prec /)                                               ! Qsq in MeV^2
    
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
    !print*, "Big Q2 range [1.e-3;4.e-2] GeV²"
    print*, " * 0.3 < th_mu "!< 2. mrad"
    print*, " * Emu > 70 GeV"
    print*, " * Eph > 10 MeV"
    !print*, " * -12. < th_ph < 12. mrad"
    !print*, " * d_detec = 30 m"
  
    call initflavour("mu-p", Mmu**2+Mproton**2+2*Mproton*100.e3)
  end subroutine INITUSER

  function QUANT(q1,q2,q3,q4,q5,q6,q7)
    real(kind=prec), intent(in) :: q1(4),q2(4),q3(4),q4(4),q5(4),q6(4),q7(4)
    real(kind=prec) :: ql1(4),ql2(4),ql3(4),ql4(4), ql5(4)
    real(kind=prec) :: th5_cut,Eph_cut,thmu_low,thmu_up,Emu_low
    real(kind=prec) :: th3,costh3,q3perp,Emu
    real(kind=prec) :: q5perp,th5,phi5,Eph
    real(kind=prec) :: Qsq
    real(kind=prec) :: quant(nr_q)
    
    character(len=3) :: Qsq_window
  
    call fix_mu

    ! Proton rest frame / lab frame
    ql1 = boost_rf(q2,q1) ! muon in
    ql2 = boost_rf(q2,q2) ! proton in (ungenutzt)
    ql3 = boost_rf(q2,q3) ! muon out
    ql4 = boost_rf(q2,q4) ! proton out (ungenutzt)
    ql5 = boost_rf(q2,q5) ! photon
  
    ! Cut Definitionen
    Eph_cut    = 10._prec 
    Qsq_window = 'BIG'  
    thmu_low   = 0.3e-3_prec 
    !thmu_up    = 2.e-3_prec   
    Emu_low    = 70.e3_prec
    !th5_cut    = 12.e-3_prec

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

    ! --- CUTS AUSWERTEN ---
    pass_cut = .true.

    ! Muon Schnitte
    if (th3 .lt. thmu_low) pass_cut = .false.
    !if (th3 .gt. thmu_up)  pass_cut = .false.
    if (Emu .lt. Emu_low)  pass_cut = .false.
  
    ! Q^2 Windows
    if (Qsq_window .eq. 'BIG') then
      if ((Qsq .lt. 1.e3_prec) .or. (Qsq .gt. 4.e4_prec)) pass_cut = .false.
    endif
    if (Qsq_window .eq. 'SMA') then
      if ((Qsq .lt. 5.e2_prec) .or. (Qsq .gt. 1.e3_prec)) pass_cut = .false.
    endif
    
    ! Photon Cuts
    if (Eph .gt. 1.e-5_prec) then
      if (Eph.lt. Eph_cut) pass_cut = .false.
    endif

      !if (abs(th5) .gt. th5_cut) then
      !  pass_cut = .false.
      !endif
    !endif
   

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

  end function QUANT

  subroutine USEREVENT(X, NDIM)
    integer :: ndim
    real(kind=prec) :: x(ndim)
    userweight = 1._prec
  end subroutine USEREVENT

end module user
