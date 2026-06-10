module user
  use mcmule
  implicit none

  integer, parameter :: nrq = 9 !29
  integer, parameter :: nrbins = 500
  
  real(kind=prec), parameter :: min_val(nrq) = (/ &
    0.07e-3_prec, 0._prec, -13e-3_prec, 0._prec, -pi, &    		! th3[rad], Emu[MeV], th5[rad], Eph[MeV], phi5[rad]
    -1._prec, &                                                	 	! costh3[]
    1._prec, &                                                  	! Qsq in MeV^2
    -0.2_prec,-0.2_prec /) !, &                                  		! x5[m], y5[m]
    ! -0.2_prec,-0.2_prec,-0.2_prec,-0.2_prec,-0.2_prec,-0.2_prec,-0.2_prec,-0.2_prec,-0.2_prec,-0.2_prec, & 	! x5_B1..x5_B10[m]
    ! -0.2_prec,-0.2_prec,-0.2_prec,-0.2_prec,-0.2_prec,-0.2_prec,-0.2_prec,-0.2_prec,-0.2_prec,-0.2_prec /) 	! y5_B1..y5_B10[m]

  real(kind=prec), parameter :: max_val(nrq) = (/ &
    8.005e-3_prec, 101.e3_prec,  13.e-3_prec, 101.e3_prec,  pi, &   	! th3[rad], Emu[MeV], th5[rad], Eph[MeV], phi5[rad]
    1._prec, &                                                  	! costh3[]
    10.e4_prec, &	                                             	! Qsq in MeV^2
    0.2_prec,0.2_prec /) !, &                                    	      	! x5[m], y5[m]
    ! 0.2_prec,0.2_prec,0.2_prec,0.2_prec,0.2_prec,0.2_prec,0.2_prec,0.2_prec,0.2_prec,0.2_prec, &  		! x5_B1..x5_B10[m]
    ! 0.2_prec,0.2_prec,0.2_prec,0.2_prec,0.2_prec,0.2_prec,0.2_prec,0.2_prec,0.2_prec,0.2_prec /)  		! y5_B1..y5_B10[m]
    
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
    print*, "Big Q2 range [1.e-3;4.e-2] GeV²"
    print*, " * 0.3 < th_mu < 2.mrad"
    print*, "No bands but cuts on x5 and y5"
    !print*, " * 1.2 < th_mu < 8.0 mrad"   !4xBigger
    !print*, " * 0.075 < th_mu < 0.5 mrad" !4xSmaller
    print*, " * Emu > 70 GeV"
    print*, " * -12. < th_ph < 12. mrad"
    print*, " * d_detector=30m"
    print*, " * Eph > 200 MeV"
  
    call initflavour("mu-p", Mmu**2+Mproton**2+2*Mproton*100.e3)
  end subroutine INITUSER

  function QUANT(q1,q2,q3,q4,q5,q6,q7)
    real(kind=prec), intent(in) :: q1(4),q2(4),q3(4),q4(4),q5(4),q6(4),q7(4)
    real(kind=prec) :: ql1(4),ql2(4),ql3(4),ql4(4), ql5(4)
    real(kind=prec) :: th5_cut,Eph_cut,thmu_low,thmu_up,Emu_low
    real(kind=prec) :: th3,costh3,q3perp,Emu
    real(kind=prec) :: q5perp,th5,phi5,Eph
    real(kind=prec) :: Qsq,Qsq_low,Qsq_up
    real(kind=prec) :: d_detec,x5,y5,x5_low,x5_up,y5_low,y5_up
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
  
    ! Cut Definitionen
    Qsq_window = 'BIG'  	![1.e-3;4.e-2] GeV²
    !Qsq_window = 'SMA'  
    thmu_window = 'NOR'	!0.3 < th_mu < 2.mrad
    !thmu_window = 'BIG'	!1.2 < th_mu < 8.0 mrad
    !thmu_window = 'SMA'		!0.075 < th_mu < 0.5 mrad
    Emu_low    = 70.e3_prec
    
    th5_cut    = 12.e-3_prec
    Eph_cut    = 200._prec 
    
    
    x5_low = -0.191_prec
    x5_up = 0.191_prec
    y5_low = -0.191_prec
    y5_up = 0.191_prec
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

    ! Muon Windows
    select case (thmu_window)
      case ('NOR')
        thmu_low = 0.3e-3_prec
        thmu_up  = 2.e-3_prec
      case ('BIG')
        thmu_low = 1.2e-3_prec
        thmu_up  = 8.e-3_prec
      case ('SMA')
        thmu_low = 0.075e-3_prec
        thmu_up  = 0.5e-3_prec
    end select

    ! Q^2 Windows
    select case (Qsq_window)
      case ('BIG')
        Qsq_low = 1.e3_prec
        Qsq_up  = 4.e4_prec
      case ('SMA')
        Qsq_low = 5.e2_prec
        Qsq_up  = 1.e3_prec
    end select

    if ((Qsq.lt.Qsq_low) .or. (Qsq.gt.Qsq_up)) pass_cut = .false.
    if ((th3.lt.thmu_low) .or. (th3.gt.thmu_up)) pass_cut = .false.
    if (Emu.lt.Emu_low) pass_cut = .false.
    
    ! Photon Cuts
    if (Eph.gt.Eph_cut) then
      if (abs(th5).gt.th5_cut) pass_cut = .false.
    endif
    
    if ((x5.lt.x5_low).or.(x5.gt.x5_up)) pass_cut = .false.
    !if ((y5.lt.y5_low).or.(y5.gt.y5_up)) pass_cut = .false.
   

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
    
    ! last_hist_nr = 9
    ! n_bands = 10 !must be an even number

    ! ! Banded slices
    ! bin_width = 0.0382_prec !ECal2 with 10x cells with 38.2 mm x 38.2 mm ->active area x&y: [-19.1;19.1]
    ! band_min = -(n_bands/2.0_prec * bin_width)
    
    ! ! Y slices (x5)
    ! offset_y = last_hist_nr
    ! do i=1,n_bands
    !  write(str_i,'(I0)') i
    !  names(offset_y+i) = 'x5_B'//trim(str_i)
    !  pass_cut(offset_y+i) = (ql5(4) .gt. Eph_cut) .and. &
    !                        (band_min + (i-1)*bin_width .le. y5) .and. &
    !                        (y5 .lt. band_min + i*bin_width)
    !  quant(offset_y+i) = x5
    ! end do
    
    ! ! X slices (y5)
    ! offset_x = offset_y + n_bands
    ! do i=1,n_bands
    !  write(str_i,'(I0)') i
    !  names(offset_x+i) = 'y5_B'//trim(str_i)
    !  pass_cut(offset_x+i) = (ql5(4) .gt. Eph_cut) .and. &
    !                          (band_min + (i-1)*bin_width .le. x5) .and. &
    !                          (x5 .lt. band_min + i*bin_width)
    !  quant(offset_x+i) = y5
    ! end do

  end function QUANT

  subroutine USEREVENT(X, NDIM)
    integer :: ndim
    real(kind=prec) :: x(ndim)
    userweight = 1._prec
  end subroutine USEREVENT

end module user
