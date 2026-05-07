                 !!!!!!!!!!!!!!!!!!!!!
                     MODULE  USER
                 !!!!!!!!!!!!!!!!!!!!!

  use mcmule
  implicit none

  integer, parameter :: nrq = 36
  integer, parameter :: nrbins = 500
  
  real(kind=prec), parameter :: min_val(nrq) = (/ &
    0.295e-3_prec, 95.e3_prec, -13e-3_prec, 50._prec, -pi, &	! th3[rad], Emu[MeV], th5[rad], Eph[MeV], phi5[rad]
    0.295e-3_prec,  95.e3_prec, -12e-3_prec, 50._prec, -pi, &  ! th3_cms[rad], Emu_cms[MeV], th5_cms[rad], Eph_cms[MeV], phi5_cms[rad]
    -0.191_prec,-0.191_prec, &                                  ! x5[m], y5[m]
    0._prec,0._prec, &						! ql5(2),ql5(1)
    0.999998_prec, & 						!costh3[]
    0._prec,&                                                	!Qsq in MeV^2
    -0.191_prec,-0.191_prec,-0.191_prec,-0.191_prec,-0.191_prec,-0.191_prec,-0.191_prec,-0.191_prec,-0.191_prec,-0.191_prec, & ! x5_B1..x5_B10[m]
    -0.191_prec,-0.191_prec,-0.191_prec,-0.191_prec,-0.191_prec,-0.191_prec,-0.191_prec,-0.191_prec,-0.191_prec,-0.191_prec /) ! y5_B1..y5_B10[m]

  real(kind=prec), parameter :: max_val(nrq) = (/ &
    2.005e-3_prec, 101.e3_prec,  13.e-3_prec, 101.e3_prec,  pi, & 	! th3[rad], Emu[MeV], th5[rad], Eph[MeV], phi5[rad]
    2.005e-3_prec, 101.e3_prec,  12.e-3_prec, 101.e3_prec,  pi, &  	! th3_cms[rad], Emu_cms[MeV], th5_cms[rad], Eph_cms[MeV], phi5_cms[rad]
    0.191_prec,0.191_prec, &                                    	! x5[m], y5[m]
    650._prec,650._prec, &						! ql5(2),ql5(1)
    0.999999_prec, &							!costh3[]
    5.e4_prec,&                                                     	!Qsq in MeV^2
    0.191_prec,0.191_prec,0.191_prec,0.191_prec,0.191_prec,0.191_prec,0.191_prec,0.191_prec,0.191_prec,0.191_prec, &  ! x5_B1..x5_B10[m]
    0.191_prec,0.191_prec,0.191_prec,0.191_prec,0.191_prec,0.191_prec,0.191_prec,0.191_prec,0.191_prec,0.191_prec /)  ! y5_B1..y5_B10[m]
    
  integer :: userdim = 0
  integer :: namesLen=12
  integer :: filenamesuffixLen=10
  integer :: nq=nrq
  integer :: nbins=nrbins
  integer :: bin_kind = 0       !! 0 for d\sigma/dQ; +1 for Q d\sigma/dQ;

  contains

  SUBROUTINE FIX_MU
  musq = mM**2
  END SUBROUTINE FIX_MU

  SUBROUTINE INITUSER
  print*, "Welcome to Mary's McMule userfile <3"
  print*, "Full Q2 range [0.;0.05] GeV²]"
  print*, " * 0.3 < th_mu < 2. mrad"
  print*, " * Emu > 70 GeV"
  print*, " * Eph > 50MeV"
  print*, " * -12. < th_ph < 12. mrad"
  print*, " * d_detec = 30 m"
  
  call initflavour("mu-p", Mmu**2+Mproton**2+2*Mproton*100.e3)
  END SUBROUTINE

  FUNCTION QUANT(q1,q2,q3,q4,q5,q6,q7)
  real(kind=prec), intent(in) :: q1(4),q2(4),q3(4),q4(4),q5(4),q6(4),q7(4)
  real(kind=prec) :: ql1(4),ql2(4),ql3(4),ql4(4), ql5(4),ql6(4),ql7(4)
  real(kind=prec) :: th3,costh3,q3perp,q5perp,th5,Emu,Eph,Eph_cut
  real(kind=prec) :: phi5,Qsq
  real(kind=prec) :: q3perp_cms,th3_cms,Emu_cms,Eph_cms,th5_cms,q5perp_cms,phi5_cms
  real(kind=prec) :: d_detec,x5,y5
  real(kind=prec) :: quant(nrq)
  integer :: n_bands,i,offset_x,offset_y,last_hist_nr,nr_bandhists
  real(kind=prec) :: band_min, band_max, bin_width
  character(len=3) :: str_i
  real(kind=prec) :: thmu_low, thmu_up,Emu_low

  call fix_mu

  ! proton rest frame / lab frame
  ql1 = boost_rf(q2,q1) !muon in
  ql2 = boost_rf(q2,q2) !proton in
  ql3 = boost_rf(q2,q3) !muon out
  ql4 = boost_rf(q2,q4) !proton out
  ql5 = boost_rf(q2,q5) !photon
  ql6 = boost_rf(q2,q6)
  ql7 = boost_rf(q2,q7)
  
  ! Cuts
  Eph_cut = 200._prec
  !thmu_low = 1.35e-3
  !thmu_up = 1.65e-3
  thmu_low = 0.3e-3
  thmu_up = 2.e-3
  Emu_low = 70.e3
  
  d_detec = 30._prec
  ! proton rest frame / lab frame
  q3perp = sqrt(ql3(1)**2 + ql3(2)**2)
  th3 = atan2(q3perp, ql3(3))
  Emu = ql3(4)
  
  ! cms frame
  q3perp_cms = sqrt(q3(1)**2 + q3(2)**2)
  th3_cms = atan2(q3perp_cms,q3(3))
  Emu_cms = q3(4)

  pass_cut = .true.

  ! Muon cuts  
  if(th3 .lt. thmu_low) pass_cut = .false.
  if(th3 .gt. thmu_up) pass_cut = .false.
  if(Emu .lt. Emu_low) pass_cut = .false.

  costh3 = cos(th3)
  
  Qsq = - ( (ql1(4)-ql3(4))**2 - (ql1(1)-ql3(1))**2 - (ql1(2)-ql3(2))**2 - (ql1(3)-ql3(3))**2 )
  if ((Qsq.lt.0.0).and.(Qsq.gt.5.e4)) pass_cut = .false.
  
  Eph = ql5(4)
  q5perp = sqrt(ql5(1)**2 + ql5(2)**2)
  th5 = atan2(q5perp, ql5(3)) !always >0
  phi5 = atan2(ql5(2),ql5(1))
  Eph_cms = q5(4)
  q5perp_cms = sqrt(q5(1)**2 + q5(2)**2)
  th5_cms = atan2(q5perp_cms, q5(3))
  phi5_cms = atan2(q5(2),q5(1))
  x5 = d_detec*tan(th5)*cos(phi5)
  y5 = d_detec*tan(th5)*sin(phi5) !>0 if phi5>0 and <0 if phi5<0
  if (Eph .gt. Eph_cut) then
    if (abs(th5) .gt. 12.e-3) pass_cut = .false.
  end if          ! <-- closes the Eph_cut check

  if(.not.all(pass_cut)) return
  !Information for x,y bands
  last_hist_nr = 16
  n_bands = 10 !must be an even number	
  nr_bandhists = 2*n_bands

  ! Lab values
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

  ! CMS values
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
  
  names(13) = "ql5(2)"
  quant(13) = ql5(2)
  names(14) = "ql5(1)"
  quant(14) = ql5(1)
  
  names(15) = "costh3"
  quant(15) = costh3
  
  names(16) = "Qsq"
  quant(16) = Qsq

  ! Banded slices
  bin_width = 0.0382_prec !ECal2 with 10x cells with 38.2 mm x 38.2 mm ->active area x&y: [-19.1;19.1]
  band_min = -(n_bands/2.0_prec * bin_width)
  band_max = n_bands/2.0_prec * bin_width
  ! Y slices (x5)
  offset_y = last_hist_nr
  do i=1,n_bands
    write(str_i,'(I0)') i
    names(offset_y+i) = 'x5_B'//trim(str_i)
    pass_cut(offset_y+i) = (ql5(4) .gt. Eph_cut) .and. &
                          (band_min + (i-1)*bin_width .le. y5) .and. &
                          (y5 .lt. band_min + i*bin_width)
    quant(offset_y+i) = x5
  end do
  ! X slices (y5)
  offset_x = offset_y + n_bands
  do i=1,n_bands
    write(str_i,'(I0)') i
    names(offset_x+i) = 'y5_B'//trim(str_i)
    pass_cut(offset_x+i) = (ql5(4) .gt. Eph_cut) .and. &
                            (band_min + (i-1)*bin_width .le. x5) .and. &
                            (x5 .lt. band_min + i*bin_width)
    quant(offset_x+i) = y5
  end do

  END FUNCTION QUANT

  SUBROUTINE USEREVENT(X, NDIM)
  integer :: ndim
  real(kind=prec) :: x(ndim)
  userweight = 1._prec
  END SUBROUTINE USEREVENT

                 !!!!!!!!!!!!!!!!!!!!!!!
                     END MODULE  USER
                 !!!!!!!!!!!!!!!!!!!!!!!
