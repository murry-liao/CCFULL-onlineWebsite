!'''
!@File    :   CC_python.py Ec_version solved by modified Numerove or Discrete Basis method 
!@Time    :   2025/06/22 
!@Author  :   Liao ZeHong, Kouichi Hagino
!@Version :   1.0
!@Contact :   liaozh26@mail2.sysu.edu.cn
!
!gfortran ./ccfull.f90
!'''


Module ccfull_initialization_mod
    Implicit None   
  !  """Principal global variables
  !  P                   - penetrability 
  !  Sigma               - fusion cross section, unit mb
  !  Spin                - mean angular momentum
  !    
  !  Apro, Zpro, Rpro    - atomic #, proton # and radius of the projectile
  !  Atar, Ztar, Rtar    - those of the target
  !  ReduceMass          - reduced mass, unit MeV/c**2
  !  Amu                 - nucleon mass, unit MeV/c**2
  !  E                   - bombarding energy in the center of mass frame, unit MeV
  !    
  !  V0,R0,A0            - depth, range, and dissuseness parameters of uncoupled 
  !                        nuclear potential, which is assumed to be a Woods-Saxon form
  !    
  !  IVIBROTT (IVIBROTP) - option for intrinsic degree of freedom 
  !                        = -1; no excitation (inert)
  !                        =  0; vibrational coupling 
  !                        =  1; rotational coupling
  !
  !  Ntar (Npro)         - the number of levels to be included    
  !
  !  BetaTn (BetaPn)     - i am not very clear about it (Liao)
  !  BetaT (BetaP)       - defeormation parameter
  !  OmegaT (OmegaP)     - excitation energy of the oscillator
  !  LambdaT (LambdaP)   - multipolarity
  !  NphononT (NphononP) - the number of phonon to be included
  !  Beta2T (Beta2P)     - static quadrupole deformation parameter
  !  Beta4T (Beta4P)     - static hexadecapole deformation parameter
  !  Beta6T              - static target hexacontatetrapole deformation parameter
  !  BetaT2n             - is the same as BetaTn but for the second mode
  !  BetaT2              - is the same as BETAT but for the second mode
  !  OmegaT2             - is the same as OmegaT but for the second mode
  !  LambdaT2            - is the same as LambdaT but for the second mode
  !  NphononT2           - is the same as NphononT but for the second mode
  !
  !  E2T (E2P)           - excitation energy of 2+ state in a rotational band
  !  NrotT (NrotP)       - the number of levels in the rotational band to be included 
  !                        (up to I^pi=2*NROT+ states are included)
  !
  !  L                   - angular momentum of the relative motion
  !  CPOT                - coupling matrix
  !  """
  Real(8), parameter :: Hbar = 197.329d0
  Real(8), parameter :: Pi   = 3.141592653d0
  Real(8), parameter :: Amu  = 938.0
  Real(8), parameter :: e2   = 1.44

  integer, parameter :: Nlevelmax = 30
  integer, Public ::  Apro, Zpro, Atar, Ztar, L_i
  Real(8), Public ::  R0P, R0T, Rpro, Rtar, ReduceMass, h2m
  Real(8), Public ::  OmegaT, BetaT, OmegaT2, BetaT2, E2T, Beta2T, Beta4T, Beta6T
  Real(8), Public ::  OmegaP, BetaP, E2P, Beta2P, Beta4P
  Real(8), Public ::  V0, R0, A0, Emin, Emax, dE, R_max, R_min, dr
  Real(8), Public ::  BetaTn, BetaPn, BetaT2n, BetaP2n
  Real(8), Public ::  R_barrier, V_barrier, curv, R_bottom, V_bottom
  Real(8), Public ::  R_barrier_l, V_barrier_l, curv_l, R_bottom_l, V_bottom_l
  Real(8), Public ::  t
  Real,    Public ::  Ftr, Qtrans
  Real,    Public ::  E
  Integer, Public ::  IVIBROTP, IVIBROTT, LambdaT, NphononT, LambdaT2, NphononT2
  Integer, Public ::  NrotT, Ntar, LambdaP, NphononP, NrotP, Npro
  Integer, Public ::  Ntrans, iq
  Integer, Public ::  Nlevel, R_iterat, Nx, Nx_cc, Nx1, Nx2
  Real(8), Allocatable, Public :: Pot_para1(:), Pot_para2(:), Pot_para3(:)
  Real(8), Allocatable, Public :: Ev(:), eps(:)
  Real(8), Allocatable, Public :: omeahv(:), omeahv2(:), omeahvp(:), erott(:), erotp(:)
  Real(8), Allocatable, Public :: betnahv(:,:), betcahv(:,:), betnahv2(:,:), betcahv2(:,:)
  Real(8), Allocatable, Public :: betnahvp(:,:), betcahvp(:,:), bett(:,:), betp(:,:)
  Real(8), Allocatable, Public :: Evec(:,:), CPOT0(:,:), CPOTH(:,:), CPOT(:,:,:)
  Real(8), Allocatable, Public :: H(:,:), unit(:,:)
  ! Snapshot interface for the reduced-basis emulator.  The flattened vector
  ! uses the DBM ordering idx = radial_index*Nlevel + channel, with radial
  ! indices 0:R_iterat+1.  Its first Nx_cc entries are the DBM interior vector.
  Complex(8), Allocatable, Public :: wf(:)
  Real(8), Allocatable, Public :: Numerov_radius(:)
  Real(8), Allocatable, Public :: Echannel(:)
  Real(8), Public :: Numerov_snapshot_penetrability = 0.d0
  Real(8), Public :: Numerov_snapshot_consistency = 0.d0
  Logical, Public :: Numerov_snapshot_ready = .False.

  Integer, Allocatable, Public :: imutual(:,:)
  character(len=2), dimension(0:109) :: Element

  data Element / &
       'O ','H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
       'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ', &
       'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
       'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In', &
       'Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm', &
       'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re', &
       'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra', &
       'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md', &
       'No','Lr','XX','X1','X2','X3','X4','04' /

Contains

  Subroutine Read_Input()
    Character(len=512) :: input_line
    Integer :: input_status

    Open(10, File='ccfull.inp', Status = 'Unknown')
    Open(20, File='sigma.dat',  Status = 'Unknown')
    write(20,'(4A15)') 'E',  'sigma',  '<L>',  'P (L = 0)'
    Open(30, File='Output.dat', Status = 'Unknown')
    Read(10,*)Apro, Zpro, Atar, Ztar
    Read(10,*)R0P, IVIBROTP, R0T, IVIBROTT

    Rpro = R0P * Apro ** (1.d0 / 3.d0)
    Rtar = R0T * Atar ** (1.d0 / 3.d0)
    ReduceMass = Apro * Atar * 1.d0 / (Apro + Atar) * Amu

    h2m = Hbar**2 / (2 * ReduceMass)
    
    If (IVIBROTT == 0) Then

      Read(10,*)OmegaT, BetaT, LambdaT, NphononT
      Read(10,*)OmegaT2, BetaT2, LambdaT2, NphononT2
      Ntar = NphononT

    Else If (IVIBROTT == 1) Then

      Read(10,'(A)') input_line
      Read(input_line,*,Iostat=input_status) E2T, Beta2T, Beta4T, Beta6T, NrotT
      If (input_status /= 0) Then
        ! Backward compatibility with the original four-field rotor record.
        Read(input_line,*) E2T, Beta2T, Beta4T, NrotT
        Beta6T = 0.0d0
      End If
      Read(10,*)
      Ntar = NrotT
    Else
      Read(10,*)
      Read(10,*)
      Ntar = 0
    End If

    If (IVIBROTP == 0) Then

      Read(10,*)OmegaP, BetaP, LambdaP, NphononP
      Npro = NphononP

    Else If (IVIBROTP == 1) Then

      Read(10,*)E2P, Beta2P, Beta4P, NrotP
      Npro = NrotP
    Else
      Read(10,*)
      Npro = 0
    End If
    Nlevel=(Npro + 1)*(Ntar + 1)

    Read(10,*)Qtrans, Ftr
    Read(10,*)V0, R0, A0
    Read(10,*)Emin, Emax, dE
    Read(10,*)R_max, dr
  End Subroutine Read_Input

  Subroutine Output_information()
      Character(len=1) :: ans


      print '(A, I0 ,A, A,I0 ,A)', 'System: ', Int(Apro), Element(Zpro),'+', Int(Atar), Element(Ztar)

      print '(A, F8.3, A, F8.3, A, F8.3, A, F8.3, A)', &
            "Simulation range E from ", Emin, " MeV To ", Emax, " MeV, dE = ", dE, " MeV"

      print '(A, F8.4, A, F8.4, A, F8.3, A, F8.3, A)', &
            "Simulation range R from ", R_min, " fm To ", R_max, " fm, dR = ", dr, " fm"

      print '(A, F8.3, A, F8.3, A, F8.3, A, F8.3, A)', &
            "Potential parameters: V0= ", V0, "(MeV), a= ", A0, "(fm), r0= ", R0, "(fm), R= ",&
             R0*(Apro**(1.0/3.0) + Atar**(1.0/3.0)), "(fm)"

      print '(A, F8.4, A)', "Coulomb barrier position :", R_barrier, "fm"

      print '(A, F8.4, A)', "Coulomb barrier energy   :", V_barrier, "MeV"

      print '(A, F8.4, A)', "Coulomb barrier curv     :", curv, "MeV"

      print '(A, F8.4, A)', "Coulomb bottom position  :", R_bottom, "fm"

      print '(A, F10.4, A)', "Coulomb bottom energy    :", V_bottom, "MeV"

      write(30,'(A, I0 ,A, A,I0 ,A)') 'System: ', Int(Apro), Element(Zpro),'+', Int(Atar), Element(Ztar)

      write(30,'(A, F8.3, A, F8.3, A, F8.3, A, F8.3, A)') &
            "Simulation range E from ", Emin, " MeV To ", Emax, " MeV, dE = ", dE, " MeV"

      write(30,'(A, F8.4, A, F8.4, A, F8.3, A, F8.3, A)') &
            "Simulation range R from ", R_min, " fm To ", R_max, " fm, dR = ", dr, " fm"

      write(30,'(A, F8.3, A, F8.3, A, F8.3, A, F8.3, A)') &
            "Potential parameters: V0= ", V0, "(MeV), a= ", A0, "(fm), r0= ", R0, "(fm), R= ",&
             R0*(Apro**(1.0/3.0) + Atar**(1.0/3.0)), "(fm)"

      write(30,'(A, F8.4, A)') "Coulomb barrier position :", R_barrier, "fm"

      write(30,'(A, F8.4, A)') "Coulomb barrier energy   :", V_barrier, "MeV"

      write(30,'(A, F8.4, A)') "Coulomb barrier curv     :", curv, "MeV"

      write(30,'(A, F8.4, A)') "Coulomb bottom position  :", R_bottom, "fm"

      write(30,'(A, F10.4, A)') "Coulomb bottom energy    :", V_bottom, "MeV"

      if (Ntar /= 0) then
          if (IVIBROTT == 0) then
              write(*,'(A, F6.3, A, F6.3, A, I0, A, I0)') "Phonon Excitation in the targ.: beta=",&
                         BetaT, ", omega=", OmegaT, " (MeV), Lambda=", LambdaT, ", Nph=", NphononT
              BetaTn = BetaT
              Write(*,*)' Different beta_N from beta_C for this mode(n/y)?'
              Read(*,*)ans
              IF(ans .Eq. 'Y' .Or. ans .Eq. 'y') Then
                Write(*,*)'beta_N=?'
                Read(*,*)BetaTn
              End If
              write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the targ.: beta_N=",&
                                                 BetaTn, ", beta_C=", BetaT, ", r0=", R0T, " (fm),"
              write(*,'(A, F6.3, A, I0, A, I0)') "                              omega=", OmegaT, "&
                                                 (MeV), Lambda=", LambdaT, ", Nph=", NphononT
              write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the targ.: beta_N=",&
                                                   BetaTn, ", beta_C=", BetaT, ", r0=", R0T, " (fm),"
              write(30,'(A, F6.3, A, I0, A, I0)') "                              omega=", OmegaT, &
                                                      " (MeV), Lambda=", LambdaT, ", Nph=", NphononT
          else if (IVIBROTT == 1) then
              write(*,'(A, F6.3, A, F6.3, A, F6.3, A, F6.3, A)') &
                    "Rotational Excitation in the targ.: beta2=", Beta2T, ", beta4=", Beta4T, &
                    ", beta6=", Beta6T, ", r0=", R0T, " (fm),"
              write(*,'(A, F6.3, A, I0)') "                                   E2=", E2T, " (MeV), Nrot=", NrotT
              write(30,'(A, F6.3, A, F6.3, A, F6.3, A, F6.3, A)') &
                    "Rotational Excitation in the targ.: beta2=", Beta2T, ", beta4=", Beta4T, &
                    ", beta6=", Beta6T, ", r0=", R0T, " (fm),"
              write(30,'(A, F6.3, A, I0)') "                                   E2=", E2T, " (MeV), Nrot=", NrotT
          end if
      end if

      if (NphononT2 /= 0) then
          write(*,'(A, F6.3, A, F6.3, A, I0, A, I0)') "Phonon Excitation in the targ.: beta=", BetaT2, ", omega=",&
                                                         OmegaT2, " (MeV), Lambda=", LambdaT2, ", Nph=", NphononT2
          BetaT2n = BetaT2
          Write(*,*)' Different beta_N from beta_C for this mode(n/y)?'
          Read(*,*)ans
          IF(ans .Eq. 'Y' .Or. ans .Eq. 'y') Then
            Write(*,*)'beta_N=?'
            Read(*,*)BetaT2n
          End If
          write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the targ.: beta_N=", BetaT2n,&
                                                     ", beta_C=", BetaT2, ", r0=", R0T, " (fm),"
          write(*,'(A, F6.3, A, I0, A, I0)') "                              omega=", OmegaT2, " (MeV), Lambda=",&
                                             LambdaT2, ", Nph=", NphononT2
          write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the targ.: beta_N=", BetaT2n,&
                                                     ", beta_C=", BetaT2, ", r0=", R0T, " (fm),"
          write(30,'(A, F6.3, A, I0, A, I0)') "                              omega=", OmegaT2, " (MeV), Lambda=",&
                                               LambdaT2, ", Nph=", NphononT2
          Call Mutual
      End if

      print *, "------------------------------"
      print *, "Mode of excitation for projectile"

      if (Npro /= 0) then
          if (IVIBROTP == 0) then
              write(*,'(A, F6.3, A, F6.3, A, I0, A, I0)') "Phonon Excitation in the proj.: beta=", BetaP, ", omega=",&
                                                           OmegaP, " (MeV), Lambda=", LambdaP, ", Nph=", NphononP
              BetaPn = BetaP
              Write(*,*)' Different beta_N from beta_C for this mode(n/y)?'
              Read(*,*)ans
              IF(ans .Eq. 'Y' .Or. ans .Eq. 'y') Then
                Write(*,*)'beta_N=?'
                Read(*,*)BetaPn
              End If
              write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the proj.: beta_N=", BetaPn, ", beta_C=",&
                                                         BetaP, ", r0=", R0P, " (fm),"
              write(*,'(A, F6.3, A, I0, A, I0)') "                              omega=", OmegaP, " (MeV), Lambda=",&
                                                                                       LambdaP, ", Nph=", NphononP
              write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the proj.: beta_N=", BetaPn,&
                                                         ", beta_C=", BetaP, ", r0=", R0P, " (fm),"
              write(30,'(A, F6.3, A, I0, A, I0)') "                              omega=", OmegaP, " (MeV), Lambda=",&
                                                   LambdaP, ", Nph=", NphononP
          else if (IVIBROTP == 1) then
              write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Rotational Excitation in the proj.: beta2=", Beta2P,&
                                                         ", beta4=", Beta4P, ", r0=", R0P, " (fm),"
              write(*,'(A, F6.3, A, I0)') "                                   E2=", E2P, " (MeV), Nrot=", NrotP
              write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Rotational Excitation in the proj.: beta2=", Beta2P,&
                                                         ", beta4=", Beta4P, ", r0=", R0P, " (fm),"
              write(30,'(A, F6.3, A, I0)') "                                   E2=", E2P, " (MeV), Nrot=", NrotP
          end if
      end if
      Call Anharmonicity
      Call grotation

      print *, "------------------------------"
      print *, "Transfer coupled mode"

      if (Ntrans /= 0) then
          write(*,'(A, F6.3, A, F6.3, A)') "Transfer channel: Strength= ", Ftr, ", Q = ", Qtrans, " MeV"
          write(30,'(A, F6.3, A, F6.3, A)') "Transfer channel: Strength= ", Ftr, ", Q = ", Qtrans, " MeV"
      end if
      print *, "------------------------------"
  End Subroutine Output_information


  Subroutine Initialize_CCFull()

    Implicit None
    Allocate(imutual(0:Nlevelmax,0:Nlevelmax))
    Allocate(Ev(Nlevelmax))
    Allocate(Evec(Nlevelmax,Nlevelmax))
    Allocate(betnahv(0:Nlevelmax,0:Nlevelmax))
    Allocate(betcahv(0:Nlevelmax,0:Nlevelmax))
    Allocate(omeahv(0:Nlevelmax))
    Allocate(betnahv2(0:Nlevelmax,0:Nlevelmax))
    Allocate(betcahv2(0:Nlevelmax,0:Nlevelmax))
    Allocate(omeahv2(0:Nlevelmax+1))
    Allocate(betnahvp(0:Nlevelmax,0:Nlevelmax))
    Allocate(betcahvp(0:Nlevelmax,0:Nlevelmax))
    Allocate(omeahvp(0:Nlevelmax))
    Allocate(bett(0:Nlevelmax,0:Nlevelmax))
    Allocate(erott(0:Nlevelmax+1))
    Allocate(betp(0:Nlevelmax,0:Nlevelmax))
    Allocate(erotp(0:Nlevelmax))
    Allocate(eps(0:Nlevelmax))

    Call Read_Input()

    Call PotShape(0, R_barrier, V_barrier, curv, R_bottom, V_bottom)
    R_iterat = int((R_max - R_bottom) / dr)
    R_min = R_bottom
    R_max = R_bottom + dr * R_iterat
    t = Hbar**2 / (2 * ReduceMass * dr**2)
    
    Allocate(CPOT(Nlevelmax,Nlevelmax,0:R_iterat+2))
    Allocate(CPOT0(Nlevelmax,Nlevelmax))
    Allocate(CPOTH(Nlevelmax,Nlevelmax))
    Allocate(Echannel(Nlevelmax))

    Call Output_information()

    ! Keep the same vector dimensions and ordering as the DBM implementation.
    Nx = R_iterat
    Nx_cc = Nlevel * Nx
    Nx1 = Nx_cc + Nlevel
    Nx2 = Nx_cc + 2 * Nlevel

  End Subroutine Initialize_CCFull


End Module ccfull_initialization_mod


Subroutine PotShape(l, rb, vb, curv_out, rmin, vmin)
  Use ccfull_initialization_mod
  
  Implicit None
  Integer :: n, i
  Integer, Intent(In) :: l
  Real(8), Intent(Out) :: rb, vb, curv_out, rmin, vmin
  Real(8) :: r, u0, u1, ra, rb_tmp, tolk, u, ddv00
  Real(8) :: dV_dr, V
  external dV_dr, V
  r = 50.5D0
  u0 = dV_dr(r, l)

  Do
    r = r - 1.0D0
    If (r < 0.0D0) Then
      rb = -5.0D0
      vb = -5.0D0
      curv_out = -1.0D0
      rmin = -1.0D0
      vmin = -1.0D0
      Return
    End If
    u1 = dV_dr(r, l)
    If (u0 * u1 <= 0.0D0) Exit
    u0 = u1
  End Do

  ra = r + 1.0D0
  rb_tmp = r
  tolk = 1.0D-6
  n = Int(Log10(Abs(rb_tmp - ra) / tolk) / Log10(2.0D0) + 0.5D0)

  Do i = 1, n
    r = (ra + rb_tmp) / 2.0D0
    u = dV_dr(r, l)
    If (u0 * u < 0.0D0) Then
      rb_tmp = r
    Else
      ra = r
    End If
  End Do

  rb = r
  ddv00 = (dV_dr(rb + 1.0D-5, l) - dV_dr(rb - 1.0D-5, l)) / (2.0D-5)

  If (ddv00 > 0.0D0) Then
    Print *, "Error: Second derivative positive at barrier top."
    Stop
  End If

  curv_out = Hbar * Sqrt(Abs(ddv00) / ReduceMass)
  vb = V(rb, l)

  ra = rb - 0.5D0
  u0 = dV_dr(ra, l)
  rb_tmp = 0.5D0
  u1 = dV_dr(rb_tmp, l)
  n = Int(Log10(Abs(rb - ra) / tolk) / Log10(2.0D0) + 0.5D0)

  Do i = 1, n
    r = (ra + rb_tmp) / 2.0D0
    u = dV_dr(r, l)
    If (u0 * u < 0.0D0) Then
      rb_tmp = r
    Else
      ra = r
    End If
  End Do

  rmin = r
  vmin = V(rmin, l)

End Subroutine PotShape

Real(8) Function Vn(r)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Real(8) :: R12
  R12 = R0 * (Apro**(1.0D0/3.0D0) + Atar**(1.0D0/3.0D0))
  Vn = -V0 / (1.0D0 + Exp((r - R12) / A0))
End Function Vn

Real(8) Function W(r)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Real(8) :: rw
  rw = 1.0D0 * (Apro**(1.0D0/3.0D0) + Atar**(1.0D0/3.0D0))
  W = 50.0D0 / (1.0D0 + Exp((r - rw) / 0.3D0))
End Function W

Real(8) Function dVn_dr(r)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Real(8) :: R12, exp_term
  R12 = R0 * (Apro**(1.0D0/3.0D0) + Atar**(1.0D0/3.0D0))
  exp_term = Exp((r - R12) / A0)
  dVn_dr = V0 / A0 * exp_term / (1.0D0 + exp_term)**2
End Function dVn_dr

Real(8) Function dVcent_dr(r, l)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Integer, Intent(In) :: l
  dVcent_dr = -2.0D0 * l * (l + 1.0D0) * Hbar**2 / (2.0D0 * ReduceMass * r**3)
End Function dVcent_dr

Real(8) Function dVc_dr(r)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  dVc_dr = -Zpro * Ztar * Hbar / 137.0D0 / r**2
End Function dVc_dr

Real(8) Function dV_dr(r, l)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Integer, Intent(In) :: l
  Real(8)             :: dVc_dr, dVn_dr, dVcent_dr
  external dVc_dr, dVn_dr, dVcent_dr
  dV_dr = dVn_dr(r) + dVc_dr(r) + dVcent_dr(r, l)
End Function dV_dr

Real(8) Function Vc(r)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Vc = Zpro * Ztar * Hbar / 137.0D0 / r
End Function Vc

Real(8) Function Vr(r, l)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Integer, Intent(In) :: l
  Vr = h2m * l * (l + 1.0D0) / r**2
End Function Vr

Real(8) Function V(r, l)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Integer, Intent(In) :: l
  Real(8)             :: Vn, Vc, Vr
  external Vn, Vc, Vr
  V = Vn(r) + Vc(r) + Vr(r, l)
End Function V

Real(8) Function VnCC(r, Xt)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r, Xt
  Real(8) :: R12
  R12 = R0 * (Apro**(1.0D0/3.0D0) + Atar**(1.0D0/3.0D0))
  VnCC = -V0 / (1.0D0 + Exp((r - R12 - Xt) / A0))
End Function VnCC

Real(8) Function Fct(r)
  Use ccfull_initialization_mod
  IMPLICIT NONE
  Real(8), Intent(In) :: r
  Real(8) ::  result, Lambda

  Lambda = LambdaT
  If (r > Rtar) Then
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r)**Lambda
  Else
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar)**Lambda
  End If

  result = result * BetaT / Sqrt(4.D0 * Pi)

  Fct = result
  Return
End Function Fct

Real(8) Function Fct2(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer, Parameter :: Lambda = 2

  If (r > Rtar) Then
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r) ** Lambda
  Else
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar) ** Lambda
  End If

  Fct2 = result
  Return
End Function Fct2

Real(8) Function Fct3(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer, Parameter :: Lambda = 3

  If (r > Rtar) Then
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r) ** Lambda
  Else
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar) ** Lambda
  End If

  result = result * (Beta4T + 7.D0 * Beta2T**2 / 7.D0 / Sqrt(Pi))

  Fct3 = result
  Return
End Function Fct3

Real(8) Function Fct4(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer, Parameter :: Lambda = 4

  If (r > Rtar) Then
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r) ** Lambda
  Else
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar) ** Lambda
  End If

  result = result * (Beta4T + 9.D0 * Beta2T**2 / 7.D0 / Sqrt(Pi))

  Fct4 = result
  Return
End Function Fct4

Real(8) Function Fct6(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer, Parameter :: Lambda = 6

  If (r > Rtar) Then
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r) ** Lambda
  Else
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar) ** Lambda
  End If

  ! Linear Coulomb Y_60 coupling.  The nuclear beta6 coupling is kept to
  ! all orders through the deformed Woods-Saxon eigenchannel construction
  ! in coupled_matrix0/VnCC, consistently with the Morton et al. analysis.
  Fct6 = result * Beta6T
  Return
End Function Fct6

Real(8) Function Fcp(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result, Lambda

  Lambda = LambdaP

  If (r > Rpro) Then
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rpro / r) ** Lambda
  Else
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rpro * (r / Rpro) ** Lambda
  End If

  result = result * BetaP / SQRT(4.D0 * Pi)

  Fcp = result
  Return
End Function Fcp

Real(8) Function Fcp2(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = 2

  If (r > Rpro) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rpro / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rpro * (r / Rpro)**Lambda
  End If

  Fcp2 = result
  Return
End Function Fcp2

Real(8) Function Fcp3(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = 3
  result = 0.D0

  If (r > Rpro) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rpro / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rpro * (r / Rpro)**Lambda
  End If

  result = result * (Beta4P + 7.D0 * Beta2P**2 / 7.D0 / Sqrt(Pi))

  Fcp3 = result
  Return
End Function Fcp3

Real(8) Function Fcp4(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = 4
  result = 0.D0

  If (r > Rpro) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rpro / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rpro * (r / Rpro)**Lambda
  End If

  result = result * (Beta4P + 9.D0 * Beta2P**2 / 7.D0 / Sqrt(pi))

  Fcp4 = result
  Return
End Function Fcp4

Real(8) Function Fctt(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = LambdaT2
  result = 0.D0

  If (r > Rtar) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar)**Lambda
  End If

  result = result * BetaT2 / SQRT(4.D0 * Pi)

  Fctt = result
  Return
End Function Fctt

Real(8) Function Fct2v(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = 2
  result = 0.D0

  If (r > Rtar) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar)**Lambda
  End If

  Fct2v = result
  Return
End Function Fct2v

Real(8) Function Fcp2v(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = 2
  result = 0.D0

  If (r > Rpro) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rpro / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rpro * (r / Rpro)**Lambda
  End If

  Fcp2v = result
  Return
End Function Fcp2v

Real(8) Function Ftrans(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result, dVn_dr
  external dVn_dr

  result = 0.D0
  result = Ftr * dVn_dr(r)

  Ftrans = result
  Return
End Function Ftrans

Real(8) Function CG(j1, m1, j2, m2, j3, m3)
    ! Clebsch-Gordan coefficient <j1 m1 j2 m2 | j3 m3>
    implicit real*8 (a-h, o-z)
    external fact

    if (m1 + m2 .ne. m3) then
        cg = 0.d0
        return
    end if

    if (j3 .lt. abs(j1 - j2)) then
        cg = 0.d0
        return
    end if

    if (j3 .gt. j1 + j2) then
        cg = 0.d0
        return
    end if

    ka = j1 + j2 - j3
    kb = j3 + j1 - j2
    kc = j2 + j3 - j1
    kd = j1 + j2 + j3 + 1

    del = sqrt(fact(ka) * fact(kb) * fact(kc) / fact(kd))

    cg = 0.d0
    do n = 0, max(j1 + j2 - j3, j1 - m1, j2 + m2)
        ka1 = j1 + j2 - j3 - n
        if (ka1 .lt. 0.d0) cycle

        ka2 = j3 - j2 + m1 + n
        if (ka2 .lt. 0.d0) cycle

        ka3 = j3 - j1 - m2 + n
        if (ka3 .lt. 0.d0) cycle

        ka4 = j1 - m1 - n
        if (ka4 .lt. 0.d0) cycle

        ka5 = j2 + m2 - n
        if (ka5 .lt. 0.d0) cycle

        cg = cg + (-1.d0)**n / &
                 (fact(n) * fact(ka1) * fact(ka2) * fact(ka3) * fact(ka4) * fact(ka5))
    end do

    cg = cg * sqrt(fact(j1 + m1) * fact(j1 - m1))
    cg = cg * sqrt(fact(j2 + m2) * fact(j2 - m2))
    cg = cg * sqrt(fact(j3 + m3) * fact(j3 - m3))

    cg = cg * sqrt(2.d0 * j3 + 1.d0) * del

    return
End Function CG

Subroutine mdiag(a,n,d,v)

	implicit real*8(a-h,o-z)
	parameter(nmax=100)
  parameter (nlevelmax=30)
  dimension a(nlevelmax,nlevelmax)
	dimension b(nmax),z(nmax),d(nlevelmax),v(nlevelmax,nlevelmax)

	do 12 ip=1,n
	do 11 iq=1,n
	v(ip,iq)=0.d0
 11	continue
	v(ip,ip)=1.d0
 12	continue

	do 13 ip=1,n
	b(ip)=a(ip,ip)
	d(ip)=b(ip)
	z(ip)=0.d0
 13	continue

	nrot=0

	do 24 i=1,50
	sm=0.d0

	do 15 ip=1,n-1
	do 14 iq=ip+1,n
	sm=sm+abs(a(ip,iq))
 14	continue
 15     continue

	if(sm.eq.0.d0) return

	if(i.lt.4) then
	tresh=0.2d0*sm/n**2
	else
	tresh=0.d0
	endif

	do 22 ip=1,n-1
	do 21 iq=ip+1,n
	g=100.d0*abs(a(ip,iq))
	
	if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq)))) then
	   a(ip,iq)=0.d0
	elseif(abs(a(ip,iq)).gt.tresh) then
	   h=d(iq)-d(ip)

	   if(abs(h)+g.eq.abs(h)) then
	     t=a(ip,iq)/h
	   else
	     theta=0.5d0*h/a(ip,iq)
	     t=1.d0/(abs(theta)+sqrt(1.d0+theta**2))
	     if(theta.lt.0.d0) t=-t
	   endif

	c=1.d0/sqrt(1.d0+t**2)
	s=t*c
	tau=s/(1.d0+c)
	h=t*a(ip,iq)
	z(ip)=z(ip)-h
	z(iq)=z(iq)+h
	d(ip)=d(ip)-h
	d(iq)=d(iq)+h
	a(ip,iq)=0.d0

	do 16 j=1,ip-1
	g=a(j,ip)
	h=a(j,iq)
	a(j,ip)=g-s*(h+g*tau)
	a(j,iq)=h+s*(g-h*tau)
 16	continue

	do 17 j=ip+1,iq-1
	g=a(ip,j)
	h=a(j,iq)
	a(ip,j)=g-s*(h+g*tau)
	a(j,iq)=h+s*(g-h*tau)
 17	continue

	do 18 j=iq+1,n
	g=a(ip,j)
	h=a(iq,j)
	a(ip,j)=g-s*(h+g*tau)
	a(iq,j)=h+s*(g-h*tau)
 18	continue

	do 19 j=1,n
	g=v(j,ip)
	h=v(j,iq)
	v(j,ip)=g-s*(h+g*tau)
	v(j,iq)=h+s*(g-h*tau)
 19	continue
	
	nrot=nrot+1
	endif

 21	continue
 22	continue

	do 23 ip=1,n
	b(ip)=b(ip)+z(ip)
	d(ip)=b(ip)
	z(ip)=0.d0

 23	continue
 24	continue

	return

End Subroutine  mdiag

SUBROUTINE matinv(nmax, c, d)
  IMPLICIT NONE
  
  ! Argument declarations
  INTEGER, INTENT(IN) :: nmax
  COMPLEX(8), DIMENSION(30,30), INTENT(INOUT) :: c, d
  
  ! Parameter and local variable declarations
  COMPLEX(8) :: u(30), v(30), t, a, b
  REAL(8) :: deter
  INTEGER :: m, n, j, k, l
  
  ! Initialize
  deter = 1.0_8
  d = (0.0_8, 0.0_8)
  DO m = 1, nmax
      d(m, m) = (1.0_8, 0.0_8)
  END DO
  
  ! Main inversion loop
  DO n = 1, nmax
      t = c(n, n)
      
      ! Check for zero pivot element
      IF (ABS(t) < 1.0E-10_8) THEN
          j = n
          DO
              j = j + 1
              IF (j > nmax) THEN
                  PRINT *, 'matrix not invertible'
                  RETURN
              END IF
              
              t = c(n, j)
              deter = -deter
              IF (ABS(t) > 1.0E-10_8) EXIT
          END DO
          
          ! Swap rows
          DO k = 1, nmax
              u(k) = c(n, k)
              v(k) = d(n, k)
              c(n, k) = c(j, k)
              d(n, k) = d(j, k)
              c(j, k) = u(k)
              d(j, k) = v(k)
          END DO
      END IF
      
      ! Eliminate column
      DO k = 1, nmax
          IF (k == n) CYCLE
          
          a = c(k, n) / c(n, n)
          DO l = 1, nmax
              c(k, l) = c(k, l) - a * c(n, l)
              d(k, l) = d(k, l) - a * d(n, l)
          END DO
      END DO
      
      ! Normalize row
      b = c(n, n)
      deter = b * deter
      DO m = 1, nmax
          c(n, m) = c(n, m) / b
          d(n, m) = d(n, m) / b
      END DO
  END DO
  
END SUBROUTINE matinv


Real(8) Function fact(n)
    Implicit none
    Integer n, i

    If (n .lt. 0) Then
        fact = 0.d0
        Return
    End If

    If (n .eq. 0) Then
        fact = 1.d0
        Return
    End If

    fact = 1.d0
    Do i = 1, n
        fact = fact * i * 1.d0
    End Do

    Return
End Function fact

Subroutine grotation()
    Use ccfull_initialization_mod
    Implicit None
    Integer :: it, jt, ip, jp, i
    Character(len=1) :: ans
    Real(8) :: erott0, erotp0

    Do it = 0, Ntar
      Do jt = 0, Ntar
        bett(it,jt) = 0.0D0
        If (it == jt - 1 .Or. it == jt + 1 .Or. it == jt) bett(it,jt) = Beta2T
      End Do
      erott(it) = 2.0D0 * it * (2.0D0 * it + 1.0D0) / 6.0D0 * E2T
    End Do

    If (IVIBROTT == 1 .And. Ntar > 1) Then
      Print *, ''
      Print *, 'Generalised E2 couplings for the target rotor (n/y)?'
      Read(*,'(A)') ans
      If (ans == 'y' .Or. ans == 'Y') Then
        Print *, 'Beta2_targ=', Beta2T
        Do it = 1, Ntar - 1
          jt = it + 1
          Print *, 'Modify beta2 for transition from', 2*it, '+ to', 2*jt, '+ ? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            Print *, 'BETA2=?'
            Read(*,*) bett(it,jt)
            bett(jt,it) = bett(it,jt)
            Print *, 'Beta2 set to: ', bett(it,jt)
          End If
        End Do

        Print *, 'Reorientation terms:'
        Do it = 1, Ntar
          jt = it
          Print *, 'Modify beta2 for transition from', 2*it, '+ to', 2*jt, '+ ? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            Print *, 'BETA2=?'
            Read(*,*) bett(it,jt)
            bett(jt,it) = bett(it,jt)
          End If
        End Do

        Print *, 'Excitation energies:'
        Do i = 2, Ntar
          Print *, 'Energy of the', 2*i, '+ state for pure rotor =', erott(i)
          Print *, 'Modify this energy? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            erott0 = erott(i)
            Print *, 'Energy=?'
            Read(*,*) erott(i)
            Print *, 'Modified from', erott0, ' to ', erott(i)
          End If
        End Do
      End If
    End If

    Do ip = 0, Npro
      Do jp = 0, Npro
        betp(ip,jp) = 0.0D0
        If (ip == jp - 1 .Or. ip == jp + 1 .Or. ip == jp) betp(ip,jp) = Beta2P
      End Do
      erotp(ip) = 2.0D0 * ip * (2.0D0 * ip + 1.0D0) / 6.0D0 * E2P
    End Do

    If (IVIBROTP == 1 .And. Npro > 1) Then
      Print *, 'Generalised E2 couplings for the projectile rotor (n/y)?'
      Read(*,'(A)') ans
      If (ans == 'y' .Or. ans == 'Y') Then
        Print *, 'Beta2_proj=', Beta2P
        Do ip = 1, Npro - 1
          jp = ip + 1
          Print *, 'Modify beta2 for transition from', 2*ip, '+ to', 2*jp, '+ ? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            Print *, 'BETA2=?'
            Read(*,*) betp(ip,jp)
            betp(jp,ip) = betp(ip,jp)
          End If
        End Do

        Print *, 'Reorientation terms:'
        Do ip = 1, Npro
          jp = ip
          Print *, 'Modify beta2 for transition from', 2*ip, '+ to', 2*jp, '+ ? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            Print *, 'BETA2=?'
            Read(*,*) betp(ip,jp)
            betp(jp,ip) = betp(ip,jp)
          End If
        End Do

        Print *, 'Excitation energies:'
        Do i = 2, Npro
          Print *, 'Energy of the', 2*i, '+ state for pure rotor =', erotp(i)
          Print *, 'Modify this energy? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            erotp0 = erotp(i)
            Print *, 'Energy=?'
            Read(*,*) erotp(i)
            Print *, 'Modified from', erotp0, ' to ', erotp(i)
          End If
        End Do
      End If
    End If

End Subroutine grotation

Subroutine Mutual()
    Use ccfull_initialization_mod
    Implicit None
    Integer :: i, j
    Character(len=1) :: ans, sign1, sign2
    Integer :: imut

    Do i = 0, Nlevelmax
      Do j = 0, Nlevelmax
        imutual(i,j) = 0
      End Do
    End Do

    Print *, ''
    Print *, 'Mutual excitations in the *target* nucleus'
    Print *, 'Include the mutual excitations (y/n)?'
    Read(*,'(A)') ans

    If (ans == 'n' .Or. ans == 'N') Then
      imut = 0
      Nlevel = (Ntar + NphononT2 + 1) * (Npro + 1)
      Do i = 0, Ntar
        Do j = 0, NphononT2
          If (i == 0 .Or. j == 0) imutual(i,j) = 1
        End Do
      End Do
      Return
    End If

    Print *, 'All the possible mutual excitation channels (n/y)?'
    Read(*,'(A)') ans

    If (ans == 'y' .Or. ans == 'Y') Then
      imut = 1
      Nlevel = (Ntar + 1) * (NphononT2 + 1) * (Npro + 1)
      Do i = 0, Ntar
        Do j = 0, NphononT2
          imutual(i,j) = 1
        End Do
      End Do
      Return
    End If

    imut = 2
    If (Mod(LambdaT, 2) == 0) Then
      sign1 = '+'
    Else
      sign1 = '-'
    End If
    If (Mod(LambdaT2, 2) == 0) Then
      sign2 = '+'
    Else
      sign2 = '-'
    End If

    Do i = 0, Ntar
      Do j = 0, NphononT2
        If (i == 0 .Or. j == 0) Then
          imutual(i,j) = 1
        Else
          Print *, 'Include (', LambdaT, sign1, '^', i, ',', LambdaT2, sign2, '^', j, ') state ? (y/n)'
          Read(*,'(A)') ans
          If (ans == 'n' .Or. ans == 'N') Then
            imutual(i,j) = 0
          Else
            imutual(i,j) = 1
          End If
        End If
      End Do
    End Do

    Print *, 'Excited states in the target to be included:'
    Nlevel = 0

    Do i = 0, Ntar
      Do j = 0, NphononT2
        If (imutual(i,j) == 0) Cycle
        Print *, '(', LambdaT, sign1, '^', i, ',', LambdaT2, sign2, '^', j, ') state'
        Nlevel = Nlevel + 1
      End Do
    End Do

    Nlevel = Nlevel * (Npro + 1)

End Subroutine Mutual

Subroutine Anharmonicity()
  Use ccfull_initialization_mod
  Implicit None

  Character(1) :: ans, sign
  Integer :: it, jt, i

  ! Target: AHV Couplings (1st mode)
  Do it = 0, Ntar
    Do jt = 0, Ntar
      betnahv(it, jt) = 0.0
      betcahv(it, jt) = 0.0

      If (it == jt - 1) Then
        betnahv(it, jt) = BetaTn * Sqrt(Real(jt))
        betcahv(it, jt) = BetaT  * Sqrt(Real(jt))
      Else If (jt == it - 1) Then
        betnahv(it, jt) = BetaTn * Sqrt(Real(it))
        betcahv(it, jt) = BetaT  * Sqrt(Real(it))
      End If
    End Do
    omeahv(it) = it * OmegaT
  End Do

  If (IVIBROTT == 0 .And. Ntar > 1) Then
    Write(*, '(A)') 'AHV couplings for the first mode in the target phonon (n/y)?'
    Read(*, '(A1)') ans
    If (ans == 'Y' .Or. ans == 'y') Then
      Write(*, '(A)') '**** AHV Couplings in the target (the 1st mode)'
      If ((-1.)**LambdaT ==  1) sign = '+'
      If ((-1.)**LambdaT == -1) sign = '-'

      Do it = 0, Ntar
        Do jt = 0, Ntar
          If (it > jt)Cycle
          If ((it == 0 .And. jt <= 1)) Cycle
          Write(*, *) ' '
          Write(*, '(A,I2,2A,I2,A,I2,2A,I2,A)') 'Transition from the', LambdaT, sign, '^', it,&
                                                 ' to the', LambdaT, sign, '^', jt, 'state:'
          Write(*, '(A,2G15.5)') '   beta_N and beta_C in the HO limit=', betnahv(it, jt), betcahv(it, jt)
          Write(*, '(2A)') '    Modify these beta_N a/o beta_C (n/y)?'
          Read(*, '(A1)') ans
          If (ans == 'Y' .Or. ans == 'y') Then
            Write(*, *) '   beta_N and beta_C =?'
            Read(*, *) betnahv(it, jt), betcahv(it, jt)
            Write(*, '(G15.5, G15.5)') betnahv(it, jt), betcahv(it, jt)
          End If
        End Do
      End Do

      Write(*, *) 'Excitation energy for the first mode:'
      Do i = 2, Ntar
        Write(*, '(A,I2,A,G15.5)') '    Energy of the', i, '-phonon state in the HO=', i * OmegaT
        Write(*, *) '   Modify this energy(n/y)?'
        Read(*, '(A1)') ans
        If (ans == 'Y' .Or. ans == 'y') Then
          Write(*, *) '   Energy=?'
          Read(*, *) omeahv(i)
          Write(*, '(A,I2,A,G15.5,A,G12.5,A)') 'Energy of the', i, '-phonon state=', omeahv(i), '(HO: ', i * OmegaT, ')'
        End If
      End Do
    End If
  End If


  ! Target: AHV Couplings (2st mode)
  Do it = 0, NphononT2
    Do jt = 0, NphononT2
      betnahv2(it, jt) = 0.0
      betcahv2(it, jt) = 0.0

      If (it == jt - 1) Then
        betnahv2(it, jt) = BetaT2n * Sqrt(Real(jt))
        betcahv2(it, jt) = BetaT2  * Sqrt(Real(jt))
      Else If (jt == it - 1) Then
        betnahv(it, jt) = BetaT2n * Sqrt(Real(it))
        betcahv(it, jt) = BetaT2  * Sqrt(Real(it))
      End If
    End Do
    omeahv2(it) = it * OmegaT2
  End Do

  If (NphononT2 > 1) Then
    Write(*, '(A)') 'AHV couplings for the second mode in the target phonon (n/y)?'
    Read(*, '(A1)') ans
    If (ans == 'Y' .Or. ans == 'y') Then
      Write(*, '(A)') '**** AHV Couplings in the target (the 1st mode)'
      If ((-1.)**LambdaT2 ==  1) sign = '+'
      If ((-1.)**LambdaT2 == -1) sign = '-'

      Do it = 0, NphononT2
        Do jt = 0, NphononT2
          If (it > jt)Cycle
          If ((it == 0 .And. jt <= 1)) Cycle
          Write(*, *) ' '
          Write(*, '(A,I2,2A,I2,A,I2,2A,I2,A)') 'Transition from the', LambdaT2, sign, '^', it,&
                                                     ' to the', LambdaT2, sign, '^', jt, 'state:'
          Write(*, '(A,2G15.5)') '   beta_N and beta_C in the HO limit=', betnahv2(it, jt), betcahv2(it, jt)
          Write(*, '(2A)') '    Modify these beta_N a/o beta_C (n/y)?'
          Read(*, '(A1)') ans
          If (ans == 'Y' .Or. ans == 'y') Then
            Write(*, *) '   beta_N and beta_C =?'
            Read(*, *) betnahv2(it, jt), betcahv2(it, jt)
            Write(*, '(G15.5, G15.5)') betnahv2(it, jt), betcahv2(it, jt)
          End If
        End Do
      End Do

      Write(*, *) 'Excitation energy for the second mode:'
      Do i = 2, Ntar
        Write(*, '(A,I2,A,G15.5)') '    Energy of the', i, '-phonon state in the HO=', i * OmegaT2
        Write(*, *) '   Modify this energy(n/y)?'
        Read(*, '(A1)') ans
        If (ans == 'Y' .Or. ans == 'y') Then
          Write(*, *) '   Energy=?'
          Read(*, *) omeahv2(i)
          Write(*, '(A,I2,A,G15.5,A,G12.5,A)') 'Energy of the', i, '-phonon state=', omeahv2(i), '(HO: ', i * OmegaT2, ')'
        End If
      End Do
    End If
  End If


  Do it = 0, Npro
    Do jt = 0, Npro
      betnahvp(it, jt) = 0.0
      betcahvp(it, jt) = 0.0

      If (it == jt - 1) Then
        betnahvp(it, jt) = BetaPn * Sqrt(Real(jt))
        betcahvp(it, jt) = BetaP  * Sqrt(Real(jt))
      Else If (jt == it - 1) Then
        betnahvp(it, jt) = BetaPn * Sqrt(Real(it))
        betcahvp(it, jt) = BetaP  * Sqrt(Real(it))
      End If
    End Do
    omeahvp(it) = it * OmegaP
  End Do

  If (IVIBROTP == 0 .And. Npro > 1) Then
      Write(*, '(A)') 'AHV couplings for the first mode in the projectile phonon (n/y)?'
      Read(*, '(A1)') ans
      If (ans == 'Y' .Or. ans == 'y') Then
        Write(*, '(A)') '**** AHV Couplings in the projectile'
        If ((-1.)**LambdaP ==  1) sign = '+'
        If ((-1.)**LambdaP == -1) sign = '-'

        Do it = 0, Npro
          Do jt = 0, Npro
            If (it > jt)Cycle
            If ((it == 0 .And. jt <= 1)) Cycle
            Write(*, *) ' '
            Write(*, '(A,I2,2A,I2,A,I2,2A,I2,A)') 'Transition from the', LambdaP, sign, '^', it, ' to the', LambdaP, sign,&
                                                   '^', jt, 'state:'
            Write(*, '(A,2G15.5)') '   beta_N and beta_C in the HO limit=', betnahvp(it, jt), betcahvp(it, jt)
            Write(*, '(2A)') '    Modify these beta_N a/o beta_C (n/y)?'
            Read(*, '(A1)') ans
            If (ans == 'Y' .Or. ans == 'y') Then
              Write(*, *) '   beta_N and beta_C =?'
              Read(*, *) betnahvp(it, jt), betcahvp(it, jt)
              Write(*, '(G15.5, G15.5)') betnahvp(it, jt), betcahvp(it, jt)
            End If
          End Do
        End Do

        Write(*, *) 'Excitation energy for the first mode:'
        Do i = 2, Npro
          Write(*, '(A,I2,A,G15.5)') '    Energy of the', i, '-phonon state in the HO=', i * OmegaP
          Write(*, *) '   Modify this energy(n/y)?'
          Read(*, '(A1)') ans
          If (ans == 'Y' .Or. ans == 'y') Then
            Write(*, *) '   Energy=?'
            Read(*, *) omeahvp(i)
            Write(*, '(A,I2,A,G15.5,A,G12.5,A)') 'Energy of the', i, '-phonon state=', omeahvp(i), '(HO: ', i * OmegaP, ')'
          End If
        End Do
      End If
    End If

End Subroutine Anharmonicity

Subroutine coupled_matrix0()
  Use ccfull_initialization_mod
  implicit none
  real(8) :: A(Nlevelmax, Nlevelmax)
  Real(8) :: C, CG
  integer :: i, j, ip, it, it2, jp, jt, jt2
  external CG

  !initialization of A matrix
  Do i = 1, Nlevelmax
    Do j = 1, Nlevelmax
     A(i,j) = 0.0d0
    End Do
  End Do

  If((Nlevel - Ntrans) == 1)Return
  i = 0
  Do ip = 0, Npro
    Do it = 0, Ntar
      Do it2 = 0, NphononT2

        if (NphononT2 .Ne. 0)Then
          if(imutual(it, it2) == 0)Cycle
        End if

        i = i + 1 
        j = 0
        Do jp = 0, Npro
          Do jt = 0, Ntar
            Do jt2 = 0, NphononT2

              if (NphononT2 .Ne. 0)Then
                if(imutual(jt, jt2) == 0)Cycle
              End if

              j = j + 1 

              If (i > j)Then
                A(i,j) = A(j,i)
                Cycle
              End If

              C = 0.0d0

              if (ip == jp .and. it2 == jt2) then
                  if (IVIBROTT == 0) then
                      C = Rtar * betnahv(it, jt) / sqrt(4.0d0 * pi)

                  else
                      C = Rtar * bett(it, jt) * sqrt((2 * 2 * it + 1) * 5.0d0 * (2 * 2 * jt + 1) / (4.0d0 * pi)) &
                          * CG(2 * it, 0, 2, 0, 2 * jt, 0)**2 / (2.0d0 * 2 * jt + 1)

                      C = C + Rtar * Beta4T * sqrt((2 * 2 * it + 1) * 9.0d0 * (2 * 2 * jt + 1) / (4.0d0 * pi)) &
                          * CG(2 * it, 0, 4, 0, 2 * jt, 0)**2 / (2.0d0 * 2 * jt + 1)

                      C = C + Rtar * Beta6T * sqrt((2 * 2 * it + 1) * 13.0d0 * (2 * 2 * jt + 1) / (4.0d0 * pi)) &
                          * CG(2 * it, 0, 6, 0, 2 * jt, 0)**2 / (2.0d0 * 2 * jt + 1)
                  end if
              end if

              if (ip == jp .and. it == jt) then
                  C = C + Rtar * betnahv2(it2, jt2) / sqrt(4.0d0 * pi)
              end if

              if (it == jt .and. it2 == jt2) then
                  if (IVIBROTP == 0) then
                      C = C + Rpro * betnahvp(ip, jp) / sqrt(4.0d0 * pi)
                  else
                      C = C + Rpro * betp(ip, jp) * sqrt((2 * 2 * ip + 1) * 5.0d0 * (2 * 2 * jp + 1) / (4.0d0 * pi)) &
                          * CG(2 * ip, 0, 2, 0, 2 * jp, 0)**2 / (2.0d0 * 2 * jp + 1)

                      C = C + Rpro * Beta4P * sqrt((2 * 2 * ip + 1) * 9.0d0 * (2 * 2 * jp + 1) / (4.0d0 * pi)) &
                          * CG(2 * ip, 0, 4, 0, 2 * jp, 0)**2 / (2.0d0 * 2 * jp + 1)
                  end if
              end if

              A(i,j) = C


            End Do
          End Do
        End Do

      End Do
    End Do
  End Do

  Call mdiag(A, Nlevel - Ntrans, Ev, Evec)

  Return

End Subroutine coupled_matrix0

Subroutine coupled_matrix(r, cpot_matrix)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  real(8), intent(out) :: cpot_matrix(Nlevelmax, Nlevelmax)
  Real(8) :: C, A, C0
  Real(8) :: Fct, Fct2, Fct3, Fct4, Fct6, Fct2v, Fctt
  Real(8) :: Ftrans, Fcp, Fcp2, Fcp3, Fcp4, Fcp2v
  Real(8) :: CG, Vn, VnCC
  integer :: i, j, k, ip, it, it2, jp, jt, jt2
  external CG, Fct, Fct2, Fct3, Fct4, Fct6, Fct2v, Fctt
  external Ftrans, Fcp, Fcp2, Fcp3, Fcp4, Fcp2v, Vn, VnCC
  
  Do i = 1, nlevelmax
    eps(i) = 0.
  End Do

  Do i = 1, Nlevelmax
    Do j = 1, Nlevelmax
     cpot_matrix(i,j) = 0.0d0
    End Do
  End Do

  If((Nlevel - Ntrans) == 1)Return

  i = 0
  Do ip = 0, Npro
    Do it = 0, Ntar
      Do it2 = 0, NphononT2

        if (NphononT2 .Ne. 0)Then
          if(imutual(it, it2) == 0.)Cycle
        End if

        i = i + 1 
        j = 0
        Do jp = 0, Npro
          Do jt = 0, Ntar
            Do jt2 = 0, NphononT2

            if (NphononT2 .Ne. 0)Then
              if(imutual(jt, jt2) == 0.)Cycle
            End if

            j = j + 1 

            If (i > j)Then

              cpot_matrix(i,j) = cpot_matrix(j,i)
              Cycle
            End If

            C = 0.0d0
            ! nuclear coupling
            Do k = 1, Nlevel-Ntrans
              C = C + VnCC(r, Ev(k))*Evec(i, k)*Evec(j, k)

            End Do

            ! ---- Target contribution ----
            if (ip == jp .and. it2 == jt2) then
                if (IVIBROTT == 0) then
                    if (it /= jt) then
                        A = betcahv(it, jt) * Fct(r)
                        if (BetaT /= 0.d0) A = A / BetaT
                        
                    else
                        A = betcahv(it, jt) * Fct2v(r) / sqrt(4.d0 * pi)

                    end if
                    C = C + A
                else
                    C = C + sqrt((2 * 2 * it + 1) * 5.d0 * (2 * 2 * jt + 1) / (4.d0 * pi)) &
                          * CG(2 * it, 0, 2, 0, 2 * jt, 0)**2 / (2.d0 * 2 * jt + 1) * Fct2(r) &
                          * (bett(it, jt) + 2.d0 * sqrt(5.d0 / pi) * bett(it, jt)**2 / 7.d0)

                    C = C + sqrt((2 * 2 * it + 1) * 9.d0 * (2 * 2 * jt + 1) / (4.d0 * pi)) &
                          * CG(2 * it, 0, 4, 0, 2 * jt, 0)**2 / (2.d0 * 2 * jt + 1) * Fct4(r)

                    C = C + sqrt((2 * 2 * it + 1) * 13.d0 * (2 * 2 * jt + 1) / (4.d0 * pi)) &
                          * CG(2 * it, 0, 6, 0, 2 * jt, 0)**2 / (2.d0 * 2 * jt + 1) * Fct6(r)
                end if
            end if

            ! ---- Projectile contribution ----
            if (it == jt .and. it2 == jt2) then
                if (IVIBROTP == 0) then
                    if (ip /= jp) then
                        A = betcahvp(ip, jp) * Fcp(r)
                        if (BetaP /= 0.d0) A = A / BetaP
                    else
                        A = betcahvp(ip, jp) * Fcp2v(r) / sqrt(4.d0 * pi)
                    end if
                    C = C + A
                else
                    C = C + sqrt((2 * 2 * ip + 1) * 5.d0 * (2 * 2 * jp + 1) / (4.d0 * pi)) &
                          * CG(2 * ip, 0, 2, 0, 2 * jp, 0)**2 / (2.d0 * 2 * jp + 1) * Fcp2(r) &
                          * (betp(ip, jp) + 2.d0 * sqrt(5.d0 / pi) * betp(ip, jp)**2 / 7.d0)

                    C = C + sqrt((2 * 2 * ip + 1) * 9.d0 * (2 * 2 * jp + 1) / (4.d0 * pi)) &
                          * CG(2 * ip, 0, 4, 0, 2 * jp, 0)**2 / (2.d0 * 2 * jp + 1) * Fcp4(r)
                end if
            end if

            ! ---- BetaT2 part (vibrational excitation of target) ----
            if (ip == jp .and. it == jt) then
                if (it2 /= jt2) then
                    A = betcahv2(it2, jt2) * Fctt(r)
                    if (BetaT2 /= 0.d0) A = A / BetaT2
                else
                    A = betcahv2(it2, jt2) * Fct2v(r) / sqrt(4.d0 * pi)
                end if
                C = C + A
            end if

            ! ---- Excitation energy shift ----
            If (it == jt .and. ip == jp .and. it2 == jt2) then
                if (IVIBROTT == 0) then
                    C = C + omeahv(it)
                    eps(i) = eps(i) + omeahv(it)
                else
                    C = C + erott(it)
                    eps(i) = eps(i) + erott(it)
                end if

                if (IVIBROTP == 0) then
                    C = C + omeahvp(ip)
                    eps(i) = eps(i) + omeahvp(ip)
                else
                    C = C + erotp(ip)
                    eps(i) = eps(i) + erotp(ip)
                end if

                C = C + omeahv2(it2)
                eps(i) = eps(i) + omeahv2(it2)
            End if

            cpot_matrix(i, j) = C


            End Do
          End Do
        End Do

      End Do
    End Do
  End Do

  C0 = cpot_matrix(1,1)
  Do i = 1, Nlevel-Ntrans
    cpot_matrix(i,i) = cpot_matrix(i,i) - Vn(r)
    !cpot_matrix(i,i) = cpot_matrix(i,i) - C0
  End Do

  ! transfer coupling
  If (Ntrans == 1) Then
      cpot_matrix(1, Nlevel) = Ftrans(r)
      cpot_matrix(Nlevel, 1) = Ftrans(r)
      cpot_matrix(Nlevel, Nlevel) = cpot_matrix(1, 1) - Qtrans + (Zpro + iq) * (Ztar - iq) &
                                    / r * Hbar / 137.0 - Zpro * Ztar / r * Hbar / 137.0
  End If

  Return

End Subroutine coupled_matrix

Subroutine rkutta00(psi0, phi0, psi1)
  Use ccfull_initialization_mod
  Implicit None
  complex*16, Intent(In)  ::  psi0(Nlevelmax), phi0(Nlevelmax)
  complex*16, intent(out) ::  psi1(Nlevelmax)
  complex*16 :: ai
  complex*16, dimension(Nlevelmax) :: ak1, ak2, ak3, ak4
  complex*16, dimension(Nlevelmax) :: bk1, bk2, bk3, bk4

  Real(8) ::  fac, r, rh, rpp, V
  Integer ::  j1, i0, ic, is
  external V

  ai = (0.d0, 1.d0)
  fac = dr * (2.d0 * ReduceMass / Hbar**2)

  r  = R_min
  rh = R_min + dr/2.d0
  rpp= R_min + dr
  Do j1 = 1, Nlevel
      ak1(j1) = 0.0D0
      ak2(j1) = 0.0D0
      ak3(j1) = 0.0D0
      ak4(j1) = 0.0D0
      bk1(j1) = 0.0D0
      bk2(j1) = 0.0D0
      bk3(j1) = 0.0D0
      bk4(j1) = 0.0D0
  End Do

  Do i0 = 1, Nlevel
      Do ic = 1, Nlevel
        ak1(i0) = ak1(i0) + fac * CPOT(i0, ic, 0) * psi0(ic)
      End Do
      ak1(i0) = ak1(i0) - fac * (E - V(r, L_i)) * psi0(i0)
      bk1(i0) = dr * phi0(i0)
  End Do

  Do i0 = 1, Nlevel
      Do ic = 1, Nlevel
        ak2(i0) = ak2(i0) + fac * CPOTH(i0, ic) * (psi0(ic) + 1.0D0 / 2.0D0 * bk1(ic))
      End Do
      ak2(i0) = ak2(i0) - fac * (E - V(rh, L_i)) * (psi0(i0) + 1.0D0 / 2.0D0 * bk1(i0))
      bk2(i0) = dr * (phi0(i0) + 1.0D0 / 2.0D0 * ak1(i0))
  End Do

  Do i0 = 1, Nlevel
      Do ic = 1, Nlevel
        ak3(i0) = ak3(i0) + fac * CPOTH(i0, ic) * (psi0(ic) + 1.0D0 / 2.0D0 * bk2(ic))
      End Do
      ak3(i0) = ak3(i0) - fac * (E - V(rh, L_i)) * (psi0(i0) + 1.0D0 / 2.0D0 * bk2(i0))
      bk3(i0) = dr * (phi0(i0) + 1.0D0 / 2.0D0 * ak2(i0))
  End Do

  Do i0 = 1, Nlevel
      Do ic = 1, Nlevel
        ak4(i0) = ak4(i0) + fac * CPOT(i0, ic, 1) * (psi0(ic) + bk3(ic))
      End Do
      ak4(i0) = ak4(i0) - fac * (E - V(rpp, L_i)) * (psi0(i0) + bk3(i0))
      bk4(i0) = dr * (phi0(i0) + ak3(i0))
  End Do

  Do is = 1, nlevel
      psi1(is) = psi0(is) + (1.0D0 / 6.0D0) * (bk1(is) + 2.0D0 * bk2(is) + 2.0D0 * bk3(is) + bk4(is))
  End Do

  Return

End Subroutine rkutta00

Subroutine stabilize(xi1, xi, aa, ir1)
  Use ccfull_initialization_mod
  Implicit None
  Integer, Intent(In)  ::ir1
  complex*16, Intent(Inout)  ::  aa(Nlevelmax, Nlevelmax)
  complex*16, Intent(Inout)  ::  xi(Nlevelmax, Nlevelmax)
  complex*16, Intent(Inout)  ::  xi1(Nlevelmax, Nlevelmax)
  complex*16 :: ai
  complex*16, dimension(Nlevelmax, Nlevelmax)  :: psi
  complex*16, dimension(Nlevelmax, Nlevelmax)  :: cc, cin
  complex*16, dimension(Nlevelmax, Nlevelmax)  :: aa0
  complex*16, dimension(Nlevelmax, Nlevelmax)  :: xid, xid1
  Real(8) ::  r, V, fac
  Integer ::  i, j, k, i0, ic, ich
  external  V

  Do i = 1, nlevel
      Do j = 1, nlevel
        aa0(i, j)  = aa(i, j)
        xid(i, j)  = xi(i, j)
        xid1(i, j) = xi1(i, j)
      End Do
  End Do

  ai = (0.d0, 1.d0)
  fac = dr**2 * (2.d0 * ReduceMass / Hbar**2)
  r = R_min + ir1*dr

  Do i0 = 1, Nlevel
    Do ic = 1, Nlevel
        cc(i0, ic) = -fac / 12.0D0 * CPOT(i0, ic, ir1)
        If (i0 == ic) Then
          cc(i0, ic) = cc(i0, ic) - fac / 12.0D0 * (V(r, L_i) - E) + 1.0D0
        End If
    End Do
  End Do

  Call matinv(Nlevel,cc,cin)

  Do ich = 1, nlevel
      Do i0 = 1, nlevel
        psi(i0, ich) = 0.0D0
        Do ic = 1, nlevel
            psi(i0, ich) = psi(i0, ich) + cin(i0, ic) * xi(ic, ich)
        End Do
      End Do
  End Do

  Do i = 1, nlevel
      Do j = 1, nlevel
        aa(i, j) = 0.0D0
        Do k = 1, nlevel
            aa(i, j) = aa(i, j) + psi(i, k) * aa0(k, j)
        End Do
      End Do
  End Do

  call matinv(Nlevel,psi,cin)

  Do i = 1, Nlevel
      Do j = 1, Nlevel
        xi(i, j)  = 0.0D0
        xi1(i, j) = 0.0D0
        Do k = 1, Nlevel
            xi(i, j)  = xi(i, j) + xid(i, k)  * cin(k, j)
            xi1(i, j) = xi1(i, j) + xid1(i, k) * cin(k, j)
        End Do
      End Do
  End Do

  Return
End Subroutine stabilize

Subroutine Numerov_core(P, capture_snapshot)
  Use ccfull_initialization_mod
  Implicit None
  real(8), intent(out) :: P
  logical, intent(in) :: capture_snapshot
  complex*16, dimension(Nlevelmax) :: psi, psi0, psi1
  complex*16, dimension(Nlevelmax, Nlevelmax) :: xi, xi0, xi1
  complex*16, dimension(Nlevelmax) :: phi0
  complex*16, dimension(Nlevelmax, Nlevelmax) :: bb, bin
  complex*16, dimension(Nlevelmax, Nlevelmax) :: bb2
  complex*16, dimension(Nlevelmax, Nlevelmax) :: cc, cin
  complex*16, dimension(Nlevelmax, Nlevelmax) :: dd0, dd1
  complex*16, dimension(Nlevelmax) :: dd
  real(8), dimension(0:200) :: fcw, gcw, fpcw, gpcw
  real(8), dimension(0:200) :: sigmad
  integer, dimension(0:200) :: iexp
  real(8), dimension(Nlevelmax) :: ech, ech2
  complex*16 :: k, kk, k2
  complex*16 :: ai
  complex*16 :: cwup0, cwdown0, cwup1, cwdown1
  complex*16 :: dummy, dummy2
  complex*16, dimension(Nlevelmax, Nlevelmax) :: xi1d, xi0d
  complex*16, dimension(Nlevelmax, Nlevelmax) :: aa, bb0, bb20
  complex*16, allocatable :: fundamental(:,:,:)
  Integer ::  i, j, lc, io, io2, i0, ir, ii, ik, ic, ich, ibarrier, idx
  Integer ::  j1, j2
  Real(8) ::  V, fac, r, r1, r2
  Real(8) ::  ak, ec, eta, rho
  external V

  If (capture_snapshot) Then
    Numerov_snapshot_ready = .False.
    Allocate(fundamental(Nlevel,Nlevel,0:R_iterat+1))
    fundamental = (0.d0,0.d0)
  End If

  ai = (0.d0, 1.d0)
  fac = dr**2 * (2.d0 * ReduceMass / Hbar**2)
  ibarrier = (R_barrier - R_min)/dr
  Do lc = 0, 200
    fcw(lc)   = 0.0d0
    gcw(lc)   = 0.0d0
    fpcw(lc)  = 0.0d0
    gpcw(lc)  = 0.0d0
    sigmad(lc) = 0.0d0
    iexp(lc)   = 0
  End Do

  Do i = 1, Nlevelmax
      Do j = 1, Nlevelmax
          bb(i,j)  = 0.0d0
          bin(i,j) = 0.0d0
          cc(i,j)  = 0.0d0
          cin(i,j) = 0.0d0
          aa(i,j)  = 0.0d0
      End Do
      dd(i) = 0.0d0
  End Do

  Do i = 1, Nlevel
      aa(i,i) = 1.0d0
  End Do


  Do io = 1, Nlevel

    Do j1 = 1, Nlevel
        psi(j1)  = 0.0d0
        psi0(j1) = 0.0d0
        psi1(j1) = 0.0d0
        phi0(j1) = 0.0d0
    End Do

    If (io == 1) Then
        Do io2 = 1, Nlevel
            ech(io2) = E - V(R_min, L_i) - CPOT(io2, io2, 0)
            ech2(io2) = E - CPOT(io2, io2, R_iterat)
            Echannel(io2) = ech(io2)

        End Do
    End If

    If (ech(io) > 0.d0) Then
        k = sqrt(2.d0 * ReduceMass / Hbar**2 * ech(io))
        psi0(io) = exp(-ai * k * R_min)
        phi0(io) = -ai * k * psi0(io)

    Else
        k = sqrt(2.d0 * ReduceMass / hbar**2 * abs(ech(io)))
        psi0(io) = exp(k * R_min)
        phi0(io) = k * psi0(io)
    End If


    call rkutta00(psi0,phi0,psi1)

    do i0 = 1, Nlevel
        xi0(i0, io) = (1.d0 - fac/12.d0 * (V(R_min, L_i) - e)) * psi0(i0)
        xi1(i0, io) = (1.d0 - fac/12.d0 * (V(R_min+dr, L_i) - e)) * psi1(i0)
        do ic = 1, Nlevel
            xi0(i0, io) = xi0(i0, io) - fac/12.d0 * CPOT(i0, ic, 0) * psi0(ic)
            xi1(i0, io) = xi1(i0, io) - fac/12.d0 * CPOT(i0, ic, 1) * psi1(ic)
        end do
    end do
  End Do

  If (capture_snapshot) Then
    Call CaptureNumerovFundamental(0,xi0,aa,fundamental)
    Call CaptureNumerovFundamental(1,xi1,aa,fundamental)
  End If

  !----------------------------------------------  iterations start
  Do ir = 2, R_iterat + 1
    r = R_min + dr * ir
    !r0 = R_min + dr * (ir - 2)
    r1 = R_min + dr * (ir - 1)

    ! Standard Numerov in the xi = (1 - T) u representation, T = dr**2 Q / 12:
    !   xi(n+1) = [12 (1 - T(n))**(-1) - 10] xi(n) - xi(n-1).
    ! This is the stencil reduced by the CCFULL-CC-ROM emulator.  Hagino's
    ! modified Numerov used (sqrt(3) + dr**2 Q / sqrt(12))**2 - 1, which
    ! truncates (1 - T)**(-1) after second order in T.
    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            dd0(i0, ic) = -fac / 12.d0 * CPOT(i0, ic, ir - 1)
            If (i0 == ic) Then
                dd0(i0, ic) = dd0(i0, ic) + 1.d0 - fac / 12.d0 * (V(r1, L_i) - E)
            End If
        End Do
    End Do

    Call matinv(Nlevel, dd0, dd1)
    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            dd1(i0, ic) = 12.d0 * dd1(i0, ic)
            If (i0 == ic) Then
                dd1(i0, ic) = dd1(i0, ic) - 10.d0
            End If
        End Do
    End Do

    Do ich = 1, Nlevel
        Do i0 = 1, Nlevel
            xi(i0, ich) = -xi0(i0, ich)
            Do ic = 1, nlevel
                xi(i0, ich) = xi(i0, ich) + dd1(i0, ic) * xi1(ic, ich)
            End Do
        End Do
    End Do

    If (capture_snapshot) Call CaptureNumerovFundamental(ir,xi,aa,fundamental)

    If(ir == R_iterat+1)Exit
    If(ir == ibarrier) Call stabilize(xi1,xi,aa,ir)
    

    Do ich = 1, nlevel
        Do i0 = 1, nlevel
            xi0(i0, ich) = xi1(i0, ich)
            xi1(i0, ich) = xi(i0, ich)
        End Do
    End Do


  End Do


  !--------------------------------------------------------------
  !  matching to the coulomb wave function at rmax

  Do io = 1, Nlevel

    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            cc(i0, ic) = -fac / 12.d0 * CPOT(i0, ic, R_iterat - 1)
            If (i0 == ic) Then
                cc(i0, ic) = cc(i0, ic) - fac / 12.d0 * (V(R_max - dr, L_i) - E) + 1.d0
            End If
        End Do
    End Do


    Call matinv(Nlevel,cc,cin)
    Do i0 = 1, Nlevel
        psi0(i0) = 0.d0
        Do ic = 1, Nlevel
            psi0(i0) = psi0(i0) + cin(i0, ic) * xi0(ic, io)
        End Do
    End Do


    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            cc(i0, ic) = -fac / 12.d0 * CPOT(i0, ic, R_iterat + 1)
            If (i0 == ic) Then
                cc(i0, ic) = cc(i0, ic) - fac / 12.d0 * (V(R_max + dr, L_i) - E) + 1.d0
            End If
        End Do
    End Do

    Call matinv(Nlevel,cc,cin)

    Do i0 = 1, Nlevel
        psi(i0) = 0.d0
        Do ic = 1, Nlevel
            psi(i0) = psi(i0) + cin(i0, ic) * xi(ic, io)
        End Do
    End Do


    Do ii = 1, Nlevel
        ! coulomb wave function
        ec = E - CPOT(ii, ii, R_iterat - 1)

        If (ec < 0.d0) Then
            r1 = R_max - dr
            r2 = R_max + dr
            ak = sqrt(2.d0 * ReduceMass * abs(ec) / Hbar**2)
            
            bb0(ii, io) = (exp(-ak * r2) * psi0(ii) - exp(-ak * r1) * psi(ii)) &
                          / (exp(ak * (r1 - r2)) - exp(-ak * (r1 - r2)))
            bb20(ii, io) = -(exp(ak * r2) * psi0(ii) - exp(ak * r1) * psi(ii)) &
                          / (exp(ak * (r1 - r2)) - exp(-ak * (r1 - r2)))
        Else
            rho = sqrt(2.d0 * ReduceMass * ec) / hbar * (R_max - dr)
            eta = (Zpro * Ztar / 137.d0) * Sqrt(ReduceMass / (2.d0 * ec))


            !Call dfcoul(eta, rho, fcw, fpcw, gcw, gpcw, sigmad, L_i, iexp)
            Call myCOULFG(L_i*1.d0, eta, rho, fcw(L_i), fpcw(L_i), gcw(L_i), gpcw(L_i), iexp(L_i))
  
 
            cwup0 = (gcw(L_i) + ai * fcw(L_i))
            cwdown0 = gcw(L_i) - ai * fcw(L_i)

            ec = E - CPOT(ii, ii, R_iterat + 1)
            rho = Sqrt(2.d0 * ReduceMass * ec) / Hbar * (R_max + dr)
            eta = (Zpro * Ztar / 137.d0) * Sqrt(ReduceMass / (2.d0 * ec))
            !Call dfcoul(eta, rho, fcw, fpcw, gcw, gpcw, sigmad, L_i, iexp)
            Call myCOULFG(L_i*1.d0, eta, rho, fcw(L_i), fpcw(L_i), gcw(L_i), gpcw(L_i), iexp(L_i))
            cwup1 = (gcw(L_i) + ai * fcw(L_i))
            cwdown1 = gcw(L_i) - ai * fcw(L_i)

            bb0(ii, io) = (cwup0 * psi(ii) - cwup1 * psi0(ii)) &
                          / (cwup0 * cwdown1 - cwup1 * cwdown0)
            bb20(ii, io) = (cwdown1 * psi0(ii) - cwdown0 * psi(ii)) &
                          / (cwup0 * cwdown1 - cwup1 * cwdown0)
        End If
    End Do

  End Do
  !===============================================================
  !                                        penetration probability

  Do i = 1, Nlevel
      Do j = 1, Nlevel
        bb(i,j) = 0.d0
        bb2(i,j) = 0.d0
        Do j2 = 1, nlevel
            bb(i,j) = bb(i,j) + bb0(i,j2) * aa(j2,j)
            bb2(i,j) = bb2(i,j) + bb20(i,j2) * aa(j2,j)
        End Do
      End Do
  End Do

  Call matinv(Nlevel,bb,bin)

  P = 0.d0

  Do io = 1, Nlevel
      If (ech(io) .lt. 0.d0) Cycle
      k = sqrt((2.d0 * ReduceMass / Hbar**2 * ech(io)))
      kk = sqrt(2.d0 * ReduceMass / Hbar**2 * e)
      p = p + (abs(bin(io,1)))**2 * k / kk
  End Do


  ! Reconstruct the physical entrance-channel scattering solution.  The
  ! fundamental matrices were stored in the original inner-boundary basis;
  ! bin(:,1) normalizes them to one incoming wave in the entrance channel.
  If (capture_snapshot) Then
    If (Allocated(wf)) Deallocate(wf)
    If (Allocated(Numerov_radius)) Deallocate(Numerov_radius)
    Allocate(wf(Nx2),Numerov_radius(0:R_iterat+1))
    wf = (0.d0,0.d0)
    Do ir=0,R_iterat+1
      Numerov_radius(ir)=R_min+dr*ir
      Do ich=1,Nlevel
        idx=ir*Nlevel+ich
        Do io=1,Nlevel
          wf(idx)=wf(idx)+fundamental(ich,io,ir)*bin(io,1)
        End Do
      End Do
    End Do

    ! At the IWBC radius, an open-channel basis wave has unit modulus.  This
    ! gives a direct, independent check that the exported snapshot retains P.
    Numerov_snapshot_penetrability=0.d0
    kk=sqrt(2.d0*ReduceMass/Hbar**2*E)
    Do ich=1,Nlevel
      If (ech(ich) > 0.d0) Then
        k=sqrt(2.d0*ReduceMass/Hbar**2*ech(ich))
        Numerov_snapshot_penetrability=Numerov_snapshot_penetrability+ &
          abs(wf(ich))**2*k/kk
      End If
    End Do
    Numerov_snapshot_consistency=abs(Numerov_snapshot_penetrability-P)
    Numerov_snapshot_ready=.True.
    Deallocate(fundamental)
  End If

  Return

End Subroutine Numerov_core


Subroutine CaptureNumerovFundamental(ir,xi_modified,aa,fundamental)
  Use ccfull_initialization_mod
  Implicit None
  Integer, Intent(In) :: ir
  Complex(8), Intent(In) :: xi_modified(Nlevelmax,Nlevelmax)
  Complex(8), Intent(In) :: aa(Nlevelmax,Nlevelmax)
  Complex(8), Intent(Inout) :: fundamental(Nlevel,Nlevel,0:R_iterat+1)
  Complex(8) :: cc(Nlevelmax,Nlevelmax),cin(Nlevelmax,Nlevelmax)
  Complex(8) :: transformed(Nlevelmax,Nlevelmax)
  Integer :: i,j,k
  Real(8) :: fac,r,V
  External V

  fac=dr**2*(2.d0*ReduceMass/Hbar**2)
  r=R_min+dr*ir
  cc=(0.d0,0.d0)
  transformed=(0.d0,0.d0)
  Do i=1,Nlevel
    Do j=1,Nlevel
      cc(i,j)=-fac/12.d0*CPOT(i,j,ir)
      If (i == j) cc(i,j)=cc(i,j)-fac/12.d0*(V(r,L_i)-E)+1.d0
      Do k=1,Nlevel
        transformed(i,j)=transformed(i,j)+xi_modified(i,k)*aa(k,j)
      End Do
    End Do
  End Do
  Call matinv(Nlevel,cc,cin)
  Do i=1,Nlevel
    Do j=1,Nlevel
      fundamental(i,j,ir)=(0.d0,0.d0)
      Do k=1,Nlevel
        fundamental(i,j,ir)=fundamental(i,j,ir)+cin(i,k)*transformed(k,j)
      End Do
    End Do
  End Do
End Subroutine CaptureNumerovFundamental


Subroutine Numerov(P)
  Implicit None
  Real(8), Intent(Out) :: P
  Call Numerov_core(P,.False.)
End Subroutine Numerov


Subroutine NumerovSnapshot(P)
  Implicit None
  Real(8), Intent(Out) :: P
  Call Numerov_core(P,.True.)
End Subroutine NumerovSnapshot


Subroutine WriteNumerovSnapshot(filename)
  Use ccfull_initialization_mod
  Implicit None
  Character(Len=*), Intent(In) :: filename
  Integer :: unit_number,ir,ich,idx
  If (.Not. Numerov_snapshot_ready) Error Stop 'Numerov snapshot is not ready'
  Open(Newunit=unit_number,File=filename,Status='Replace',Action='Write')
  Write(unit_number,'(A)') '# ccfull_numerov_snapshot_v1'
  Write(unit_number,'(A,I0)') '# nlevel=',Nlevel
  Write(unit_number,'(A,I0)') '# radial_points=',R_iterat+2
  Write(unit_number,'(A,I0)') '# nx_cc=',Nx_cc
  Write(unit_number,'(A,I0)') '# nx2=',Nx2
  Write(unit_number,'(A,ES24.16)') '# energy_MeV=',Real(E,8)
  Write(unit_number,'(A,I0)') '# angular_momentum=',L_i
  Write(unit_number,'(A,ES24.16)') '# penetrability=',Numerov_snapshot_penetrability
  Write(unit_number,'(A,ES24.16)') '# consistency_abs=',Numerov_snapshot_consistency
  Write(unit_number,'(A)') 'radial_index,r_fm,channel,real,imag'
  Do ir=0,R_iterat+1
    Do ich=1,Nlevel
      idx=ir*Nlevel+ich
      Write(unit_number,'(I0,",",ES24.16,",",I0,",",ES24.16,",",ES24.16)') &
        ir,Numerov_radius(ir),ich,Real(wf(idx)),Aimag(wf(idx))
    End Do
  End Do
  Close(unit_number)
End Subroutine WriteNumerovSnapshot

Subroutine myCOULFG(angL,eta,rho,myf,mydf,myg,mydg,iv)
  double precision  angL,eta, rho, myf, mydf,myg,mydg
  integer iv,lmax 
  double precision XLMIN, XLMAX
  double precision XX, ETA1
  integer  MODE1, KFN, IFAIL
  double precision , allocatable, dimension(:):: FC,GC,FCP,GCP
  intent(in) angL
  intent(in) eta
  intent(in) rho
  intent(out) myf
  intent(out) mydf
  intent(out) myg
  intent(out) mydg
  intent(out) iv     
  XLMIN=angL
  XLMAX=angL
  lmax=IDINT(XLMAX)+1
  allocate (FC(lmax),GC(lmax),FCP(lmax),GCP(lmax))
  MODE1=1
  KFN=0 
  XX=rho
  ETA1=eta
  call COULFG(XX,ETA1,XLMIN,XLMAX,FC,GC,FCP,GCP,MODE1,KFN,IFAIL)
  myf=FC(lmax)
  mydf=FCP(lmax)
  myg=GC(lmax)
  mydg=GCP(lmax) 
  iv=IFAIL
  return 
End Subroutine

SUBROUTINE COULFG(XX, ETA1, XLMIN, XLMAX, FC, GC, FCP, GCP, MODE1, KFN, IFAIL)
    IMPLICIT NONE
    ! Argument type declarations
    DOUBLE PRECISION, INTENT(IN) :: XX, ETA1, XLMIN, XLMAX
    DOUBLE PRECISION, DIMENSION(*), INTENT(OUT) :: FC, GC, FCP, GCP
    INTEGER, INTENT(IN) :: MODE1, KFN
    INTEGER, INTENT(OUT) :: IFAIL

    ! Local variable declarations
    DOUBLE PRECISION :: ACCUR, ACC, ACC4, ACCH, ETA, GJWKB, PACCQ, X, XLM, E2MM1, DELL, XLL, XI, FCL, PK, PX, F, D, C
    DOUBLE PRECISION :: PK1, EK, RK2, TK, DF, FPL, XL, RL, EL, SL, FCL1, TA, WI, P, Q, AR, AI, BR, BI, DR, DI, DP, DQ
    DOUBLE PRECISION :: A, B, GAM, W, ALPHA, BETA, FCM, GCL, GPL, FJWKB, GCL1
    INTEGER :: MODE, NFP, NPQ, IEXP, M1, LXTRA, L1, LP, MAXL, L, LP1
    LOGICAL :: ETANE0, XLTURN

    ! Common block
    COMMON /STEE/ PACCQ, NFP, NPQ, IEXP, M1

    ! Constants
    DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, TEN2 = 1.0D2, ABORT = 2.0D4
    DOUBLE PRECISION, PARAMETER :: HALF = 0.5D0, TM30 = 1.0D-30, BIG = 1.0D+100
    DOUBLE PRECISION, PARAMETER :: RT2DPI = 0.7978845608286535587989211986876373D0

    ! Initialize accuracy and other variables
    ACCUR = 1.0D-16
    MODE = 1
    IF (MODE1 == 2 .OR. MODE1 == 3) MODE = MODE1
    IFAIL = 0
    IEXP = 1
    NPQ = 0
    ETA = ETA1
    GJWKB = ZERO
    PACCQ = ONE
    IF (KFN /= 0) ETA = ZERO
    ETANE0 = (ETA /= ZERO)
    ACC = ACCUR * 10.0D0
    ACC4 = ACC * TEN2 * TEN2
    ACCH = SQRT(ACC)

    ! Check for small XX
    IF (XX <= ACCH) THEN
        IFAIL = -1
        RETURN
    END IF

    X = XX
    XLM = XLMIN
    IF (KFN == 2) XLM = XLM - HALF
    IF (XLM <= -ONE .OR. XLMAX < XLMIN) THEN
        IFAIL = -2
        RETURN
    END IF

    E2MM1 = ETA*ETA + XLM*XLM + XLM
    XLTURN = (X*(X - TWO*ETA) < XLM*XLM + XLM)
    DELL = XLMAX - XLMIN + ACC
    LXTRA = INT(DELL)
    XLL = XLM + DBLE(LXTRA)

    M1 = MAX(INT(XLMIN + ACC), 0) + 1
    L1 = M1 + LXTRA

    XI = ONE / X
    FCL = ONE
    PK = XLL + ONE
    PX = PK + ABORT
    F = ETA/PK + PK*XI
    IF (ABS(F) < TM30) F = TM30
    D = ZERO
    C = F

    DO
        PK1 = PK + ONE
        EK = ETA / PK
        RK2 = ONE + EK*EK
        TK = (PK + PK1)*(XI + EK/PK1)
        D = TK - RK2 * D
        C = TK - RK2 / C
        IF (ABS(C) < TM30) C = TM30
        IF (ABS(D) < TM30) D = TM30
        D = ONE / D
        DF = D * C
        F = F * DF
        IF (D < ZERO) FCL = -FCL
        PK = PK1
        IF (PK > PX) THEN
            IFAIL = -3
            RETURN
        END IF
        IF (ABS(DF-ONE) < ACC) EXIT
    END DO

    NFP = PK - XLL - 1
    IF (LXTRA /= 0) THEN
        FCL = FCL / BIG
        FPL = FCL * F
        IF (MODE == 1) FCP(L1) = FPL
        FC(L1) = FCL
        XL = XLL
        RL = ONE
        EL = ZERO
        DO LP = 1, LXTRA
            IF (ETANE0) EL = ETA / XL
            IF (ETANE0) RL = SQRT(ONE + EL*EL)
            SL = EL + XL*XI
            L = L1 - LP
            FCL1 = (FCL * SL + FPL) / RL
            FPL = FCL1 * SL - FCL * RL
            FCL = FCL1
            FC(L) = FCL
            IF (MODE == 1) FCP(L) = FPL
            IF (MODE /= 3 .AND. ETANE0) GC(L+1) = RL
            IF (ABS(FCL) > BIG) THEN
                DO LP1 = L, M1+LXTRA
                    IF (MODE == 1) FCP(LP1) = FCP(LP1) * 1.0D-20
                    FC(LP1) = FC(LP1) * 1.0D-20
                END DO
                FCL = FC(L)
                FPL = FPL * 1.0D-20
            END IF
            XL = XL - ONE
        END DO
        IF (FCL == ZERO) FCL = ACC
        F = FPL / FCL
    END IF

    IF (XLTURN) CALL JWKB(X, ETA, MAX(XLM, ZERO), FJWKB, GJWKB, IEXP)
    IF (IEXP > 1 .OR. GJWKB > ONE/(ACCH*TEN2)) THEN
        W = FJWKB
        GAM = GJWKB * W
        P = F
        Q = ONE
    ELSE
        XLTURN = .FALSE.
        TA = TWO * ABORT
        PK = ZERO
        WI = ETA + ETA
        P = ZERO
        Q = ONE - ETA*XI
        AR = -E2MM1
        AI = ETA
        BR = TWO*(X - ETA)
        BI = TWO
        DR = BR/(BR*BR + BI*BI)
        DI = -BI/(BR*BR + BI*BI)
        DP = -XI*(AR*DI + AI*DR)
        DQ = XI*(AR*DR - AI*DI)
        DO
            P = P + DP
            Q = Q + DQ
            PK = PK + TWO
            AR = AR + PK
            AI = AI + WI
            BI = BI + TWO
            D = AR*DR - AI*DI + BR
            DI = AI*DR + AR*DI + BI
            C = ONE/(D*D + DI*DI)
            DR = C * D
            DI = -C * DI
            A = BR*DR - BI*DI - ONE
            B = BI*DR + BR*DI
            C = DP*A - DQ*B
            DQ = DP*B + DQ*A
            DP = C
            IF (PK > TA) THEN
                IFAIL = -4
                RETURN
            END IF
            IF (ABS(DP) + ABS(DQ) < (ABS(P) + ABS(Q)) * ACC) EXIT
        END DO
        NPQ = INT(PK / TWO)
        PACCQ = HALF * ACC / MIN(ABS(Q), ONE)
        IF (ABS(P) > ABS(Q)) PACCQ = PACCQ * ABS(P)
        GAM = (F - P) / Q
        IF (Q <= ACC4 * ABS(P)) THEN
            IFAIL = -5
            RETURN
        END IF
        W = ONE / SQRT((F - P)*GAM + Q)
    END IF

    ALPHA = ZERO
    IF (KFN == 1) ALPHA = XI
    IF (KFN == 2) ALPHA = XI * HALF
    BETA = ONE
    IF (KFN == 1) BETA = XI
    IF (KFN == 2) BETA = SQRT(XI) * RT2DPI
    FCM = SIGN(W, FCL) * BETA
    FC(M1) = FCM
    IF (MODE /= 3) THEN
        IF (.NOT. XLTURN) THEN
            GCL = FCM * GAM
        ELSE
            GCL = GJWKB * BETA
        END IF
        IF (KFN /= 0) GCL = -GCL
        GC(M1) = GCL
        GPL = GCL*(P - Q/GAM) - ALPHA*GCL
        IF (MODE /= 2) THEN
            GCP(M1) = GPL
            FCP(M1) = FCM*(F - ALPHA)
        END IF
    END IF

    IF (LXTRA == 0) RETURN

    W = BETA * W / ABS(FCL)
    MAXL = L1 - 1
    DO L = M1, MAXL
        IF (MODE /= 3) THEN
            XL = XL + ONE
            IF (ETANE0) EL = ETA / XL
            IF (ETANE0) RL = GC(L+1)
            SL = EL + XL*XI
            GCL1 = ((SL - ALPHA)*GCL - GPL) / RL
            IF (ABS(GCL1) > BIG) THEN
                IFAIL = L1 - L
                RETURN
            END IF
            GPL = RL*GCL - (SL + ALPHA)*GCL1
            GCL = GCL1
            GC(L+1) = GCL1
            IF (MODE /= 2) THEN
                GCP(L+1) = GPL
                FCP(L+1) = W*(FCP(L+1) - ALPHA*FC(L+1))
            END IF
        END IF
        FC(L+1) = W * FC(L+1)
    END DO

END SUBROUTINE COULFG

subroutine JWKB(XX, ETA1, XL, FJWKB, GJWKB, IEXP)
    implicit none
    ! ===== 输入输出参数 =====
    real(8), intent(in)  :: XX, ETA1, XL
    real(8), intent(out) :: FJWKB, GJWKB
    integer, intent(out) :: IEXP

    ! ===== 内部变量 =====
    real(8) :: ZERO, HALF, ONE, SIX, TEN
    real(8) :: DZ, RL35
    real(8) :: aloge, X, ETA, GH2, XLL1, HLL, HL, SL, RL2, GH
    real(8) :: PHI, PHI10

    ! ===== 常量赋值 =====
    ZERO = 0.0d0
    HALF = 0.5d0
    ONE  = 1.0d0
    SIX  = 6.0d0
    TEN  = 10.0d0
    DZ   = 0.0d0
    RL35 = 35.0d0

    aloge = log(TEN)

    ! ===== 核心计算 =====
    X   = XX
    ETA = ETA1

    GH2  = X * (2.0d0*ETA - X)
    XLL1 = max(XL*XL + XL, DZ)

    if (GH2 + XLL1 <= ZERO) then
        FJWKB = 0.0d0
        GJWKB = 0.0d0
        IEXP  = 0
        return
    end if

    HLL = XLL1 + SIX / RL35
    HL  = sqrt(HLL)
    SL  = ETA/HL + HL/X
    RL2 = ONE + ETA*ETA/HLL
    GH  = sqrt(GH2 + HLL) / X

    PHI = X*GH - HALF * ( HL*log((GH+SL)**2 / RL2) - log(GH) )
    if (ETA /= ZERO) PHI = PHI - ETA*atan2(X*GH, X - ETA)

    PHI10 = -PHI * aloge
    IEXP  = int(PHI10)

    if (IEXP > 70) then
        GJWKB = TEN**(PHI10 - dble(IEXP))
    else
        GJWKB = exp(-PHI)
        IEXP  = 0
    end if

    FJWKB = HALF / (GH * GJWKB)

end subroutine JWKB

Subroutine WHIT(HETA, R, XK, E, LL, F, FD, IE)
    implicit none
    ! ====== 输入输出参数 ======
    integer, intent(in) :: LL
    integer, intent(inout) :: IE
    real(8), intent(in) :: HETA, R, XK, E
    real(8), intent(out) :: F(LL+1), FD(LL+1)

    ! ====== 内部变量 ======
    integer :: L, LP1, LM, LMP1, IS, H
    integer :: I, M, M1, M2, JS, IFEQL
    real(8) :: fpmax, fpminl
    real(8) :: EE, AK, ETA, RHO, RHOA, PJE, A, B, C, D
    real(8) :: T(12), S(7)

    fpmax = 1.0d290
    L = LL + 1

    EE = -1.0d0
    AK = XK
    ETA = HETA
    LP1 = L + 1
    RHO = AK * R

    S(:) = 0.0d0

    ! ====== 设置 LM ======
    if (L <= 50) then
        LM = 60
    else
        LM = L + 10
    end if
    LMP1 = LM + 1
    IS = 7

    PJE = 30.0d0 * RHO + 1.0d0
    H = max(int(PJE), 4)
    H = RHO / H

    RHOA = 10.0d0 * (ETA + 1.0d0)
    if (RHOA <= RHO) then
        IFEQL = 1
        RHOA = RHO
    else
        IFEQL = 0
    end if

    PJE = RHOA / H + 0.5d0
    RHOA = H * int(PJE)

    if (IFEQL == 1) then
        if (RHOA <= RHO + 1.5d0 * H) then
            RHOA = RHO + 2.0d0 * H
        end if
    end if

    if (EE < 0.0d0) then
        ! ----------- 级数展开部分 ----------
        C = 1.0d0 / RHOA
        A = 1.0d0
        B = 1.0d0 - C * ETA
        F(1) = A
        FD(1) = B
        do M = 1, 26
            D = 0.5d0 * (ETA + dble(M - 1)) * (ETA + dble(M)) * C / dble(M)
            A = -A * D
            B = -B * D - A * C
            F(1) = F(1) + A
            FD(1) = FD(1) + B
        end do

        A = -ETA * log(2.0d0 * RHOA) - RHOA
        fpminl = -log(fpmax)
        if (IE == 0 .and. A < fpminl) IE = int(fpminl - A)
        A = exp(A + dble(IE))

        F(1) = A * F(1)
        FD(1) = A * FD(1) * (-1.0d0 - 2.0d0 * ETA / RHOA)

        if (IFEQL == 1) then
            S(IS) = F(1)
            if (IS == 7) then
                IS = 6
                RHOA = RHOA + H
                ! 再做一次展开
                C = 1.0d0 / RHOA
                A = 1.0d0
                B = 1.0d0 - C * ETA
                F(1) = A
                FD(1) = B
                do M = 1, 26
                    D = 0.5d0 * (ETA + dble(M - 1)) * (ETA + dble(M)) * C / dble(M)
                    A = -A * D
                    B = -B * D - A * C
                    F(1) = F(1) + A
                    FD(1) = FD(1) + B
                end do
                A = -ETA * log(2.0d0 * RHOA) - RHOA
                fpminl = -log(fpmax)
                if (IE == 0 .and. A < fpminl) IE = int(fpminl - A)
                A = exp(A + dble(IE))
                F(1) = A * F(1)
                FD(1) = A * FD(1) * (-1.0d0 - 2.0d0 * ETA / RHOA)
            end if
        else
            F(1) = T(1)
            FD(1) = T(2)
        end if

    else
        stop "WHIT1: EE >= 0 not implemented"
    end if

    ! ----------- 递推关系 ----------
    C = 1.0d0 / RHO
    do M = 1, L-1
        A = ETA / dble(M)
        B = A + C * dble(M)
        F(M+1) = (B * F(M) - FD(M)) / (A + 1.0d0)
        FD(M+1) = (A - 1.0d0) * F(M) - B * F(M+1)
    end do

    do M = 1, L
        FD(M) = AK * FD(M)
    end do

    return
End Subroutine

Program CCFull_Main
  Use ccfull_initialization_mod
  Implicit None
  Integer ::  i, j, ir, ie, il
  Integer ::  Iestep
  Real(8) ::  r, rh, sigma, spin, s0, P_0, P, V
  external V
  Call Initialize_CCFull()
  Write(*,*)'Nlevle = ', Nlevel
  Call coupled_matrix0

  Do ir = 0, R_iterat + 1
    r = R_min + dr * ir
    Call coupled_matrix(r, CPOT0)
    Do  i = 1, Nlevel
      Do  j = 1, Nlevel
        CPOT(i, j, ir) = CPOT0(i, j)
      End Do
    End Do
  End Do
  

  rh = R_min + dr / 2.0
  Call coupled_matrix(rh, CPOTH)


  Iestep = int((Emax - Emin)/dE)
  Do ie = 0, Iestep

    E = Emin + dE * ie
    sigma = 0.
    spin  = 0.
    P_0 = 0.

    Do il = 0, 200
      L_i = il
      s0 = sigma

      Call PotShape(L_i, R_barrier_l, V_barrier_l, curv_l, R_bottom_l, V_bottom_l)
      If( V(R_bottom_l, L_i) > E .or. R_bottom_l < 0 )Then
          P = 0.
          exit  
      End If

      ! The standalone physical model retains the original observable path.
      ! Reduced-basis training calls NumerovSnapshot directly from its driver.
      Call Numerov(P)
      
      sigma = sigma + (2.0D0 * L_i + 1.0D0) * P * pi * Hbar**2 / (2.0D0 * ReduceMass * E) * 10.0D0
      spin  = spin  + (2.0D0 * L_i + 1.0D0) * P * pi * Hbar**2 / (2.0D0 * ReduceMass * E) * 10.0D0 * L_i

      If (L_i == 0) Then
          P_0 = P
      End If

      If(sigma - s0 < s0*1.d-4)Exit 
      !Write(*,*)E, L_i, P

    End Do
    !stop

    if (sigma < 0.01D0) then
        write(*,'(A,F8.3,A,E15.5,A,F8.5,A,E15.5)') ' E = ', E, ' MeV,   Sigma = ', sigma, &
                                        ' mb, <L> = ', spin / sigma, ' hbar, <P_0> = ', P_0
        write(20,'(F15.5,E15.5,F15.5,E15.5)') E,  sigma,  spin / sigma,  P_0
    else
        write(*,'(A,F8.3,A,F15.5,A,F8.5,A,E15.5)') ' E = ', E, ' MeV,   Sigma = ', sigma,&
                                       ' mb, <L> = ', spin / sigma, ' hbar, <P_0> = ', P_0
        write(20,'(4F15.5)') E,  sigma,  spin / sigma,  P_0
    end if

  End Do

  Print *, "Program completed successfully."

End Program CCFull_Main
