!***********************************************************************
!
!    Copyright (c) 2016, Lawrence Livermore National Security, LLC.
!                        Produced at the Lawrence Livermore National
!                        Laboratory.
!                        Written by Nicolas Schunck, schunck1@llnl.gov
!
!    LLNL-CODE-728299 All rights reserved.
!    LLNL-CODE-573953 All rights reserved.
!
!    Copyright 2017, R. Navarro Perez, N. Schunck, R. Lasseri, C. Zhang,
!                    J. Sarich
!    Copyright 2012, M.V. Stoitsov, N. Schunck, M. Kortelainen, H.A. Nam,
!                    N. Michel, J. Sarich, S. Wild
!    Copyright 2005, M.V. Stoitsov, J. Dobaczewski, W. Nazarewicz, P.Ring
!
!    This file is part of HFBTHO.
!
!    HFBTHO is free software: you can redistribute it and/or modify it
!    under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    HFBTHO is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with HFBTHO. If not, see <http://www.gnu.org/licenses/>.
!
!    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
!    LICENSE
!
!    Our Preamble Notice
!
!      A. This notice is required to be provided under our contract
!         with the U.S. Department of Energy (DOE). This work was
!         produced at the Lawrence Livermore National Laboratory under
!         Contract No. DE-AC52-07NA27344 with the DOE.
!      B. Neither the United States Government nor Lawrence Livermore
!         National Security, LLC nor any of their employees, makes any
!         warranty, express or implied, or assumes any liability or
!         responsibility for the accuracy, completeness, or usefulness
!         of any information, apparatus, product, or process disclosed,
!         or represents that its use would not infringe privately-owned
!         rights.
!      C. Also, reference herein to any specific commercial products,
!         process, or services by trade name, trademark, manufacturer
!         or otherwise does not necessarily constitute or imply its
!         endorsement, recommendation, or favoring by the United States
!         Government or Lawrence Livermore National Security, LLC. The
!         views and opinions of authors expressed herein do not
!         necessarily state or reflect those of the United States
!         Government or Lawrence Livermore National Security, LLC, and
!         shall not be used for advertising or product endorsement
!         purposes.
!
!    The precise terms and conditions for copying, distribution and
!    modification are contained in the file COPYING.
!
!***********************************************************************

! ==================================================================== !
!                                                                      !
!                        DFT SOLVER PACKAGE                            !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------
!> This module provides the main \theCode DFT solver. It includes
!> routines for the calculation and diagonalization of the HFB matrix;
!> the definition of densities on the quadrature mesh and in configuration
!> space; the self-consistent loop the calculation of expectation values of
!> observables. It also includes a package of a few routines to perform
!> particle number projection.
!>
!>  @author
!>    Mario Stoitsov, Nicolas Schunck, Markus kortelainen, Rodrigo Navarro Perez
!----------------------------------------------------------------------
!  Subroutines: - HFBTHO_DFT_SOLVER
!               - heading
!               - thodefh
!               - preparer
!               - coordinateLST
!               - iter
!               - hfbdiag
!               - ALambda
!               - Canonical
!               - resu
!               - initialize_HFBTHO_NAMELIST
!               - read_HFBTHO_NAMELIST
!               - check_consistency
!               - initialize_HFBTHO_SOLVER
!               - base0
!               - base
!               - gaupol
!               - start
!               - printRHO
!               - gfv
!               - sdiag
!               - nucleus
!               - stab
!               - ord
!               - tracesln
!               - tracesln_qp
!               - densitln
!               - coulom1
!               - coulom
!               - coulom_test
!               - HartreeDir
!               - optHFBTHO
!               - DENSIT
!               - field
!               - gamdel
!               - gamdel_gogny
!               - expectpj
!               - densitpj
!               - coulompj
!               - broyden_min
!               - expect
!               - Constraint_or_not
!               - getLagrange
!               - requested_blocked_level
!----------------------------------------------------------------------
Module HFBTHO_solver

  Use HFBTHO_utilities
  Use HFBTHO
  Use HFBTHO_multipole_moments
  Use HFBTHO_collective
  Use HFBTHO_fission_fragments
  Use HFBTHO_gogny
  Use HFBTHO_io

  Implicit None

Contains

  !-----------------------------------------------------------------------------------------------
  !
  !   Axially-deformed configurational constrained and/or unconstrained Hartree-Fock-Bogoliubov
  !   calculations with Skyrme-like functionals and delta pairing using the Harmonic-Oscillator
  !   (HO), and/or Transformed HO (THO) basis with or without reflection symmetry imposed, with
  !   or without the Lipkin-Nogami procedure. The solver can handle all Skyrme-like functionals,
  !   DME-functionals, Fayans-functionals, calculate infinite nuclear matter properties, finite
  !   nuclei (even-even, odd-even, odd-odd), and neutron drops, isoscalar and isovector monopo-
  !   le FAM QRPA calculations for spherical and deformed nuclei.
  !
  !   All necessary input variables contain the suffix _INI. Below, the complete list of these
  !   variables with some example values:
  !
  !   ======== hfbtho_NAMELIST.dat
  !   n00_INI=20; npr1_INI=70; npr2_INI=50;  kindhfb_INI=-1; inin_INI=-1
  !   b0_INI=2.234776; q_INI=0.0; cdef_INI=0.0; cqad_INI=0.5; skyrme_INI='SLY4'; nkblo_INI=0
  !   ILST_INI=0; keypj_INI=1; iproj_INI=0; npr1pj_INI=0;
  !   icou_INI=2; IDEBUG_INI=0; npr2pj_INI=0;
  !   Parity_INI=.False.; epsi_INI=0.00001_pr; MAX_ITER_INI=101
  !   Add_Pairing_INI=.False.; DO_FITT_INI=.False.; Print_PTHO_Namelist_INI=.True.
  !
  !   ======== from read_UNEDF_NAMELIST
  !   DMEORDER=-1; DMELDA=0; use_TMR_pairing=0
  !   HBZERO=20.73553000000000;    E2CHARG=1.4399784085965135; CRHO(0)=-933.3423749999999;  CRHO(1)=830.0524855000001;
  !   CDRHO(0)=861.0625000000000;  CDRHO(1)=-1064.2732500000;  CTAU(0)=57.12868750000000;   CTAU(1)=24.65673650000000;
  !   CRDR(0)=-76.99620312499999;  CRDR(1)=15.65713512500000;  CRDJ(0)=-92.25000000000000;  CRDJ(1)=-30.7500000000000;
  !   CJ(0)=17.20961150000000;     CJ(1)=64.57581250000000;    CPV0(0)=-258.2000000000000;  CPV0(1)=-258.2000000000000;
  !   CPV1(0)=0.5000000000000000;  CPV1(1)=0.500000000000000;  SIGMA=0.1666666666666667;    CEXPAR=1.000000000000000;
  !   E_NM=-15.97214914144462;     K_NM=229.9009644826037;     SMASS_NM =1.439546988976078; RHO_NM =0.1595387567117334;
  !   ASS_NM =32.00430281505202;   LASS_NM=45.96175148046161;  VMASS_NM =1.249838547196253;
  !   MPI=0.6995945261023822;      GA=1.290000000000000;       FPI=0.4683223517486062;      C1=-0.1598130000000000;
  !   C3 =-0.6708200000000;        C4 =0.6708200000000000;     CD =-2.062000000000000;      CE=-0.6250000000000;
  !   LAMBDAX =3.547896604156107;  USE_INM=.false.;            USE_CM_COR =.true.;          USE_DME3N_TERMS=.true.;
  !   USE_J2TERMS =.true.;         USE_CHARGE_DENSITY=.false.; PRINT_NAMELIST=.true.;
  !
  ! Memo:
  !  -  inin_INI switches scratch unconstrained (inin=1,2,3) or constrained (inin=100,200,300)
  !     calculations. Unconstrained mode begins with a small number of constrained iterations.
  !  -  inin_INI switches unconstrained (inin=-1,-2,-3) or constrained (inin=-100,-200,-300)
  !     calculations from a previous solution if the latter exists. If not, the solver sets
  !     inin=Abs(inin) and resumes from scratch.
  !  -  The same holds for odd nuclei. If even-even solution for the odd nucleus does not exists
  !     it is calculated first.
  !  -  Print_Screen=T/F for n00_INI=+/-: output is generated and written in thoout.dat file
  !	only if n00_INI>0. For n00_INI<0, the number of shells is set to abs(n00_INI) but all
  !     output is supressed.
  !
  !  -  At the end of the solution, the solver provides all results in the arrays
  !
  !                  nucname,ereslbl(1:2),eres(1:ierest)
  !
  !     which contain:
  !
  !        LBL,BLKN,BLKZ,Jsi,JININ,A,N',Z,Efn,Efp,JEtot,Jbett,Jbetn,Jbetp,JQt,JQn,JQp
  !        JpEn,JpEp,JpDn,JpDp,JAsn,JAsp,Jrt,Jrn,Jrp,Jrc,Jht,Jhn,Jhp,Jqht,Jqhn,Jqhp,
  !        JKINt,JKINn,JKINp,JSO,JCDIR,JCEX,JDisn,JDisp,JV2Mn,JV2Mp,JILST,JKIND,JL,
  !        JECMPAV1,JECMPAV2,JECMPAV3,JA,JN,JZ,ITER,UEtot,Ubett,Ubetn,Ubetp,UQt,UQn,
  !        UQp,Uln,Ulp,UpEn,UpEp,UpDn,UpDp,UAsn,UAsp,Urt,Urn,Urp,Urc,Uht,Uhn,Uhp,
  !        Uqht,Uqhn,Uqhp,UKINT,UKINN,UKINP,USO,UCDIR,UCEX,UDisn,UDisp,UV2Mn,UV2Mp,
  !        UECMT,UECMN,UECMP,UROTT,UROTN,UROTP,USQUJT,USQUJN,USQUJP,UCRANT,UCRANN,
  !        UCRANP,UERIGT,UERIGN,UERIGP,EHFBLN,EHFB,LNbet,LNben,LNbep,LNQt,LNQn,LNQp,
  !        LNpEn,LNpEp,LNpDn,LNpDp,LNrt,LNrn,LNrp,LNrC,LNam2n,LNam2p,LNe2n,LNe2p,
  !        BlEqpN,BlDEqpN,BlOvrN,BlEqpZ,BlDEqpZ,BlOvrZ
  !
  !  -  Standard inputs from               'hfbtho_NAMELIST.dat'
  !  -  USer-defined functional read from  'UNEDF_NAMELIST.DAT'
  !  -  Final solution stored to files     '*.hel' and/or '*.tel'
  !  -  Output is written to files         'thoout.dat', 'thodef.dat' and 'hodef.dat'
  !  -  Output files *.dat may exist as    'thoout','thores','hodenp',
  !                                      'thodenp','thodef','thoene',
  !                                      'thoprc','dat0.1.2.3.4'
  !  -  External accuracy pr/ipr always set in UNEDF module
  !
  !  -  n00    Number of oscillator shells
  !             n00>0 prints to thoout.dat & screen
  !             n00<0 no print at all
  !             n00=0 program stops (NB!)
  !  -  b0     Oscillator Basis parameter b0>0 (If b0<0 it takes a default value)
  !  -  beta0  Value of Basis deformation parameter
  !  -  AN     Number of neutrons N
  !  -  AZ     Number of protons Z
  !  -  FTST'  Fayance forces label
  !  -  kind   Kind of calculations 1: noLN, -1:LN
  !  -  inin   Unconstraint Iterations from peviouse solution
  !             -1: (from spherical *.hel or *.tel file)
  !             -2: (from prolate   *.hel or *.tel file)
  !             -3: (from oblate    *.hel or *.tel file)
  !            Unconstraint Iterations fom scratch with
  !            preliminary constraint at deformation Cbeta, (i):
  !              1: Spherical scratch
  !              2: Prolate scratch
  !              3: Oblate scratch
  !            Constrained calculation (icstr) at Cbeta, see (i)
  !              100, 200, 300 fom scratch
  !             -100,-200,-300 fom previouse solution
  !  -  blNeutrons: a group responsible for blocking a particular neutron level
  !     The group consists of 5 numbers, e.g., for 7-[ 3, 0, 3]: 7 -1 3 0 3
  !      k1  2 \times \Omega
  !         =0: the whole group (k) is disregarded (n0 blocking)
  !         >0: blocking in N+1 nucleus
  !         <0: blocking in N-1 nucleus
  !      k2  parity (+1 or -1); NB! when k2=0, the ground state walker is applied
  !      k3,k4,k5  Nilson quantum numbers
  !  -  blProtons: exactly the same as (k) but for protons
  !
  !-----------------------------------------------------------------------------------------------
  Subroutine HFBTHO_DFT_SOLVER
    Use HFBTHO_THO
    Implicit None
    Integer(ipr) :: iw,ib,j,i,it,l,maxi0,icstr0,iterMax,icons,il,kickoff,iexit,ncons_eff
    Real(pr)     :: epsi0,qq,f,f1,f2,f3,r,g,g1
    !-------------------------------------------------------------
    ! Initializing all according to *_INI values
    !-------------------------------------------------------------
    Call initialize_HFBTHO_SOLVER
    If(ierror_flag.Ne.0) Return
#if(DO_MASSTABLE==1 || DRIP_LINES==1 || DO_PES==1)
    If(lout.le.lfile) Open(lfile,file='thoout'//row_string//'.dat',status='unknown')
#else
#ifdef SUPRESS_OUTPUT
    If(lout.le.lfile) Open(lfile,file='thoout.dat',status='unknown')
#else
    If(lout.lt.lfile) Open(lfile,file='thoout.dat',status='unknown')
#endif
#endif
    Call Constraint_or_not(inin_INI,inin,icstr)
    If(ierror_flag.Ne.0) Return
    !-------------------------------------------------------------------------
    ! Loop recalculating eventually the even-even solution for an odd nucleus
    !-------------------------------------------------------------------------
    irestart=0; iexit=0
    matrix_elements_calculated = .false.
    Do
       n00=Abs(n00_INI);  b0=b0_INI;           q=q_INI;           iLST=iLST_INI;
       maxi=MAX_ITER_INI; npr(1)=npr_INI(1);   npr(2)=npr_INI(2); npr(3)=npr(1)+npr(2);
       skyrme=skyrme_INI; kindhfb=kindhfb_INI
       keypj=keypj_INI;   iproj=iproj_INI;     npr1pj=npr1pj_INI; npr2pj=npr2pj_INI;
       nkblo=nkblo_INI;
       basis_HFODD=basis_HFODD_INI
       !-------------------------------------------------------------
       ! Define the set of constraints
       !-------------------------------------------------------------
       numberCons=0; kickoff=0
       Do l=1,lambdaMax
          If(lambda_active(l).Gt.0) numberCons = numberCons + 1
          If(lambda_active(l).Lt.0) kickoff = kickoff + 1
       End Do
       ! Add constraint on the neck
       If(neck_constraints) Then
          numberCons = numberCons + 1
          neckLag = zero
          !kickoff = kickoff + 1
       End If
       !
       ncons_eff = numberCons + kickoff
       If(.Not.Allocated(multLag)) Allocate(multLag(1:lambdaMax)); multLag=zero
       If(.Not.Allocated(multLambda)) Allocate(multLambda(1:ncons_eff)); multLambda=0
       If(.Not.Allocated(multRequested)) Allocate(multRequested(0:lambdaMax)); multRequested=zero
       !
       icons=0
       Do l=1,lambdaMax
          If(lambda_active(l).Gt.0) Then
             icons=icons+1
             multLambda(icons)=lambda_values(l)
          End If
          multRequested(l) = expectation_values(l)
       End Do
       If(neck_constraints) Then
          icons=icons+1
          multLambda(icons)=0
       End If
       !-------------------------------------------------------------
       ! Blocking
       !-------------------------------------------------------------
       Do it=1,2
          If(nkblo(it,1).Ne.0) Then
             If(nkblo(it,1).Gt.0) Then
               ! particle state
               npr(it)=npr(it)+1
               iparenti(it)=-1
             Else
               ! hole state
               npr(it)=npr(it)-1
               iparenti(it)=+1
             End If
             nkblo(it,1)=Abs(nkblo(it,1))
             If(nkblo(it,2).Eq.0) Then
               ! ground state walker
               keyblo(it)=keyblo(it)+1 !nkblo(it,1)
               If(keyblo(it).Eq.blomax(it)) irestart=0
             Else
               irestart=0
             End If
          End If
       End Do
       !-------------------------------------------------------------
       ! HFB+HO calculations
       !-------------------------------------------------------------
       If(ILST.Le.0) Then
          icacou=0; icahartree=0
          Call preparer(.True.)
          If(ierror_flag.Ne.0) Return
          ! Reading input file: if not possible, start from scratch
          Call inout(1,iexit)
          If(iexit>0) Then
             Call start
          Else
             Call gamdel(.false.,.false.)
          End If
          If(ierror_flag.Ne.0) Return
          !-------------------------------------------------------------
          ! Preliminary constrained calculations
          !-------------------------------------------------------------
          If(kickoff.Gt.0) Then
             icstr0=icstr; epsi0=epsi; ! remember accuracy
             icstr=1                   ! constraint true
             epsi=1.0_pr               ! small accuracy
             if(is_nedf) then
                iterMax = maxi; maxi = 25
             else
                iterMax = maxi; maxi = 10
             endif
             numberCons=0
             Do l=1,lambdaMax
                If(Abs(lambda_active(l)).Gt.0) Then
                   numberCons=numberCons+1
                   multLambda(numberCons)=lambda_values(l)
                End If
             End Do
             Do iw=lout,lfile
                If(Parity) Then
                   Write(iw,'(/,a,i3,a,i2,a,/)') '  ### INITIAL STAGE(constrained calculations, reflection symmetry used)'
                Else
                   Write(iw,'(/,a,i3,a,i2,a,/)') '  ### INITIAL STAGE(constrained calculations, no reflection symmetry used)'
                End If
             End Do
             Call iter(.True.)     ! small constraint iterations
             If(ierror_flag.Ne.0) Return
             ! For the next phase, use only true constraints
             icstr=icstr0; epsi=epsi0
             maxi = iterMax
             numberCons=0
             Do l=1,lambdaMax
                If(lambda_active(l).Gt.0) Then
                   numberCons=numberCons+1
                   multLambda(numberCons)=lambda_values(l)
                End If
             End Do
          End If
          !-------------------------------------------------------------
          ! REGULAR HFB+HO ITERATIONS
          !-------------------------------------------------------------
          Do iw=lout,lfile
             If(Parity) Then
                Write(iw,'(/,a,i3,a,i2,a,/)')    '  ### REGULAR STAGE (reflection symmetry imposed)'
             Else
                Write(iw,'(/,a,i3,a,i2,a,/)')    '  ### REGULAR STAGE (no reflection symmetry imposed)'
             End If
          End Do
          Call iter(.True.)
          If(ierror_flag.Ne.0) Return
          Call resu(1)
          If(ierror_flag.Ne.0) Return
       End If
       !! write LST function on disk
       !Open(unit=66,file='LST.dat',status='unknown')
       !Write(66,'("#",15X,"R",20X,"f(R)",20X,"f^(1)",20X,"f^(2)",20X,"f^(3)")')
       !Do il=1,170
       !   qq=Real(il-1)/10.0_pr*1.0_pr
       !   If(il.Eq.1) Call thofun(0,g,f,f1,f2,f3,g1,.True.,.True.)
       !   Call thofun(1,qq,f,f1,f2,f3,r,.False.,.True.)
       !   Write(66,'(6E24.10)') r,qq,f1,f2,f3
       !End Do
       !Close(66)
       !-------------------------------------------------------------
       ! HFB+THO calculations from HFB+HO
       !-------------------------------------------------------------
       If(ILST.Lt.0) Then
          ILST1=1; icacou=0; icahartree=0
          Call coordinateLST(.False.) ! THO basis
          If(ierror_flag.Ne.0) Return
          Call densit                 ! THO densities
          If(ierror_flag.Ne.0) Return
          Call field                  ! Nuclear fields
          If(ierror_flag.Ne.0) Return
          Call iter(.True.)           ! HFB+THO iterations
          If(ierror_flag.Ne.0) Return
          Call resu(1)                ! print/record results
          If(ierror_flag.Ne.0) Return
       End If
       !-------------------------------------------------------------
       ! HFB+THO calculations from *.tel
       !-------------------------------------------------------------
       If(ILST.Gt.0) Then
          If(inin.Gt.0) Then
             ierror_flag=ierror_flag+1
             ierror_info(ierror_flag)= ' Forbidden iLST>0, inin>0 '
             Return
          End If
          icacou=0; icahartree=0
          Call preparer(.True.)
          If(ierror_flag.Ne.0) Return
          ! Reading input file: if not possible, start from scratch
          Call inout(1,iexit)
          If(iexit>0) Then
             Call start
          Else
             Call gamdel(.false.,.false.)
          End If
          If(ierror_flag.Ne.0) Return
          Call iter(.True.)           ! HFB+THO iterations
          If(ierror_flag.Ne.0) Return
          Call resu(1)                ! print/record results
          If(ierror_flag.Ne.0) Return
       End If
       !-------------------------------------------------------------
       ! Go for the requested blocking state in a case of odd nuclei
       ! if restarted due to corrupted/missing previous solution
       !-------------------------------------------------------------
       inin=-Abs(inin)
       If(irestart.Eq.0) Exit
    End Do
    !
  End Subroutine HFBTHO_DFT_SOLVER
  !=======================================================================
  !> Print heading to screen 'lout' and to tape thoout.dat 'lfile'
  !=======================================================================
  Subroutine heading
#if(USE_OPENMP==1)
    Use omp_lib
#endif
#if(USE_MPI==1)
    use mpi
#endif
    Implicit None
#if(USE_MPI==1)
    Integer(ipr) :: rank,ierr
#endif
    Integer(ipr) :: iw,idt(8),numThreads,idThread
    Character(len=12) rcl(3)
    Character(len=50) today
    !
#if(USE_OPENMP==1)
!$OMP PARALLEL SHARED(iw,lout,lfile) PRIVATE(numThreads,idThread)
    numThreads = omp_get_num_threads()
    idThread = omp_get_thread_num()
#if(DO_MASSTABLE==1 || DO_PES==1)
    If (idThread .Eq. 0.and.irow.le.1) Then
#elif(DRIP_LINES==1)
    If (idThread .Eq. 0.and.irow.le.1.and.calc_counter.eq.1) Then
#elif(USE_MPI==1)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    If (idThread .Eq. 0 .and. rank .eq. 0) Then
#else
    If (idThread .Eq. 0) Then
#endif
        Do iw=lout,lfile
           Write(iw,'("Multi-threading framework with OpenMP:",i2," threads/task")') numThreads
        End Do
    End If
!$OMP END PARALLEL
#endif

    Call Date_and_time(rcl(1),rcl(2),rcl(3),idt)
    Write(today,'(a,i2,a,i2,a,i4,a,i2,a,i2,a)')'(',idt(2),'/',idt(3),'/',idt(1),', ',idt(5),':',idt(6),')'
    Do iw=lout,lfile
       Write(iw,*)
       Write(iw,'("  ==========================================")')
       Write(iw,'("           FORTRAN 95 CODE (KIND=",i2,") ")') pr
       Write(iw,'("               Version: ",a)') Version
       Write(iw,'("  ==========================================")')
       Write(iw,'("       AXIALLY DEFORMED CONFIGURATIONAL     ")')
       Write(iw,'("     HARTREE-FOCK-BOGOLIUBOV CALCULATIONS   ")')
       Write(iw,'("                     WITH                   ")')
       Write(iw,'("      SKYRME+DELTA PAIRING OR GOGNY EDFs    ")')
       Write(iw,'("                  USING THE                 ")')
       Write(iw,'("             HARMONIC-OSCILLATOR            ")')
       Write(iw,'("                    AND/OR                  ")')
       Write(iw,'("        TRANSFORMED HARMONIC-OSCILLATOR     ")')
       Write(iw,'("                    BASIS                   ")')
       Write(iw,'("                     ---                    ")')
       Write(iw,'("                v1.66  (2005)               ")')
       Write(iw,'("    Stoitsov, Dobaczewski, Nazarewicz, Ring ")')
       Write(iw,'("                v2.00d (2012)               ")')
       Write(iw,'("        Stoitsov, Schunck, Kortelainen      ")')
       Write(iw,'("                v3.00  (2016)               ")')
       Write(iw,'("            Schunck, Navarro Perez          ")')
       Write(iw,'("  ==========================================")')
       Write(iw,'("    Nucleus: ",a," (A=",i4,", N=",i3,", Z=",i3,")")') nucname,npr(1)+npr(2),npr(1),npr(2)
       If(Parity) Then
          Write(iw,'("       Reflection Symmetry Imposed       ")')
       Else
          Write(iw,'("      No Reflection Symmetry Imposed     ")')
       End If
       Write(iw,'("            ",a)') today
       Write(iw,'("  ==========================================")')
       Write(iw,*)
       Write(iw,*)
    End Do
  End Subroutine heading
  !=======================================================================
  !> Print labels to hodef.dat or/and thodef.dat files
  !=======================================================================
  Subroutine thodefh(iw1)
    Implicit None
    Integer(ipr) :: iw1
    hlabels(1)='LBL';       hlabels(11)='JEtot';    hlabels(21)='JpDp';
    hlabels(2)='BLKN';      hlabels(12)='Jbett';    hlabels(22)='JAsn';
    hlabels(3)='BLKZ';      hlabels(13)='Jbetn';    hlabels(23)='JAsp';
    hlabels(4)='Jsi';       hlabels(14)='Jbetp';    hlabels(24)='Jrt';
    hlabels(5)='JININ';     hlabels(15)='JQt';      hlabels(25)='Jrn';
    hlabels(6)='A';         hlabels(16)='JQn';      hlabels(26)='Jrp';
    hlabels(7)='N';         hlabels(17)='JQp';      hlabels(27)='Jrc';
    hlabels(8)='Z';         hlabels(18)='JpEn';     hlabels(28)='Jht';
    hlabels(9)='Efn';       hlabels(19)='JpEp';     hlabels(29)='Jhn';
    hlabels(10)='Efp';      hlabels(20)='JpDn';     hlabels(30)='Jhp';
    !
    hlabels(31)='Jqht';    hlabels(41)='JDisp';    hlabels(51)='JN';
    hlabels(32)='Jqhn';    hlabels(42)='JV2Mn';    hlabels(52)='JZ';
    hlabels(33)='Jqhp';    hlabels(43)='JV2Mp';    hlabels(53)='ITER';
    hlabels(34)='JKINt';   hlabels(44)='JILST';    hlabels(54)='UEtot';
    hlabels(35)='JKINn';   hlabels(45)='JKIND';    hlabels(55)='Ubett';
    hlabels(36)='JKINp';   hlabels(46)='JL';       hlabels(56)='Ubetn';
    hlabels(37)='JSO';     hlabels(47)='JECMPAV1'; hlabels(57)='Ubetp';
    hlabels(38)='JCDIR';   hlabels(48)='JECMPAV2'; hlabels(58)='UQt';
    hlabels(39)='JCEX';    hlabels(49)='JECMPAV3'; hlabels(59)='UQn';
    hlabels(40)='JDisn';   hlabels(50)='JA';       hlabels(60)='UQp';
    !
    hlabels(61)='Uln';     hlabels(71)='Urp';      hlabels(81)='UKINP';
    hlabels(62)='Ulp';     hlabels(72)='Urc';      hlabels(82)='USO';
    hlabels(63)='UpEn';    hlabels(73)='Uht';      hlabels(83)='UCDIR';
    hlabels(64)='UpEp';    hlabels(74)='Uhn';      hlabels(84)='UCEX';
    hlabels(65)='UpDn';    hlabels(75)='Uhp';      hlabels(85)='UDisn';
    hlabels(66)='UpDp';    hlabels(76)='Uqht';     hlabels(86)='UDisp';
    hlabels(67)='UAsn';    hlabels(77)='Uqhn';     hlabels(87)='UV2Mn';
    hlabels(68)='UAsp';    hlabels(78)='Uqhp';     hlabels(88)='UV2Mp';
    hlabels(69)='Urt';     hlabels(79)='UKINT';    hlabels(89)='UECMT';
    hlabels(70)='Urn';     hlabels(80)='UKINN';    hlabels(90)='UECMN';
    !
    hlabels(91)='UECMP';   hlabels(101)='UERIGT';  hlabels(111)='LNQp';
    hlabels(92)='UROTT';   hlabels(102)='UERIGN';  hlabels(112)='LNpEn';
    hlabels(93)='UROTN';   hlabels(103)='UERIGP';  hlabels(113)='LNpEp';
    hlabels(94)='UROTP';   hlabels(104)='EHFBLN';  hlabels(114)='LNpDn';
    hlabels(95)='USQUJT';  hlabels(105)='EHFB';    hlabels(115)='LNpDp';
    hlabels(96)='USQUJN';  hlabels(106)='LNbet';   hlabels(116)='LNrt';
    hlabels(97)='USQUJP';  hlabels(107)='LNben';   hlabels(117)='LNrn';
    hlabels(98)='UCRANT';  hlabels(108)='LNbep';   hlabels(118)='LNrp';
    hlabels(99)='UCRANN';  hlabels(109)='LNQt';    hlabels(119)='LNrC';
    hlabels(100)='UCRANP'; hlabels(110)='LNQn';    hlabels(120)='LNam2n';
    !
    hlabels(121)='LNam2p';
    hlabels(122)='LNe2n';
    hlabels(123)='LNe2p';
    hlabels(124)='BlEqpN';
    hlabels(125)='BlDEqpN';
    hlabels(126)='BlOvrN';
    hlabels(127)='BlEqpZ';
    hlabels(128)='BlDEqpZ';
    hlabels(129)='BlOvrZ';
    !
    Write(iw1,'((1x,a,2x),6x,660(a,2x))') hlabels
    !
    ! HELP
    !Do i=1,129
    ! Write(iw1,'(1x,i3,a,a)',advance='NO') i,':',trim(hlabels(i))
    !End Do
    ! 1:LBL  2:BLKN  3:BLKZ  4:Jsi  5:JININ  6:A  7:N  8:Z  9:Efn 10:Efp
    ! 11:JEtot 12:Jbett 13:Jbetn 14:Jbetp 15:JQt 16:JQn 17:JQp 18:JpEn 19:JpEp 20:JpDn
    ! 21:JpDp 22:JAsn 23:JAsp 24:Jrt 25:Jrn 26:Jrp 27:Jrc 28:Jht 29:Jhn 30:Jhp
    ! 31:Jqht 32:Jqhn 33:Jqhp 34:JKINt 35:JKINn 36:JKINp 37:JSO 38:JCDIR 39:JCEX 40:JDisn
    ! 41:JDisp 42:JV2Mn 43:JV2Mp 44:JILST 45:JKIND 46:JL 47:JECMPAV1 48:JECMPAV2 49:JECMPAV3 50:JA
    ! 51:JN 52:JZ 53:ITER 54:UEtot 55:Ubett 56:Ubetn 57:Ubetp 58:UQt 59:UQn 60:UQp
    ! 61:Uln 62:Ulp 63:UpEn 64:UpEp 65:UpDn 66:UpDp 67:UAsn 68:UAsp 69:Urt 70:Urn
    ! 71:Urp 72:Urc 73:Uht 74:Uhn 75:Uhp 76:Uqht 77:Uqhn 78:Uqhp 79:UKINT 80:UKINN
    ! 81:UKINP 82:USO 83:UCDIR 84:UCEX 85:UDisn 86:UDisp 87:UV2Mn 88:UV2Mp 89:UECMT 90:UECMN
    ! 91:UECMP 92:UROTT 93:UROTN 94:UROTP 95:USQUJT 96:USQUJN 97:USQUJP 98:UCRANT 99:UCRANN 100:UCRANP
    ! 101:UERIGT 102:UERIGN 103:UERIGP 104:EHFBLN 105:EHFB 106:LNbet 107:LNben 108:LNbep 109:LNQt 110:LNQn
    ! 111:LNQp 112:LNpEn 113:LNpEp 114:LNpDn 115:LNpDp 116:LNrt 117:LNrn 118:LNrp 119:LNrC 120:LNam2n
    ! 121:LNam2p 122:LNe2n 123:LNe2p 124:BlEqpN 125:BlDEqpN 126:BlOvrN 127:BlEqpZ 128:BlDEqpZ 129:BlOvrZ
  End Subroutine thodefh
  !=======================================================================
  !> Allocates arrays depending on the nunber of shells/states and on
  !> the number of Gauss points
  !=======================================================================
  Subroutine thoalloc
    Implicit None
    Integer :: ier,ib,ND
    !
    ! number of int.points
    If(Parity) Then
       ngh=ngh_INI; ngl=ngl_INI; nleg=nleg_INI     !Yesp
    Else
       ngh=2*ngh_INI; ngl=ngl_INI; nleg=nleg_INI   !Nop
    End If
    !
    !nbx=2*n00+1                   ! maximal number of k-blocks
    !ntx=(n00+1)*(n00+2)*(n00+3)/6 ! max.num. p/n levels
    !nzx=n00                       ! maximal nz-quantum number
    !nrx=n00/2+1                   ! maximal nr-quantum number
    !nlx=n00                       ! maximal ml-quantum number
    !ndx=(n00+2)*(n00+2)/4         ! maximal dim. of one k-block
    !nhhdim=number of nonzero HH matrix elements
    !
    nzrlx=(nzx+1)*(nrx+1)*(nlx+1)   ! phy(:,:,nzrlx)
    nghl=ngh*ngl                    ! nghl=ngh*ngl
    nqx=ndx*ndx; nb2x=nbx+nbx; ndx2=ndx+ndx
    ilnqx=ilpj*nqx; ilnghl=ilpj*nghl
    nhfbx=ndx+ndx; nhfbqx=nhfbx*nhfbx; nkx=ntx; ndxs=ndx*(ndx+1)/2
    !-----------------------------------------
    !Arrays depending on gauss points
    !-----------------------------------------
    If(Allocated(xleg)) Deallocate(xleg,wleg)
    If(nleg.Gt.0) Allocate(xleg(nleg),wleg(nleg))
    If(Allocated(xh)) Deallocate(xh,wh,xl,sxl,wl,vc ,vhbn,vn,vrn,vzn,vdn,vsn,dvn,vhbp,vp,vrp,vzp,vdp,vsp,dvp,  &
                                 vSZFIn,vSFIZn,vSRFIn,vSFIRn,vSZFIp,vSFIZp,vSRFIp,vSFIRp, &
                                 fl,fli,fh,fd,fp1,fp2,fp3,fp4,fp5,fp6, fs1,fs2,fs3,fs4,fs5,fs6, &
                                 wdcor,wdcori,cou,vDHartree,vhart00,vhart01,vhart11)
    Allocate(xh(ngh),wh(ngh),xl(ngl),sxl(ngl),wl(ngl),vc(nghl,nghl))
    Allocate(vhbn(nghl),vn(nghl),vrn(nghl),vzn(nghl),vdn(nghl),vsn(nghl),dvn(nghl),  &
             vhbp(nghl),vp(nghl),vrp(nghl),vzp(nghl),vdp(nghl),vsp(nghl),dvp(nghl),  &
             vSZFIn(nghl),vSFIZn(nghl),vSRFIn(nghl),vSFIRn(nghl),  &
             vSZFIp(nghl),vSFIZp(nghl),vSRFIp(nghl),vSFIRp(nghl))
    Allocate(fl(nghl),fli(nghl),fh(nghl),fd(nghl),fp1(nghl),fp2(nghl),fp3(nghl),  &
             fp4(nghl),fp5(nghl),fp6(nghl),fs1(nghl),fs2(nghl),fs3(nghl),fs4(nghl),  &
             fs5(nghl),fs6(nghl),wdcor(nghl),wdcori(nghl),cou(nghl),vDHartree(nghl,2), &
             vhart00(nghl,nghl),vhart01(nghl,nghl),vhart11(nghl,nghl))
    If(Allocated(aka)) Deallocate(aka,ro,tau,dro,dj,NABLAR,NABLAZ,SZFI,SFIZ,SRFI,SFIR)
    Allocate(aka(nghl,2),ro(nghl,2),tau(nghl,2),dro(nghl,2),dj(nghl,2),  &
             SZFI(nghl,2),SFIZ(nghl,2),SRFI(nghl,2),SFIR(nghl,2),NABLAR(nghl,2),NABLAZ(nghl,2))
    If(Allocated(qfield)) Deallocate(qfield)
    Allocate(qfield(nghl,lambdaMax+1)) ! constraining fields: lambdaMax multipoles + neck
    If(Allocated(MEFFn)) Deallocate(MEFFn,MEFFp)
    Allocate(MEFFn(nghl),MEFFp(nghl))
    If(Allocated(geff_inv)) Deallocate(geff_inv)
    Allocate(geff_inv(nghl,2))
    !-----------------------------------------
    ! Arrays depending on configurations
    !-----------------------------------------
    If(Allocated(rk)) Deallocate(rk,ak,qh,qh1,ql,ql1,nz,nr,nl,ns,npar,id  &
         ,ia,ikb,ipb,ka,kd,tb,txb,numax,ek,dk,vk,vk1,uk,vkmax,ddc,ddc1,hfb1,lcanon)
    Allocate(rk(nqx,nb2x),ak(nqx,nb2x),qh(0:nzx,1:ngh+1)  &
         ,qh1(0:nzx,1:ngh+1),ql(0:nrx,0:nlx,1:ngl+1),ql1(0:nrx,0:nlx,1:ngl+1)  &
         ,nz(ntx),nr(ntx),nl(ntx),ns(ntx),npar(ntx),id(nbx),ia(nbx),ikb(nbx),lcanon(0:nbx,2)  &
         ,ipb(nbx),ka(nbx,2),kd(nbx,2),tb(ntx),txb(nbx),numax(0:nkx,2)  &
         ,ek(nkx,2),dk(nkx,2),vk(nkx,2),vk1(nkx,2),uk(nkx,2),vkmax(nkx,2)  &
         ,ddc(ndx,nkx,2),ddc1(ndx,nkx,2),hfb1(nhfbx,2))
    If(finite_range.or.coulomb_gaussian) Then
       If(Allocated(nrr)) Deallocate(nrr,nll,nss,noo,nzz,nzzx)
       Allocate(nrr(ntx),nll(ntx),nss(ntx),noo(ntx),nzz(ntx,ntx),nzzx(ntx))
    End If
    !-----------------------------------------
    ! HFB Arrays
    !-----------------------------------------
    If(Allocated(erhfb)) Deallocate(erhfb,drhfb,erhfb1,drhfb1)
    Allocate(erhfb(nkx),drhfb(nkx),erhfb1(nkx),drhfb1(nkx))
    If(Allocated(hfb)) Deallocate(hfb,zhfb,evvk,hfbcan,evvkcan)
    Allocate(hfb(ndx2,ndx2),zhfb(ndx2),evvk(ndx2),hfbcan(ndx,ndx),evvkcan(ndx))
    If(Allocated(AN)) Deallocate(AN,ANk,PFIU,PFID,FIU,FID,FIUR,FIDR,FIUD2N,FIDD2N,FIUZ,FIDZ)
    Allocate(AN(nqx),ANk(nqx),PFIU(ndx),PFID(ndx),FIU(ndx),FID(ndx)  &
         ,FIUR(ndx),FIDR(ndx),FIUD2N(ndx),FIDD2N(ndx),FIUZ(ndx),FIDZ(ndx))
    !-----------------------------------------
    ! Optimal LAPACK storage
    !-----------------------------------------
    If(Allocated(alwork)) Deallocate(alwork,lwork)
#if(SWITCH_ESSL==0)
    ialwork=1+6*ndx+2*ndx**2; ilwork=3+5*ndx;
    Allocate(alwork(ialwork),lwork(ilwork));alwork = 0.0; lwork = 1
#else
    ialwork=0; ilwork=5*ndx;
    Allocate(alwork(1),lwork(ilwork));alwork = 0.0; lwork = 0
#endif
    !ialwork=1; ilwork=1;
    !If(Allocated(alwork)) Deallocate(alwork,lwork)
    !Allocate(alwork(ialwork),lwork(ilwork))
    !ier=0; Call DSYEVD('V','L',ndx2,hfb,ndx2,evvk,ALWORK,-1,LWORK,-1,ier)
    !If(ier.Ne.0) Then
    !   ierror_flag=ierror_flag+1
    !   ierror_info(ierror_flag)='STOP: FATAL ERROR CONDITION IN DSYEVD'
    !   Return
    !End If
    !ialwork=Int(alwork(1)); ilwork=lwork(1)
    !If(Allocated(alwork)) Deallocate(alwork,lwork)
    !Allocate(alwork(ialwork),lwork(ilwork))
    !-----------------------------------------
    ! Eqp, U,V
    !-----------------------------------------
    If(Allocated(RVqpN)) Deallocate(RVqpN,RVqpP,RUqpN,RUqpP,REqpN,REqpP)
    Allocate(RVqpN(nuv),RVqpP(nuv),RUqpN(nuv),RUqpP(nuv),REqpN(nqp),REqpP(nqp))
    If(Allocated(KpwiP)) Deallocate(KpwiP,KpwiN,KqpN,KqpP)
    Allocate(KpwiN(nqp),KpwiP(nqp),KqpN(nqp),KqpP(nqp))
    If(Allocated(fn_T)) Deallocate(fn_T,fp_T)
    Allocate(fn_T(nqp),fp_T(nqp))
    !-----------------------------------------
    ! PNP ARRAYS: CONF. AND GAUGE ANGLE
    !-----------------------------------------
    If(Allocated(exp1iphy))Deallocate(ropj,taupj,dropj,djpj,akapj,coupj,pjk  &
         ,SZFIpj,SFIZpj,SRFIpj,SFIRpj,epj,cpj,ypj,rpj,ddepj,phypj,sinphy  &
         ,exp1iphy,exp2iphy,exp1iphym,exp2iphym)
    Allocate(ropj(nghl,ilpj,2),taupj(nghl,ilpj,2),dropj(nghl,ilpj,2)  &
         ,djpj(nghl,ilpj,2),akapj(nghl,ilpj,2),coupj(nghl,ilpj),pjk(ilpj,2)  &
         ,SZFIpj(nghl,ilpj,2),SFIZpj(nghl,ilpj,2),SRFIpj(nghl,ilpj,2)  &
         ,SFIRpj(nghl,ilpj,2),epj(ilpj,2),cpj(nkx,ilpj,2),ypj(nkx,ilpj,2)  &
         ,rpj(nkx,ilpj,2),ddepj(nqx,ilpj,nb2x),phypj(ilpj),sinphy(ilpj),  &
         exp1iphy(ilpj),exp2iphy(ilpj),exp1iphym(ilpj),exp2iphym(ilpj))
    !-----------------------------------------
    ! FIELDS INITIALIZATION (NB! optimize)
    !-----------------------------------------
    ro=zero;     tau=zero;    dro=zero;    dj=zero;  aka=zero; rk=zero
    vn=zero;     vsn=zero;    vhbn=zero;   vrn=zero; vzn=zero; vdn=zero;
    vp=zero;     vsp=zero;    vhbp=zero;   vrp=zero; vzp=zero; vdp=zero;
    dvn=zero;    dvp=zero;
    vSFIZn=zero; vSZFIn=zero; vSFIRn=zero; vSRFIn=zero;  vDHartree=zero;
    vSFIZp=zero; vSZFIp=zero; vSFIRp=zero; vSRFIp=zero;
    ! Jason
    If(Allocated(allhfb)) Then
       Do ib=1,oldnb
          Deallocate(allhfb(ib)%arr,allevvk(ib)%arr,allalwork(ib)%arr,alllwork(ib)%arr)
       End Do
       Deallocate (allhfb,allevvk,allalwork,alllwork)
       Deallocate (allIALWORK,allILWORK,allISUPPZ)
    End If
    If (Allocated(allibro)) Deallocate(allibro)
    !
  End Subroutine thoalloc
  !=======================================================================
  !> Setup routine: initializes most variables - including NAMELISTS; set
  !> up the basis and the quadrature grid; computes matrix elements of the
  !> Gogny force (if finite range); prints summary information of the run
  !=======================================================================
  Subroutine preparer(lpr)
    Use HFBTHO_gauss
    Implicit None
    Logical :: lpr
    Integer(ipr) :: iw,l,icount
    Real(pr) :: amas_base
    !
    If(n00.Eq.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)=' STOP: No more nuclei pass to the solver'
       Return
    End If
    !-----------------------------------------
    ! select the symbol of the nucleus
    !-----------------------------------------
    Call nucleus(1,npr(2),nucname)
    If(ierror_flag.Ne.0) Return
    !-----------------------------------------
    ! print headings to screen/'thoout.dat'
    !-----------------------------------------
    If(lpr) Then
       Call heading
       Call print_functional_parameters()
       Do iw=lout,lfile
          If(ierror_flag.Ne.0) Return
          If(Print_HFBTHO_Namelist) Then
             Write(iw,'(100(2x,a,f15.8))')
             Write(iw,'(100(2x,a,f15.8))') 'NAMELIST CONTENT (copy/past to hfbtho_NAMELIST.dat and modify)'
             Write(iw,'(100(2x,a,f15.8))') '-------------------------------------------------------------'
             Write(iw,HFBTHO_GENERAL)
             Write(iw,HFBTHO_INITIAL)
             Write(iw,HFBTHO_ITERATIONS)
             Write(iw,HFBTHO_FUNCTIONAL)
             Write(iw,HFBTHO_PAIRING)
             Write(iw,HFBTHO_CONSTRAINTS)
             Write(iw,HFBTHO_BLOCKING)
             Write(iw,HFBTHO_PROJECTION)
             Write(iw,HFBTHO_FEATURES)
             Write(iw,HFBTHO_TDDFT)
             Write(iw,HFBTHO_NECK)
             Write(iw,HFBTHO_TEMPERATURE)
             Write(iw,HFBTHO_DEBUG)
        End If
       End Do
    End If
    !-----------------------------------------
    ! particle number as real variable
    !-----------------------------------------
    tz(1)=Real(npr(1),Kind=pr); tz(2)=Real(npr(2),Kind=pr)
    If(fragment_properties) Then
       tz=tz_fragments(1:2)
    End If
    amas=tz(1)+tz(2)
    ! In case of blocking, force the same oscillator length for the core
    ! even-even nucleus
    amas_base = amas
    If(.Not.odd_noBlock) Then
       If(nkblo_INI(1,1).Gt.0.And.npr(1).Eq.2*(npr(1)/2)) amas_base = amas_base+one
       If(nkblo_INI(1,1).lt.0.And.npr(1).Eq.2*(npr(1)/2)) amas_base = amas_base-one
    End If
    If(.Not.odd_noBlock) Then
       If(nkblo_INI(2,1).Gt.0.And.npr(2).Eq.2*(npr(2)/2)) amas_base = amas_base+one
       If(nkblo_INI(2,1).lt.0.And.npr(2).Eq.2*(npr(2)/2)) amas_base = amas_base-one
    End If
    drhoi=zero
    !-----------------------------------------
    ! default combinations
    !-----------------------------------------
    chargee2=e2charg
    coex=-chargee2*(three/pi)**p13; cex=-0.750_pr*coex
    !--------------------------------------------------------
    ! hbzero from forces [hqc**2/(two*amu)] EOedit for SV-min
    !--------------------------------------------------------
    if(.not.hb0_charge_dependent) then
       hbzeron = hbzero
       hbzerop = hbzero
    !   hb0=hbzero; If (use_cm_cor) hb0=hb0*(one-one/amas)
    !   hb0n=hbzeron; If (use_cm_cor) hb0n=hb0n*(one-one/amas)
    !   hb0p=hbzerop; If (use_cm_cor) hb0p=hb0p*(one-one/amas)
    !else if () then
      !STUFF
    else
    !   hbzero = 0.5_pr*(hbzeron + hbzerop)
       hb0=hbzero; If (use_cm_cor) hb0=hb0*(one-one/amas)
       hb0n=hbzeron; If (use_cm_cor) hb0n=hb0n*(one-one/amas)
       hb0p=hbzerop; If (use_cm_cor) hb0p=hb0p*(one-one/amas)
    endif
    hb0=hbzero; If (use_cm_cor) hb0=hb0*(one-one/amas)
    hb0n=hbzeron; If (use_cm_cor) hb0n=hb0n*(one-one/amas)
    hb0p=hbzerop; If (use_cm_cor) hb0p=hb0p*(one-one/amas)
    !-----------------------------------------
    ! basis parameter q
    !-----------------------------------------
    beta0=q; q=Exp((3.0_pr*Sqrt(5.0_pr/(16.0_pr*pi)))*beta0)
    !-----------------------------------------
    ! basis parameters b0,bp,bz
    !-----------------------------------------
    If(b0.Le.zero) Then
       ! define oscillator frequency from default with empirical factor 1.2,
       ! and set length accordingly
       r00=r0*amas_base**p13; r02=r00**2; r04=r02**2
       hom=41.0_pr*amas_base**(-p13)*r0
       b0=Sqrt(two*hbzero/hom)
    Else
       ! define oscillator frequency from user-defined length, and set default
       ! empirical factor accordingly
       hom=hqc**2/(amn*b0**2)
       r0=(hom/41.0_pr)*amas**(p13)
       r00=r0*amas**p13; r02=r00**2; r04=r02**2
    End If
    bp=b0*q**(-one/6.0_pr); bz=b0*q**(one/3.0_pr); bpp=bp*bp
    !-----------------------------------------
    ! constraint in terms of beta
    !-----------------------------------------
    ty20=Sqrt(5.0_pr/pi)*hom/b0**2/two
    !-----------------------------------------
    ! projection: number of grid points
    !-----------------------------------------
    keypj=Max(1,keypj); ilpj=keypj; ilpj2=ilpj**2;
    If(iproj.Eq.0) Then
       npr1pj=npr(1); npr2pj=npr(2)
    Else
       npr1pj=npr(1)+npr1pj; npr2pj=npr(2)+npr2pj
    End If
    !-----------------------------------------
    ! blocking window
    !-----------------------------------------
    !pwiblo=Min(Max(25.0_pr/Sqrt(Real(npr(1)+npr(2),Kind=pr)),2.0_pr),8.0_pr)
    pwiblo=1.0_pr
    !-----------------------------------------
    ! THO
    !-----------------------------------------
    ass=zero; iasswrong=0
    !-----------------------------------------
    ! iterations
    !-----------------------------------------
    etot=zero; varmas=zero; rms=zero; ept=-two; del=one; alast=-seven; siold=one
    varmasNZ=zero; pjmassNZ=zero; ass=zero; skass=zero
    !---------------------------------------------------------
    ! statistics to screen('lout')/file('lfile')
    !---------------------------------------------------------
    If(lpr) Then
       Do iw=lout,lfile
          Write(iw,*)
          Write(iw,'(a)')             '  ---------------------------------------'
          Write(iw,'(a)')             '        Characteristics of the run       '
          Write(iw,'(a)')             '  ---------------------------------------'
          Write(iw,'(a,i5)')          '  Output file ................: ',lfile
          Write(iw,'(a,2x,a2,i4)')    '  Nucleus ....................: ',nucname,npr(1)+npr(2)
          Write(iw,'(a,i5)')          '  Number of HO shells ........: ',n00
          Write(iw,'(a,f20.14)')      '  HO length b0 (fm) ..........: ',b0
          Write(iw,'(a,f8.3,a,f8.3)') '  Basis deformation ..........:  beta0=',beta0,' q=',q
          Write(iw,'(a,5(1x,e15.8))') '  HO: b0,1/b0,bp,bz,q ........: ',b0,one/b0,bp,bz,q
          Write(iw,'(a,3(1x,e15.8))') '  h**2/(2m_n), cmc, e**2 .....: ',hbzeron,hb0n,chargee2
          Write(iw,'(a,3(1x,e15.8))') '  h**2/(2m_p), cmc, e**2 .....: ',hbzerop,hb0p,chargee2
          Write(iw,'(a,2(1X,e15.8))') '  hom=f*41.0_pr*A^{-1/3}, f...: ',hom,r0
          If(automatic_basis) Then
             Write(iw,'(a)')       '    Adjusted basis is ........:   ON'
          End If
          If(iLST.Eq.0)  Then         ! HFB+HO case only
             iLST1=0
             Write(iw,'(a)')          '  THO basis is ...............:  OFF'
          Else                        ! HFB+THO case
             Write(iw,'(a)')          '  THO basis is ...............:   ON'
             If(iLST.Gt.0) Then       ! HFB+THO  only
                iLST1=1
                If(inin.Gt.0) Then
                   ierror_flag=ierror_flag+1
                   ierror_info(ierror_flag)=' Stop: Forbidden iLST>0, inin>0 combination.'
                   Return
                End If
                Write(iw,'(a)')       '    THO parameters from tholst.wel'
             Else                     ! HFB+THO after HFB+HO
                iLST1=0
                Write(iw,'(a)')       '    HFB+THO after a HFB+HO run    '
             End If
          End If
          Write(iw,'(a,i5)')          '  Maximal number of iterations: ',maxi
          Write(iw,'(a,f6.3)')        '  Initial mixing parameter ...: ',xmix
          If(inin.Eq.1)  Then
             Write(iw,'(a)')          '  Initial w.f. ...............:  from spherical scratch'
          End If
          If(inin.Eq.2)  Then
             Write(iw,'(a)')          '  Initial w.f. ...............:  from prolate scratch'
          End If
          If(inin.Eq.3)  Then
             Write(iw,'(a)')          '  Initial w.f. ...............:  from oblate scratch'
          End If
          If(inin.Lt.0) Then
             Write(iw,'(a)')          '  Initial wave functions from :  tape'
          End If
          Write(iw,'(a,3x,a)')        '  Energy functional ..........: ',skyrme
          If(finite_range) Then
             Write(iw,'(a,f6.2,a)')   '    with a finite-range central force'
          End If
          If(icou.Eq.-4) Write(iw,'(a)')    '    direct Coulomb by sum of Gaussians, exchange Coulomb with Slater approximation'
          If(icou.Eq.-3) Write(iw,'(a)')    '    direct, exchange and pairing Coulomb by sum of Gaussians'
          If(icou.Eq.-2) Write(iw,'(a)')    '    direct and exchange Coulomb by sum of Gaussians'
          If(icou.Eq.-1) Write(iw,'(a)')    '    direct Coulomb force only by sum of Gaussians'
          If(icou.Eq. 0) Write(iw,'(a)')    '    without Coulomb forces'
          If(icou.Eq. 1) Write(iw,'(a)')    '    direct Coulomb by substitution method'
          If(icou.Eq. 2) Write(iw,'(a)')    '    direct Coulomb by substitution method, exchange Coulomb with Slater approximation'
          If(kindhfb.Lt.0) Then
             Write(iw,'(a)')          '  Lipkin-Nogami procedure is .:   ON'
          Else
             Write(iw,'(a)')          '  Lipkin-Nogami procedure is .:  OFF'
          End If
          If(ilpj-1.Eq.0) Then
             Write(iw,'(a)')          '  PAV procedure is ...........:  OFF'
          Else
             Write(iw,'(a)')          '  PAV procedure is ...........:   ON'
             Write(iw,'(a,i5)')       '    Number of gauge points....: ',keypj
          End If
          If(icstr.Eq.0) Then
             Write(iw,'(a)')          '  Constraint calculation is ..:  OFF'
          Else
             Write(iw,'(a)')          '  Constraint calculation is ..:   ON'
             icount=0
             Do l=1,8
                If(Abs(lambda_active(l)).Gt.0) Then
                   icount=icount+1
                   Write(iw,'(a,i1,a,i1,a,f8.3)') '    Constraint ',icount,' .............: lambda=',l, &
                                                  ' Ql=',multRequested(l)
                End If
             End Do
             If(neck_constraints) Then
                icount=icount+1
                Write(iw,'(a,i1,a,a,f8.3)') '    Neck       ',icount,' .............: lambda=0', &
                                                  ' Ql=',neckRequested
             End If
          End If
          If(keyblo(1).Ne.0) Then
             Write(iw,'(a)')          '  Neutron blocking is ........:   ON'
          End If
          If(keyblo(2).Ne.0) Then
             Write(iw,'(a)')          '  Proton blocking is .........:   ON'
          End If
          If(switch_on_temperature) Then
             Write(iw,'(a,f6.2,a)')   '  Temperature T ..............: ',temper,' MeV'
          Else
             Write(iw,'(a,f6.2)')     '  Temperature T ..............:   0.00 MeV'
          End If
          If(pairing_regularization) Then
             Write(iw,'(a,f6.2,a)')   '  Pairing regularization is ..:   ON'
          Else
             Write(iw,'(a,f6.2,a)')   '  Pairing regularization is ..:  OFF'
          End If
          If(collective_inertia) Then
             Write(iw,'(a,f6.2,a)')   '  Collective inertia is ......:   ON'
          Else
             Write(iw,'(a,f6.2,a)')   '  Collective inertia is ......:  OFF'
          End If
          If(fission_fragments) Then
             Write(iw,'(a,f6.2,a)')   '  Fission fragments are ......:   ON'
          Else
             Write(iw,'(a,f6.2,a)')   '  Fission fragments are ......:  OFF'
          End If
          Write(iw,'(a,i3)')          '  Restart indicator ..........: ',inin
          If(nbroyden.Eq.0) Then
             Write(iw,'(a,i3)')       '  Linear mixing ..............: ',nbroyden
          Else
             Write(iw,'(a,i3)')       '  Broyden mixing (#iterations): ',nbroyden
          End If
       End Do
    End If
    !-----------------------------------------
    ! BASIS, GAUSS POINTS, HOWF
    !-----------------------------------------
    Call gfv                      ! factorials
    If(ierror_flag.Ne.0) Return
    Call base0(lpr)               ! basis space (calculate configurational space)
    If(ierror_flag.Ne.0) Return
    Call thoalloc                 ! global allocation
    If(ierror_flag.Ne.0) Return
    Call gausspoints              ! GAUSS mesh points
    If(ierror_flag.Ne.0) Return
    Call recompute_coulomb_expansion ! expansion of the Coulomb potential on Gaussians
    If(ierror_flag.Ne.0) Return
    Call base(lpr)                ! oscillator configurations (set up quantum numbers)
    If(ierror_flag.Ne.0) Return
    Call gaupol(lpr)              ! basis wf at gauss mesh points
    If(ierror_flag.Ne.0) Return
    If(finite_range.or.coulomb_gaussian) Then
       Call gogny_matrix_elements
    endif
    If(ierror_flag.Ne.0) Return
    !
  End Subroutine preparer
  !====================================================================
  !> Defines and stores the weights and nodes of integration in the
  !> dimensionless coordinates \f$ \xi\f$ and \f$ \eta \f$. Defines the
  !> scaling transformation for the THO basis.
  !====================================================================
  Subroutine coordinateLST(lpr)
    Use HFBTHO_THO, Only: f01234
    Implicit None
    Logical :: lpr
    Integer(ipr) :: i,il,ih
    If(iLST1.Eq.0) Then
       ! HO-basis
       Do il=1,ngl
          Do ih=1,ngh
             i=ih+(il-1)*ngh
             fh(i)=bz*xh(ih)
             fl(i)=bp*Sqrt(xl(il))
             wdcor(i)=pi*wh(ih)*wl(il)*bz*bp*bp
             wdcori(i)=one/wdcor(i)
          End Do
       End Do
    Else
       ! THO basis
       Call f01234(.False.)
       If(ierror_flag.Ne.0) Return
    End If
    !
    Call optHFBTHO                ! optimal HO/THO combinations
    If(ierror_flag.Ne.0) Return
    !
  End Subroutine coordinateLST
  !====================================================================
  !> Self-consistent loop
  !====================================================================
  Subroutine iter(lpr)
    Implicit None
    Logical :: lpr
    Real(pr)       :: time5,assprn,delln(2)
    Real(pr), Save :: time
    Integer(ipr)   :: iw,it,ite
    integer(ipr)   :: wct1,wct2,countrate,countmax
    !---------------------------------------------------
    ! print to screen('lout')/thoout.dat('lfile')
    !---------------------------------------------------
    Do iw=lout,lfile
       If(iLST.Eq.0) Then
          Write(iw,'(a,f7.3,4(a,i3),a)')  &
               '  |HFB+HO> iterations(b0=',b0,', Nsh=',n00,  &
               ', inin=',inin,', N=',npr(1),', Z=',npr(2),')...'
       Else
          If(iLST1.Eq.0.Or.iasswrong(3).Ne.0) Then
             If(iasswrong(3).Ne.0) Then
                Write(iw,'(a,f7.3,a,i3,a)')  &
                     '  |HFB+THO substituted by HFB+HO> iterations (b0=',  &
                     b0,', Nsh=',n00,')...'
             Else
                Write(iw,'(a,f7.3,a)')'  towards |hfb+tho> iterations...'
                Write(iw,'(a,f7.3,a)')
                Write(iw,'(a,f7.3,a,i3,a)')  &
                     '  |Preliminary HFB+HO> iterations (b0=',b0,', Nsh=',n00,')...'
             End If
          Else
             If(itass.Eq.1) Then
                Write(iw,'(2(a,f7.3),a,i3,a)')  &
                     '  |HFB+THO> iterations(b0=',b0,', neutron density decay=',  &
                     decay,', Nsh=',n00,')...'
             Else
                Write(iw,'(2(a,f7.3),a,i3,a)')  &
                     '  |HFB+THO> iterations(b0=',b0,', proton density decay=',  &
                     decay,', Nsh=',n00,')...'
             End If
          End If
       End If
       Write(iw,1)
       Write(iw,'(20(a))') '  i','          si ','    mix ','  beta',    &
      &    '      Etot ','      A ','      rn','      rp ','        En', &
      &    '      Dn','      Ep','      Dp','        Ln  ','    Lp ',    &
      &    '    time  time(Gog.)'
       Write(iw,1)
1      Format(2x,130('-'))
    End Do
    !---------------------------------------------------------------------
    ! main hfb iteration loop
    !---------------------------------------------------------------------
    iError_in_HO=0; iError_in_THO=0; time=0.0_pr; time5=0.0_pr
    Do ite=1,maxi
       Call system_clock(wct1,countrate,countmax)
       !
       iiter=ite
       !
       If (lpr.Or.iiter.Eq.1) Then
          assprn=ass(1); If(assprn.Gt.ass(2)) assprn=-ass(2) ! protons come with '-'
          delLN=del+frdel; If(kindhfb.Lt.0) delLN=del+frdel+ala2         ! LN case
          ! during iterations print
          Do iw=lout,lfile
             !If(Max(Abs(drhoi(1)),Abs(drhoi(2))).Gt.1.0e-10_pr) Then
             !   Write(*,*) '  WARNING! Int(Dro)=',Max(Abs(drhoi(1)),Abs(drhoi(2)))
             !End If
             Write(iw,2) iiter,bbroyden,si,xmix,bet,etot,varmas,rms(1),rms(2),ept(1)+frept(1),delLN(1), &
                         ept(2)+frept(2),delLN(2),alast(1),alast(2),time,wct_gogny
          End Do
       End If
       !-------------------------------------------------
       ! HFBDIAG
       !-------------------------------------------------
       Do it=itmin,itmax
          Call hfbdiag(it,0)   ! hfb diagonalization with minimal canonical
          If(ierror_flag.Ne.0) Return
       End Do
       !-------------------------------------------------
       ! EXPECT, DENSIT, COULOMB, FIELD, GAMDEL
       !-------------------------------------------------
       Call expect(.False.)    ! expectation values
       If (numberCons.Gt.0) Call getLagrange(ite)   ! new Lagrange parameters for constraints
       If(ierror_flag.Ne.0) Return
       Call field              ! new fields
       If(ierror_flag.Ne.0) Return
       Call gamdel(.false.,.true.)    ! hf-matrix
       If(ierror_flag.Ne.0) Return
       !-------------------------------------------------
       ! Dumping control (old linear mixing)
       !-------------------------------------------------
       xmix0=0.1 !original 0.1
       If(si.Lt.siold) Then
          xmix=Min(xmax,xmix * 1.130_pr);  !old value 1.13
       Else
          xmix=xmix0
       End If
       siold=si
       !-------------------------------------------------
       ! time per iteration
       !-------------------------------------------------
       call system_clock(wct2,countrate,countmax)
       time=(wct2-wct1)/real(countrate,kind=pr)
       time5=time5+time
       !-------------------------------------------------
       ! Solution is OK within the iteration limit
       !-------------------------------------------------
       If(iiter.Ge.2.And.si.Lt.epsi) Then
          If(iLST1.Eq.0) Then
             iError_in_HO=0
          Else
             iError_in_THO=0
          End If
          ! iteration interrupted print
          If(.Not.lpr) Then
             delLN=del+frdel; If(kindhfb.Lt.0) delLN=del+frdel+ala2
             Do iw=lout,lfile
                Write(iw,3) iiter,bbroyden,si,xmix,bet,etot,varmas,rms(1),rms(2),ept(1)+frept(1),delLN(1), &
                     ept(2)+frept(2),delLN(2),alast(1),alast(2),time !Max(Abs(drhoi(1)),Abs(drhoi(2)))
                Write(iw,'(a,f8.3,a)') '  Total CPU time=',time5/60.0_pr,' minutes'
             End Do
          End If
          ! converged print
          Do iw=lout,lfile
             Write(iw,4) iiter,si,iError_in_HO,iError_in_THO
             Write(iw,'(a,f8.3,a)') '  Total CPU time=',time5/60.0_pr,' minutes'
          End Do
          iiter=iiter+1
          Return
       End If
       !-------------------------------------------------
       ! Slow convergence and lambda >0 (stop iterations)
       !-------------------------------------------------
       If(iiter.Ge.1000.And.(alast(1).Gt.zero.Or.alast(2).Gt.zero)) Exit
       !
    End Do    ! ite
    iiter=iiter+1
    !-------------------------------------------------
    ! Solution interrupted due to iterations limit
    !-------------------------------------------------
    If(iLST1.Eq.0) Then
       iError_in_HO=-1
    Else
       iError_in_THO=-1
    End If
    delLN=del; If(kindhfb.Lt.0) delLN=del+ala2
    ! iterations limit print
    Do iw=lout,lfile
       Write(iw,2) iiter,bbroyden,si,xmix,bet,etot,varmas,rms(1),rms(2),ept(1),delLN(1), &
            ept(2),delLN(2),alast(1),alast(2),Max(Abs(drhoi(1)),Abs(drhoi(2)))
       Write(iw,5) iiter,si,iError_in_HO,iError_in_THO
       Write(iw,'(a,f8.3,a)') '  Total CPU time=',time5/60.0_pr,' minutes'
    End Do
    !-------------------------------------------------
2   Format(i4,a,1x,f12.8,f5.2,f7.3,f13.6,1x,f6.1,2(f8.3),' | ',4(f8.3),' | ',4(f8.3))
!2   Format(i4,a,1x,f12.8,f5.2,f7.3,f13.6,1x,f6.1,2(f8.3),' | ',4(f8.3),' | ',20(f8.3))
3   Format(2x,130('-'),/,'  *   iteration interrupted after',i4,' steps   si=',f17.10,' ho=',i3,' tho=',i3,/,2x,130('-'))
4   Format(2x,130('-'),/,'  *   iteration converged   after',i4,' steps   si=',f17.10,' ho=',i3,' tho=',i3,/,2x,130('-'))
5   Format(2x,130('-'),/,'  *   iterations limit interrupt after',i4,' steps   si=',f17.10,' ho=',i3,' tho=',i3,/,2x,130('-'))
    !-------------------------------------------------
  End Subroutine iter
  !====================================================================
  !> Block-diagonalization of the HFB matrix, loop over the Fermi level and
  !> definition of the density matrix and pairing tensor in configuration
  !> space.
  !====================================================================
  Subroutine hfbdiag(it,icanon)
#if(USE_OPENMP==1)
    Use omp_lib
#endif
    Implicit None
    Logical :: lpr_pwi,norm_to_improve
    Character(Len=1) :: char1,char2,char3
    Integer(ipr) :: iw,it,i0,icanon,ibiblo,ier,i,j,k,k0,kl,lc,ib,nd,  &
                    nhfb,n1,n2,kaib,m,ndk,nd1,nd2,kdib,k1,k2,id1,id2, &
                    n12,n21,ntz,nhhph,nhhpp,ibro,ibroib,i_uv,i_eqp,jj,&
                    tid,IL,IU,NUMFOU,jlwork,jalwork,ldw,ldi,ii
    Real(pr) :: al,al2,emin,hla,dla,pn,eqpe,ela,enb,enb1,ekb, &
                s1,s2,s3,alnorm,sitest,fac1,fac2,fT,exponent, &
                VL,VU,ABSTOL,buffer
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:),f_T(:)
    Integer(ipr), Allocatable :: ISUPPZ(:),lwork_p(:)
    Real(pr), Allocatable :: alwork_p(:),eigenv(:),eigenf(:,:),hfbmat(:,:)
    Real(pr), External :: DLAMCH
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('hfbdiag',0)
    !
    If(it.Eq.1) Then
       EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN; f_T=>fn_T
    Else
       EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP; f_T=>fp_T
    End If
    KpwiPo=0; KqpPo=0; f_T=zero
    !
    nhhph=(it-1)*nhhdim; nhhpp=(it+1)*nhhdim
    If(.Not. Allocated(allhfb)) Then
       oldnb = nb ! for destroying data structures in next computation
       Allocate(allhfb((nb)),allevvk((nb)),allalwork((nb)),alllwork((nb)),allIALWORK(nb),allILWORK(nb))
       Allocate(allISUPPZ(nb))
       Do ib=1,nb
          nhfb=2*id(ib)
          Allocate(allhfb(ib)%arr(1:nhfb,1:nhfb))
          Allocate(allevvk(ib)%arr(1:nhfb))
          !jalwork=1+6*nhfb+2*nhfb**2; allIALWORK(ib)=jalwork  ! DSYEVD
          !jlwork=3+5*nhfb; allILWORK(ib)=jlwork               ! DSYEVD
#if(SWITCH_ESSL==0)
          jalwork=26*nhfb; allIALWORK(ib)=jalwork  ! DSYEVR
          jlwork=10*nhfb; allILWORK(ib)=jlwork     ! DSYEVR
#else
          jalwork=8*nhfb; allIALWORK(ib)=jalwork   ! DSYEVX
          jlwork=5*nhfb; allILWORK(ib)=jlwork      ! DSYEVX
#endif
          Allocate(allalwork(ib)%arr(1:jalwork))
          Allocate(alllwork(ib)%arr(1:jlwork))
          Allocate(allISUPPZ(ib)%arr(1:2*nhfb)) ! DSYEVR
       End Do
    End If
    !
    ABSTOL=2.0_pr*DLAMCH('S')
    !
    If (.Not. Allocated(allibro)) Then
       Allocate(allibro(1:NB))
       allibro(1)=0
       Do ib=2,NB
          allibro(ib) = allibro(ib-1) + (ID(ib-1)*(ID(ib-1)+1)/2)
       End Do
    End If
    !
    !------------------------------------------------------------------
    ! Loop the internal normalization
    !------------------------------------------------------------------
    !sitest=Max(Min(0.10_pr,si*0.010_pr),0.000010_pr)
    sitest=Min(0.10_pr,si*0.010_pr)
    norm_to_improve=.True.; inner(it)=-1; sumnz(it)=one
    Do While(norm_to_improve)
       !
       inner(it)=inner(it)+1
       !
       If(Abs(sumnz(it)).Lt.sitest.Or.inner(it).Eq.20) norm_to_improve=.False.
       !
       sumnz(it)=zero; entropy(it)=zero; v2min(it)=one; Dispersion(it)=zero
       !
       kl=0; emin=1000.0_pr; al=ala(it)
       !
       ! blocking
       If(iparenti(it).Eq.0) blomax(it)=0
       blo123d(it)=0; blok1k2d(it)=0; blocanon(it)=0;
       ibiblo=bloblo(keyblo(it),it)
       !------------------------------------------------------------------
       ! Runs over blocks
       !------------------------------------------------------------------
       i_uv=0; i_eqp=0
       lc=0; lcanon(0,it)=0; klmax=0
!$OMP Parallel Default(None) &
!$OMP& SHARED(nb,id,ia,it,nbx,allibro,brin,allhfb,allevvk, &
!$OMP&        allALWORK,allLWORK,allIALWORK,allILWORK,nhhph,nhhpp,al, &
!$OMP&        zhfb,ndx2,allISUPPZ,ABSTOL) &
!$OMP& PRIVATE(ib,nd,nhfb,i0,m,ibro,n1,nd1,n2,nd2,hla,dla,ier,tid,char1,char2, &
!$OMP&         NUMFOU,IL,IU,VL,VU,eigenf,eigenv,hfbmat,ISUPPZ,alwork_p,lwork_p, &
!$OMP&         ldw,ldi,char3,jalwork,jlwork)
#if(USE_OPENMP==1)
       tid = OMP_GET_THREAD_NUM()
#endif
!$OMP DO SCHEDULE(DYNAMIC)
       Do ib=1,nb
          nd=id(ib); nhfb=nd+nd; i0=ia(ib); m=ib+(it-1)*nbx; ibro=allibro(ib)
#if(SWITCH_ESSL==0)
          jalwork=26*nhfb; jlwork=10*nhfb
#else
          jalwork=8*nhfb; jlwork=5*nhfb
#endif
          allhfb(ib)%arr(1:nhfb,1:nhfb)=0.0_pr; allevvk(ib)%arr(1:nhfb)=0.0_pr
          allALWORK(ib)%arr(1:jalwork)=0.0_pr; allLWORK(ib)%arr(1:jlwork)=0;
          allISUPPZ(ib)%arr(1:2*nhfb)=0
          !------------------------------------------------------------------
          !  hfb-matrix
          !------------------------------------------------------------------
          Allocate(hfbmat(nhfb,nhfb))
          Do n1=1,nd
             nd1=n1+nd
             Do n2=1,n1
                nd2=n2+nd; ibro=ibro+1
                hla=brin(nhhph+ibro); dla=brin(nhhpp+ibro)
                hfbmat(n1,n2)=hla;    hfbmat(nd2,n1)=dla
                hfbmat(nd1,n2)=dla;   hfbmat(nd1,nd2)=-hla
             End Do
             hfbmat(n1,n1)  =hfbmat(n1,n1)  -al
             hfbmat(nd1,nd1)=hfbmat(nd1,nd1)+al
          End Do
          char1='V'; char2='I'; char3='L'; NUMFOU=0
          VL=0.0_pr ;VU=0.0_pr; IL=1; IU=nhfb; ldw=allIALWORK(ib); ldi=allILWORK(ib)
          Allocate(eigenv(nhfb)); eigenv(1:nhfb)=0.0_pr
          Allocate(eigenf(nhfb,nhfb)); eigenf(1:nhfb,1:nhfb)=0.0_pr
          Allocate(alwork_p(ldw)); Allocate(lwork_p(ldi))
#if(SWITCH_ESSL==0)
          ier=0; Allocate(ISUPPZ(2*nhfb))
          Call DSYEVR(char1,char2,char3,nhfb,hfbmat,nhfb,VL,VU,IL,IU,ABSTOL,NUMFOU,    &
                      eigenv,eigenf,nhfb,ISUPPZ,alwork_p,ldw,lwork_p,ldi,ier)
#else
          ier=0; Allocate(ISUPPZ(nhfb))
          Call DSYEVX('V','A','L',nhfb,hfbmat,nhfb,VL,VU,IL,IU,ABSTOL,NUMFOU, &
                      eigenv,eigenf,nhfb,alwork_p,ldw,lwork_p,ISUPPZ,ier)
#endif
          allevvk(ib)%arr(1:nhfb) = eigenv(1:nhfb)
          allhfb(ib)%arr(1:nhfb,1:nhfb) = eigenf(1:nhfb,1:nhfb)
          If(ier.NE.0) Then
             Write(6,*)'The algorithm failed to compute eigenvalues.'
#if(USE_OPENMP==1)
             Write(6,*)'I am',tid,' and I am working on array ',ib,ier
#endif
          End If
          Deallocate(eigenf,eigenv,hfbmat,ISUPPZ,alwork_p,lwork_p)
       End Do ! ib
!$OMP End Do
!$OMP End Parallel
       Do ib=1,NB
          nd=id(ib); nhfb=nd+nd; i0=ia(ib); m=ib+(it-1)*nbx; ibro=allibro(ib)
          !------------------------------------------------------------------
          ! Blocking
          !------------------------------------------------------------------
          ! external blocking
          If(iiter.Eq.1.And.inner(it).Eq.0) Then
             If(iparenti(it).Ne.0.And.keyblo(it).Eq.0) Then
                ! eventually charging
                !   keyblo(it)=1
                !   bloblo(keyblo(it),it)=ib
                !   blo123(keyblo(it),it)=requested level (k0)
                Call requested_blocked_level(ib,it)
                If(ierror_flag.Ne.0) Return
                ibiblo=bloblo(keyblo(it),it)
             End If
          End If
          ! general blocking
          k0=0
          If(ibiblo.Eq.ib) Then
             If(iiter.Eq.1.And.inner(it).Eq.0) Then
                ! blocked level as in the even-even nucleus
                k0=blo123(keyblo(it),it); ndk=k0+nd
                Do n2=1,nd
                   nd2=n2+nd
                   hfb1(n2,it)=allhfb(ib)%arr(n2,ndk)    !U
                   hfb1(nd2,it)=allhfb(ib)%arr(nd2,ndk)  !V
                End Do
                ! number of states in the block to be tested
                blocross(it)=Min(blomax(it)+10,nd)
             End If
             ! overlap between new and old blocked levels
             s3=zero
             Do n1=1,blocross(it)
                ndk=n1+nd; s1=zero
                Do n2=1,nd
                   nd2=n2+nd
                   s1=s1+Abs(hfb1(nd2,it)*allhfb(ib)%arr(nd2,ndk)) !VV
                   s1=s1+Abs(hfb1(n2,it)*allhfb(ib)%arr(n2,ndk))   !UU
                End Do
                If(s1.Gt.s3) Then
                   s3=s1; k0=n1
                End If
             End Do
             blo123d(it)=k0
             If(.Not.norm_to_improve) Then
                ! find maximal HO component
                ndk=k0+nd
                s1=zero
                Do n1=1,nd
                   nd1=n1+nd
                   hfb1(n1,it)=allhfb(ib)%arr(n1,ndk); hfb1(nd1,it)=allhfb(ib)%arr(nd1,ndk)
                   s2=Max(s1,Abs(allhfb(ib)%arr(n1,ndk)),Abs(allhfb(ib)%arr(nd1,ndk)))
                   If(s2.Gt.s1) Then
                      s1=s2; i=n1+i0  ! labels in k[k1,k2] numbering
                   End If
                End Do
                ! print blocked state
                Do iw=lout,lfile
                   Write(iw,'(4x,a,2(a,i3),2x,3(a,1x,f12.8,1x),(i3,a,i3,1x),a)')  &
                        protn(it),' Blocking: block=',ib,  &
                        ' state=',k0,  &
                        ' Eqp=',allevvk(ib)%arr(k0+nd),  &
                        ' Dqpe=',allevvk(ib)%arr(k0+nd)-eqpmin(it),  &
                        ' Ovlp=',s3  &
                        , keyblo(it),'/',blomax(it)  &
                        , tb(i)
                End Do
                ! ieresbl=6, 'BLKN','BLKZ'
                ereslbl(it)=tb(i)
                If(it.Eq.1) Then
                   ! 'BlEqpN','BlDEqpN','BlOvrN'
                   eresbl(1)=allevvk(ib)%arr(k0+nd); eresbl(2)=allevvk(ib)%arr(k0+nd)-eqpmin(it); eresbl(3)=s1
                Else
                   ! 'BlEqpZ','BlDEqpZ','BlOvrZ'
                   eresbl(4)=allevvk(ib)%arr(k0+nd); eresbl(5)=allevvk(ib)%arr(k0+nd)-eqpmin(it); eresbl(6)=s1
                End If
             End If
          End If
          !------------------------------------------------------------------
          ! Run over all qp states k in the block
          !------------------------------------------------------------------
          kaib=kl
          Do k=1,nd
             ndk=k+nd
             ! referent spectra
             pn=zero
             Do i=1,nd
                hla=allhfb(ib)%arr(i+nd,ndk)**2; pn=pn+hla
             End Do
             ! Blocking
             If(k.Eq.k0) Then
                n1=k0+nd
                Do i=1,nd
                   hla=allhfb(ib)%arr(i+nd,n1)**2; dla=allhfb(ib)%arr(i,n1)**2; pn=pn-half*(hla-dla)
                End Do
             End If
             eqpe=allevvk(ib)%arr(nd+k); ela=eqpe*(one-two*pn)
             enb=ela+al;                 ekb=Sqrt(Abs(eqpe**2-ela**2))
             i_eqp=i_eqp+1
             !------------------------------------------------------------------
             ! cut-off condition: energy pwi + Fermi cut-off function
             !------------------------------------------------------------------
                                                              exponent=Huge(1.0_pr)
             If(Abs(100.0_pr*(enb-pwi)).Lt.Log(Huge(1.0_pr))) exponent=Exp(100.0_pr*(enb-pwi))
             If(basis_HFODD) Then
                lpr_pwi=enb.Le.pwi !jacek sharp cut off for hfodd
                !lpr_pwi=enb.Le.pwi.Or.Abs(one/(one+exponent)).Gt.cutoff_tol
             Else
                lpr_pwi=enb.Le.pwi.Or.Abs(one/(one+exponent)).Gt.cutoff_tol
             End If
             if(finite_range) then
                lpr_pwi = .true.
             endif
             !------------------------------------------------------------------
             ! Remember the whole qp solution
             !------------------------------------------------------------------
             If(.Not.norm_to_improve) Then
                EqpPo(i_eqp)=eqpe                       ! Eqp_k
                If(lpr_pwi) KqpPo(kl+1)=i_eqp           ! below pwi otherwise zero
                If(lpr_pwi) KpwiPo(kl+1)=i_uv           ! below pwi otherwise zero
                Do n2=1,nd
                   nd2=n2+nd; i_uv=i_uv+1
                   UqpPo(i_uv)=allhfb(ib)%arr(n2,ndk)   ! U_ak
                   VqpPo(i_uv)=allhfb(ib)%arr(nd2,ndk)  ! V_ak
                End Do
             End If
             !------------------------------------------------------------------
             ! Define Fermi-Dirac occupations
             !------------------------------------------------------------------
             fT=zero
             If(switch_on_temperature.And.temper.Gt.1.D-12) Then
                fT = half*(one-Tanh(half*eqpe/temper))
                ! factor two comes from K>0 states only
                buffer = zero
                If(fT.Gt.zero.And.fT.Lt.one) Then
                   buffer = two*fT*Log(fT) + two*(one-fT)*Log(one-fT)
                End If
                entropy(it) = entropy(it) - buffer
                f_T(i_eqp) = fT
             End If
             !------------------------------------------------------------------
             ! Pairing window
             !------------------------------------------------------------------
             If(lpr_pwi) Then
                kl=kl+1                                      !number of active states
                If(k0.Eq.k) blok1k2d(it)=kl                  !blocking: dynamic #: k[k1,k2] numbering
                If((eqpe.Le.emin).And.(pn.Gt.0.0001_pr)) Then   !to avoid unocc at magic numbers
                   emin=eqpe; alnorm=pn                      !min qpe and its occupation
                End If
                erhfb(kl)=enb; drhfb(kl)=ekb; uk(kl,it)=pn     !ref.s.p. energies, deltas, occupancies
                sumnz(it)=sumnz(it)+two*pn+two*(one-two*pn)*fT !internal normalization
             End If
          End Do ! End k
          !
          If(.Not.norm_to_improve) Then
             !
             !------------------------------------------------------------------
             !  Density matrices (computed only when norm_to_improve = False,
             !  i.e., when the particle number is conserved within 10^-5)
             !------------------------------------------------------------------
             kdib=kl-kaib; ka(ib,it)=kaib; kd(ib,it)=kdib
             k1=kaib+1; k2=kaib+kdib
             eqpe=0.0_pr
             Do n2=1,nd
                Do n1=n2,nd
                   s1=zero; s2=zero
                   If(k1.Le.k2) Then
                      Do k=k1,k2
                         ! temperature
                         fac1 = one; fac2 = zero
                         If(switch_on_temperature) Then
                            ii=KqpPo(k); fac1=one-f_T(ii); fac2=f_T(ii)
                         End If
                         nd1=KpwiPo(k)+n1; nd2=KpwiPo(k)+n2
                         s1=s1+VqpPo(nd1)*fac1*VqpPo(nd2)+UqpPo(nd1)*fac2*UqpPo(nd2)
                         s2=s2+UqpPo(nd1)*fac1*VqpPo(nd2)+VqpPo(nd1)*fac2*UqpPo(nd2) &
                              +VqpPo(nd2)*fac1*UqpPo(nd1)+UqpPo(nd2)*fac2*VqpPo(nd1)
                      End Do
                      s1=two*s1; s2=half*s2               ! two:due to m-projection, half:due to symmetrization
                      ! blocking
                      If(ibiblo.Eq.ib) Then
                         i=blok1k2d(it); id1=KpwiPo(i)+n1; id2=KpwiPo(i)+n2
                         s1=s1-VqpPo(id1)*VqpPo(id2)+UqpPo(id1)*UqpPo(id2)
                         s2=s2-half*(UqpPo(id1)*VqpPo(id2)+VqpPo(id1)*UqpPo(id2))
                      End If
                   End If
                   n12=n1+(n2-1)*nd; n21=n2+(n1-1)*nd
                   rk(n12,m)=s1; rk(n21,m)=s1              !  V V'
                   ak(n12,m)=-s2; ak(n21,m)=-s2            !- U V', ak=half*(pairing density)
                   hfbcan(n1,n2)=s1; allhfb(ib)%arr(n1,n2)=s1
                End Do !n1
             End Do !n2
             !------------------------------------------------------------------
             ! Canonical basis
             !------------------------------------------------------------------
             If(k1.Le.k2) Then
                Call Canonical(it,icanon,k2,k1,nd,i0,lc,ib,ibiblo,m,ibro)
                If(ierror_flag.Ne.0) Return
             End If
             lcanon(ib,it)=lc
             !
          End If
          !
       End Do !ib
       !
       If(kl.Eq.0) Then
          ierror_flag=ierror_flag+1
          ierror_info(ierror_flag)=' STOP: kl=zero, no states below pwi!!!'
          Return
       End If
       If(iparenti(it).Ne.0.And.ibiblo.Eq.0) Then
          ierror_flag=ierror_flag+1
          ierror_info(ierror_flag)='STOP: No blocking candidate found!!!'
          Return
       End If
       eqpmin(it)=emin; klmax(it)=kl; sumnz(it)=sumnz(it)-tz(it)
       !------------------------------------------------------------------
       ! Lambda search
       !------------------------------------------------------------------
       Call ALambda(al,it,kl)
       If(ierror_flag.Ne.0) Return
       If(keyblo(it).Eq.0) Then
          ala(it)=al
       Else
          ala(it)=ala(it)+0.50_pr*(al-ala(it))
       End If
       ! NB! 'alast' instead of 'al' at small pairing
       alast(it)=al
       If(Abs(ept(it)+frept(it)).Lt.0.0001_pr.And.(.Not.switch_on_temperature)) Then
          ntz=Int(tz(it)+0.1_pr); ntz=ntz/2
          Do k=1,kl
             drhfb(k)=erhfb(k)
          End Do
          Call ord(kl,drhfb)
          alast(it)=drhfb(ntz)  !last bound s.p. energy
       End If
       !------------------------------------------------------------------
       ! THO asymptotic decay
       !------------------------------------------------------------------
       ! density asymptotic decay \rho(r)->Exp(-ass(it)*r)
       ! ass(it)=2*Sqrt((E_min-\lambda)/((A-1)/A)*hbar**2/(2*m)))
       al2=zero
       If(kindhfb.Lt.0) Then
          al2=al+two*ala2(it)*(one-two*alnorm) ! al=al+two*ala2(it)
       End If
       if(it.eq.1) then
          al2=(emin-al2)/hb0n
       else
          al2=(emin-al2)/hb0p
       endif
       ! wrong asymptotic
       iasswrong(it)=0; If(al2.Le.zero) iasswrong(it)=1; ass(it)=two*Sqrt(Abs(al2))
       !
    End Do ! While(norm_to_improve)
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('hfbdiag',1)
    !
  End Subroutine hfbdiag
  !=======================================================================
  !> Determination of the Fermi energy based on the value of the particle
  !> number and the current density matrix
  !=======================================================================
  Subroutine ALambda(al,it,kl)
    Implicit None
    Integer(ipr) :: it,i,k,kl,icze,lit,ntz,iw
    Real(pr), Save :: fm7=1.0e-7_pr,fm10=1.0e-10_pr
    Real(pr) :: al,vh,xinf,xsup,esup,ez,dez,dfz,dvh,y,a,b,einf,absez,sn
    Real(pr) :: fT,dfT
    Real(pr), Pointer :: f_T(:)
    !-------------------------------------------------
    ! Fermi-Dirac occupations
    !-------------------------------------------------
    If(switch_on_temperature) Then
       If(it.Eq.1) Then
          f_T=>fn_T
       Else
          f_T=>fp_T
       End If
    End If
    !-------------------------------------------------
    ! Chemical potential without pairing
    !-------------------------------------------------
    If(CpV0(it-1).Eq.zero.and.trim(skyrme).ne.'D1S') Then
       ntz=Int(tz(it)+0.1_pr); ntz=ntz/2
       Do k=1,kl
          drhfb(k)=erhfb(k)
       End Do
       Call ord(kl,drhfb)
       If (ntz.Lt.kl) Then
          al=half*(drhfb(ntz)+drhfb(ntz+1))
       Else
          al=drhfb(ntz)+0.001_pr
       End If
       Return
    End If
    !-------------------------------------------------
    ! Chemical potential with pairing
    !-------------------------------------------------
    xinf=-1000.0_pr; xsup=1000.0_pr; esup=one; icze=0
    Do lit=1,500
       sn=zero;dez=zero;dfz=zero
       Do i=1,kl
          vh=zero; dvh=zero; fT=zero; dfT=zero
          y=erhfb(i)-al; a=y*y+drhfb(i)**2; b=Sqrt(a)
          !
          If(switch_on_temperature.And.temper.Gt.1.e-12_pr) Then
             fT =half*(one-Tanh(half*b/temper))
             dfT=y/b/temper*fT*(one-fT)
             f_T(i)=fT
          Else
             fT =zero
             dfT=zero
          End If
          !
          If(b.Gt.zero)  vh=half*(one-y/b)
          !
          !If(b.Lt.fm7.And.icze.Eq.1) vh=-einf/(esup-einf) !no pairing
          If(vh.Lt.1.e-12_pr) vh = zero
          If((vh-one).Gt.1.e-12_pr)  vh = one
          If(b.Gt.zero) dvh=half*drhfb(i)**2/(a*b)         ! D[ez,al](i)
          ! blocking
          If(i.Eq.blok1k2d(it)) Then
             vh=half; dvh=zero
          End If
          sn=sn+two*vh+two*(one-two*vh)*fT
          dez=dez+two*(one-two*fT)*dvh
          dfz=dfz+two*(one-two*vh)*dfT   ! D[ez,al]
       End Do
       ez=sn-tz(it); absez=Abs(ez)/tz(it)
       dez=dez+dfz
       !-------------------------------------------------
       ! Correcting bounds
       !-------------------------------------------------
       If(ez.Lt.zero) Then
          xinf=Max(xinf,al); einf=ez
       Else
          xsup=Min(xsup,al); esup=ez
       End If
       If(lit.Eq.1) Then
          If(absez.Le.0.10_pr) Then
             al=al-ez
          Else
             al=al-0.10_pr*Sign(one,ez)
          End If
       Else
          al=al-ez/(dez+1.e-20_pr)                      ! Newton method
       End If
       If(xsup-xinf.Lt.fm7) icze=1                      ! low/upp close
       If(al.Lt.xinf.Or.al.Gt.xsup) al=half*(xinf+xsup) ! mean upp/low
       If(absez.Le.fm10) Return
    End Do
    !-------------------------------------------------
    ! Low accuracy warning
    !-------------------------------------------------
    Do iw=lout,lfile
      Write(iw,'(a,2(e12.5,2x),a,2(2x,f8.4),a,i2)') ' Low accuracy=',sn,ez,' for N,Z=',tz,' it=',it
    End Do
  End Subroutine Alambda
  !=========================================================================================
  !> Determination of the canonical basis by diagonalization of the density matrix
  !=========================================================================================
  Subroutine Canonical(it,icanon,k2,k1,nd,i0,lc,ib,ibiblo,m,ibroib)
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: it,i0,icanon,ibiblo,i,iw,k,kk,lc,ib,nd,n1,n2,m,nd1,k1,k2,n12,ier
    Integer(ipr) :: nhhph,nhhpp,ibro,ibroib,il,iu,NUMFOU
    Integer(ipr), Allocatable :: ifail(:)
    Real(pr) :: s1,s2,vx,h1,d1,h2,d2,ddn1,ddn2,vl,vu,ABSTOL
    Real(pr), Allocatable :: hh(:,:),de(:,:),ewavef(:,:)
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
    Integer(ipr), Allocatable :: ISUPPZ(:)
    Real(pr), Allocatable :: eigenv(:),eigenf(:,:)
    Real(pr), External :: DLAMCH
    !
    If(it.Eq.1) Then
       EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN
    Else
       EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP
    End If
    !
#if(SWITCH_ESSL==0)
    If(Allocated(alwork)) Deallocate(alwork,lwork)
    ialwork=26*nd; ilwork=10*nd
    Allocate(ALWORK(ialwork),LWORK(ilwork)); ALWORK = 0.0_pr; LWORK = 1
#else
    If(Allocated(alwork)) Deallocate(alwork,lwork)
    ialwork=8*nd; ilwork=5*nd;
    Allocate(ALWORK(ialwork),LWORK(ilwork)); ALWORK = 0.0_pr; LWORK = 0
#endif
    !
    ABSTOL=2.0_pr*DLAMCH('S')
    !
    If(Abs(ept(it)).Lt.0.0001_pr.And.(.Not.switch_on_temperature)) Then
       !------------------------------------------------------
       ! No pairing => just taking the HF states
       !------------------------------------------------------
       Do k=1,nd
          kk=k1+k-1; lc=lc+1                              ! total number of the canonical states
          ddc(1:nd,lc,it)=zero; vk(lc,it)=zero            ! zeros: nd could be larger then k2-k1+1
          If(kk.Gt.k2) Cycle
          vx=zero
          Do i=1,nd
             h1=VqpPo(KpwiPo(kk)+i)**2; vx=vx+h1
          End Do
          If (vx.Le.zero) vx=zero                         ! roundoff errors
          If (vx.Ge.one ) vx=one
          Do i=1,nd
             If(vx.Ge.half) Then
                ddc(i,lc,it)=VqpPo(KpwiPo(kk)+i)          ! (ph) s.p. orbitals in conf.space
             Else
                ddc(i,lc,it)=UqpPo(KpwiPo(kk)+i)          ! (ph) s.p. orbitals in conf.space
             End If
          End Do
          Dispersion(it)=Dispersion(it)+four*vx*(one-vx)  ! internal P/N Dispersion
          If(Abs(vx-half).Le.v2min(it)) Then
             v2min(it)=Abs(vx-half); v2minv(it)=vx        ! divergent condition
             lcc=lc
          End If
          vk(lc,it)=vx                                    ! (ph) s.p. occupations v^2
          !------------------------------------------------------
          ! RESU only
          !------------------------------------------------------
          If(icanon.Ne.0) Then
             ek(lc,it)=EqpPo(KqpPo(kk))*(one-two*vx)+ala(it)      ! (ph) s.p. energies
             dk(lc,it)=zero                                       ! (ph) s.p. deltas
          End If
       End Do !k
    Else
       !------------------------------------------------------
       ! Pairing => calculate canonical basis
       !------------------------------------------------------
#if(SWITCH_ESSL==0)
       VL=0.0_pr; VU=0.0_pr; IL=1; IU=nd; NUMFOU=0
       Allocate(ISUPPZ(2*nd))
       Allocate(eigenv(nd)); eigenv(1:nd)=0.0_pr
       Allocate(eigenf(nd,nd)); eigenf(1:nd,1:nd)=0.0_pr
       ier=0; Call DSYEVR('V','A','L',nd,hfbcan,ndx,VL,VU,IL,IU,ABSTOL,NUMFOU, &
                          eigenv,eigenf,nd,ISUPPZ,ALWORK,ialwork,LWORK,ilwork,ier)
       evvkcan(1:nd) = eigenv(1:nd)
       hfbcan(1:nd,1:nd) = eigenf(1:nd,1:nd)
       Deallocate(eigenv,eigenf,ISUPPZ)
       !ier=0; Call DSYEVD('V','L',nd,hfbcan,ndx,evvkcan,ALWORK,ialwork,LWORK,ilwork,ier)
#else
       vl=0.0_pr; vu=0.0_pr; il=1; iu=1; m=0; abstol=0.0_pr
       Allocate(ifail(nd),ewavef(ndx,ndx))
       ier=0; Call DSYEVX('V','A','L',nd,hfbcan,ndx,vl,vu,il,iu,abstol,m, &
                          evvkcan,ewavef,ndx,ALWORK,ialwork,LWORK,ifail,ier)
       hfbcan(1:ndx,1:ndx) = ewavef(1:ndx,1:ndx)
       Deallocate(ifail,ewavef)
#endif
       ! bug in LAPACK
       If(ier.Gt.0) Then
          Do iw=lout,lfile
             Write(iw,*) 'FATAL ERROR CONDITION IN CANONICAL DSYEVR, ier=',ier,'(RECOVERED)'
          End Do
          Do n2=1,nd
             Do n1=n2,nd
                vx=allhfb(ib)%arr(n1,n2)
                hfbcan(n2,n1)=vx; hfbcan(n1,n2)=vx
             End Do
          End Do
          Call sdiag(ndx,nd,hfbcan,evvkcan,hfbcan,zhfb,+1)
       End If
       !------------------------------------------------------
       ! Eigenvalues and wavefunctions
       !------------------------------------------------------
       Do k=1,nd
          lc=lc+1                                        ! total number of the canonical states
          Do i=1,nd
             ddc(i,lc,it)=hfbcan(i,k)                    ! (ph) canon orbitals in conf.space
          End Do
          vx=evvkcan(k)*half
          If (vx.Le.zero) vx=zero                        ! roundoff errors
          If (vx.Ge.one ) vx=one
          ! blocking
          If(ibiblo.Eq.ib.And.vx.Gt.0.49_pr.And.vx.Le.0.51_pr) blocanon(it)=lc
          Dispersion(it)=Dispersion(it)+four*vx*(one-vx) ! internal P/N Dispersion
          If(Abs(vx-half).Le.v2min(it)) Then
             v2min(it)=Abs(vx-half); v2minv(it)=vx       ! divergent condition
             lcc=lc
          End If
          vk(lc,it)=vx                                   ! (ph) canon occupations v^2
          !------------------------------------------------------
          ! RESU only
          !------------------------------------------------------
          If(icanon.Ne.0) Then
             ! canonical energies and deltas (no physical meaning in PNP)
             nhhph=(it-1)*nhhdim; nhhpp=(it+1)*nhhdim
             Allocate(hh(nd,nd),de(nd,nd))
             ibro=ibroib
             Do n1=1,nd
                Do n2=1,n1
                   ibro=ibro+1
                   vx=brin(nhhph+ibro); hh(n2,n1)=vx; hh(n1,n2)=vx
                   vx=brin(nhhpp+ibro); de(n2,n1)=vx; de(n1,n2)=vx
                End Do
             End Do
             h1=zero; d1=zero
             Do n2=1,nd
                h2=zero; d2=zero
                Do n1=1,nd
                   ddn1=hfbcan(n1,k)
                   h2=h2+ddn1*hh(n1,n2)
                   d2=d2+ddn1*de(n1,n2)
                End Do
                ddn2=hfbcan(n2,k)
                h1=h1+h2*ddn2
                d1=d1+d2*ddn2
             End Do
             ek(lc,it)=h1                              ! (ph) canon s.p. energies
             dk(lc,it)=d1                              ! (ph) canon s.p. deltas
             Deallocate(hh,de)
          End If
          !
       End Do !k
    End If
    !------------------------------------------------------
    ! RESU only
    !------------------------------------------------------
    If(icanon.Ne.0) Then
       !------------------------------------------------------
       ! Find maximal HO components of all qp states
       !------------------------------------------------------
       Do k=k1,k2
          s1=zero
          Do n1=1,nd
             nd1=nd+n1
             s2=Max(s1,Abs(VqpPo(KpwiPo(k)+n1)),Abs(UqpPo(KpwiPo(k)+n1)))
             If(s2.Gt.s1) Then
                s1=s2
                vkmax(k,it)=s1                       ! maximal overlap
                numax(k,it)=n1+i0                    ! its number in k[k1,k2] numbering
             End If
          End Do
       End Do
       !------------------------------------------------------
       ! Searching for possible blocking candidates
       !------------------------------------------------------
       If(iparenti(it).Eq.0) Then
          n1=0
          Do k=k1,k2
             n1=n1+1
              ! Search within |(1-2*N)*Eqpe| lover than 'pwiblo'
              ! The levels number n1 is 1,2,3,... for the given block ([123] numbering)
             If(Abs(EqpPo(KqpPo(k))-eqpmin(it)).Le.pwiblo) Then
                blomax(it)=blomax(it)+1                         ! blocked state #, maximel # of block candidates
                If(blomax(it).Gt.bloall) Then
                   ierror_flag=ierror_flag+1
                   ierror_info(ierror_flag)='Too many blocking candidates! Increase bloall and run again'
                   Return
                End If
                bloblo(blomax(it),it)=ib                        ! block where to block
                blo123(blomax(it),it)=n1                        ! state # [123] numbering
                blok1k2(blomax(it),it)=k                        ! state # k[k1,k2] numbering
                bloqpdif(blomax(it),it)=Abs(EqpPo(KqpPo(k))-eqpmin(it))
             End If
          End Do
       End If
    End If
    !
  End Subroutine Canonical
  !=======================================================================
  !> Displays results at convergence: single particle energies, densities, fields
  !=======================================================================
  Subroutine resu(irecord)
    Implicit None
    Integer(ipr) :: it,iw,ib,im,m,nd,k,k0,k1,k2,j,n,imax,nhfb,irecord,iexit
    Real(pr)     :: sum,eqpe,pn,ela,enb,ek0,vk0,ekk,delb,ovmax,s,uuvv,  &
         dk0,skk,summ(4),vvs,vvc,enjacek
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
    !
    iexit=0
    !
    !--------------------------------------------
    ! last HFB run for full canon.calculations
    !--------------------------------------------
    Do it=itmin,itmax
       Call hfbdiag(it,1)       ! hfb with maximal canonical
       If(ierror_flag.Ne.0) Return
    End Do
    !PAV
    Call expect(.False.)        ! expectation values
    !
    !Call coulom_test
    !
    If(ierror_flag.Ne.0) Return
    Call field                  ! new fields
    If(ierror_flag.Ne.0) Return
    Call gamdel(.false.,.true.)        ! hf-matrix
    If(ierror_flag.Ne.0) Return
    ! inout(2): HFB matrices if regular nucleus is even-even (not applicable for non-integer
    ! fission fragments)
    If( (write_hel.And.fragment_properties) .Or.&
        (write_hel.And. .Not.fragment_properties .And. .Not.odd_noBlock .And. &
        (npr(1).Eq.2*(npr(1)/2).And.npr(2).Eq.2*(npr(2)/2))) .Or. &
        (write_hel.And. .Not.fragment_properties .And. odd_noBlock) ) Then
       Call inout(2,iexit)
       If(ierror_flag.Ne.0) Return
    End If
    !--------------------------------------------
    ! Printing densities and fields
    !--------------------------------------------
    ! Call printLST(ro(1,1),ro(1,2))      !need fix
    ! Call printLST(tau(1,1),tau(1,2))
    ! Call printLST(dro(1,1),dro(1,2))
    ! Call printLST(dj (1,1),dj (1,2))
    ! printing of fields
    ! Call printLST(vhb(1,1),vhb(1,2))
    ! Call printLST(v (1,1),v (1,2))
    ! Call printLST(vs(1,1),vs(1,2))
    !--------------------------------------------
    ! Printing quasiparticle states
    !--------------------------------------------
    Do it=itmin,itmax
       If(it.Eq.1) Then
          EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN
       Else
          EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP
       End If
       !
       If(Print_Screen) Then
          iw=lfile
          Write(iw,200) tit(it)
          Write(iw,*) ' eqp(k) -> q.p. energy '
          Write(iw,*) ' e(k)   -> referent s.p. energy '
          Write(iw,*) ' p(k)   -> occ.probability '
          Write(iw,*) ' del(k) -> referent s.p. gap '
          Write(iw,*) ' fermi energy alast=',alast(it)
          Write(iw,'(a,a)')  &
               '  #k  block#    eqp(k)     e(k)       (1-2N)E      decay        p(k)',  &
               '        del(k)    overl      labels'
200     Format(//,' #quasiparticle energies ',a,/,1x,32('-'))
       End If
       sum=zero
       Do ib=1,nb
          nd=id(ib); im=ia(ib); m=ib+(it-1)*nbx; nhfb=nd+nd
          k1=ka(ib,it)+1
          k2=ka(ib,it)+kd(ib,it)
          If(k1.Le.k2) Then
             Do k=k1,k2                                 ! print active states only
                pn=uk(k,it)                             ! qp probabilities
                j=k
                If(pn.Gt.-1.d-14) Then                   ! print If signIficant pn
                   ! main oscillator component
                   ovmax=vkmax(k,it)                    ! maximal overlap
                   imax=numax(k,it)                     ! its number
                   ! printing
                   eqpe=EqpPo(KqpPo(k))                 ! qp energies
                   if(it.eq.1) then
                      skk=two*Sqrt(Abs(eqpe-ala(it))/hb0n)  ! qp decay
                   else
                      skk=two*Sqrt(Abs(eqpe-ala(it))/hb0p)  ! qp decay
                   endif
                   ela=eqpe*(one-two*pn)
                   enb=ela+ala(it)                      ! ref. s.p. energies
                   delb=Sqrt(Abs(eqpe**2-ela**2))       ! ref. s.p. delta
                   sum=sum+two*pn                       ! particle number
                   If(Print_Screen) Then
                      iw=lfile
                      Write(iw,201) k,ib,eqpe,enb,(one-two*pn)*eqpe,skk,pn,delb,ovmax,tb(imax)
201                 Format(i4,2x,i3,1x,f12.6,f12.6,f12.6,f12.6,2x,f12.8,  &
                           2(2x,f7.4),' ',a13)
                   End If
                End If
             End Do
          End If
       End Do !ib
       !--------------------------------------------
       ! Printing canonical single particle states
       !--------------------------------------------
       If(Print_Screen) Then
          iw=lfile
          Write(iw,'(a,i4,a,i4)')  &
               '#all active are ',j,' q.p. states out of ',nt
          Write(iw,'(a,f6.1)') '#since the cut off is pwi=',pwi
          Write(iw,'(3a,f6.1)')'#check: number of ',tit(it),'=',sum
          Write(iw,100) tit(it)
          Write(iw,*) ' labels -> {2*omega}{parity}[nn=nz+2*nr+nl,nz,nl]'
          Write(iw,*) ' cqpe   -> canonical q.p. energies'
          Write(iw,*) ' ce     -> canonical s.p. energies'
          Write(iw,*) ' fermi energy=',alast(it)
          Write(iw,*) ' average cdelt=',del(it)
          Write(iw,'(a,a)')'  k0      ceqp        ce         v*v',  &
               '       u*v        cdel     overl      labels'
100     Format(//,' #canonical s.p. energies ',a,/,1x,33('-'),//)
       End If
       k0=0
       summ=zero; enjacek=zero
       Do ib=1,nb
          nd=id(ib); im=ia(ib)
          k1=ka(ib,it)+1; k2=ka(ib,it)+kd(ib,it)
          If(k1.Le.k2) Then
             Do k=1,nd
                k0=k0+1
                ! for Lipkin Nogami
                vvs=two*Sqrt(vk(k0,it))*Sqrt(one-vk(k0,it))    !2vu
                vvc=two*vk(k0,it)-one                          !2v^2-1
                summ(1)=summ(1)+vvs**2
                summ(2)=summ(2)+vvs**2*vvc
                summ(3)=summ(3)+vvs**4
                summ(4)=summ(4)+(vvs*vvc)**2
                ! search for main oscillator component
                ovmax=zero
                Do n=1,nd
                   s=Abs(ddc(n,k0,it))                         !canon orbitals in conf.space
                   If (s.Ge.ovmax) Then
                      ovmax=s; imax=n
                   End If
                End Do
                ! printing
                ek0=ek(k0,it)                                  !canon s.p. energies
                enjacek=enjacek+ek0*vk(k0,it)
                If(ek0.Lt.pwi) Then                            !print up to 'pwi'
                   vk0=vk(k0,it)                               !canon occupations v^2
                   If(vk0.Gt.-1.d-4) Then                       !print If signIficant v^2
                      dk0=-dk(k0,it)                           !canon s.p. deltas
                      ekk=Sqrt((ek0-ala(it))**2+dk(k0,it)**2)  !resulting cqpe
                      uuvv=Sqrt(Abs(vk0*(one-vk0)))            !resulting u*v
                      If(Print_Screen) Then
                         iw=lfile
                         Write(iw,101) k0,ekk+ala(it),ek0,vk0,uuvv,dk0,ovmax,tb(im+imax)
101                    Format(i4,2f12.6,2(1x,f12.8),2(2x,f7.4),' ',a13)
                      End If
                   End If
                End If
             End Do !k0
          End If
       End Do !ib
       !--------------------------------------------
       ! Lipkin-Nogami
       !--------------------------------------------
       ssln(1,it)=summ(1)
       ssln(2,it)=summ(2)
       ssln(3,it)=summ(4)*summ(1)-summ(2)**2+summ(1)**3/4.0_pr-half*summ(3)*summ(1)
       If(Print_Screen) Then
          iw=lfile
          Write(iw,*) ' Sum canonical e_v*V^2_k=',two*enjacek
       End If
    End Do !it
    !--------------------------------------------
    ! To thoout.dat, thodef.dat and hodef.dat
    !--------------------------------------------
    If(irecord.Ne.0) Then
       iappend=1
       Call expect(.True.)    !print  & record HFB+PAV results
       If(ierror_flag.Ne.0) Return
       iappend=0
    Else
       Call expect(.True.)    !print HFB+PAV results
       If(ierror_flag.Ne.0) Return
    End If
    !
  End Subroutine resu
  !=======================================================================
  !> Initializes all NAMELISTS
  !=======================================================================
  Subroutine initialize_HFBTHO_NAMELIST
    Implicit None
    ! HFBTHO_GENERAL
    number_of_shells    = 10
    oscillator_length   = -one
    basis_deformation   = zero
    proton_number       = 24
    neutron_number      = 26
    type_of_calculation = 1
    ! HFBTHO_INITIAL
    beta2_deformation = zero
    beta2_deformation = zero
    ! HFBTHO_ITERATIONS
    number_iterations = 100
    accuracy          = 1.D-5
    restart_file      = -1
    ! HFBTHO_FUNCTIONAL
    functional          = 'SLY4'
    add_initial_pairing = .False.
    type_of_coulomb     = 2
    include_3N_force    = .False.
    ! HFBTHO_PAIRING
    user_pairing    = .False.
    vpair_n         = -300.0_pr
    vpair_p         = -300.0_pr
    pairing_cutoff  =   60.0_pr
    pairing_feature =    0.5_pr
    ! HFBTHO_CONSTRAINTS
    lambda_values       = (/ 0, 0, 0, 0, 0, 0, 0, 0 /)
    lambda_active       = (/ 0, 0, 0, 0, 0, 0, 0, 0 /)
    expectation_values  = (/ 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr /)
    ! HFBTHO_BLOCKING
    proton_blocking  = (/ 0, 0, 0, 0, 0 /)
    neutron_blocking = (/ 0, 0, 0, 0, 0 /)
    ! HFBTHO_PROJECTION
    switch_to_THO    = 0
    projection_is_on = 0
    gauge_points     = 1
    delta_Z          = 0
    delta_N          = 0
    ! HFBTHO_TEMPERATURE
    set_temperature = .False.
    temperature     = zero
    ! HFBTHO_FEATURES
    collective_inertia     = .False.
    fission_fragments      = .False.
    pairing_regularization = .False.
    automatic_basis        = .False.
    ! HFBTHO_TDDFT
    filter                 = .False.
    fragment_properties    = .False.
    real_Z                 = 24.0_pr
    real_N                 = 26.0_pr
    ! HFBTHO_NECK
    set_neck_constrain     = .False.
    neck_value             = 0.5_pr
    ! HFBTHO_DEBUG
    number_Gauss        =  40
    number_Laguerre     =  40
    number_Legendre     =  80
    compatibility_HFODD = .False.
    number_states       = 500
    force_parity        = .True.
    print_time          = 0
    !
  End Subroutine initialize_HFBTHO_NAMELIST
  !=======================================================================
  !> Read all NAMELISTS from hfbtho_NAMELIST.dat
  !=======================================================================
  Subroutine read_HFBTHO_NAMELIST
    Implicit None
    Integer(ipr) :: ios,lnamelist=16
    Open(lnamelist,file='hfbtho_NAMELIST.dat',DELIM='APOSTROPHE') ! 'QUOTE'
    !
    ierror_flag = 0
    !
    ! General input data
    Read(UNIT=lnamelist,NML=HFBTHO_GENERAL,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_GENERAL read'
       Return
    End If
    !
    ! Deformations of the initial WS solution
    Read(UNIT=lnamelist,NML=HFBTHO_INITIAL,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_INITIAL read'
       Return
    End If
    !
    ! Iterations
    Read(UNIT=lnamelist,NML=HFBTHO_ITERATIONS,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_ITERATIONS read'
       Return
    End If
    !
    ! Type of functional
    Read(UNIT=lnamelist,NML=HFBTHO_FUNCTIONAL,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_FUNCTIONAL read'
       Return
    End If
    !
    ! Characteristics of pairing
    Read(UNIT=lnamelist,NML=HFBTHO_PAIRING,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_PAIRING read'
       Return
    End If
    !
    ! Constraints
    Read(UNIT=lnamelist,NML=HFBTHO_CONSTRAINTS,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_CONSTRAINTS read'
       Return
    End If
    !
    ! Blocking
    Read(UNIT=lnamelist,NML=HFBTHO_BLOCKING,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_BLOCKING read'
       Return
    End If
    !
    ! Particle number projection
    Read(UNIT=lnamelist,NML=HFBTHO_PROJECTION,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_PROJECTION read'
       Return
    End If
    !
    ! Finite temperature
    Read(UNIT=lnamelist,NML=HFBTHO_TEMPERATURE,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_TEMPERATURE read'
       Return
    End If
    !
    ! Various features of the calculation
    Read(UNIT=lnamelist,NML=HFBTHO_FEATURES,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_FEATURES read'
       Return
    End If
    !
    ! Interface with TDDFT codes
    Read(UNIT=lnamelist,NML=HFBTHO_TDDFT,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_TDDFT read'
       Return
    End If
    !
    ! Constraint on the neck
    Read(UNIT=lnamelist,NML=HFBTHO_NECK,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_NECK read'
       Return
    End If
    !
    ! Debug
    Read(UNIT=lnamelist,NML=HFBTHO_DEBUG,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_DEBUG read'
       Return
    End If
    !
    Close(lnamelist)
    !
  End Subroutine read_HFBTHO_NAMELIST
  !=======================================================================
  !> Performs some basis consistency checks of input data
  !=======================================================================
  Subroutine check_consistency
    Implicit None
    Integer(ipr) :: counter, i
    Real(pr) :: A, preset_inin(3)
    !Character(30), Dimension(:) :: preset_forces(41)
    Character(30), Dimension(:) :: preset_forces(42)   !EOedit for SV-min
    !
    If((abs(n00_INI).Lt.1).Or.(abs(n00_INI).GT.50)) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_of_shells = ",i6," out-of-bounds: [1,50]")') &
             n00_INI
       Return
    End If
    !
    If((npr_INI(1).Lt.1).Or.(npr_INI(2).Lt.1)) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("Z = ",i6," N = ",i6," out-of-bounds: (Z,N)>1")') &
             npr_INI(2),npr_INI(1)
       Return
    End If
    !
    If(Abs(kindhfb_INI).Ne.1) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("type_of_calculation = ",i6," unrecognized: (-1,1)")') &
             kindhfb_INI
       Return
    End If
    !
    If(epsi_INI.Lt.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("accuracy = ",e24.12," out-of-bounds: >0")') &
             epsi_INI
       Return
    End If
    !
    preset_inin( 1) = 1
    preset_inin( 2) = 2
    preset_inin( 3) = 3
    !
    counter=0
    Do i=1, 3
       If(Abs(inin_INI).Eq.preset_inin(i)) Then
          counter=1
          Exit
       End If
    End Do
    !
    If(counter.Eq.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("restart_file = ",i6," unrecognized: see list in publi")') &
             inin_INI
       Return
    End If
    !
    preset_forces( 1) = 'SIII'
    preset_forces( 2) = 'SKM*'
    preset_forces( 3) = 'SKP'
    preset_forces( 4) = 'SLY4'
    preset_forces( 5) = 'SLY5'
    preset_forces( 6) = 'SLY6'
    preset_forces( 7) = 'SLY7'
    preset_forces( 8) = 'SKI3'
    preset_forces( 9) = 'SKO'
    preset_forces(10) = 'SKX'
    preset_forces(11) = 'UNE0'
    preset_forces(12) = 'UNE1'
    preset_forces(13) = 'UNE2'
    preset_forces(14) = 'N0LO'
    preset_forces(15) = 'N1LO'
    preset_forces(16) = 'N2LO'
    preset_forces(17) = 'FITS'
    preset_forces(18) = 'D1'
    preset_forces(19) = 'D1p'
    preset_forces(20) = 'D1S'
    preset_forces(21) = 'D1N'
    preset_forces(22) = 'T0X0'
    preset_forces(23) = 'DME_LO'
    preset_forces(24) = 'DME_NLO'
    preset_forces(25) = 'DME_N2LO'
    preset_forces(26) = 'DME_NLOD'
    preset_forces(27) = 'DME_N2LOD'
    preset_forces(28) = 'REG_LO'
    preset_forces(29) = 'REG_NLO'
    preset_forces(30) = 'REG_N2LO'
    preset_forces(31) = 'REG_NLOD'
    preset_forces(32) = 'REG_N2LOD'
    preset_forces(33) = 'NEDF'
    preset_forces(34) = 'SeaLL1'
    preset_forces(35) = 'NEDF1'
    preset_forces(36) = 'NEDF2'
    preset_forces(37) = 'NEDF3'
    preset_forces(38) = 'NEDF4'
    preset_forces(39) = 'NEDF5'
    preset_forces(40) = 'HFB1'
    preset_forces(41) = 'SKM*mod'
    preset_forces(42) = 'SV-min'  !EOedit for SV-min
    !
    counter=0
    !Do i=1, 41
    Do i=1, 42   !EOedit for SV-min
       If(Trim(skyrme_INI).Eq.Trim(preset_forces(i))) Then
          counter=1
          Exit
       End If
    End Do
    ! Functional must be in preset list
    If(counter.Eq.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("functional = ",a30," unrecognized: see list in publi")') &
             skyrme_INI
       Return
    End If
    ! Pairing cut-off must be positive
    If(pwi_INI.Lt.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("pairing_cutoff = ",i4," out-of-bounds: >=0")') &
             pwi_INI
       Return
    End If
    ! Pairing cut-off must be positive
    If(cpv1_INI.Lt.0.0.Or.cpv1_INI.Gt.1.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("pairing_feature = ",i4," out-of-bounds: [0.0,1.0]")') &
             cpv1_INI
       Return
    End If
    ! Options for Coulomb: -3, -2, -1, 0, 1, 2
    If(icou_INI.Lt.-4.Or.icou_INI.Gt.2) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("type_of_coulomb = ",i4," unrecognized: (-4,-3,-2,-1,0,1,2)")') &
             icou_INI
       Return
    End If
    ! Choices of basis (HO or THO): -1, 0, 1
    If(Abs(iLST_INI).Gt.1) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("switch_to_THO = ",i4," unrecognized: (-1,0,1)")') &
             iLST_INI
       Return
    End If
    ! At least one gauge point if projection is required
    If(keypj_INI.Le.0.And.iproj_INI.Ne.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("gauge_points = ",i4," out-of-bounds: >=0")') &
             keypj_INI
       Return
    End If
    ! Number of protons must be greater than 0 for projection
    If((npr_INI(1)+npr1pj_INI).Lt.1.And.iproj_INI.Ne.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("delta_N = ",i4," out-of-bounds: N+dN>=1")') &
             npr1pj_INI
       Return
    End If
    ! Number of neutrons must be greater than 0 for projection
    If((npr_INI(2)+npr2pj_INI).Lt.1.And.iproj_INI.Ne.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("delta_Z = ",i4," out-of-bounds: Z+dZ>=1")') &
             npr2pj_INI
       Return
    End If
    ! Temperature must be positive
    If(temper.Lt.zero) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("temperature = ",i4," out-of-bounds: T>=0")') &
             temper
       Return
    End If
    ! Number of Gauss-Laguerre integration points between 0 and 100
    If(ngh_INI.Lt.1.Or.ngh_INI.Gt.100) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_Gauss = ",i4," out-of-bounds: [1,100]")') &
             ngh_INI
       Return
    End If
    ! Number of Gauss-Hermite integration points between 0 and 100
    If(ngl_INI.Lt.1.Or.ngl_INI.Gt.100) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_Laguerre = ",i4," out-of-bounds: [1,100]")') &
             ngl_INI
       Return
    End If
    ! Number of Gauss-Legendre integration points lower than 100
    If(nleg_INI.Gt.100) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_Legendre = ",i4," out-of-bounds: [-infty,100]")') &
             nleg_INI
       Return
    End If
    ! Number of Gauss-Legendre integration points between 1 and 100 for PNP
    If((nleg_INI.Lt.1.Or.nleg_INI.Gt.100).And.iproj_INI.Ne.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_Legendre = ",i4," out-of-bounds: [1,100]")') &
             nleg_INI
       Return
    End If
    ! Number of basis states must be greater than 0
    If(nstate_INI.Lt.1.And.basis_HFODD_INI) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_states = ",i4," out-of-bounds: >0")') &
             nstate_INI
       Return
    End If
    !
  End Subroutine check_consistency
  !=======================================================================
  !> Initializes the solver based on user-defined input data. In particular
  !> defines the parameters of the functional
  !=======================================================================
  Subroutine initialize_HFBTHO_SOLVER
    !---------------------------------------------------------------------
    ! default parameters
    !---------------------------------------------------------------------
    Implicit None
    Real(pr) :: A
    !------------------------------------
    ! tapes
    !------------------------------------
    lwin=41; lwou=42;  lwel=52; lres=57; lin=3
    !------------------------------------
    ! From Namelist or default values
    !------------------------------------
    nstate                = nstate_INI
    epsi                  = epsi_INI                  ! stop criteria
    Add_Pairing           = Add_Pairing_INI           ! add pairing starting from file
    icou                  = icou_INI                  ! Coulomb flag, see routine
    DO_FITT               = DO_FITT_INI               ! calculates quantities for reg.optimization
    IDEBUG                = IDEBUG_INI                ! debug
    Parity                = Parity_INI                ! reflection symmetry
    Print_HFBTHO_Namelist = Print_HFBTHO_Namelist_INI ! Print Namelist
    !---------------------------------------------------------------------
    ! Pairing set by user
    !---------------------------------------------------------------------
    rho_c=0.160_pr
    If(set_pairing) Then
       CpV0(0)=V0n_INI
       CpV0(1)=V0p_INI
       CpV1(0)=cpv1_INI
       CpV1(1)=cpv1_INI
       pwi=pwi_INI
    Else
       pwi=60.0_pr
    End If
    !------------------------------------
    ! output control
    !------------------------------------
    If(n00_INI.Gt.0) Then
       Print_Screen=.True.
#if(DO_MASSTABLE==1 || DRIP_LINES==1 || DO_PES==1)
       if(mpi_size.eq.1) then
          lfile=lout+1                              ! condensed output to screen & full thoout.dat file
       else
          lout = lfile                              ! full output to thoout.dat file only
       endif
#else
       lfile=lout+1                                 ! condensed output to screen & full thoout.dat file
#endif
    Else
       Print_Screen=.False.
#if(DO_MASSTABLE==1 || DRIP_LINES==1 || DO_PES==1)
       lout = lfile                                 ! condensed output to thoout.dat file only
#else
#ifdef SUPRESS_OUTPUT
       lout = lfile
#else
       lfile=lout-1                                 ! no output to screen & thoout.dat
#endif
#endif
    End If
    !------------------------------------
    ! Pi
    !------------------------------------
    PI=four*Atan(one)
    !------------------------------------
    ! blocking
    !------------------------------------
    bloblo=0; blo123=0; blok1k2=0;  keyblo=0
    blomax=0; nkblo=0;  iparenti=0; irestart=0
    blocanon=0;         eqpmin=zero
    !------------------------------------
    ! buffers
    !------------------------------------
    eres=zero;  eresu=zero;  eresl=zero;
    eresj=zero; eresbl=zero; ereslbl=' 00[00,00,00]'
    !------------------------------------
    ! def parameters
    !------------------------------------
    ffdef3=Sqrt(five/(four*pi))/two
    ffdef4=Sqrt(117.0_pr)/(four*pi) !MCcomment 117=6^2+9^2, don't know the meaning of this.
    ffdef5=Sqrt(nine/(four*pi))/eight
    ffdef6=Sqrt(five*pi)/three
    ffdef7=Sqrt(pi)/four
    !------------------------------------
    ! former linear mixing
    !------------------------------------
    xmix0=0.3_pr             ! lowest mixing parameter  (redefined later)
    xmix =0.3_pr             ! initial mixing parameter (changes every iteration)
    xmax =1.0_pr             ! mario
    !------------------------------------
    ! misc (redefined later)
    !------------------------------------
    rehfbcan=0.0_pr; depnp=0.0_pr; ala2=0.00_pr
    ept=-2.0_pr; del=1.0_pr; ala=-7.0_pr
    ala1(1)=-14.6851_pr; ala1(2)=-3.7522_pr; si=1.0_pr
    iqrpa=0; icacou=0;  icacoupj=0; icahartree=0; iasswrong=0
    iError_in_HO=0;  iError_in_THO=0
    ECMHFB=0.0_pr; ECMPAV=0.0_pr
    If(use_full_cm_cor) Then
       !A = npr_INI(1) + npr_INI(2)
       A = tz(1) + tz(2)
       facECM = A/(A-1.0_pr)
    End If
    entropy(:)=zero
    !------------------------------------
    ! Saxon-Woods: von koepf und ring, z.phys. (1991)
    !------------------------------------
    v0ws=-71.28_pr; akv=0.4616_pr; r0v=1.2334_pr; av=0.6150_pr
    vso=11.1175_pr; rso=1.1443_pr; aso=0.6476_pr
    !------------------------------------
    ! fixed text
    !------------------------------------
    tp(1)='+'; tp(2)='-'; tis(1)='n'; tis(2)='p';
    tit(1)='neutrons'; tit(2)='protons '
    tl(0)='s'; tl(1)='p'; tl(2)='d'; tl(3)='f'; tl(4)='g'
    tl(5)='h'; tl(6)='i'; tl(7)='j'; tl(8)='k'; tl(9)='l'
    tl(10)='m'; tl(11)='n'; tl(12)='o'; tl(13)='p'; tl(14)='q'
    tl(15)='r'; tl(16)='s'; tl(17)='t'; tl(18)='u'; tl(19)='v'; tl(20)='w'
    !------------------------------------
    ! fixed parity sign
    !------------------------------------
    tpar(1)=+1; tpar(2)=-1;
    !------------------------------------
    ! physical constants
    !------------------------------------
    amn=938.90590_pr
    amu=931.4940130_pr; r0=1.20_pr
    alphi=137.036020_pr; hqc=197.328910_pr
    !------------------------------------
    ! Coulomb
    !------------------------------------
    !chargee2=hqc/alphi
    !chargee2=1.43997841_pr
    !-----------------------------------
    ! set the loops over particle types
    !-----------------------------------
    itmin=1 ; itmax = 2;
    If(npr_INI(1).Eq.0) itmin = 2
    If(npr_INI(2).Eq.0) itmax = 1
    !-----------------------------------
    ! error flag and info
    !-----------------------------------
    ierror_flag=0
    ierror_info(ierror_flag)='No errors in the solver!'
    !
    Call set_functional_parameters(skyrme_INI,.False.)
    !-----------------------------------
    ! set multipole moments units
    !-----------------------------------
    Call moments_setUnits
    !
  End Subroutine initialize_HFBTHO_SOLVER
  !=======================================================================
  !> Adjust deformation beta2 and oscillator frequency of the basis and
  !> deformations beta2,beta4 of the initial Woods-Saxon potential
  !=======================================================================
  Subroutine adjust_basis(q2val)
    Implicit None
    Real(pr), Intent(IN) :: q2val
    Real(pr) :: OMEGA0,hbarc,mass_neut
    ! Loose fit based on fission of 240Pu, see PRC 90, 054305 (2014)
    If(q2val.Gt.30.0) Then
       OMEGA0=0.1*EXP(-0.02_pr*Abs(q2val))*Abs(q2val)+6.5_pr
    Else
       OMEGA0=8.1464_pr
    End If
    hbarc=197.328910_pr; mass_neut=938.90590_pr
    oscillator_length=hbarc/Sqrt(OMEGA0*mass_neut)
    basis_deformation=0.05*Sqrt(q2val)
    beta2_deformation=basis_deformation
    beta4_deformation=beta2_deformation/100.0_pr
  End Subroutine adjust_basis
  !=======================================================================
  !> Counts and orders basis states in cylindrical coordinates
  !=======================================================================
  Subroutine base0(lpr)
    Implicit None
    Logical :: lpr
    Integer(ipr) :: iw,k,nre,nze,ke,la,le,ip,ir,iz,il,is,Iall,ilauf,jlauf,ib,nd
    Integer(ipr) :: NOSCIL
    Real(pr), Allocatable :: e(:)
    Real(pr) :: hbz,hbp,ee
    !
    If(n00.Gt.n00max) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='STOP: too large n00 versus n00max'
       Return
    End If
    !-----------------------------------------------
    ! MAXIMUM NUMBER OF THE HO SHELLS (n00,NOSCIL)
    ! (7,120),(8,165),(9,220),(10,286),(11,364)
    ! (12,455),(14,680),(16,969),(18,1330),(20,1771)
    !-----------------------------------------------
    NOSCIL=(n00+1)*(n00+2)*(n00+3)/6
    !-----------------------------------------------
    ! count all states for n00max
    !-----------------------------------------------
    nze=n00max; nre=n00max/2; ke=n00max
    If(basis_HFODD) Then
       nze=n00; nre=n00/2; ke=n00
    End If
    Iall=0;
    Do k=1,ke+1
       la=k-1; le=min0(ke,k)
       Do ip=1,2
          Do ir=0,nre
             Do iz=0,nze
                Do il=la,le
                   Do is=+1,-1,-2
                      If (iz+2*ir+il.Gt.n00max)    Cycle
                      If (il+(is+1)/2.Ne.k)        Cycle
                      If (Mod(iz+il,2).Ne.ip-1)    Cycle
                      Iall=Iall+1
                   End Do
                End Do
             End Do
          End Do
       End Do
    End Do
    !-----------------------------------------------
    ! charge all energies for n00max
    !-----------------------------------------------
    Allocate(e(Iall))
    hbz=two*hbzero/bz**2; hbp=two*hbzero/bp**2;
    Iall=0;
    Do k=1,ke+1
       la=k-1; le=min0(ke,k)
       Do ip=1,2
          Do ir=0,nre
             Do iz=0,nze
                Do il=la,le
                   Do is=+1,-1,-2
                      If (iz+2*ir+il.Gt.n00max)    Cycle
                      If (il+(is+1)/2.Ne.k)        Cycle
                      If (Mod(iz+il,2).Ne.ip-1)    Cycle
                      Iall=Iall+1
                      e(Iall)=hbz*(Real(iz,Kind=pr)+half) &
                         +hbp*(two*Real(ir,Kind=pr)+Real(il,Kind=pr)+one)
                   End Do
                End Do
             End Do
          End Do
       End Do
    End Do
    !-----------------------------------------------
    ! sort energies and derive base cut-off energy
    !-----------------------------------------------
    Call ord(Iall,e);
    If(Iall.Gt.NOSCIL) Then
       EBASECUT=E(NOSCIL)+1.0D-5
    Else
       EBASECUT=E(Iall)+1.0D-5
    End If
    If(basis_HFODD.And.nstate.Le.NOSCIL) EBASECUT=E(nstate)+1.0D-5
    Deallocate(e)
    !-----------------------------------------------
    ! calculate the actual states
    !-----------------------------------------------
    nze=n00max; nre=n00max/2; ke=n00max
    If(basis_HFODD) Then
       nze=n00; nre=n00/2; ke=n00
    EndIf
    ib=0; ilauf=0; ndx=0; nzx=0; nrx=0; nlx=0; nqp=0; nuv=0; nnx=0; nrlx=0; nox=0
    ! loop over k-quantum number
    Do k=1,ke+1
       la=k-1; le=min0(ke,k)
       ! loop over parity
       If(.Not.Parity) jlauf=ilauf !Nop
       Do ip=1,2
          If(Parity) jlauf=ilauf !Yesp
          Do ir=0,nre
             Do iz=0,nze
                Do il=la,le
                   Do is=+1,-1,-2
                      If (iz+2*ir+il.Gt.n00max)    Cycle
                      If (il+(is+1)/2.Ne.k)        Cycle
                      If (Mod(iz+il,2).Ne.(ip-1))  Cycle
                      ee=hbz*(Real(iz,Kind=pr)+half)&
                    +hbp*(two*Real(ir,Kind=pr)+Real(il,Kind=pr)+one)
                      If(ee.Lt.EBASECUT) Then
                         ilauf=ilauf+1
                         nzx=Max(nzx,iz); nrx=Max(nrx,ir); nlx=Max(nlx,il)
                         nnx=Max(nnx,iz+2*ir+il); nrlx=Max(nrlx,2*ir+il)
                         nox=max(nox,k-1)
                      End If
                   End Do
                End Do
             End Do
          End Do
          If(Parity) Then                !Yesp
             If (ilauf.Gt.jlauf) Then
                ib=ib+1
                nd=ilauf-jlauf
                ndx=Max(ndx,nd)
                nqp=nqp+nd; nuv=nuv+nd*nd
             End If
          End If
       End Do
       If(.Not.Parity) Then              !Nop
          If(ilauf.Gt.jlauf) Then
             ib=ib+1
             nd=ilauf-jlauf
             ndx=Max(ndx,nd)
             nqp=nqp+nd; nuv=nuv+nd*nd
          End If
       End If
    End Do
    nbx=ib; ntx=ilauf
    !-----------------------------------------------
    ! print statistics
    !-----------------------------------------------
    If(lpr) Then
       Do iw=lout,lfile
          Write(iw,*)
          Write(iw,'(a)')  '  ---------------------------------------'
          Write(iw,'(a)')  '        Harmonic Oscillator Basis        '
          Write(iw,'(a)')  '  ---------------------------------------'
          Write(iw,'(a,2(i6,2x),a)') '  NUV, NQP:                      ',nuv,nqp
          Write(iw,'(a,2(i6,2x),a)') '  Comparison with bookkeeping spherical basis:'
          Write(iw,'(a,2(i6,2x),a)') '  n00:                           ',n00,n00,  &
               'Maximal number of shells'
          Write(iw,'(a,2(i6,2x),a)') '  nbx, 2*n00+1:                  ',nbx,2*n00+1,  &
               'Maximal number of K-blocks'
          Write(iw,'(a,2(i6,2x),a)') '  ntx, (n00+1)*(n00+2)*(n00+3)/6 ',ntx,(n00+1)*(n00+2)*(n00+3)/6,  &
               'Max.num. p/n levels'
          Write(iw,'(a,2(i6,2x),a)') '  nzx, n00:                      ',nzx,n00,  &
               'Maximal nz-quantum number'
          Write(iw,'(a,2(i6,2x),a)') '  nrx, n00/2  :                  ',nrx,n00/2,  &
               'Maximal nr-quantum number'
          Write(iw,'(a,2(i6,2x),a)') '  nlx, n00:                      ',nlx,n00,  &
               'Maximal ml-quantum number'
          Write(iw,'(a,2(i6,2x),a)') '  ndx, (n00+2)*(n00+2)/4:        ',ndx,(n00+2)*(n00+2)/4,  &
               'Maximal dim. of one k-block'
          Write(iw,*)
       End Do
    End If
    !
  End Subroutine base0
  !=======================================================================
  !> Define quantum numbers for basis states and computes matrix sizes
  !=======================================================================
  Subroutine base(lpr)
    Implicit None
    Logical :: lpr
    Integer(ipr) :: nze,nre,ke,ib,ilauf,jlauf,nom,nnm,  &
         k,la,le,ip,ir,iz,il,is,nn,ND,IBX,N1,N2,iw, &
         klauf
    Real(pr) :: hbz,hbp,ee
    !
    hbz=two*hbzero/bz**2; hbp=two*hbzero/bp**2;
    !
    nze=n00max; nre=n00max/2; ke=n00max
    If(basis_HFODD) Then
       nze=n00; nre=n00/2; ke=n00
    End If
    ib=0; ilauf=0; nzm=0; nrm=0; nlm=0; nom=0; nnm=0; jlauf = 0
    If(finite_range.or.coulomb_gaussian) Then
       If(Allocated(ib_zrls)) Deallocate(i_zrls,ib_zrls)
       If(.not.allocated(ib_zrls)) Then
          Allocate( i_zrls(0:nze,0:nre,0:n00max,0:1))
          Allocate(ib_zrls(0:nze,0:nre,0:n00max,0:1))
       End If
       i_zrls = 0
       ib_zrls = 0
       !-----------------------------------------------
       ! loop over k-quantum number with z=0
       !-----------------------------------------------
       Do k=1,n00max+1
          la=k-1; le=min0(n00max,k)
          klauf = ilauf
          Do ir=0,nre
             Do il=la,le
                Do is=+1,-1,-2
                   jlauf = 0
                   do ip = 1,2
                      do iz = 0,nze
                         If (iz+2*ir+il.Gt.n00max) Cycle
                         If (il+(is+1)/2.Ne.k)     Cycle
                         If (Mod(iz+il,2).Ne.ip-1) Cycle
                         ee=hbz*(Real(iz,Kind=pr)+half) &
                           +hbp*(two*Real(ir,Kind=pr)+Real(il,Kind=pr)+one)
                         If(ee.Lt.EBASECUT) Then
                            jlauf=jlauf+1
                            nzz(ilauf+1,jlauf) = iz
                         End If
                      enddo
                   enddo
                   if(jlauf.ne.0) then
                      ilauf = ilauf + 1
                      nrr(ilauf)=ir; nll(ilauf)=il; nss(ilauf)=is;
                      noo(ilauf) = k; nzzx(ilauf) = jlauf
                   endif
                End Do
             End Do
          End Do
       End Do ! end k
       nttx=ilauf
    End If
    !
    nze=n00max; nre=n00max/2; ke=n00max
    If(basis_HFODD) Then
       nze=n00; nre=n00/2; ke=n00
    End If
    ib=0; ilauf=0; nzm=0; nrm=0; nlm=0; nom=0; nnm=0
    !-----------------------------------------------
    ! loop over k-quantum number
    !-----------------------------------------------
    Do k=1,ke+1
       la=k-1; le=min0(ke,k)
       ! loop over parity
       If(.Not.Parity) jlauf=ilauf !Nop
       Do ip=1,2
          If(Parity) jlauf=ilauf   !Yesp
          Do ir=0,nre
             Do iz=0,nze
                Do il=la,le
                   Do is=+1,-1,-2
                      If (iz+2*ir+il.Gt.n00max) Cycle
                      If (il+(is+1)/2.Ne.k)     Cycle
                      If (Mod(iz+il,2).Ne.ip-1) Cycle
                      ee=hbz*(Real(iz,Kind=pr)+half)&
                    +hbp*(two*Real(ir,Kind=pr)+Real(il,Kind=pr)+one)
                      If(ee.Lt.EBASECUT) Then
                         ilauf=ilauf+1
                         If (ilauf.Gt.ntx) Then
                            ierror_flag=ierror_flag+1
                            ierror_info(ierror_flag)='STOP: in base: ntx too small'
                            Return
                         End If
                         nz(ilauf)=iz; nr(ilauf)=ir; nl(ilauf)=il; ns(ilauf)=is; npar(ilauf)=ip
                         If(finite_range.or.coulomb_gaussian) Then
                            i_zrls(iz,ir,il,(is+1)/2) = ilauf
                            ib_zrls(iz,ir,il,(is+1)/2) = ib+1
                         End If
                         nn =iz+2*ir+il
                         Write(tb(ilauf),100) 2*k-1,tp(ip),nn,iz,il
100                      Format(i2,a1,'[',i2,',',i2,',',i2,']')
                         Do iw=lout,lfile
                            If(lpr.And.IDEBUG.Gt.10) &
                               Write(iw,'(i4,a,i2,a,i2,a,i2,a,i2,a,i2,a,2x,a,1x,a,f14.8)')  &
                                 ilauf,'   nn=',nn,'   nz=',iz,'   nr=',ir,  &
                                 '   ml=',il,'  ms=',is,' /2',tb(ilauf),'e=',ee
                         End Do
                         nzm=Max(nzm,iz); nrm=Max(nrm,ir); nlm=Max(nlm,il)
                         nom=Max(nom,2*k-1); nnm=Max(nnm,iz+2*ir+il)
                      End If
                   End Do
                End Do
             End Do
          End Do
          !-----------------------------------------------
          ! Block memory
          !-----------------------------------------------
          If(Parity) Then                !Yesp
             If (ilauf.Gt.jlauf) Then
                ib=ib+1
                ia(ib)=jlauf; id(ib)=ilauf-jlauf
                ikb(ib)=k; ipb(ib)=ip
                Write(txb(ib),'(i3,a,i2,a,a1)') ib,'. block:  k=',k+k-1,'/2',tp(ip)
                !ir=(ib+1)/2
                !Write(*,*)  ib,2*k-1,'2*Omega=',2*ir - 1
                Do iw=lout,lfile
                   If(lpr.And.IDEBUG.Gt.10) Write(iw,'(/,a,i3,a,a1)')'  For the above block:  k=',k+k-1,'/2',tp(ip)
                End Do
             End If
             If(id(ib).Eq.0) Then
                ierror_flag=ierror_flag+1
                ierror_info(ierror_flag)='STOP: in base Block Memory(1)'
                Return
             End If
          End If
       End Do ! end of ip
       !-----------------------------------------------
       ! Block memory
       !-----------------------------------------------
       If(.Not.Parity) Then               !Nop
          If (ilauf.Gt.jlauf) Then
             ib=ib+1
             ia(ib)=jlauf; id(ib)=ilauf-jlauf
             nn = nz(ilauf)+2*nr(ilauf)+nl(ilauf); ip = 2 - Mod(nn,2)
             ikb(ib)=k; ipb(ib)=ip
             Write(txb(ib),'(i3,a,i2,a,a1)') ib,'. block:  k=',k+k-1,'/2',tp(ip)
             Do iw=lout,lfile
                If(lpr.And.IDEBUG.Gt.10) Write(iw,'(/,a,i3,a,a1)')'  For the above block:  k=',k+k-1,'/2',tp(ip)
             End Do
          End If
          If(id(ib).Eq.0) Then
             ierror_flag=ierror_flag+1
             ierror_info(ierror_flag)='STOP: in base Block Memory(2)'
             Return
          End If
       End If
    End Do ! end k
    nb=ib;  nt=ilauf
    !-----------------------------------------------
    ! broyden/linear mixing (storage)
    !-----------------------------------------------
    nhhdim=0
    Do ib=1,NB
       ND=ID(ib)
       Do N1=1,ND
          Do N2=1,N1
             nhhdim=nhhdim+1
          End Do
       End Do
    End Do
    nhhdim2=2*nhhdim; nhhdim3=3*nhhdim; nhhdim4=4*nhhdim
    If(Allocated(brin)) Deallocate(brin,brout)
    If(neck_constraints) Then
       Allocate(brin(nhhdim4+lambdaMax+1),brout(nhhdim4+lambdaMax+1))
    Else
       Allocate(brin(nhhdim4+lambdaMax),brout(nhhdim4+lambdaMax))
    End If
    !-----------------------------------------------
    ! Print statistics
    !-----------------------------------------------
    If(lpr) Then
       Do iw=lout,lfile
          Write(iw,'(a,i4)')   '  Actual basis used'
          Write(iw,'(a,i4)')   '  Number of blocks: nb .......: ',nb
          Write(iw,'(a,i4)')   '  Number of levels: nt .......: ',nt
          Write(iw,'(a,i4)')   '  Maximal 2*omega : nom ......: ',nom
          Write(iw,'(a,i4)')   '  Maximal nz:       nzm ......: ',nzm
          Write(iw,'(a,i4)')   '  Maximal nr:       nrm ......: ',nrm
          Write(iw,'(a,i4)')   '  Maximal ml:       nlm ......: ',nlm
          Write(iw,'(a,i4)')   '  Maximal N=nz+2*nr+nl .......: ',nnm
          Write(iw,'(a,i4)')   '  2 x biggest block dim. .....: ',ndx2
          Write(iw,'(a,i8)')   '  Non-zero elements of h .....: ',nhhdim
          Write(iw,'(a,i8)')   '  Number of Broyden elements .: ',nhhdim4
          Write(iw,'(a,i4)')
       End Do
    End If
    If(nzm.Ge.n00max.Or.(nom-1)/2.Eq.n00max) Then
       Write(*,*) 'nzm=',nzm,'  (nom-1)/2=',(nom-1)/2,'  n00max=',n00max
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='STOP: Please increase n00max to have correct basis'
    End If
  End Subroutine base
  !=======================================================================
  !> Subroutine gaupol() computes the harmonic oscillator wave functions
  !> and their first derivatives in cylindrical coordinates \f$ (\rho,\theta,z )\f$. The full HO
  !> wave function reads
  !>  \f[
  !>     \psi_{n_{r}\Lambda n_{z}}(\rho, \theta, z) =
  !>     \psi_{n_{r}}^{|\Lambda|}(\eta) \psi_{n_{z}}(\xi) \frac{e^{i\Lambda\theta}}{\sqrt{2\pi}}
  !>                  \chi_{\Sigma}(\sigma),
  !>  \f]
  !>  where the dimensionless coordinates are \f$ \eta = \beta_{\perp}^{2}\rho^2 \f$,
  !>  \f$ \xi = \beta_{z}z \f$, \f$ \Lambda \f$ is the eigenvalue of \f$ \hat{\ell}_{z} \f$,
  !>  \f$ \Sigma \f$ the eigenvalue of \f$ \hat{s}_{z} \f$ and \f$ \sigma \f$ is a number
  !>  that labels the spin eigenstate in coordinate space, by opposition to \f$ \chi_{\Sigma} \f$,
  !>  which refers to a vector in spin space. Under the action of time-reversal, the HO
  !>  basis states turn into \f$ \psi_{n_{r}-\Lambda n_{z}}(\rho, \theta, z) \f$, where
  !>  the spatial part proportional to \f$ \rho, z \f$ is unchanged.
  !=======================================================================
  Subroutine gaupol(lpr)
    Implicit None
    Logical :: lpr
    Real(pr) :: w0,z,x,s,s0,s1,w00,w4pii,dsq,d1,d2,d3,d4,hs0,hs1
    Integer(ipr) :: ih,il,iw,ix,n,l,n1,n2
    !-----------------------------------------------
    !> - For the z-coordinate, the HO wave function reads
    !>     \f[
    !>        \psi_{n_{z}}(z) = \mathcal{N}_{n_{z}}\beta_{z}^{1/2} e^{-\xi^{2}/2} H_{n_{z}}(\xi),
    !>     \f]
    !>   with the normalization factor
    !>     \f[
    !>        \mathcal{N}_{n_{z}}
    !>      = \left( \frac{1}{\sqrt{\pi}2^{n_{z}}n_{z}!} \right)^{1/2}.
    !>     \f]
    !>   The first derivative of the HO wave function is
    !>     \f[
    !>        \frac{\partial\psi_{n_{z}}}{\partial z} = \beta_{z}
    !>            \left[ \sqrt{2n_{z}} \psi_{n_{z}-1}(\xi) - \xi \psi_{n_{z}}(\xi) \right].
    !>     \f]
    !>   In the code, variable \f$ \mathtt{qh} \f$ contains \f$ (1/\beta^{1/2}_{z})\psi_{n_{z}} \f$
    !>   on the nodes of the Gauss-Hermite integration mesh; variable \f$ \mathtt{qh1} \f$
    !>   contains the value of \f$ (1/\beta^{3/2}_{z})\partial\psi_{n_{z}}/\partial z \f$ on that same
    !>   mesh. Note that this is consistent with a dimension analysis: \f$ \psi_{n_{z}}(z) \f$
    !>   has dimension \f$ [L]^{-1/2} \f$ so that the integral over space of the square
    !>   \f$ |\psi_{n_{z}}(z)|^2 \f$ is dimensionless and can be interpreted as a probability.
    !>   Similarly, the dimension of the derivative \f$ \partial\psi_{n_{z}}/\partial z \f$ is
    !>   \f$ [L]^{-3/2} \f$. In the code, objects \f$ \mathtt{qh} \f$ and \f$ \mathtt{qh1} \f$
    !>   are thus dimensionless.
    !-----------------------------------------------
    w4pii=pi**(-0.250_pr)
    Do ih=1,ngh
       z=xh(ih); w0=w4pii*Exp(-half*z*z)
       w0 =w0*Sqrt(wh(ih))
       qh(0,ih)=w0;       qh(1,ih)=sq(2)*w0*z
       qh1(0,ih)=-w0*z;   qh1(1,ih)=sq(2)*w0*(one-z*z)
       Do n=2,nzm
          qh(n,ih)=sqi(n)*(sq(2)*z*qh(n-1,ih)-sq(n-1)*qh(n-2,ih))
          qh1(n,ih)=sq(n+n)*qh(n-1,ih)-z*qh(n,ih)
       End Do
    End Do
    !-----------------------------------------------
    !> - For the \f$\rho\f$-coordinate, the HO wave function reads
    !>     \f[
    !>        \psi_{n_{r}}^{\Lambda}(\eta) =
    !>        \mathcal{N}_{n_{r}}^{\Lambda} \beta_{\perp}\sqrt{2}
    !>              \eta^{|\Lambda|/2} e^{-\eta/2} L_{n_{r}}^{|\Lambda|}(\eta),
    !>     \f]
    !>   with the normalization factor
    !>     \f[
    !>        \mathcal{N}_{n_{r}}^{\Lambda}
    !>      = \left( \frac{n_{r}}{ (n_{r} + |\Lambda|)!} \right)^{1/2}.
    !>     \f]
    !>   The first derivative of the HO wave function is
    !>     \f[
    !>        \frac{\partial\psi_{n_{r}}^{\Lambda}}{\partial\rho} = \frac{\beta_{\perp}}{\sqrt{\eta}}
    !>            \left[ (2n_{r} + |\Lambda| - \eta)\psi_{n_{r}}^{\Lambda}(\eta)
    !>                   -2\sqrt{n_{r}(n_{r}+|\Lambda|)}\psi_{n_{r}-1}^{\Lambda}(\eta) \right].
    !>     \f]
    !>   In the code, variable \f$ \mathtt{ql} \f$ contains the value of
    !>   \f$ (1/\beta_{\perp}\sqrt{2})\psi_{n_{r}}^{\Lambda} \f$ on the nodes of the Gauss-Laguerre
    !>   integration mesh; variable \f$ \mathtt{ql1} \f$ contains the value of
    !>   \f$ (\sqrt{\eta}/\beta_{\perp}^{2}\sqrt{2}) \partial\psi_{n_{r}}^{\Lambda}/\partial \rho \f$
    !>   on that same mesh.
    !-----------------------------------------------
    Do il=1,ngl
       x=xl(il); w00=sq(2)*Exp(-half*x)
       Do l=0,nlm
          w0=w00*Sqrt(half*wl(il)*x**l)
          ql(0,l,il)=wfi(l)*w0;         ql(1,l,il)=(l+1-x)*wfi(l+1)*w0
          ql1(0,l,il)=(l-x)*wfi(l)*w0;  ql1(1,l,il)=(Real(l*l+l,Kind=pr) &
                                                  -x*Real(l+l+3,Kind=pr)+x*x)*wfi(l+1)*w0
          Do n=2,nrm
             dsq=sq(n)*sq(n+l); d1=Real(n+n+l-1,Kind=pr)-x
             d2=sq(n-1)*sq(n-1+l); d3=n+n+l-x; d4=two*dsq
             ql(n,l,il)=(d1*ql(n-1,l,il)-d2*ql(n-2,l,il))/dsq
             ql1(n,l,il)=d3*ql(n,l,il)-d4*ql(n-1,l,il)
          End Do
       End Do
    End Do
    !-----------------------------------------------
    ! Test accuracy for Hermite orthonormalization
    !-----------------------------------------------
    hs0=zero; hs1=two
    Do n1=0,nzm
       Do n2=0,n1
          If (Mod(n1-n2,2).Eq.0) Then
             s=zero
             Do ih=1,ngh
                s=s+qh(n1,ih)*qh(n2,ih)
             End Do
             If(n1.Ne.n2) Then
                hs0=Max(s,hs0)
             Else
                hs1=Min(s,hs1)
             End If
          End If
       End Do
    End Do
    !-----------------------------------------------
    ! Test accuracy for Laguerre orthonormalization
    !-----------------------------------------------
    s0=zero; s1=two
    Do l=0,nlm
       Do n1=0,nrm
          Do n2=0,n1
             s=zero
             Do il=1,ngl
                s=s+ql(n1,l,il)*ql(n2,l,il)
             End Do
             If(n1.Ne.n2) Then
                s0=Max(s,s0)
             Else
                s1=Min(s,s1)
             End If
          End Do
       End Do
    End Do
    !-----------------------------------------------
    ! print accuracy
    !-----------------------------------------------
    If(lpr) Then
       Do iw=lout,lfile
          Write(iw,'(a)')  '  ---------------------------------------'
          Write(iw,'(a)')  '            Integration Meshes           '
          Write(iw,'(a)')  '  ---------------------------------------'
          Write(iw,'(a,i3)')  &
               '  Number of Gauss-Hermite mesh points ngh ....: ',ngh
          Write(iw,'(a,i3)')  &
               '  Number of Gauss-Laguerre mesh points ngl ...: ',ngl
          Write(iw,'(a,i3)')  &
               '  Number of Gauss-Legendre mesh points nleg ..: ',nleg
          Write(iw,'(a)') &
               '  Integration boundaries'
          Write(iw,'(2(a,f12.8))')  &
               '    Hermite  - from xh(1)  =',xh(1),  ' to xh(ngh)   =',xh(ngh)
          Write(iw,'(2(a,f12.8))')  &
               '    Laguerre - From xl(1)  =',xl(1),  ' to xl(ngl)   =',xl(ngl)
          If(nleg.Gt.0) Then
             Write(iw,'(2(a,f12.8))')  &
                  '    Legendre - From xleg(1)=',xleg(1),' to xleg(nleg)=',xleg(nleg)
          End If
          Write(iw,*)  &
               ' Max.dev.in:     Orthogonality            Normalization'
          Write(iw,*) ' Hermite  ',hs0,Abs(one-hs1)
          Write(iw,*) ' Laguerre ',s0,Abs(one-s1)
       End Do
    End If
    !-----------------------------------------------
    ! debug
    !-----------------------------------------------
    If (lpr.And.IDEBUG.Gt.20) Then
       ix=3
       Do iw=lout,lfile
          Write(iw,*) ' nz    qh(nz,ih=1,...)'
          Do n=0,nzm
             Write(iw,'(i4,3f15.8)') n,(qh(n,ih),ih=1,ix)
             Write(iw,'(i4,3f15.8)') n,(qh1(n,ih),ih=1,ix)
             Write(iw,*) ' '
          End Do
          Do l=0,nlm
             Write(iw,*) ' nr ml    ql(nr,l,il=1,...)'
             Do n=0,nrm
                Write(iw,'(i4,i3,3f15.8)') n,l,(ql(n,l,il),il=1,ix)
                Write(iw,'(i4,i3,3f15.8)') n,l,(ql1(n,l,il),il=1,ix)
                Write(iw,*) ' '
             End Do
          End Do
       End Do
       !-----------------------------------------------
       ! Test for Hermite polynomials normalization
       !-----------------------------------------------
       Do n1=0,nzm
          Do n2=0,n1
             If (Mod(n1-n2,2).Eq.0) Then
                s=zero
                Do ih=1,ngh
                   s=s+qh(n1,ih)*qh(n2,ih)
                End Do
                Do iw=lout,lfile
                   Write(iw,100) n1,n2,s
                End Do
100             Format(' Gauss-Hermite: n1=',i3,'  n2=',i3,f20.8)
             End If
          End Do
       End Do
       !-----------------------------------------------
       ! Test for Laguerre polynomials normalization
       !-----------------------------------------------
       Do l=0,nlm
          Do n1=0,nrm
             Do n2=0,n1
                s=zero
                Do il=1,ngl
                   s=s+ql(n1,l,il)*ql(n2,l,il)
                End Do
                Do iw=lout,lfile
                   Write(iw,101) l,n1,n2,s
101                Format(' Gauss Laguerre: l='  &
                        ,i2,' n1=',i3,'  n2=',i3,f20.8)
                End Do
             End Do
          End Do
       End Do
    End If
    !
    Call coordinateLST(.False.)  ! coordinate LST
    !
  End Subroutine gaupol
  !=======================================================================
  !> Scratch initialization of the self-consistent loop by diagonalization
  !> of the Woods-Saxon Hamiltonian
  !=======================================================================
  Subroutine start
    !---------------------------------------------------------------------
    ! initializes scratch Saxon-Woods potentials
    !---------------------------------------------------------------------
    Implicit None
    Integer(ipr) :: iw,i,ih,il,ihl,it,ita
    Real(pr) :: zb(ngh),rrb(ngl),rb(ngl),rav,rao,vpws,vls,betas,gamma,fac, &
                facb,zz,rr,r,ctet,cphi,p2,p20,p22,s,u,w,f,rc,c,beta00,     &
                b2_ws,b4_ws,pleg2,pleg4
    !----------------------------------------------------------------------------
    ! Re-initializing all again since scratch calculation
    !----------------------------------------------------------------------------
    Call initialize_HFBTHO_SOLVER
    If(ierror_flag.Ne.0) Return
    Call Constraint_or_not(inin_INI,inin,icstr)
    If(ierror_flag.Ne.0) Return
    Do it=itmin,itmax
       If(npr(it).Ne.2*(npr(it)/2) .And. .Not.fragment_properties .And. .Not.odd_noBlock) Then
          irestart=irestart+1; npr(it)=npr_INI(it)
       End If
       If(npr(it).Ne.2*(npr(it)/2) .And. .Not.fragment_properties .And. odd_noBlock) Then
          irestart=irestart+1; iparenti(it)=0
       End If
    End Do
    npr(3)=npr(1)+npr(2)
    If(irestart.Ne.0) Then
       ! odd nucleus requested but no even-even solution, recalculate the even-even nucleus from scratch
       If(.Not.odd_noBlock) Then
          Do iw=lout,lfile
             Write(iw,'(1x,a,2i4)')
             Write(iw,'(1x,a,2i4)') ' Initialization for the even-even core (N,Z)=: ',npr(1:2)
          End Do
       Else
          Do iw=lout,lfile
             Write(iw,'(1x,a,2i4)')
             Write(iw,'(1x,a,2i4)') ' Initialization without blocking for (N,Z)=: ',npr(1:2)
          End Do
       End If
    Else
       ! scratch for the even-even nucleus requested
       Do iw=lout,lfile
          Write(iw,'(1x,a,2i4)')
          Write(iw,'(a,a,3i4)')    '  Scratch initialization for the nucleus: ',nucname,npr(1:2)
          Write(iw,'(1x,a,2i4)')
       End Do
    End If
    n00=Abs(n00_INI);  b0=b0_INI;           q=q_INI; iLST=iLST_INI
    maxi=MAX_ITER_INI; inin=inin_INI;
    skyrme=skyrme_INI; kindhfb=kindhfb_INI
    iproj=iproj_INI;   npr1pj=npr1pj_INI;   npr2pj=npr2pj_INI;
    icacou=0; icahartree=0
    !
    Call preparer(.False.)
    !
    If(ierror_flag.Ne.0) Return
    inin=Abs(inin)          ! positive even if inin_INI is not
    If(Abs(b2_0).Gt.1.5_pr) b2_0=1.5_pr ! Avoid crazy initial points (quadrupole deformation)
    If(Abs(b4_0).Gt.1.0_pr) b4_0=1.0_pr ! Avoid crazy initial points (hexadecapole deformation)
    !-----------------------------------
    ! Saxon-Woods potentials
    !-----------------------------------
    Do iw=lout,lfile
       Write(iw,'(/,a)') '  Initial potentials of Saxon-Woods shape '
    End Do
    beta00=bet     ! wf to requested deformation
    Do iw=lout,lfile
       Write(iw,'(a,2f14.8)') '  v0ws   =',v0ws
       Write(iw,'(a,2f14.8)') '  kappa  =',akv
       Write(iw,'(a,2f14.8)') '  vs0    =',vso
       Write(iw,'(a,2f14.8)') '  r0     =',r0v
       Write(iw,'(a,2f14.8)') '  a      =',av
       Write(iw,'(a,2f14.8)') '  r0-so  =',rso
       Write(iw,'(a,2f14.8)') '  a-so   =',aso
       Write(iw,'(a,f14.8)')  '  b2_ws  =',b2_0
       Write(iw,'(a,f14.8)')  '  b4_ws  =',b4_0
    End Do
    !-----------------------------------
    ! Densities
    !-----------------------------------
    Do it=itmin,itmax
       ita=3-it; rav=r0v(it)*amas**p13; rao=rso(it)*amas**p13
       vpws=v0ws*(one-akv*(npr(it)-npr(ita))/amas)
       vls=half*(hqc/amu)**2*vpws*vso(it)
       ! Deformations of the surface
       b2_ws = b2_0 * Sqrt(5.0_pr/(4.0_pr*pi))
       b4_ws = b4_0 * Sqrt(9.0_pr/(4.0_pr*pi))
       ! Volume conservation condition
       !gamma=zero
       !fac= one+betas*Cos( gamma*pi/180.0_pr)
       !fac=(one+betas*Cos((gamma+120.0_pr)*pi/180.0_pr))*fac
       !fac=(one+betas*Cos((gamma-120.0_pr)*pi/180.0_pr))*fac
       !fac=fac**(-p13)
       fac = two + (143.0_pr*Sqrt(five)*b2_ws**3 + 1287.0_pr*b2_ws**2*b4_ws + 390.0_pr*Sqrt(five)*b2_ws*b4_ws**2 &
            + 243.0_pr*b4_ws**3)/(2002.0_pr*Pi**1.5_pr) + (three*(b2_ws**2 + b4_ws**2))/(two*Pi)
       fac=(two/fac)**(p13)
       ! z,r-coordinates in fm
       zb=xh*bz; rrb=xl*bp*bp; rb=Sqrt(rrb)
       Do ih=1,ngh
          zz=zb(ih)**2
          Do il=1,ngl
             rr=rrb(il)+zz; r=Sqrt(rr); ctet=zz/rr
             ! Deformed surface
             !p20=3.0_pr*ctet-one; p22=Sqrt(3.0_pr)*cphi
             !p2=p20*Cos(gamma*pi/180.0_pr)+p22*Sin(gamma*pi/180.0_pr)
             pleg2 = half*(three*ctet - one)
             pleg4 = (35.0_pr*ctet**2 - 30.0_pr*ctet + three)/eight
             facb=fac*(one + b2_ws*pleg2 + b4_ws*pleg4)
             ! Woods-Saxon potential
             u= vpws/( one+Exp( (r-rav*facb) / av(it) ))
             w=-vls /( one+Exp( (r-rao*facb) / aso(it)))
             ihl=ih+(il-1)*ngh
             If(it.Eq.1) Then
                vhbn(ihl)=hb0n; vn(ihl)=u; vsn(ihl)=w;
                vrn(ihl)=zero; vzn(ihl)=zero; vdn(ihl)=zero;
                vSFIZn(ihl)=zero; vSZFIn(ihl)=zero;
                vSFIRn(ihl)=zero; vSRFIn(ihl)=zero;
             Else
                vhbp(ihl)=hb0p; vp(ihl)=u; vsp(ihl)=w;
                vrp(ihl)=zero; vzp(ihl)=zero; vdp(ihl)=zero;
                vSFIZp(ihl)=zero; vSZFIp(ihl)=zero;
                vSFIRp(ihl)=zero; vSRFIp(ihl)=zero;
             End If
             ro(ihl,it)=u
             aka(ihl,it)=5.0d-3*Exp((r-rav*facb)/2.0_pr)
          End Do
       End Do
       !s=npr(it)/Sum(ro(:,it))
       s=tz(it)/Sum(ro(:,it))
       Do il=1,ngl
          Do ih=1,ngh
             ihl=ih+(il-1)*ngh
             f=s/(pi*wh(ih)*wl(il)* bz*bp*bp); ro(ihl,it)=f*ro(ihl,it)
          End Do
       End Do
       !-----------------------------------
       ! pairing
       !-----------------------------------
       Do il=1,nghl
          If(it.Eq.1) Then
             dvn(il)=-100.0_pr*aka(il,it)
          Else
             dvp(il)=-100.0_pr*aka(il,it)
          End If
       End Do
    End Do
    !-----------------------------------
    ! coulomb
    !-----------------------------------
    If(icou.le.0) Then
       cou=zero
    Else
       rc=r0v(2)*amas**p13
       Do il=1,ngl
          Do ih=1,ngh
             r=Sqrt(zb(ih)**2+rrb(il))
             If (r.Lt.rc) Then
                c=half*(3/rc-r*r/(rc**3))
             Else
                c=one/r
             End If
             !cou(ih+(il-1)*ngh)=c*npr(2)/alphi
             cou(ih+(il-1)*ngh)=c*tz(2)/alphi
          End Do
       End Do
    End If
    !-----------------------------------
    ! initial ph+pp matrix elements
    !-----------------------------------
    ak=0.1_pr; rk=0.1_pr ! initial density matrix elements (improve later)
    brin=zero             ! initial matrix elements to zero
    iiter=0               ! iteration number iiter to zero
    Call gamdel(.true.,.true.)
    !
  End Subroutine start
  !=======================================================================
  !> Display the density matrix \f$ \rho \f$ and pairing tensor \f$ \kappa \f$
  !> as a function of the (dimensionless) radius \f$ r = \sqrt{\eta^2 + \xi^2}\f$
  !=======================================================================
  Subroutine printRHO
    Implicit None
    Integer(ipr), Save :: ihli,Ifle
    Ifle=76+iLST1
#if(DO_MASSTABLE==1 || DRIP_LINES==1 || DO_PES==1)
    If(Ifle.Eq.76) Open(Ifle,file='ho_den'//row_string//'.dat',status='unknown')
    If(Ifle.Eq.77) Open(Ifle,file='tho_den'//row_string//'.dat',status='unknown')
#else
    If(Ifle.Eq.76) Open(Ifle,file='ho_den.dat',status='unknown')
    If(Ifle.Eq.77) Open(Ifle,file='tho_den.dat',status='unknown')
#endif
    Write(Ifle,*) 'r  denN  denP  akaN  akaP '
    Do ihli=1,nghl
       Write(Ifle,'(12(1x,e16.8))') Sqrt(fh(ihli)**2+fl(ihli)**2)  &
            ,ro(ihli,1),ro(ihli,2),aka(ihli,1),aka(ihli,2)
    End Do
    Close(Ifle)
  End Subroutine printRHO
  !=======================================================================
  !> Calculates sign, Sqrt, factorials, etc. of integers and half integers
  !>   - \f$ \mathtt{iv(n) } = (-1)^n    \f$
  !>   - \f$ \mathtt{sq(n) } = \sqrt{n}, \f$
  !>   - \f$ \mathtt{sqi(n)} = 1/\sqrt{n}\f$
  !>   - \f$ \mathtt{fak(n)} = n!        \f$
  !>   - \f$ \mathtt{wf(n) } = \sqrt{n!} \f$
  !>   - \f$ \mathtt{wfi(n)} = 1/\sqrt{n!}\f$
  !=======================================================================
  Subroutine gfv
    Implicit None
    Integer(ipr) :: i,igfv
    Parameter(igfv=170)               !maximal number for GFV
    If(Allocated(iv)) Deallocate(iv,fak,fi,sq,sqi,wf,wfi)
    Allocate(iv(-igfv:igfv),fak(0:igfv),fi(0:igfv),sq(0:igfv),sqi(0:igfv))
    Allocate(wf(0:igfv),wfi(0:igfv))
    iv(0)=1; sq(0)=zero; sqi(0)=1.0d30
    fak(0)=one; fi(0)=one; wf(0)=one; wfi(0)=one
    Do i=1,igfv
       iv(i)=-iv(i-1)
       iv(-i) = iv(i)
       sq(i)=Sqrt(Real(i,Kind=pr)); sqi(i)=one/sq(i)
       fak(i)=Real(i,Kind=pr)*fak(i-1); fi(i)=one/fak(i)
       wf(i)=sq(i)*wf(i-1); wfi(i)=one/wf(i)
    End Do
  End Subroutine gfv
  !=======================================================================
  !> Diagonalization of a real, symemtric matrix (backup routine if LAPACK fails)
  !=======================================================================
  Subroutine sdiag(nmax,n,a,d,x,e,is)
    !---------------------------------------------------------------------
    ! A   matrix to be diagonalized
    ! D   eigenvalues,  X   eigenvectors, E   auxiliary field
    ! IS=1  eigenvalues are ordered (major component of X is positive)
    ! 0  eigenvalues are not ordered
    !---------------------------------------------------------------------
    Use HFBTHO_utilities, Only: pr,ipr
    Implicit None
    Integer(ipr), Save :: i,j,j1,k,l,im
    Integer(ipr)       :: n,nmax,is
    Real(pr), Save :: f,g,h,hi,s,p,b,r,pra,c
    Real(pr) :: a(nmax,nmax),x(nmax,nmax),e(n),d(n)
    Real(pr), Save :: tol=1.0D-32,eps=9.0D-12,one=1.0_pr,zero=0.0_pr
    !
    If (n.Le.1) Then
       d(1)=a(1,1); x(1,1)=one
       Return
    End If
    Do i=1,n
       Do j=1,i
          x(i,j)=a(i,j)
       End Do
    End Do
    ! householder-reduktion
    i=n
15 Continue
    If (i.Ge.2) Then
       l=i-2
       f=x(i,i-1); g=f; h=zero
       If (l.Gt.0) Then
          Do k=1,l
             h=h+x(i,k)*x(i,k)
          End Do
       End If
       s=h+f*f
       If (s.Lt.tol) Then
          h=zero
          Go To 100
       End If
       If (h.Gt.zero) Then
          l=l+1; g=Sqrt(s)
          If (f.Ge.zero) g=-g
          h=s-f*g; hi=one/h; x(i,i-1)=f-g; f=zero
          If (l.Gt.0) Then
             Do j=1,l
                x(j,i)=x(i,j)*hi
                s=zero
                Do k=1,j
                   s=s+x(j,k)*x(i,k)
                End Do
                j1=j+1
                If (l.Ge.j1) Then
                   Do k=j1,l
                      s=s+x(k,j)*x(i,k)
                   End Do
                End If
                e(j)=s*hi; f=f+s*x(j,i)
             End Do
          End If
          f=f*hi*0.50_pr
          If (l.Gt.0) Then
             Do j=1,l
                s=x(i,j); e(j)=e(j)-f*s; p=e(j)
                Do  k=1,j
                   x(j,k)=x(j,k)-s*e(k)-x(i,k)*p
                End Do
             End Do
          End If
       End If
100  Continue
       d(i)=h; e(i-1)=g; i=i-1
       Go To 15
       ! Bereitstellen der Transformationmatrix
    End If
    d(1)=zero; e(n)=zero; b=zero; f=zero
    Do i=1,n
       l=i-1
       If (d(i).Eq.0.) Go To 221
       If (l.Gt.0) Then
          Do J=1,L
             s=zero
             Do k=1,l
                s=s+x(i,k)*x(k,j)
             End Do
             Do k=1,l
                x(k,j)=x(k,j)-s*x(k,i)
             End Do
          End Do
       End If
221  Continue
       d(i)=x(i,i)
       x(i,i)=one
       If (l.Gt.0) Then
          Do j=1,l
             x(i,j)=zero; x(j,i)=zero
          End Do
       End If
    End Do
    ! Diagonalisieren der Tri-Diagonal-Matrix
    Do l=1,n
       h=eps*(Abs(d(l))+ Abs(e(l)))
       If (h.Gt.b) b=h
       ! Test fuer Splitting
       Do  j=l,n
          If (Abs(e(j)).Le.b) Exit
       End Do
       ! test fuer konvergenz
       If (j.Eq.l) Go To 300
340  p=(d(l+1)-d(l))/(2.0_pr*e(l))
       r=Sqrt(p*p+one); pra=p+r
       If (p.Lt.zero) pra=p-r
       h=d(l)-e(l)/pra
       Do i=l,n
          d(i)=d(i)-h
       End Do
       f=f+h
       ! QR-transformation
       p=d(j); c=one; s=zero; i=j
360  i=i-1
       If (i.Lt.l) Go To 362
       g=c*e(i); h=c*p
       If ( Abs(p).Ge.Abs(e(i))) Then
          c=e(i)/p
          r=Sqrt(c*c+one); e(i+1)=s*p*r; s=c/r; c=one/r
          Go To 365
       End If
       c=p/e(i)
       r=Sqrt(c*c+one); e(i+1)=s*e(i)*r; s=one/r; c=c/r
365  p=c*d(i)-s*g
       d(i+1)=h+s*(c*g+s*d(i))
       Do k=1,n
          h=x(k,i+1); x(k,i+1)=x(k,i)*s+h*c
          x(k,i)=x(k,i)*c-h*s
       End Do
       Go To 360
362  e(l)=s*p
       d(l)=c*p
       If ( Abs(e(l)).Gt.b) Go To 340
       ! konvergenz
300  d(l)=d(l)+f
    End Do
    If (is.Eq.0) Return
    ! ordnen der eigenwerte
    Do i=1,n
       k=i; p=d(i); j1=i+1
       If (j1.Le.n) Then
          Do j=j1,n
             If (d(j).Ge.p) Cycle
             k=j; p=d(j)
          End Do
          If (k.Eq.i) Cycle
          d(k)=d(i); d(i)=p
          Do j=1,n
             p=x(j,i); x(j,i)=x(j,k)
             x(j,k)=p
          End Do
       End If
    End Do
    ! signum
    Do  k=1,n
       s=zero
       Do i=1,n
          h=Abs(x(i,k))
          If (h.Gt.s) Then
             s=h; im=i
          End If
       End Do
       If (x(im,k).Lt.zero) Then
          Do i=1,n
             x(i,k)=-x(i,k)
          End Do
       End If
    End Do
  End Subroutine sdiag
  !=======================================================================
  !> Determines the symbol of the current element (is=1) or the proton
  !> number of a given element (is=2)
  !=======================================================================
  Subroutine nucleus(is,npr2,te)
    Use HFBTHO_utilities, Only: pr,ipr
    Use HFBTHO, Only: ierror_flag,ierror_info
    Implicit None
    Integer(ipr) :: is,npr2,np
    Integer(ipr) :: maxz
    Parameter (maxz=133)
    Character(2) te
    Character(2*maxz+2) t
    T(  1: 40)=' n HHeLiBe B C N O FNeNaMgAlSi P ScLaR K'
    T( 41: 80)='CaScTi VCrMnFeCoNiCuZnGaGeASSeBrKrRbSr Y'
    T( 81:120)='ZrNbMoTcRoRhPdAgCdInSnSbTe IXeCsBaLaCePr'
    T(121:160)='NdPmSmEuGdTbDyHoErTmYbLuHfTa WReOsIrPtAu'
    T(161:200)='HgTlPbBiPoAtRnFrRaAcThPa UNpPuAmCmBkCfEs'
    T(201:240)='FmMdNoLrRfDbSgBhHsMtDsRgCnNhFlMcLvTsOgXx'
    T(241:265)='XxXxXxXxXxXxXxXxXxXxXxXxX'
    If (is.Eq.1) Then
       If (npr2.Lt.0.Or.npr2.Gt.maxz) Then
          ierror_flag=ierror_flag+1
          ierror_info(ierror_flag)='STOP: in nucleus npr2 is wrong:'
          Return
       End If
       te=t(2*npr2+1:2*npr2+2)
       Return
    Else
       Do np=0,maxz
          If (te.Eq.t(2*np+1:2*np+2)) Then
             npr2=np
             Return
          End If
       End Do
    End If
    ierror_flag=ierror_flag+1
    ierror_info(ierror_flag)='STOP: in nucleus the nucleus is unknown!'
  End Subroutine nucleus
  !=======================================================================
  !> Given 'Z' returns mass number 'A' on the stability line
  !=======================================================================
  Subroutine stab(npr2,npr3)
    Implicit None
    Integer(ipr) :: npr2,npr3
    Real(pr), Save :: sn,sz,dsn,c,c5
    c=0.0060_pr; c5=5.0_pr*c/3.0_pr; sz=npr2; sn=npr2
    Do While(Abs(dsn).Lt.1.0d-5)
       dsn=sz-sn+c*(sn+sz)**(5.0_pr/3.0_pr)
       dsn=dsn/(-1.0_pr+c5*(sn+sz)**(2.0_pr/3.0_pr))
       sn=sn-dsn
    End Do
    npr3=Int(sn)
  End Subroutine stab
  !=======================================================================
  !> Orders a vector 'e' of size 'n' by ascending values
  !=======================================================================
  Subroutine ord(n,e)
    Implicit None
    Integer(ipr) :: n,i,k,j
    Real(pr), Save :: p
    Real(pr) :: e(n)
    Do i=1,n
       k=i; p=e(i)
       If (i.Lt.n) Then
          Do j=i+1,n
             If (e(j).Lt.p) Then
                k=j; p=e(j)
             End If
          End Do
          If (k.Ne.i) Then
             e(k)=e(i); e(i)=p
          End If
       End If
    End Do
  End Subroutine ord
  !=======================================================================
  !> Computes the Lipkin-Nogami parameter \f$ \lambda_2 \f$ for a monopole
  !> pairing force of the type
  !>    \f[
  !>         \bar{v}_{ijkl} = -G \delta_{j\bar{i}}\delta_{l\bar{k}}
  !>                           \text{sign}(j)\text{sign}(l),
  !>    \f]
  !> The calculation is performed in the canonical basis
  !>    \f[
  !>        \lambda_{2} = \frac{G}{4}\frac{\displaystyle
  !>          \sum_{i>0} u_{i}v_{i}^{3}\sum_{i>0} u_{i}^{3}v_{i}
  !>          - \sum_{i>0} (u_{i}v_{i})^{4} }{ \displaystyle
  !>          \left( \sum_{i>0} u_{i}^{2}v_{i}^{2} \right)^{2}
  !>          - \sum_{i>0} (u_{i}v_{i})^{4} }.
  !>    \f]
  !=======================================================================
  Subroutine tracesln
    Implicit None
    Integer(ipr) :: iw,it,ib,k1,k2,kkk,k
    Real(pr) :: AAV,SNtor,SDtor
    Real(pr) :: S_U1V1,S_U1V3,S_U2V2,S_U3V1,S_U4V4
    Real(pr) :: U_ACTU,U_ACTU2,U_ACTU3,U_ACTU4
    Real(pr) :: V_ACTU,V_ACTU2,V_ACTU3,V_ACTU4
    !
    etr=zero
    Do it=itmin,itmax
       S_U1V1=ZERO; S_U1V3=ZERO; S_U2V2=ZERO; S_U3V1=ZERO; S_U4V4=ZERO
       Do ib=1,nb
          k1=ka(ib,it)+1; k2=ka(ib,it)+kd(ib,it)
          If(k1.Le.k2) Then
             kkk=lcanon(ib-1,it)
             Do k=1,id(ib)
                kkk=kkk+1; aav=vk(kkk,it)               ! v^2
                U_ACTU=Sqrt(AAV);       U_ACTU2=U_ACTU*U_ACTU
                U_ACTU3=U_ACTU2*U_ACTU; U_ACTU4=U_ACTU2*U_ACTU2
                V_ACTU=Sqrt(ONE-AAV);   V_ACTU2=V_ACTU*V_ACTU
                V_ACTU3=V_ACTU2*V_ACTU; V_ACTU4=V_ACTU2*V_ACTU2
                S_U1V1=S_U1V1+U_ACTU  * V_ACTU
                S_U1V3=S_U1V3+U_ACTU  * V_ACTU3
                S_U2V2=S_U2V2+U_ACTU2 * V_ACTU2     !Tr r (1-r)
                S_U3V1=S_U3V1+U_ACTU3 * V_ACTU
                S_U4V4=S_U4V4+U_ACTU4 * V_ACTU4     !Tr (1-r)^2 r^2
             End Do
          End If
       End Do !ib
       SNtor=8.0_pr*(S_U3V1*S_U1V3-S_U4V4)
       SDtor=32.0_pr*(S_U2V2*S_U2V2-S_U4V4)
       Geff(it)=del(it)**2/ept(it)
       ala2(it)=-Geff(it)*(SNtor/SDtor)
       If(ala2(it).Ge.10.0_pr) ala2(it)=4.0_pr    ! ala2 goes to hell
       etr(it)=-four*ala2(it)*S_U2V2              ! to total energy
    End Do !it
    etr(3)=etr(1)+etr(2)                          !to total energy
    Do iw=lout,lfile
       Write(iw,'(26x,a,2(1x,f7.3),a,3(1x,f9.3),a,2(1x,f7.3))')  &
            '  f: ala2(n,p)=',ala2,' #eln(n,p,t)=',etr,' #del+ala2=',del+ala2
    End Do
  End Subroutine tracesln
  !=======================================================================
  !> Computes the Lipkin-Nogami parameter \f$ \lambda_2 \f$ for a monopole
  !> pairing force of the type
  !>    \f[
  !>         \bar{v}_{ijkl} = -G \delta_{j\bar{i}}\delta_{l\bar{k}}
  !>                           \text{sign}(j)\text{sign}(l),
  !>    \f]
  !> The calculation is performed in the qp basis
  !>    \f[
  !>       \lambda_{2} = \frac{G}{4}
  !>         \frac{\text{Tr}^{<} \kappa^{*}\rho\;\text{Tr}^{<} (1-\rho)\kappa
  !>         -\sum_{ij} [ \rho(1-\rho) ]_{\bar{i}\bar{j}} [\rho(1-\rho)]_{ij}
  !>         }{ \displaystyle
  !>         \left( \text{Tr}^{>} \rho(1 - \rho) \right)^{2}
  !>         - \text{Tr}^{>} \rho^{2}(1 - \rho)^{2}.}
  !>    \f]
  !=======================================================================
  Subroutine tracesln_qp
    Implicit None
    Integer(ipr) :: iw,nd,ib,i1,i2,n2,ibitnb,it,i1n2nd,i2n2nd,i1i2nd
    Real(pr) :: frit,frit2,ftit
    Real(pr) :: etr2(2),trk(2),trk1(2),SNtor(2),SDtor(2),Sum(2)
    !
    ! initialization
    etr=zero; etr2=zero; trk=zero; trk1=zero
    ! loop over the blocks
    Do ib=1,nb
       nd=id(ib)
       ! Traces for neutrons and protons
       Do i2=1,nd           ! index alpha
          Do i1=i2,nd         ! index beta.ge.alpha
             sum=zero
             Do n2=1,nd
                i1n2nd=Max(i1,n2)+(Min(i1,n2)-1)*nd
                i2n2nd=Max(i2,n2)+(Min(i2,n2)-1)*nd
                Do it=itmin,itmax
                   ibitnb=ib+(it-1)*nbx
                   Sum(it)=Sum(it)+rk(i1n2nd,ibitnb)*rk(i2n2nd,ibitnb)*p14
                End Do !it
             End Do !n2
             i1i2nd=i1+(i2-1)*nd
             Do it=itmin,itmax
                ibitnb=ib+(it-1)*nbx
                frit=rk(i1i2nd,ibitnb)*half
                ftit=ak(i1i2nd,ibitnb)
                frit2=Sum(it)
                If(i1.Eq.i2) Then
                   etr(it)=etr(it)+frit-frit**2                  ! Tr r (1-r)
                   etr2(it)=etr2(it)+(one-two*frit+frit2)*frit2  ! Tr (1-r)^2 r^2
                   trk(it)=trk(it)+frit*ftit                     ! Tr r k
                   trk1(it)=trk1(it)+ftit -frit*ftit             ! Tr k (1-r)
                Else
                   etr(it)=etr(it) -two*frit**2                  ! Tr r (1-r)
                   etr2(it)=etr2(it)+two*(-two*frit+frit2)*frit2 ! Tr (1-r)^2 r^2
                   trk(it)=trk(it)+two*frit*ftit                 ! Tr r k
                   trk1(it)=trk1(it)-two*ftit*frit               ! Tr k (1-r)
                End If
             End Do !it
          End Do !i1
       End Do !i2
    End Do !ib
    ! total traces
    Do it=itmin,itmax
       SNtor(it)=8.0_pr*(trk1(it)*trk(it)-etr2(it))
       SDtor(it)=32.0_pr*(etr(it)**2     -etr2(it))
       Geff(it)=del(it)**2/ept(it)
       ala2(it)=-( SNtor(it)/SDtor(it) )*Geff(it)
       If(ala2(it).Ge.10.0_pr) ala2(it)=4.0_pr  ! in case ala2 goes to hell
       etr(it)=-four*ala2(it)*etr(it)           ! to total energy
    End Do
    etr(3)=etr(1)+etr(2)         !to total energy
    Do iw=lout,lfile
       Write(iw,'(26x,a,2(1x,f7.3),a,3(1x,f9.3),a,2(1x,f7.3))')  &
            '  #LN: ala2(n,p)=',ala2,' #eln(n,p,t)=',etr,' #del+ala2=',del+ala2
    End Do
  End Subroutine tracesln_qp
  !=======================================================================
  !> Calculates the densities in r-space at gauss-mesh points
  !> corrected due to Lipkin-Nogami
  !=======================================================================
  Subroutine densitln
    Implicit None
    Integer(ipr) :: ih,il,ib,nd,i0,i01,i02,n1,n2,nza,nzb,nra,nrb,nla,  &
         nlb,nsa,nsb,it,ml,ihli,k,k0(2),k00(2),k1,k2
    Real(pr) :: fr(2),vvs,vvc,ssln1,ssln2,ssln3,vks
    Real(pr) :: qla,qlb,qlab,qha,qhb,qhlab,qhab,sro
    !
    k0=0; k00=0
    ro=zero
    ! loop over the blocks
    Do ib=1,nb
       k00=k0
       nd=id(ib); i0=ia(ib)
       Do n2=1,nd
          i02=i0+n2; nzb=nz(i02); nrb=nr(i02);
          nlb=nl(i02); nsb=ns(i02)
          Do n1=1,n2
             i01=i0+n1; nza=nz(i01); nra=nr(i01)
             nla=nl(i01); nsa=ns(i01)
             k0=k00
             Do it=itmin,itmax
                k1=ka(ib,it)+1
                k2=ka(ib,it)+kd(ib,it)
                fr(it)=zero
                If(k1.Le.k2) Then
                   Do k=1,nd
                      k0(it)=k0(it)+1
                      ssln1=ssln(1,it)
                      ssln2=ssln(2,it)
                      ssln3=ssln(3,it)
                      vks=vk(k0(it),it)
                      vvc=vks
                      vvs=Abs(one-vks)
                      If(vvs.Ge.1.0d-40) Then
                         vvs=two*Sqrt(vks*vvs)   !2vu
                         vvc=vks+vvs**2*p14*ssln1*((two*vks-one)*ssln1-ssln2)/ssln3
                      End If
                      fr(it)=fr(it)+two*ddc(n2,k0(it),it)*ddc(n1,k0(it),it)*vvc
                   End Do
                   If (n1.Ne.n2) Then
                      fr(it)=two*fr(it)
                   End If
                End If
             End Do
             !---diagonal in spin
             If (nsa.Eq.nsb) Then
                ml=nla
                Do il=1,ngl
                   qla=ql (nra,ml,il);    qlb=ql (nrb,ml,il)
                   qlab=qla*qlb
                   Do ih=1,ngh
                      ihli=ih+(il-1)*ngh
                      qha=qh (nza,ih); qhb=qh (nzb,ih)
                      qhab=qha*qhb
                      qhlab=qhab*qlab; sro=qhlab
                      ro(ihli,:)=ro(ihli,:)+fr(:)*sro
                   End Do   !ih
                End Do  !il
             End If
          End Do !n2
       End Do !n1
    End Do !ib
    ! set the THO weights
    Do ihli=1,nghl
       ro(ihli,:)=ro(ihli,:)*wdcori(ihli)
    End Do
  End Subroutine densitln
  !=======================================================================
  !> Coulomb field (direct part) Vautherin prescription
  !> Ref.: Phys. Rev. C 7, 296 (1973)
  !=======================================================================
  Subroutine coulom1
    Use EllipticIntegral
    Implicit None
    Integer(ipr), Save :: i,k
    Real(pr) :: zd2,rhl,y1,y2,xx1,xx2,s1,s2,e1,e2,vik,f,r,r1,r4,  &
                rr2,z,z1,zd1,x1,x2,fac1,fac2
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('coulom1',0)
    !
    If(icacou.Eq.0) Then
       !
       icacou=1
       !
       ! For parity-breaking shapes, the Coulomb potential was incorrectly
       ! calculated by assuming the two intervals [0,+\infty[ and ]-infty,0]
       ! were equivalent (see also routine coulom() below). This bug was
       ! corrected in version 200d
       If(Parity) Then
          fac1 = one;  fac2 = one
       Else
          fac1 = zero; fac2 = two
       End If
       !
       f=half*chargee2/pi
       ! See notes in subroutine coulom for explanations about some numerical
       ! factors apparently missing here.
!$OMP PARALLEL DO        &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(nghl,fl,fh,fac1,fac2,wdcor,vc,f) &
!$OMP& PRIVATE(i,r,z,r4,k,r1,z1,rr2,rhl,zd1,y1,xx1,s1,zd2,y2,xx2,s2,vik)
       Do i=1,nghl
          r=fl(i); z=fh(i)
          r4=four*r
          Do k=1,i
             r1=fl(k); z1=fh(k)
             rhl=r4*r1     ! 4 r r'
             rr2=(r+r1)**2 ! (r+r')^2
             ! z>0 part
             zd1=(z-z1)**2 ! (z-z')^2
             y1=zd1+rr2    ! d(r,z) = (r+r')^2 + (z-z')^2
             xx1=rhl/y1    ! 4 r r' / d(r,z)
             s1=Sqrt(y1)   ! sqrt(d(r,z))
             ! z<0 part
             zd2=(z+z1)**2
             y2=zd2+rr2
             xx2=rhl/y2
             s2=Sqrt(y2)
             !
             vik = f*fac2*(s1*CompleteEllipticFunction_2nd(xx1) &
                          +s2*CompleteEllipticFunction_2nd(xx2)*fac1)
             !
             vc(i,k)=vik*wdcor(k)  !wdcor=pi*wh*wl*bz*bp*bp
             vc(k,i)=vik*wdcor(i)  !wdcor=pi*wh*wl*bz*bp*bp
             !
          End Do  !k
       End Do  !i
!$OMP End Parallel Do
    End If
    ! Calculation of the coulomb field (each iteration)
    cou=zero
    Call dgemm('n','n',nghl,1,nghl,1.0_pr,vc,nghl,dro(:,2),nghl,0.0_pr,cou,nghl)
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('coulom1',1)
    !
  End Subroutine coulom1
  !=======================================================================
  !> Coulomb field (direct part), Gogny prescription
  !> Ref.: Phys. Rev. C 27, 2317 (1983)
  !=======================================================================
  Subroutine coulom
    Use bessik
    Implicit None
    Integer(ipr), Save :: i,j,k
    Real(pr), Save :: zd2,y1,y2,xx1,s1,vik,f,r,r1,fac1,fac2,rr2,z,z1,zd1,t,  &
                      bb,r2,r12,rrr,rz1,rz2,rrz1,rrz2,xx,rk1,rip1,rkp1,alpha,&
                      beta,xxx
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('coulom',0)
    !
    If(icacou.Eq.0) Then
       !
       icacou=1
       !
       ! For parity-breaking shapes, the Coulomb potential was incorrectly
       ! calculated by assuming the two intervals [0,+\infty[ and ]-infty,0]
       ! were equivalent (see also below). This bug was corrected in version
       ! 139a
       If(Parity) Then
          fac1 = one;  fac2 = one
       Else
          fac1 = zero; fac2 = two
       End If
       ! Notes:
       !   - Missing factor 2 compared to Eq. (58) CPC paper because the density
       !     ro(:,it) already contains it (see routine DENSIT) due to T-invariance
       !   - Missing factor 1/2 when applying Gauss-Legendre quadrature (from [0,1]
       !     to the proper [-1,1] interval because it will be put back in subroutine
       !     expect() and is cancelled by a factor 2 in the HF field
       !   - For conserved parity, Gauss-Hermite points are all positive, the full
       !     integral over z' is split in z'<0 and z'>0, values of z and z1 below
       !     refer to the absolute values of z' (=-z' if z'<0)
       !
       bb=50.0_pr          ! Length scale L
       beta=2.00_pr
       alpha=one/beta
       f=chargee2/Sqrt(pi) ! e^2/Sqrt(pi)
       !
!$OMP PARALLEL DO        &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(nghl,fl,fh,nleg,xleg,bb,fac1,fac2,wleg,wdcor,vc,f,alpha,beta) &
!$OMP& PRIVATE(i,r,z,k,r1,z1,rrr,rr2,zd1,zd2,rz1,rz2,rrz1,rrz2, &
!$OMP&         xx1,j,xx,y1,s1,t,y2,vik,xxx)
       Do i=1,nghl
          r = fl(i); z = fh(i)
          Do k=1,i
             !
             r1 = fl(k); z1 = fh(k)
             rrr = two*r*r1; rr2 = (r - r1)**2
             ! z>0 part
             zd1 = (z - z1)**2
             rz1 = rr2 + zd1
             ! z<0 part
             zd2 = (z + z1)**2
             rz2 = rr2 + zd2
             ! Gauss-Legendre integration over u from 0 to D
             xx1=zero
             Do j=1,nleg
                xx=(one-xleg(j)**beta)**alpha ! change of variable to 0 <= u <= 1
                xxx=(one-xleg(j)**beta)**(alpha+one)
                y1=(xleg(j)/(bb*xx))**2 ! u^2
                s1=y1*rrr               ! 2 u^2 r r'
                y2=besei0(s1)           ! I0( 2 u^2 r r' ) * exp(-2 u^2 r r')
                xx1=xx1+fac2*wleg(j)*y2*(Exp(-rz1*y1) + fac1*Exp(-rz2*y1)) / xxx
             End Do
             vik=f*xx1/bb
             !
             vc(i,k)=vik*wdcor(k)  !wdcor=pi*wh*wl*bz*bp*bp
             vc(k,i)=vik*wdcor(i)  !wdcor=pi*wh*wl*bz*bp*bp
             !
          End Do  !k
       End Do  !i
!$OMP End Parallel Do
       !
    End If

    ! Calculation of the Coulomb field
    cou=zero
    Call dgemm('n','n',nghl,1,nghl,1.0_pr,vc,nghl,ro(:,2),nghl,0.0_pr,cou,nghl)
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('coulom',1)
    !
  End Subroutine coulom
  !=======================================================================
  ! Coulomb field (direct part), Gogny prescription
  ! Ref.: Phys. Rev. C 27, 2317 (1983)
  !=======================================================================
  Subroutine coulom_test
    Use bessik
    Implicit None
    Integer(ipr), Save :: i,j,k
    Real(pr), Save :: zd2,y1,y2,xx1,s1,vik,f,r,r1,fac1,fac2,rr2,z,z1,zd1,t,  &
                      bb,r2,r12,rrr,rz1,rz2,rrz1,rrz2,xx,rk1,rip1,rkp1,alpha,&
                      beta,xxx,func
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('coulom_test',0)
    !
    !
    ! For parity-breaking shapes, the Coulomb potential was incorrectly
    ! calculated by assuming the two intervals [0,+\infty[ and ]-infty,0]
    ! were equivalent (see also below). This bug was corrected in version
    ! 139a
    If(Parity) Then
       fac1 = one;  fac2 = one
    Else
       fac1 = zero; fac2 = two
    End If
    !
    bb=5.0_pr          ! Length scale L
    beta=2.00_pr
    alpha=one/beta
    !f=chargee2/Sqrt(pi) ! e^2/Sqrt(pi)
    f=one/Sqrt(pi)       ! 1/Sqrt(pi)
    !
    Do j=1,nleg
       ! Gauss-Legendre integration over u from 0 to D
       xx=(one-xleg(j)**beta)**alpha ! change of variable to 0 <= u <= 1
       xxx=(one-xleg(j)**beta)**(alpha+one)
       !
       func=zero
       Do i=1,nghl
          r = fl(i); z = fh(i)
          Do k=1,i
             !
             r1 = fl(k); z1 = fh(k)
             rrr = two*r*r1; rr2 = (r - r1)**2
             ! z>0 part
             zd1 = (z - z1)**2
             rz1 = rr2 + zd1
             ! z<0 part
             zd2 = (z + z1)**2
             rz2 = rr2 + zd2
             y1=(xleg(j)/(bb*xx))**2 ! u^2
             s1=y1*rrr               ! 2 u^2 r r'
             y2=besei0(s1)           ! I0( 2 u^2 r r' ) * exp(-2 u^2 r r')
             xx1=fac2*wleg(j)*y2*(Exp(-rz1*y1) + fac1*Exp(-rz2*y1)) / xxx
             vik=f*xx1/bb
             !
             func=func+vik*wdcor(k)*ro(k,2)*wdcor(i)*ro(i,2)  !wdcor=pi*wh*wl*bz*bp*bp
             !
          End Do ! k
       End Do  ! i
       Write(6,'(2f30.14)') xleg(j),func
       !
    End Do  !j
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('coulom_test',1)
    !
  End Subroutine coulom_test
  !=======================================================================
  ! Hartree-field (direct part)
  !=======================================================================
  Subroutine HartreeDir
    Implicit None
    Integer(ipr) :: i,j,k
    Real(pr) :: xx1,vik00,vik01,vik11
    Real(pr) :: r,rr,rrr,r1,r2,rr1,rr2
    Real(pr) :: z,z1,zdm,zdp,rzm,rzp
    Real(pr), Allocatable :: u(:)
    If(icahartree.Eq.0) Then
       icahartree=1
       !
       If(Allocated(u)) Deallocate(u); Allocate(u(nleg))
       u=Cos(HALF*Pi*xleg)
       Do i=1,nghl
          r=fl(i); z=fh(i); rr=r*r
          Do k=1,i
             r1=fl(k); z1=fh(k); rr1=r1*r1;     rr2=two*r*r1; rrr=rr+rr1;
             zdm=(z-z1)**2;      zdp=(z+z1)**2; rzm=rrr+zdm;  rzp=rrr+zdp
             vik00=0.250_pr*Sum(wleg*( &
                  + HartreeV00(Sqrt(rzp-rr2*u)) &
                  + HartreeV00(Sqrt(rzm-rr2*u)) &
                  + HartreeV00(Sqrt(rzp+rr2*u)) &
                  + HartreeV00(Sqrt(rzm+rr2*u))))
             vik01=0.250_pr*Sum(wleg*( &
                  + HartreeV01(Sqrt(rzp-rr2*u)) &
                  + HartreeV01(Sqrt(rzm-rr2*u)) &
                  + HartreeV01(Sqrt(rzp+rr2*u)) &
                  + HartreeV01(Sqrt(rzm+rr2*u))))
             vik11=0.250_pr*Sum(wleg*( &
                  + HartreeV11(Sqrt(rzp-rr2*u)) &
                  + HartreeV11(Sqrt(rzm-rr2*u)) &
                  + HartreeV11(Sqrt(rzp+rr2*u)) &
                  + HartreeV11(Sqrt(rzm+rr2*u))))

             vhart00(i,k)=vik00*wdcor(k)         ! wdcor=pi*wh*wl*bz*bp*bp/fd
             vhart00(k,i)=vik00*wdcor(i)
             vhart01(i,k)=vik01*wdcor(k)
             vhart01(k,i)=vik01*wdcor(i)
             vhart11(i,k)=vik11*wdcor(k)
             vhart11(k,i)=vik11*wdcor(i)
          End Do  !k
       End Do  !i
       Deallocate(u)
    End If
    ! calculation of the Hartree field
    vDHartree=0.0_pr
    Do i=1,nghl
       vDHartree(:,1)=vDHartree(:,1)+vhart00(:,i)*(ro(i,1)+ro(i,2))+vhart01(:,i)*(ro(i,1)-ro(i,2))
       vDHartree(:,2)=vDHartree(:,2)+vhart11(:,i)*(ro(i,1)-ro(i,2))+vhart01(:,i)*(ro(i,1)+ro(i,2))
    End Do
  End Subroutine HartreeDir
  !=======================================================================
  !> Optimization arrays
  !>    FI2D_opt(JA,ihil) == Laplacian(r,z) HOwf
  !>    FID2D-xlamy2*FID  == Laplacian(r,z,phy) FID
  !>    FIU2D-xlapy2*FIU  == Laplacian(r,z,phy) FIU
  !=======================================================================
  Subroutine optHFBTHO
    Implicit None
    Integer(ipr) :: i,ih,il,ib,ibx,nd,nza,nra,nla,nsa
    Integer(ipr) :: ihil,laplus,im,JA,N1,N2,ndnd,n12,n21
    Real(pr)    :: qla,v2,v4,yi,y,y2,qha,qhla,xmi,u,u2,un,up,xxx
    Real(pr)    :: sml2,cnzaa,cnraa,a,b
    Real(pr)    :: FITW1,FITW2,FITW3,FITW4
    Real(pr)    :: fi1r,fi1z,fi2d,QHL1A,QH1LA,vh,vdh,vsh,hbh
    Real(pr)    :: SRFIh,SFIRh,SFIZh,SZFIh,SNABLARh,SNABLAZh
    Real(pr)    :: xlam,xlam2,xlamy,xlamy2,xlap,xlap2,xlapy,xlapy2,XLAMPY
    Real(pr)    :: bpi,bpi2,bzi,bzi2,xh2
    !
    bpi=one/bp; bpi2=bpi*bpi; bzi=one/bz; bzi2=bzi*bzi
    !
    !-----------------------------------------
    ! Allocate the optimization arrays
    !-----------------------------------------
    If(Allocated(QHLA_opt)) Deallocate(QHLA_opt,FI1R_opt,FI1Z_opt,FI2D_opt,y_opt)
    Allocate(QHLA_opt(ntx,nghl),FI1R_opt(ntx,nghl),FI1Z_opt(ntx,nghl),FI2D_opt(ntx,nghl),y_opt(nghl))
    !----------------------------------------------
    ! START BLOCKS
    !----------------------------------------------
    Do ib=1,NB
       ND=ID(ib); IM=ia(ib)
       If(Parity) Then
          LAPLUS=(ib+1)/2 !Yesp
       Else
          LAPLUS=ib       !Nop
       End If
       XLAP=LAPLUS; XLAM=XLAP-ONE; xlap2=xlap*xlap; xlam2=xlam*xlam
       !----------------------------------------------
       ! SUM OVER GAUSS INTEGRATION POINTS
       !----------------------------------------------
       Do IL=1,ngl
          v2=half/xl(il); v4=v2*v2
          Do IH=1,ngh
             ihil=ih+(il-1)*ngh; xh2=xh(ih)**2
             If(iLST1.Eq.0) Then
                ! HO-basis
                yi=Sqrt(xl(il))*bp; y=one/yi; y2=y*y
                xlamy=xlam*y; xlamy2=xlam2*y2; xlapy=xlap*y; xlapy2=xlap2*y2; XLAMPY=XLAMY+XLAPY
             Else
                ! THO-basis
                y=fli(ihil); y2=y*y; xlamy=xlam*y; u=xh(ih); u2=u*u;
                xlamy2=xlam2*y2; xlapy=xlap*y; xlapy2=xlap2*y2; XLAMPY=XLAMY+XLAPY
             End If
             y_opt(ihil)=y
             !----------------------------------------------
             ! SCAN OVER BASIS STATES
             !----------------------------------------------
             Do N1=1,ND
                JA=N1+IM; NLA=NL(JA); NRA=NR(JA); NZA=NZ(JA); NSA=NS(JA)
                SML2=NLA*NLA; CNZAA=NZA+NZA+1; CNRAA=NRA+NRA+NLA+1
                QHA=QH(NZA,IH); QLA=QL(NRA,NLA,IL); QHLA=QHA*QLA
                QHL1A=QHA*QL1(NRA,NLA,IL)*V2; QH1LA=QH1(NZA,IH)*QLA
                If(iLST1.Eq.0) Then
                   ! HO-basis
                   FI1R=(two*Sqrt(xl(il))*bpi)*QHL1A
                   FI1Z=bzi*QH1LA
                   FI2D=((xh2-CNZAA)*bzi2+four*(p14-CNRAA*V2+SML2*V4)*xl(il)*bpi2 )*QHLA
                Else
                   ! THO-basis
                   u=xh(ih); u2=u*u;
                   FI1R=FP4(IHIL)*QHLA+FP5(IHIL)*QH1LA+FP6(IHIL)*QHL1A
                   FI1Z=FP1(IHIL)*QHLA+FP2(IHIL)*QH1LA+FP3(IHIL)*QHL1A
                   FI2D=(FS1(IHIL)*QH1LA*QH1LA+FS2(IHIL)*QHL1A*QHL1A  &
                        +FOUR*FS4(IHIL)*QH1LA*QHL1A  &
                        +TWO*(FS5(IHIL)*QH1LA+FS6(IHIL)*QHL1A)*QHLA  &
                        +((U2-CNZAA)*FS1(IHIL)+(p14-CNRAA*V2+SML2*V4)*FS2(IHIL)  &
                        +FS3(IHIL))*QHLA*QHLA-TWO*(FI1R*FI1R+FI1Z*FI1Z))/(TWO*QHLA)
                End If
                QHLA_opt(JA,ihil)=QHLA; FI2D_opt(JA,ihil)=FI2D; FI1R_opt(JA,ihil)=FI1R; FI1Z_opt(JA,ihil)=FI1Z
             End Do !N1
             !
          End Do !IH
       End Do !IL
    End Do !IB
  End Subroutine optHFBTHO
  !=======================================================================
  !> Computes all Skyrme densities at Gauss quadrature points
  !=========================================================================
  Subroutine DENSIT
    Implicit None
    Integer(ipr) :: iw,nsa,nza,nra,nla,k,i,nd,il,ih,ihil,laplus,ii
    Integer(ipr) :: imen,ib,im,it,J,JJ,JA,JN,k0,k1,k2,ibiblo
    Integer(ipr) :: bb,size,ndxmax
    Parameter(ndxmax=(n00max+2)*(n00max+2)/4)
    Real(pr)     :: s,ss,sd,yi,y,y2,sml2,cnzaa,cnraa,u,u2,v2,v4
    Real(pr)     :: anik,pnik,qhla,qh1la,qhl1a,qla,qha,fi1r,fi1z,fi2d,fidd
    Real(pr)     :: xlam,xlam2,xlamy,xlamy2,xlap,xlap2,xlapy,xlapy2,XLAMPY
    Real(pr)     :: TFIU,TFID,TFIUR,TFIDR,TFIUZ,TFIDZ,TFIUD2,TFIDD2
    Real(pr)     :: TPFIU,TPFID,TPFIUR,TPFIDR,TPFIUZ,TPFIDZ,TPFIUD2,TPFIDD2
    Real(pr)     :: PIU,PIUZ,PIUR,PIUD2,PID,PIDZ,PIDR,PIDD2
    Real(pr)     :: TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,TEMP6,TEMP7,TEMP8,TEMP9,TEMP10,TEMP11,TW_T,PW_T,WGT(nghl)
    Real(pr)     :: Takaihil,Troihil,Tdjihil,Ttauihil,Tdroihil,TSRFIihil
    Real(pr)     :: TSFIRihil,TSFIZihil,TSZFIihil,TNABLARIHIL,TNABLAZIHIL
    Real(pr), Pointer :: TAKA(:),TRO(:),TDJ(:),TTAU(:),TDRO(:)
    Real(pr), Pointer :: TSRFI(:),TSFIR(:),TSFIZ(:),TSZFI(:),TNABLAR(:),TNABLAZ(:)
    Real(pr)     :: time1,time2,fk,f1k
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
#if(USE_OPENMP==1)
    Real(pr), Allocatable :: OMPTAKA(:,:,:),OMPTRO(:,:,:),OMPTDJ(:,:,:),OMPTTAU(:,:,:),OMPTDRO(:,:,:)
    Real(pr), Allocatable :: OMPTSRFI(:,:,:),OMPTSFIR(:,:,:),OMPTSFIZ(:,:,:),OMPTSZFI(:,:,:), &
                             OMPTSZIF(:,:,:),OMPTNABLAR(:,:,:),OMPTNABLAZ(:,:,:)
#else
    Real(pr) :: OMPTAKA(nghl,2),OMPTRO(nghl,2),OMPTDJ(nghl,2),OMPTTAU(nghl,2),OMPTDRO(nghl,2)
    Real(pr) :: OMPTSRFI(nghl,2),OMPTSFIR(nghl,2),OMPTSFIZ(nghl,2),OMPTSZFI(nghl,2), &
                OMPTSZIF(nghl,2),OMPTNABLAR(nghl,2),OMPTNABLAZ(nghl,2)
#endif
    !
    Real(pr) :: OMPFIU(ndxmax),OMPFID(ndxmax),OMPFIUR(ndxmax),OMPFIDR(ndxmax),OMPFIUZ(ndxmax)
    Real(pr) :: OMPFIDZ(ndxmax),OMPFIUD2N(ndxmax),OMPFIDD2N(ndxmax)
    !
    Real(pr) :: OMPPFIU(ndxmax),OMPPFID(ndxmax),OMPPFIUR(ndxmax),OMPPFIDR(ndxmax),OMPPFIUZ(ndxmax)
    Real(pr) :: OMPPFIDZ(ndxmax),OMPPFIUD2N(ndxmax),OMPPFIDD2N(ndxmax)
    !
    Real(pr) :: OMPAN(ndxmax*ndxmax),OMPANK(ndxmax*ndxmax)
    Real(pr) :: f_T(ndxmax),f1_T(ndxmax)
    Real(pr) :: dnrm2
    external dnrm2
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('densit',0)
    !
    !-----------------------------------------------
    ! ZERO N & P DENSITIES
    !-----------------------------------------------
    RO=ZERO; TAU=ZERO; DJ=ZERO; DRO=ZERO; AKA=ZERO; SZFI=ZERO; SFIZ=ZERO
    SRFI=ZERO; SFIR=ZERO; NABLAR=ZERO; NABLAZ=ZERO; VARMAS=ZERO
    !
#if(USE_OPENMP==1)
    Allocate(OMPTAKA(1:NB,1:nghl,2),OMPTRO(1:NB,1:nghl,2),OMPTDJ(1:NB,1:nghl,2),&
             OMPTTAU(1:NB,1:nghl,2),OMPTDRO(1:NB,1:nghl,2))
    Allocate(OMPTSRFI(1:NB,1:nghl,2),OMPTSFIR(1:NB,1:nghl,2),OMPTSFIZ(1:NB,1:nghl,2), &
             OMPTSZFI(1:NB,1:nghl,2),OMPTSZIF(1:NB,1:nghl,2),OMPTNABLAR(1:NB,1:nghl,2),&
             OMPTNABLAZ(1:NB,1:nghl,2))
#endif
    OMPTAKA   = ZERO; OMPTRO     = ZERO; OMPTDJ     = ZERO; OMPTTAU  = ZERO
    OMPTDRO   = ZERO; OMPTSRFI   = ZERO; OMPTSFIR   = ZERO; OMPTSFIZ = ZERO
    OMPTSZFI  = ZERO; OMPTNABLAR = ZERO; OMPTNABLAZ = ZERO;
    !
#if(USE_OPENMP==1)
!$OMP PARALLEL DO        &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(NB,npr_INI,ID,ia,nghl,NS,Parity,bloblo,keyblo,blo123d,ka,kd,KpwiN,RVqpN,RUqpN,&
!$OMP&        fn_T,KpwiP,RVqpP,RUqpP,fp_T,y_opt,QHLA_opt,FI2D_opt,FI1Z_opt,FI1R_opt,&
!$OMP&        switch_on_temperature,OMPTaka,OMPTro,OMPTdj,OMPTtau,OMPTdro,OMPTSRFI,&
!$OMP&        OMPTSFIR,OMPTSFIZ,OMPTSZFI,OMPTNABLAR,OMPTNABLAZ) &
!$OMP& PRIVATE(bb,it,ib,ND,IM,LAPLUS,XLAP,XLAM,xlap2,xlam2,ibiblo,K0,k1,k2,imen,&
!$OMP&         ompan,ompank,f_T,f1_T,J,JJ,K,ihil,y,xlamy,xlapy,XLAMPY,y2,xlamy2,&
!$OMP&         xlapy2,OMPFIU,OMPFIUZ,OMPFIUR,OMPFID,OMPFIDZ,OMPFIDR,OMPFIUD2N,OMPFIDD2N,&
!$OMP&         OMPPFIU,OMPPFIUZ,OMPPFIUR,OMPPFID,OMPPFIDZ,OMPPFIDR,OMPPFIUD2N,OMPPFIDD2N,&
!$OMP&         PIU,PIUZ,PIUR,PIUD2,PID,PIDZ,PIDR,PIDD2,JN,JA,NSA,QHLA,PNIK,FI2D,FI1Z,FI1R,&
!$OMP&         Takaihil,Troihil,Tdjihil,Ttauihil,Tdroihil,TSRFIihil,&
!$OMP&         TFIU,TFID,TFIUR,TFIDR,TFIUZ,TFIDZ,TFIUD2,TFIDD2,fk,f1k,&
!$OMP&         TPFIU,TPFID,TPFIUR,TPFIDR,TPFIUZ,TPFIDZ,TPFIUD2,TPFIDD2,&
!$OMP&         TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,TEMP6,TEMP7,TEMP8,TEMP9,TEMP10,TEMP11,TW_T,PW_T,&
!$OMP&         TSFIRihil,TSFIZihil,TSZFIihil,TNABLARIHIL,TNABLAZIHIL)
#endif
    Do bb=0,2*NB-1
       it = bb/NB + 1
       ib = Mod(bb,NB)+1
       !
       ! case of zero particle number, only flush densities
       If((npr_INI(1).Eq.0).And.(it.Eq.1)) Cycle
       If((npr_INI(2).Eq.0).And.(it.Eq.2)) Cycle
       !-----------------------------------------------
       ! SCAN OVER BLOCKS
       !-----------------------------------------------
       ND=ID(ib); IM=ia(ib)
       If(Parity) Then
          LAPLUS=(ib+1)/2 !Yesp
       Else
          LAPLUS=ib       !Nop
       End If
       XLAP=LAPLUS; XLAM=XLAP-ONE; xlap2=xlap*xlap; xlam2=xlam*xlam
       !
       ! blocking
       ibiblo=bloblo(keyblo(it),it)
       K0=0; If(ibiblo.Eq.ib) K0=blo123d(it)
       !
       !----------------------------------------------
       ! PAIRING WINDOW QP WAVE FUNCTIONS
       !----------------------------------------------
       k1=ka(ib,it)+1; k2=ka(ib,it)+kd(ib,it); imen=k2-k1+1
       If(IMEN.Gt.0) Then
          ompan=ZERO; ompank=ZERO; f_T=ZERO; f1_T=ZERO
          J=0
          If(it.Eq.1) then
             Do JJ=1,nd ! basis
                Do K=K1,K2 ! qp
                   J=J+1; I=KpwiN(K)+JJ; ompan(J)=RVqpN(I); ompank(J)=RUqpN(I)
                End Do
             End Do
             J=0
             Do K=K1,K2
                J=J+1;JJ=K !KpwiN(K)
                f_T(J)=one-fn_T(JJ);f1_T(J)=fn_T(JJ)
             End Do
          Else
             Do JJ=1,nd ! basis
                Do K=K1,K2 ! qp
                   J=J+1; I=KpwiP(K)+JJ; ompan(J)=RVqpP(I); ompank(J)=RUqpP(I)
                End Do
             End Do
             J=0
             Do K=K1,K2
                J=J+1;JJ=K !KpwiP(K)
                f_T(J)=one-fp_T(JJ);f1_T(J)=fp_T(JJ)
             End Do
          End If
          !-----------------------------------------------
          ! SCAN OVER GAUSS INTEGRATION POINTS
          !-----------------------------------------------
          Do ihil=1,nghl
             y=y_opt(ihil); xlamy =xlam*y;    xlapy =xlap*y;   XLAMPY=XLAMY+XLAPY
             y2=y*y;        xlamy2=xlam2*y2;  xlapy2=xlap2*y2
             Do K=1,IMEN
                ! V_k components
                OMPFIU(K)    = ZERO; OMPFIUZ(K)   = ZERO; OMPFIUR(K) = ZERO
                OMPFID(K)    = ZERO; OMPFIDZ(K)   = ZERO; OMPFIDR(K) = ZERO
                OMPFIUD2N(K) = ZERO; OMPFIDD2N(K) = ZERO;
                ! U_k components
                OMPPFIU(K)    = ZERO; OMPPFIUZ(K)   = ZERO; OMPPFIUR(K) = ZERO
                OMPPFID(K)    = ZERO; OMPPFIDZ(K)   = ZERO; OMPPFIDR(K) = ZERO
                OMPPFIUD2N(K) = ZERO; OMPPFIDD2N(K) = ZERO;
             End Do
             If(K0.Ne.0) Then
                PIU=ZERO;  PIUZ=ZERO; PIUR=ZERO; PIUD2=ZERO
                PID=ZERO;  PIDZ=ZERO; PIDR=ZERO; PIDD2=ZERO
             End If
             !-----------------------------------------------
             ! SUM OVER BASIS STATES
             !-----------------------------------------------
             JN=0
             Do I=1,ND
                JA=IM+I; NSA=NS(JA); JN=(I-1)*imen
                QHLA=QHLA_opt(JA,ihil); FI2D=FI2D_opt(JA,ihil)
                FI1Z=FI1Z_opt(JA,ihil); FI1R=FI1R_opt(JA,ihil)
                !-----------------------------------------------
                ! QUASIPARTICLE WF IN COORDINATE SPACE
                !-----------------------------------------------
                If (NSA.Gt.0) Then
                   ! SPIN Up
                   Call DAXPY(IMEN,-QHLA,OMPANK(JN+1),1,OMPPFIU,1)
                   ! temperature
                   If(switch_on_temperature) Then
                      Call DAXPY(IMEN,-FI2D,OMPANK(JN+1),1,OMPPFIUD2N,1)
                      Call DAXPY(IMEN,-FI1R,OMPANK(JN+1),1,OMPPFIUR,1)
                      Call DAXPY(IMEN,-FI1Z,OMPANK(JN+1),1,OMPPFIUZ,1)
                   End If
                   Call DAXPY(IMEN, QHLA,OMPAN(JN+1) ,1,OMPFIU,1)
                   Call DAXPY(IMEN, FI2D,OMPAN(JN+1) ,1,OMPFIUD2N,1)
                   Call DAXPY(IMEN, FI1R,OMPAN(JN+1) ,1,OMPFIUR,1)
                   Call DAXPY(IMEN, FI1Z,OMPAN(JN+1) ,1,OMPFIUZ,1)
                   ! blocking
                   If(K0.Ne.0) Then
                      PNIK  = OMPANK(JN+K0)
                      PIU   = PIU   + PNIK*QHLA
                      PIUD2 = PIUD2 + PNIK*FI2D
                      PIUR  = PIUR  + PNIK*FI1R
                      PIUZ  = PIUZ  + PNIK*FI1Z
                   End If
                Else
                   ! SPIN Down
                   Call DAXPY(IMEN,-QHLA,OMPANK(JN+1),1,OMPPFID,1)
                   ! temperature
                   If(switch_on_temperature) Then
                      Call DAXPY(IMEN,-FI2D,OMPANK(JN+1),1,OMPPFIDD2N,1)
                      Call DAXPY(IMEN,-FI1R,OMPANK(JN+1),1,OMPPFIDR,1)
                      Call DAXPY(IMEN,-FI1Z,OMPANK(JN+1),1,OMPPFIDZ,1)
                   End If
                   Call DAXPY(IMEN, QHLA,OMPAN(JN+1) ,1,OMPFID,1)
                   Call DAXPY(IMEN, FI2D,OMPAN(JN+1) ,1,OMPFIDD2N,1)
                   Call DAXPY(IMEN, FI1R,OMPAN(JN+1) ,1,OMPFIDR,1)
                   Call DAXPY(IMEN, FI1Z,OMPAN(JN+1) ,1,OMPFIDZ,1)
                   ! blocking
                   If(K0.Ne.0) Then
                      PNIK  = OMPANK(JN+K0)
                      PID   = PID   + PNIK*QHLA
                      PIDD2 = PIDD2 + PNIK*FI2D
                      PIDR  = PIDR  + PNIK*FI1R
                      PIDZ  = PIDZ  + PNIK*FI1Z
                   End If
                End If
             End Do ! I=1,ND
             !-----------------------------------------------
             ! DENSITIES IN COORDINATE SPACE
             !-----------------------------------------------
             Takaihil=zero;    Troihil=zero;    Tdjihil=zero;   Ttauihil=zero;  Tdroihil=zero
             TSRFIihil=zero;   TSFIRihil=zero;  TSFIZihil=zero; TSZFIihil=zero; TNABLARIHIL=zero; TNABLAZIHIL=zero
             !
             Do K=1,IMEN
                TFIU=OMPFIU(K); TFIUZ=OMPFIUZ(K); TFIUR=OMPFIUR(K); TFIUD2=OMPFIUD2N(K); TPFIU=OMPPFIU(K)
                TFID=OMPFID(K); TFIDZ=OMPFIDZ(K); TFIDR=OMPFIDR(K); TFIDD2=OMPFIDD2N(K); TPFID=OMPPFID(K)
                !
                If(switch_on_temperature) Then
                   !
                   fk=f_T(K); f1k=f1_T(K)
                   !
                   TPFIUZ=OMPPFIUZ(K); TPFIUR=OMPPFIUR(K); TPFIUD2=OMPPFIUD2N(K)
                   TPFIDZ=OMPPFIDZ(K); TPFIDR=OMPPFIDR(K); TPFIDD2=OMPPFIDD2N(K)
                   !
                   TEMP1  = (TPFIU*TFIU+TPFID*TFID)*fk-(TFIU*TPFIU+TFID*TPFID)*f1k
                            TAKAIHIL = TAKAIHIL + TEMP1
                   TEMP2  = (TFIU*TFIU+TFID*TFID)*fk+(TPFIU*TPFIU+TPFID*TPFID)*f1k
                            TROIHIL = TROIHIL + TEMP2
                   TEMP3  = (TFIUR *TFIDZ -TFIDR *TFIUZ +XLAMY*TFIU *(TFIUR -TFIDZ) -XLAPY*TFID *(TFIDR +TFIUZ)) *fk &
                          + (TPFIUR*TPFIDZ-TPFIDR*TPFIUZ+XLAMY*TPFIU*(TPFIUR-TPFIDZ)-XLAPY*TPFID*(TPFIDR+TPFIUZ))*f1k
                            TDJIHIL = TDJIHIL + TEMP3
                   !
                   TW_T=(TFIUR *TFIUR +TFIDR *TFIDR +TFIUZ *TFIUZ +TFIDZ *TFIDZ)*fk&
                       +(TPFIUR*TPFIUR+TPFIDR*TPFIDR+TPFIUZ*TPFIUZ+TPFIDZ*TPFIDZ)*f1k
                   !
                   TEMP4  = (XLAMY2*TFIU *TFIU +XLAPY2*TFID *TFID) *fk &
                          + (XLAMY2*TPFIU*TPFIU+XLAPY2*TPFID*TPFID)*f1k + TW_T
                            TTAUIHIL = TTAUIHIL + TEMP4
                   TEMP5  = (TFIU*TFIUD2+TFID*TFIDD2)*fk + (TPFIU*TPFIUD2+TPFID*TPFIDD2)*f1k + TW_T
                            TDROIHIL = TDROIHIL + TEMP5
                   TEMP6  = (TFIUR*TFID-TFIDR*TFIU)*fk + (TPFIUR*TPFID-TPFIDR*TPFIU)*f1k
                            TSRFIIHIL = TSRFIIHIL + TEMP6
                   TEMP7  = (TFIU*TFID*XLAMPY)*fk + (TPFIU*TPFID*XLAMPY)*f1k
                            TSFIRIHIL = TSFIRIHIL + TEMP7
                   TEMP8  = (XLAMY*TFIU*TFIU-XLAPY*TFID*TFID)*fk + (XLAMY*TPFIU*TPFIU-XLAPY*TPFID*TPFID)*f1k
                            TSFIZIHIL = TSFIZIHIL + TEMP8
                   TEMP9  = (TFIUZ*TFID-TFIDZ*TFIU)*fk + (TPFIUZ*TPFID-TPFIDZ*TPFIU)*f1k
                            TSZFIIHIL = TSZFIIHIL + TEMP9
                   TEMP10 = (TFIUR*TFIU+TFIDR*TFID)*fk + (TPFIUR*TPFIU+TPFIDR*TPFID)*f1k
                            TNABLARIHIL = TNABLARIHIL + TEMP10
                   TEMP11 = (TFIUZ*TFIU+TFIDZ*TFID)*fk + (TPFIUZ*TPFIU+TPFIDZ*TPFID)*f1k
                            TNABLAZIHIL = TNABLAZIHIL + TEMP11
                   !
                Else
                   !
                   TEMP1  = TPFIU*TFIU+TPFID*TFID;                  TAKAIHIL    = TAKAIHIL   + TEMP1
                   TEMP2  = TFIU*TFIU+TFID*TFID;                    TROIHIL     = TROIHIL    + TEMP2
                   TEMP3  = TFIUR*TFIDZ-TFIDR*TFIUZ  &
                           +XLAMY*TFIU*(TFIUR-TFIDZ) &
                           -XLAPY*TFID*(TFIDR+TFIUZ) ;              TDJIHIL     = TDJIHIL    + TEMP3
                   !
                   TW_T=TFIUR*TFIUR+TFIDR*TFIDR+TFIUZ*TFIUZ+TFIDZ*TFIDZ
                   !
                   TEMP4  = XLAMY2*TFIU*TFIU+XLAPY2*TFID*TFID+TW_T; TTAUIHIL    = TTAUIHIL   + TEMP4
                   TEMP5  = TFIU*TFIUD2+TFID*TFIDD2          +TW_T; TDROIHIL    = TDROIHIL   + TEMP5
                   TEMP6  = TFIUR*TFID-TFIDR*TFIU;                  TSRFIIHIL   = TSRFIIHIL  + TEMP6
                   TEMP7  = TFIU*TFID*XLAMPY;                       TSFIRIHIL   = TSFIRIHIL  + TEMP7
                   TEMP8  = XLAMY*TFIU*TFIU-XLAPY*TFID*TFID;        TSFIZIHIL   = TSFIZIHIL  + TEMP8
                   TEMP9  = TFIUZ*TFID-TFIDZ*TFIU;                  TSZFIIHIL   = TSZFIIHIL  + TEMP9
                   TEMP10 = TFIUR*TFIU+TFIDR*TFID;                  TNABLARIHIL = TNABLARIHIL+ TEMP10
                   TEMP11 = TFIUZ*TFIU+TFIDZ*TFID;                  TNABLAZIHIL = TNABLAZIHIL+ TEMP11
                   !
                End If
                !
                If(K.Ne.K0) Cycle
                !
                ! blocking
                TAKAIHIL    = TAKAIHIL    - TEMP1;                 TEMP1  = PIU*PIU+PID*PID
                TROIHIL     = TROIHIL     - HALF*(TEMP2 - TEMP1);  TEMP2  = PIUR*PIDZ-PIDR*PIUZ+XLAMY*PIU*(PIUR-PIDZ) &
                                                                                               -XLAPY*PID*(PIDR+PIUZ)
                !
                PW_T=PIUR*PIUR+PIDR*PIDR+PIUZ*PIUZ+PIDZ*PIDZ
                TDJIHIL     = TDJIHIL     - HALF*(TEMP3 - TEMP2);  TEMP3  = PW_T+XLAMY2*PIU*PIU+XLAPY2*PID*PID
                TTAUIHIL    = TTAUIHIL    - HALF*(TEMP4 - TEMP3);  TEMP4  = PW_T+PIU*PIUD2+PID*PIDD2;
                TDROIHIL    = TDROIHIL    - HALF*(TEMP5 - TEMP4);  TEMP5  = PIUR*PID-PIDR*PIU;
                TSRFIIHIL   = TSRFIIHIL   - HALF*(TEMP6 - TEMP5);  TEMP6  = PIU*PID*XLAMPY;
                TSFIRIHIL   = TSFIRIHIL   - HALF*(TEMP7 - TEMP6);  TEMP7  = XLAMY*PIU*PIU-XLAPY*PID*PID;
                TSFIZIHIL   = TSFIZIHIL   - HALF*(TEMP8 - TEMP7);  TEMP8  = PIUZ*PID-PIDZ*PIU;
                TSZFIIHIL   = TSZFIIHIL   - HALF*(TEMP9 - TEMP8);  TEMP9  = PIUR*PIU+PIDR*PID;
                TNABLARIHIL = TNABLARIHIL - HALF*(TEMP10- TEMP9);  TEMP10 = PIUZ*PIU+PIDZ*PID;
                TNABLAZIHIL = TNABLAZIHIL - HALF*(TEMP11- TEMP10)
             End Do !K
#if(USE_OPENMP==1)
             OMPTaka(ib,ihil,it)    = OMPTaka(ib,ihil,it)    + TAKAIHIL
             OMPTro(ib,ihil,it)     = OMPTro(ib,ihil,it)     + TROIHIL
             OMPTdj(ib,ihil,it)     = OMPTdj(ib,ihil,it)     + TDJIHIL
             OMPTtau(ib,ihil,it)    = OMPTtau(ib,ihil,it)    + TTAUIHIL
             OMPTdro(ib,ihil,it)    = OMPTdro(ib,ihil,it)    + TDROIHIL
             OMPTSRFI(ib,ihil,it)   = OMPTSRFI(ib,ihil,it)   + TSRFIIHIL
             OMPTSFIR(ib,ihil,it)   = OMPTSFIR(ib,ihil,it)   + TSFIRIHIL
             OMPTSFIZ(ib,ihil,it)   = OMPTSFIZ(ib,ihil,it)   + TSFIZIHIL
             OMPTSZFI(ib,ihil,it)   = OMPTSZFI(ib,ihil,it)   + TSZFIIHIL
             OMPTNABLAR(ib,IHIL,IT) = OMPTNABLAR(ib,IHIL,IT) + TNABLARIHIL
             OMPTNABLAZ(ib,IHIL,IT) = OMPTNABLAZ(ib,IHIL,IT) + TNABLAZIHIL
#else
             OMPTaka(ihil,it)    = OMPTaka(ihil,it)    + TAKAIHIL
             OMPTro(ihil,it)     = OMPTro(ihil,it)     + TROIHIL
             OMPTdj(ihil,it)     = OMPTdj(ihil,it)     + TDJIHIL
             OMPTtau(ihil,it)    = OMPTtau(ihil,it)    + TTAUIHIL
             OMPTdro(ihil,it)    = OMPTdro(ihil,it)    + TDROIHIL
             OMPTSRFI(ihil,it)   = OMPTSRFI(ihil,it)   + TSRFIIHIL
             OMPTSFIR(ihil,it)   = OMPTSFIR(ihil,it)   + TSFIRIHIL
             OMPTSFIZ(ihil,it)   = OMPTSFIZ(ihil,it)   + TSFIZIHIL
             OMPTSZFI(ihil,it)   = OMPTSZFI(ihil,it)   + TSZFIIHIL
             OMPTNABLAR(IHIL,IT) = OMPTNABLAR(IHIL,IT) + TNABLARIHIL
             OMPTNABLAZ(IHIL,IT) = OMPTNABLAZ(IHIL,IT) + TNABLAZIHIL
#endif
          End Do !ihil
       End If
    End Do !bb
#if(USE_OPENMP==1)
!$OMP End Parallel Do
!
!$OMP PARALLEL DO        &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(nghl,nb,AKA,RO,DJ,TAU,DRO,SRFI,SFIR,SFIZ,SZFI,NABLAR,NABLAZ, &
!$OMP&        OMPTaka,OMPTro,OMPTdj,OMPTtau,OMPTdro,OMPTSRFI,&
!$OMP&        OMPTSFIR,OMPTSFIZ,OMPTSZFI,OMPTNABLAR,OMPTNABLAZ)
    Do ihil = 1,nghl
       Do ib=1,nb
          AKA(ihil,1)    = AKA(ihil,1)    + OMPTaka(ib,ihil,1)
          RO(ihil,1)     = RO(ihil,1)     + OMPTro(ib,ihil,1)
          DJ(ihil,1)     = DJ(ihil,1)     + OMPTdj(ib,ihil,1)
          TAU(ihil,1)    = TAU(ihil,1)    + OMPTtau(ib,ihil,1)
          DRO(ihil,1)    = DRO(ihil,1)    + OMPTdro(ib,ihil,1)
          SRFI(ihil,1)   = SRFI(ihil,1)   + OMPTSRFI(ib,ihil,1)
          SFIR(ihil,1)   = SFIR(ihil,1)   + OMPTSFIR(ib,ihil,1)
          SFIZ(ihil,1)   = SFIZ(ihil,1)   + OMPTSFIZ(ib,ihil,1)
          SZFI(ihil,1)   = SZFI(ihil,1)   + OMPTSZFI(ib,ihil,1)
          NABLAR(ihil,1) = NABLAR(ihil,1) + OMPTNABLAR(ib,ihil,1)
          NABLAZ(ihil,1) = NABLAZ(ihil,1) + OMPTNABLAZ(ib,ihil,1)
          !
          AKA(ihil,2)    = AKA(ihil,2)    + OMPTaka(ib,ihil,2)
          RO(ihil,2)     = RO(ihil,2)     + OMPTro(ib,ihil,2)
          DJ(ihil,2)     = DJ(ihil,2)     + OMPTdj(ib,ihil,2)
          TAU(ihil,2)    = TAU(ihil,2)    + OMPTtau(ib,ihil,2)
          DRO(ihil,2)    = DRO(ihil,2)    + OMPTdro(ib,ihil,2)
          SRFI(ihil,2)   = SRFI(ihil,2)   + OMPTSRFI(ib,ihil,2)
          SFIR(ihil,2)   = SFIR(ihil,2)   + OMPTSFIR(ib,ihil,2)
          SFIZ(ihil,2)   = SFIZ(ihil,2)   + OMPTSFIZ(ib,ihil,2)
          SZFI(ihil,2)   = SZFI(ihil,2)   + OMPTSZFI(ib,ihil,2)
          NABLAR(ihil,2) = NABLAR(ihil,2) + OMPTNABLAR(ib,ihil,2)
          NABLAZ(ihil,2) = NABLAZ(ihil,2) + OMPTNABLAZ(ib,ihil,2)
       End Do
    End Do
!$OMP End Parallel Do
    Deallocate(OMPTaka,OMPTro,OMPTdj,OMPTtau,OMPTdro,OMPTSRFI,OMPTSFIR,OMPTSFIZ, &
               OMPTSZFI,OMPTNABLAR,OMPTNABLAZ)
#else
    Do it=1,2
       Do ihil = 1,nghl
          AKA(ihil,it)    = OMPTaka(ihil,it)
          RO(ihil,it)     = OMPTro(ihil,it)
          DJ(ihil,it)     = OMPTdj(ihil,it)
          TAU(ihil,it)    = OMPTtau(ihil,it)
          DRO(ihil,it)    = OMPTdro(ihil,it)
          SRFI(ihil,it)   = OMPTSRFI(ihil,it)
          SFIR(ihil,it)   = OMPTSFIR(ihil,it)
          SFIZ(ihil,it)   = OMPTSFIZ(ihil,it)
          SZFI(ihil,it)   = OMPTSZFI(ihil,it)
          NABLAR(ihil,it) = OMPTNABLAR(ihil,it)
          NABLAZ(ihil,it) = OMPTNABLAZ(ihil,it)
       End Do
    End Do
#endif
    Do it = 1,2
       TRO=>ro(:,it);         TTAU=>tau(:,it);       TDJ=>dj(:,it);     TDRO=>dro(:,it)
       TSZFI=>SZFI(:,it);     TSFIZ=>SFIZ(:,it);     TSRFI=>SRFI(:,it); TSFIR=>SFIR(:,it)
       TNABLAR=>NABLAR(:,it); TNABLAZ=>NABLAZ(:,it); TAKA=>aka(:,it)
       s=two*Sum(tro); sd=four*Sum(tdro); drhoi(it)=sd
       !Sumnz(it)=Abs(s-Real(npr(it),Kind=pr))
       Sumnz(it)=Abs(s-tz(it))
       varmas=varmas+s
       !DNFactor(it)=Real(npr(it),Kind=pr)/s
       DNFactor(it)=tz(it)/s
       !----------------------------------------------------
       ! REMOVES INT.WEIGHTS AND MULTIPLIES BY THE JACOBIAN
       !----------------------------------------------------
       !piu=two*Real(npr(it),Kind=pr)/s
       piu=two*tz(it)/s
       WGT=wdcori
       Call dscal(NGHL,piu,WGT,1)
       Tro=Tro*WGT; Ttau=Ttau*WGT; Taka=Half*Taka*WGT;
       TSRFI=TSRFI*WGT; TSFIR=TSFIR*WGT;
       TSFIZ=TSFIZ*WGT; TSZFI=TSZFI*WGT;
       Call dscal(NGHL,two,WGT,1)
       Tdro=Tdro*WGT; Tdj=Tdj*WGT
       TNABLAR=TNABLAR*WGT; TNABLAZ=TNABLAZ*WGT
       !
    End Do !it
    DNFactor(3)=DNFactor(1)+DNFactor(2)
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('densit',1)
    !----------------------------------------------------
    ! COULOMB AND HARTREE FIELDS
    !----------------------------------------------------
    If(icou.Gt.0) Then
       If(nleg.Lt.0) Then
          Call coulom1
       Else
          Call coulom
       End If
    Else
       cou=zero
    End If
    !Call HartreeDir
  End Subroutine DENSIT
  !=======================================================================
  !> Calculates fields in r-space from axially symmetric densities
  !=======================================================================
  Subroutine field
    USE UNEDF, Only: is_NEDF
    Implicit None
    Integer(ipr) :: iw,it,ita,ihli,lambda,icons
    Real(pr) :: ra,ra2,rs,rsa,rsa0,z,rrr,rear_pair
    Real(pr) :: rt,rt1,ds,da,dt,dt1,tts,tta,tt,tt1,djs,dja,djt,djt1
    Real(pr) :: rsa0A,rsa0A1,V0V1,v01a,rns,rps,rsa1,rsa12,rsa10
    Real(pr) :: rsa0An,rsa0An1,rsa0As,rsa0As1
    Real(pr) :: RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,DJ_0,DJ_1
    Real(pr) :: SZFIN,SFIZN,SRFIN,SFIRN,SZFIP,SFIZP,SRFIP,SFIRP
    Real(pr) :: SZFI_0,SFIZ_0,SRFI_0,SFIR_0,SZFI_1,SFIZ_1,SRFI_1,SFIR_1
    Real(pr) :: SNABLARN,SNABLAZN,SNABLARP,SNABLAZP
    Real(pr) :: SNABLAR_0,SNABLAZ_0,SNABLAR_1,SNABLAZ_1
    Real(pr) :: J2_0,J2_1,JabJba_0,JabJba_1
    Real(pr) :: cx,x,fac_n,fac_p,gr,Pi,ec,kc,kf,lc,a,b
    Real(pr), Dimension(0:8) :: Qval
    Real(pr),Dimension(2) :: pUr,pUt,pUNr,pUNz,pUDr,pUDj,pUFIZ,pUZFI,pUFIR,pURFI
    Real(pr),Dimension(2) :: tUr,tUt,tUNr,tUNz,tUDr,tUDj,tUFIZ,tUZFI,tUFIR,tURFI
    Real(pr), Save :: ALAMBDA=0.0_pr,AEPSI=1.0_pr,CSPR=1.0_pr
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('field',0)
    !
    Pi = 4.0_pr*Atan(1.0_pr)
    !
    ! fields
    Do ihli=1,nghl
       !
       RHO_0 =ro(ihli,1)+ro(ihli,2)     ; RHO_1 =ro(ihli,1)-ro(ihli,2)
       TAU_0 =tau(ihli,1)+tau(ihli,2)   ; TAU_1 =tau(ihli,1)-tau(ihli,2)
       DRHO_0=dro(ihli,1)+dro(ihli,2)   ; DRHO_1=dro(ihli,1)-dro(ihli,2)
       DJ_0  =dj(ihli,1)+dj(ihli,2)     ; DJ_1  =dj(ihli,1)-dj(ihli,2)
       SFIZ_0=SFIZ(ihli,1)+SFIZ(ihli,2) ; SFIZ_1=SFIZ(ihli,1)-SFIZ(ihli,2)
       SFIR_0=SFIR(ihli,1)+SFIR(ihli,2) ; SFIR_1=SFIR(ihli,1)-SFIR(ihli,2)
       SZFI_0=SZFI(ihli,1)+SZFI(ihli,2) ; SZFI_1=SZFI(ihli,1)-SZFI(ihli,2)
       SRFI_0=SRFI(ihli,1)+SRFI(ihli,2) ; SRFI_1=SRFI(ihli,1)-SRFI(ihli,2)
       SNABLAR_0=NABLAR(ihli,1)+NABLAR(ihli,2)
       SNABLAR_1=NABLAR(ihli,1)-NABLAR(ihli,2)
       SNABLAZ_0=NABLAZ(ihli,1)+NABLAZ(ihli,2)
       SNABLAZ_1=NABLAZ(ihli,1)-NABLAZ(ihli,2)
       !
       J2_0=SFIZ_0**2+SFIR_0**2+SZFI_0**2+SRFI_0**2
       J2_1=SFIZ_1**2+SFIR_1**2+SZFI_1**2+SRFI_1**2
       JabJba_0=2*(SFIZ_0*SZFI_0+SFIR_0*SRFI_0)
       JabJba_1=2*(SFIZ_1*SZFI_1+SFIR_1*SRFI_1)
       !
       tUr=zero ; tUDr=zero ; tUNr=zero ; tUNz=zero
       tUt=zero ; tUDj=zero ; tUFIZ=zero ; tUZFI=zero
       tUFIR=zero ; tURFI=zero ;
       !
       Call calculate_U_parameters(RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1, &
            (SNABLAR_0**2+SNABLAZ_0**2),(SNABLAR_1**2+SNABLAZ_1**2) )
       !
       ! FUNCTIONAL
       ! E=E+(hb0*(TAU_0+TAU_1)*HALF+hb0*(TAU_0-TAU_1)*HALF  &                         ! tau
       !+Urhotau(0,0)*RHO_0*TAU_0+Urhotau(1,0)*RHO_1*TAU_1  &                         ! rho tau
       !+Urhotau(2,0)*RHO_0*TAU_1+Urhotau(3,0)*RHO_1*TAU_0  &
       !+Urhorho(0,0)*RHO_0**2+Urhorho(1,0)*RHO_1**2  &                               ! rho^2
       !+(Urhorho(2,0)+Urhorho(3,0))*RHO_0*RHO_1  &
       !+UrhoDrho(0,0)*RHO_0*DRHO_0+UrhoDrho(1,0)*RHO_1*DRHO_1  &                     ! rho Delta rho
       !+UrhoDrho(2,0)*RHO_0*DRHO_1+UrhoDrho(3,0)*RHO_1*DRHO_0  &
       !+Unablarho(0,0)*(SNABLAR_0*SNABLAR_0+SNABLAZ_0*SNABLAZ_0)  &                  ! (nabla rho)^2
       !+Unablarho(1,0)*(SNABLAR_1*SNABLAR_1+SNABLAZ_1*SNABLAZ_1)  &
       !+(Unablarho(3,0)+Unablarho(2,0))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)  &
       !+UrhonablaJ(0,0)*RHO_0*DJ_0+UrhonablaJ(1,0)*RHO_1*DJ_1  &                     ! rho nabla J
       !+UrhonablaJ(2,0)*RHO_0*DJ_1+UrhonablaJ(3,0)*RHO_1*DJ_0  &
       !+UJnablarho(0,0)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))  &     ! J nabla rho
       !+UJnablarho(1,0)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))  &
       !+UJnablarho(2,0)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))  &
       !+UJnablarho(3,0)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))  &
       !+UJJ(0,0)*J2_0+UJJ(1,0)*J2_1  &                                               ! JJ
       !+(UJJ(3,0)+UJJ(2,0))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
       !+UJabJba(0,0)*JabJba_0 + UJabJab(1,0)*JabJba_1                                ! Jab Jba
       !+(UJabJba(3,0)+UJabJba(2,0))*(SFIZ_0*SZFI_1+SFIR_0*SRFI_1+SZFI_0*SFIZ_1+SRFI_0*SFIR_1)
       !
       ! tUr(1)=dE/d RHO_0;       tUr(2)=dE/d RHO_1
       ! tUt(1)=dE/d TAU_0;       tUt(2)=dE/d TAU_1
       ! tUDr(1)=dE/d DeltaRHO_0; tUDr(2)=dE/d DeltaRHO_1
       ! and so on ...
       !
       !TEST
       !Write(*,'(4(2x,g26.10))') UrhoDrho(0,0)-CrDr(0),CrDr(0),UrhoDrho(1,1)-CrDr(1),CrDr(1); pause
       ! Contributions in the case 'u' depends on RHO_0
       tUr(1)=tUr(1)+two*Urhorho(0,0)*RHO_0+Urhorho(0,1)*RHO_0*RHO_0+Urhorho(1,1)*RHO_1*RHO_1  &  !! rho^2
                       +(Urhorho(3,0)+Urhorho(2,0))*RHO_1+(Urhorho(3,1)+Urhorho(2,1))*RHO_0*RHO_1
       tUr(2)=tUr(2)+two*Urhorho(1,0)*RHO_1+Urhorho(0,2)*RHO_0*RHO_0+Urhorho(1,2)*RHO_1*RHO_1  &
                       +(Urhorho(3,0)+Urhorho(2,0))*RHO_0+(Urhorho(3,2)+Urhorho(2,2))*RHO_0*RHO_1
       tUr(1)=tUr(1)+vDHartree(ihli,1)
       tUr(2)=tUr(2)+vDHartree(ihli,2)
       !
       tUr(1)=tUr(1)+Urhotau(0,0)*TAU_0+Urhotau(0,1)*TAU_0*RHO_0+Urhotau(1,1)*TAU_1*RHO_1  &  !! rho tau
                    +Urhotau(2,0)*TAU_1+Urhotau(2,1)*RHO_0*TAU_1+Urhotau(3,1)*RHO_1*TAU_0
       tUt(1)=tUt(1)+Urhotau(0,0)*RHO_0+Urhotau(3,0)*RHO_1
       tUr(2)=tUr(2)+Urhotau(1,0)*TAU_1+Urhotau(1,2)*TAU_1*RHO_1+Urhotau(0,2)*TAU_0*RHO_0  &
                    +Urhotau(3,0)*TAU_0+Urhotau(3,2)*RHO_1*TAU_0+Urhotau(2,2)*RHO_0*TAU_1
       tUt(2)=tUt(2)+Urhotau(1,0)*RHO_1+Urhotau(2,0)*RHO_0
       !
       tUr(1)=tUr(1)+UrhoDrho(0,0)*DRHO_0+UrhoDrho(0,1)*RHO_0*DRHO_0+UrhoDrho(1,1)*RHO_1*DRHO_1  &  !! rho Delta rho
                    +UrhoDrho(2,0)*DRHO_1+UrhoDrho(2,1)*RHO_0*DRHO_1+UrhoDrho(3,1)*RHO_1*DRHO_0
       tUDr(1)=tUDr(1)+UrhoDrho(0,0)*RHO_0+UrhoDrho(3,0)*RHO_1
       tUr(2)=tUr(2)+UrhoDrho(1,0)*DRHO_1+UrhoDrho(1,2)*RHO_1*DRHO_1+UrhoDrho(0,2)*RHO_0*DRHO_0  &
                    +UrhoDrho(3,0)*DRHO_0+UrhoDrho(3,2)*RHO_1*DRHO_0+UrhoDrho(2,2)*RHO_0*DRHO_1
       tUDr(2)=tUDr(2)+UrhoDrho(1,0)*RHO_1+UrhoDrho(2,0)*RHO_0
       !
       tUr(1)=tUr(1)+Unablarho(0,1)*(SNABLAR_0**2+SNABLAZ_0**2)+Unablarho(1,1)*(SNABLAR_1**2+SNABLAZ_1**2)  &  !! (nabla rho)^2
                   +(Unablarho(2,1)+Unablarho(3,1))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)
       tUNr(1)=tUNr(1)+two*Unablarho(0,0)*SNABLAR_0+(Unablarho(2,0)+Unablarho(3,0))*SNABLAR_1
       tUNz(1)=tUNz(1)+two*Unablarho(0,0)*SNABLAZ_0+(Unablarho(2,0)+Unablarho(3,0))*SNABLAZ_1
       tUr(2)=tUr(2)+Unablarho(0,2)*(SNABLAR_0**2+SNABLAZ_0**2)+Unablarho(1,2)*(SNABLAR_1**2  &
            +SNABLAZ_1**2)+(Unablarho(2,2)+Unablarho(3,2))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)
       tUNr(2)=tUNr(2)+two*Unablarho(1,0)*SNABLAR_1+(Unablarho(2,0)+Unablarho(3,0))*SNABLAR_0
       tUNz(2)=tUNz(2)+two*Unablarho(1,0)*SNABLAZ_1+(Unablarho(2,0)+Unablarho(3,0))*SNABLAZ_0
       !
       tUr(1)=tUr(1)+UrhonablaJ(0,0)*DJ_0+UrhonablaJ(0,1)*DJ_0*RHO_0+UrhonablaJ(1,1)*DJ_1*RHO_1  &  !! rho nabla J
                    +UrhonablaJ(2,0)*DJ_1+UrhonablaJ(2,1)*RHO_0*DJ_1+UrhonablaJ(3,1)*RHO_1*DJ_0
       tUDj(1)=tUDj(1)+UrhonablaJ(0,0)*RHO_0+UrhonablaJ(3,0)*RHO_1
       tUr(2)=tUr(2)+UrhonablaJ(1,0)*DJ_1+UrhonablaJ(1,2)*DJ_1*RHO_1+UrhonablaJ(0,2)*DJ_0*RHO_0  &
                    +UrhonablaJ(3,0)*DJ_0+UrhonablaJ(3,2)*RHO_1*DJ_0+UrhonablaJ(2,2)*RHO_0*DJ_1
       tUDj(2)=tUDj(2)+UrhonablaJ(1,0)*RHO_1+UrhonablaJ(2,0)*RHO_0
       !
       tUr(1)=tUr(1)+UJnablarho(0,1)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))      !! J nabla rho
       tUr(1)=tUr(1)+UJnablarho(1,1)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))
       tUr(1)=tUr(1)+UJnablarho(2,1)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))
       tUr(1)=tUr(1)+UJnablarho(3,1)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))

       tUr(2)=tUr(2)+UJnablarho(0,2)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))
       tUr(2)=tUr(2)+UJnablarho(1,2)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))
       tUr(2)=tUr(2)+UJnablarho(2,2)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))
       tUr(2)=tUr(2)+UJnablarho(3,2)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))

       tUNr(1)=tUNr(1)+UJnablarho(0,0)*(SFIZ_0-SZFI_0)
       tUNr(2)=tUNr(2)+UJnablarho(1,0)*(SFIZ_1-SZFI_1)
       tUNz(1)=tUNz(1)-UJnablarho(0,0)*(SFIR_0-SRFI_0)
       tUNz(2)=tUNz(2)-UJnablarho(1,0)*(SFIR_1-SRFI_1)

       tUFIZ(1)=tUFIZ(1)+UJnablarho(0,0)*SNABLAR_0*half
       tUFIZ(2)=tUFIZ(2)+UJnablarho(1,0)*SNABLAR_1*half
       tUZFI(1)=tUZFI(1)-UJnablarho(0,0)*SNABLAR_0*half
       tUZFI(2)=tUZFI(2)-UJnablarho(1,0)*SNABLAR_1*half
       tURFI(1)=tURFI(1)+UJnablarho(0,0)*SNABLAZ_0*half
       tURFI(2)=tURFI(2)+UJnablarho(1,0)*SNABLAZ_1*half
       tUFIR(1)=tUFIR(1)-UJnablarho(0,0)*SNABLAZ_0*half
       tUFIR(2)=tUFIR(2)-UJnablarho(1,0)*SNABLAZ_1*half
       !
       !! J.J (Mario: not tested for N2LO)
       tUr(1)=tUr(1)+UJJ(0,1)*J2_0+UJJ(1,1)*J2_1  &
            +(UJJ(3,1)+UJJ(2,1))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
       tUr(2)=tUr(2)+UJJ(0,2)*J2_0+UJJ(1,2)*J2_1  &
            +(UJJ(3,2)+UJJ(2,2))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
       tUFIZ(1)=tUFIZ(1)+UJJ(0,0)*SFIZ_0+half*(UJJ(3,0)+UJJ(2,0))*SFIZ_1
       tUFIR(1)=tUFIR(1)+UJJ(0,0)*SFIR_0+half*(UJJ(3,0)+UJJ(2,0))*SFIR_1
       tUZFI(1)=tUZFI(1)+UJJ(0,0)*SZFI_0+half*(UJJ(3,0)+UJJ(2,0))*SZFI_1
       tURFI(1)=tURFI(1)+UJJ(0,0)*SRFI_0+half*(UJJ(3,0)+UJJ(2,0))*SRFI_1
       tUFIZ(2)=tUFIZ(2)+UJJ(1,0)*SFIZ_1+half*(UJJ(3,0)+UJJ(2,0))*SFIZ_0
       tUFIR(2)=tUFIR(2)+UJJ(1,0)*SFIR_1+half*(UJJ(3,0)+UJJ(2,0))*SFIR_0
       tUZFI(2)=tUZFI(2)+UJJ(1,0)*SZFI_1+half*(UJJ(3,0)+UJJ(2,0))*SZFI_0
       tURFI(2)=tURFI(2)+UJJ(1,0)*SRFI_1+half*(UJJ(3,0)+UJJ(2,0))*SRFI_0
       !! Jab.Jba
       tUr(1) = tUr(1) + UJabJba(0,1)*JabJba_0 + UJabJba(1,1)*JabJba_1 &
            +(UJabJba(3,1)+UJabJba(2,1))*(SFIZ_0*SZFI_1+SFIR_0*SRFI_1+SZFI_0*SFIZ_1+SRFI_0*SFIR_1)
       tUr(2) = tUr(2) + UJabJba(0,2)*JabJba_0 + UJabJba(1,2)*JabJba_1 &
            +(UJabJba(3,2)+UJabJba(2,2))*(SFIZ_0*SZFI_1+SFIR_0*SRFI_1+SZFI_0*SFIZ_1+SRFI_0*SFIR_1)
       tUFIZ(1) = tUFIZ(1) + UJabJba(0,0)*SZFI_0 + half*(UJabJba(3,0)+UJabJba(2,0))*SZFI_1
       tUFIR(1) = tUFIR(1) + UJabJba(0,0)*SRFI_0 + half*(UJabJba(3,0)+UJabJba(2,0))*SRFI_1
       tUZFI(1) = tUZFI(1) + UJabJba(0,0)*SFIZ_0 + half*(UJabJba(3,0)+UJabJba(2,0))*SFIZ_1
       tURFI(1) = tURFI(1) + UJabJba(0,0)*SFIR_0 + half*(UJabJba(3,0)+UJabJba(2,0))*SFIR_1
       tUFIZ(2) = tUFIZ(2) + UJabJba(1,0)*SZFI_1 + half*(UJabJba(3,0)+UJabJba(2,0))*SZFI_0
       tUFIR(2) = tUFIR(2) + UJabJba(1,0)*SRFI_1 + half*(UJabJba(3,0)+UJabJba(2,0))*SRFI_0
       tUZFI(2) = tUZFI(2) + UJabJba(1,0)*SFIZ_1 + half*(UJabJba(3,0)+UJabJba(2,0))*SFIZ_0
       tURFI(2) = tURFI(2) + UJabJba(1,0)*SFIR_1 + half*(UJabJba(3,0)+UJabJba(2,0))*SFIR_0
       !
       tUr(1)=tUr(1)+UFnonstdr(0)                                        !! other amplitudes
       tUr(2)=tUr(2)+UFnonstdr(1)
       !
       !!  External Field
       !!
       !!tUr(1)=tUr(1)+Vexternal(0,zero,fl(ihli),fh(ihli))
       !!tUr(2)=tUr(2)+Vexternal(1,zero,fl(ihli),fh(ihli))
       !
       ! Contributions in the case 'u' depends on TAU_0
       !
       tUt(1)=tUt(1)+Urhotau(0,6)*RHO_0*TAU_0  &
                    +Urhotau(1,6)*RHO_1*TAU_1+Urhotau(2,6)*RHO_0*TAU_1  &
                    +Urhotau(3,6)*RHO_1*TAU_0+Urhorho(0,6)*RHO_0**2  &
                    +Urhorho(1,6)*RHO_1**2+(Urhorho(2,6)+Urhorho(3,6))*RHO_0*RHO_1  &
                    +UrhoDrho(0,6)*RHO_0*DRHO_0+UrhoDrho(1,6)*RHO_1*DRHO_1  &
                    +UrhoDrho(2,6)*RHO_0*DRHO_1+UrhoDrho(3,6)*RHO_1*DRHO_0  &
                    +Unablarho(0,6)*(SNABLAR_0*SNABLAR_0+SNABLAZ_0*SNABLAZ_0)  &
                    +Unablarho(1,6)*(SNABLAR_1*SNABLAR_1+SNABLAZ_1*SNABLAZ_1)  &
                    +(Unablarho(2,6)+Unablarho(3,6))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)  &
                    +UrhonablaJ(0,6)*RHO_0*DJ_0+UrhonablaJ(1,6)*RHO_1*DJ_1  &
                    +UrhonablaJ(2,6)*RHO_0*DJ_1+UrhonablaJ(3,6)*RHO_1*DJ_0  &
                    +UJnablarho(0,6)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))  &
                    +UJnablarho(1,6)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))  &
                    +UJnablarho(2,6)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))  &
                    +UJnablarho(3,6)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))
       tUt(1)=tUt(1)+UJJ(0,6)*J2_0+UJJ(1,6)*J2_1  &
                   +(UJJ(2,6)+UJJ(3,6))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
       tUt(1)=tUt(1)+UJabJba(0,6)*JabJba_0+UJJ(1,6)*JabJba_1  &
            +(UJabJba(2,6)+UJabJba(3,6))*(SFIZ_0*SZFI_1+SFIR_0*SRFI_1+SZFI_0*SFIZ_1+SRFI_0*SFIR_1)
       !
       ! Contributions in the case 'u' depends on DeltaRHO_0
       !
       tUDr(1)=tUDr(1)+Urhotau(0,7)*RHO_0*TAU_0  &
                      +Urhotau(1,7)*RHO_1*TAU_1+Urhotau(2,7)*RHO_0*TAU_1  &
                      +Urhotau(3,7)*RHO_1*TAU_0+Urhorho(0,7)*RHO_0**2  &
                      +Urhorho(1,7)*RHO_1**2+(Urhorho(2,7)+Urhorho(3,7))*RHO_0*RHO_1  &
                      +UrhoDrho(0,7)*RHO_0*DRHO_0+UrhoDrho(1,7)*RHO_1*DRHO_1  &
                      +UrhoDrho(2,7)*RHO_0*DRHO_1+UrhoDrho(3,7)*RHO_1*DRHO_0  &
                      +Unablarho(0,7)*(SNABLAR_0*SNABLAR_0+SNABLAZ_0*SNABLAZ_0)  &
                      +Unablarho(1,7)*(SNABLAR_1*SNABLAR_1+SNABLAZ_1*SNABLAZ_1)  &
                      +(Unablarho(2,7)+Unablarho(3,7))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)  &
                      +UrhonablaJ(0,7)*RHO_0*DJ_0+UrhonablaJ(1,7)*RHO_1*DJ_1  &
                      +UrhonablaJ(2,7)*RHO_0*DJ_1+UrhonablaJ(3,7)*RHO_1*DJ_0  &
                      +UJnablarho(0,7)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))  &
                      +UJnablarho(1,7)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))  &
                      +UJnablarho(2,7)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))  &
                      +UJnablarho(3,7)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))
       tUDr(1)=tUDr(1)+UJJ(0,7)*J2_0+UJJ(1,7)*J2_1  &
            +(UJJ(2,7)+UJJ(3,7))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
       tUDr(1)=tUDr(1)+UJabJba(0,7)*JabJba_0+UJJ(1,7)*JabJba_1  &
            +(UJabJba(2,7)+UJabJba(3,7))*(SFIZ_0*SZFI_1+SFIR_0*SRFI_1+SZFI_0*SFIZ_1+SRFI_0*SFIR_1)
       !
       ! proton-neutron representation
       pUr(1)  =tUr(1)+tUr(2);            pUr(2)  =tUr(1)  -tUr(2)
       pUt(1)  =tUt(1)+tUt(2)+hb0n*facECM;pUt(2)  =tUt(1)  -tUt(2)+hb0p*facECM
       pUDr(1) =tUDr(1)+tUDr(2);          pUDr(2) =tUDr(1) -tUDr(2)
       pUNr(1) =tUNr(1)+tUNr(2);          pUNr(2) =tUNr(1) -tUNr(2)
       pUNz(1) =tUNz(1)+tUNz(2);          pUNz(2) =tUNz(1) -tUNz(2)
       pUDj(1) =tUDj(1)+tUDj(2);          pUDj(2) =tUDj(1) -tUDj(2)
       pUFIZ(1)=tUFIZ(1)+tUFIZ(2);        pUFIZ(2)=tUFIZ(1)-tUFIZ(2)
       pUZFI(1)=tUZFI(1)+tUZFI(2);        pUZFI(2)=tUZFI(1)-tUZFI(2)
       pUFIR(1)=tUFIR(1)+tUFIR(2);        pUFIR(2)=tUFIR(1)-tUFIR(2)
       pURFI(1)=tURFI(1)+tURFI(2);        pURFI(2)=tURFI(1)-tURFI(2)
       !
       Do it=itmin,itmax   !! loop over n  & p
          ita=3-it
          ! constraining potential
          If (numberCons.Gt.0) Then
              z=fh(ihli); rrr=fl(ihli)**2
              Call moments_valueMesh(z,rrr,Qval)
              do icons=1,numberCons
                 lambda=multLambda(icons)
                 If(lambda.Ge.1) Then
                    pUr(it)= pUr(it) - multLag(lambda)*Qval(lambda)
                 End If
                 If(lambda.Eq.0) Then
                    pUr(it)= pUr(it) - neckLag*Exp(-((z-Z_NECK*bz)/AN_VAL)**2)
                 End If
              end do
          End If
          ! Coulomb
          If(it.Eq.2) Then
             If(icou.Ge.1) pUr(it)=pUr(it)+cou(ihli)
             If(icou.Eq.2.Or.icou.Eq.-4) pUr(it)=pUr(it)+CExPar*coex*ro(ihli,it)**p13 ! Slater approximation
          End If
          rsa0=(ro(ihli,it)+ro(ihli,ita))/rho_c
          ! Pairing regularization
          If(pairing_regularization) Then
             !
             ec = pwi - ala(it) ! cut-off energy for s.p. (for comparisons with HFBRAD)
             gr = (CpV0(it-1)*(ONE-rsa0*CpV1(it-1))) ! original pairing strength
             !
             If(it.Eq.1) Then
                If(is_NEDF) Then
                   MEFFn(ihli) = 0.5_pr*(hb0n+hb0p)
                Else
                   MEFFn(ihli) = hb0n + (Ctau(0)-Ctau(1))*RHO_0 + 2.0_pr*Ctau(1)*ro(ihli,it)
                End If
                fac_n=Sqrt(one/MEFFn(ihli))
                If(pUr(1) + ec - ala(it) .lt. 0.0_pr) Then
                   kc = fac_n*Sqrt(ala(it)+ec-pUr(1))
                   lc = fac_n*Sqrt(ala(it)-ec-pUr(1))
                   kf = fac_n*Sqrt(ala(it)-pUr(1))
                   a = 0.25_pr*(kc/MEFFn(ihli)/Pi**2)*(one-0.5_pr*kf/kc*Log((kc+kf)/(kc-kf)))
                   b = 0.25_pr*(lc/MEFFn(ihli)/Pi**2)*(one-0.5_pr*kf/lc*Log((lc+kf)/(kf-lc)))
                   geff_inv(ihli,it) = one/gr - a -b
                Else
                   If(ala(it) .Gt. pUr(1)) Then
                      kc = fac_n*Sqrt(ala(it)+ec-pUr(1))
                      kf = fac_n*Sqrt(ala(it)-pUr(1))
                      geff_inv(ihli,it) = one/gr - 0.25_pr*(kc/MEFFn(ihli)/Pi**2)*(one-0.5_pr*kf/kc   &
                                                                                *Log((kc+kf)/(kc-kf)))
                   Else
                      If(pUr(1) - ec - ala(it) .lt. 0.0_pr) Then
                         kc = fac_n*Sqrt(ala(it)+ec-pUr(1))
                         kf = fac_n*Sqrt(pUr(1)-ala(it))
                         geff_inv(ihli,it) = one/gr - 0.25_pr*(kc/MEFFn(ihli)/Pi**2)*(one+kf/kc &
                                                                                    *Atan(kf/kc))
                      Else
                         geff_inv(ihli,it) = one/gr
                      End If
                   End If
                End If
                ! contribution to Delta
                dvn(ihli)=aka(ihli,it)/geff_inv(ihli,it)
                ! contribution to rearrangement
                !rear_pair =
             Else
                If(is_NEDF) Then
                   MEFFp(ihli) = 0.5_pr*(hb0n+hb0p)
                Else
                   MEFFp(ihli) = hb0p + (Ctau(0)-Ctau(1))*RHO_0 + 2.0_pr*Ctau(1)*ro(ihli,it)
                End If
                fac_p=Sqrt(one/MEFFp(ihli))
                If(pUr(2) + ec - ala(it) .lt. 0.0_pr) Then
                   kc = fac_p*Sqrt(ala(it)+ec-pUr(2))
                   lc = fac_p*Sqrt(ala(it)-ec-pUr(2))
                   kf = fac_p*Sqrt(ala(it)-pUr(2))
                   a = 0.25_pr*(kc/MEFFp(ihli)/Pi**2)*(one-0.5_pr*kf/kc*Log((kc+kf)/(kc-kf)))
                   b = 0.25_pr*(lc/MEFFp(ihli)/Pi**2)*(one-0.5_pr*kf/lc*Log((lc+kf)/(kf-lc)))
                   geff_inv(ihli,it) = one/gr - a -b
                Else
                   If(ala(it).Gt.pUr(2)) Then
                      kc = fac_p*Sqrt(ala(it)+ec-pUr(2))
                      kf = fac_p*Sqrt(ala(it)-pUr(2))
                      geff_inv(ihli,it) = one/gr - 0.25_pr*(kc/MEFFp(ihli)/Pi**2)*(one-0.5_pr*kf/kc &
                                                                                *Log((kc+kf)/(kc-kf)))
                   Else
                      If(pUr(2) - ec - ala(it) .lt. 0.0_pr) Then
                         kc = fac_p*Sqrt(ala(it)+ec-pUr(2))
                         kf = fac_p*Sqrt(pUr(2)-ala(it))
                         geff_inv(ihli,it) = one/gr - 0.25_pr*(kc/MEFFp(ihli)/Pi**2)*(one+kf/kc &
                                                                                    *Atan(kf/kc))
                      Else
                         geff_inv(ihli,it) = one/gr
                      End If
                   End If
                End If
                ! contribution to Delta
                dvp(ihli)=aka(ihli,it)/geff_inv(ihli,it)
                ! contribution to rearrangement
                !rear_pair =
             End If
          Else
              ! pairing contribution to delta dv(ihli,it)
              If(it.Eq.1) Then
                 dvn(ihli)=(CpV0(it-1)*(ONE-rsa0*CpV1(it-1)))*aka(ihli,it)
              Else
                 dvp(ihli)=(CpV0(it-1)*(ONE-rsa0*CpV1(it-1)))*aka(ihli,it)
              End If
             ! pairing contribution to rearrangement term
              If(use_TMR_pairing.eq.0) then
                 rear_pair = CpV0(it-1) *CpV1(it-1) /rho_c*aka(ihli,it)**2 &
                           + CpV0(ita-1)*CpV1(ita-1)/rho_c*aka(ihli,ita)**2
                 pUr(it)=pUr(it)-rear_pair
              endif
          End If
       End Do !it
       !
       vn(ihli)=pUr(1)       ; vp(ihli)=pUr(2)        !* RHO_ij
       vhbn(ihli)=pUt(1)     ; vhbp(ihli)=pUt(2)      !* TAU_ij
       vrn(ihli)=pUNr(1)     ; vrp(ihli)=pUNr(2)      !* NABLAr RHO__ij
       vzn(ihli)=pUNz(1)     ; vzp(ihli)=pUNz(2)      !* NABLAz RHO__ij
       vdn(ihli)=pUDr(1)     ; vdp(ihli)=pUDr(2)      !* DELTA RHO_ij
       vsn(ihli)=pUDj(1)     ; vsp(ihli)=pUDj(2)      !* NABLA . J__ij
       vSFIZn(ihli)=pUFIZ(1) ; vSFIZp(ihli)=pUFIZ(2)  !* JFIZ_ij
       vSZFIn(ihli)=pUZFI(1) ; vSZFIp(ihli)=pUZFI(2)  !* JZFI_ij
       vSFIRn(ihli)=pUFIR(1) ; vSFIRp(ihli)=pUFIR(2)  !* JFIR_ij
       vSRFIn(ihli)=pURFI(1) ; vSRFIp(ihli)=pURFI(2)  !* JRFI_ij
       !
    End Do !ihli
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('field',1)
    !
  End Subroutine field
  !=======================================================================
  !> Computes the constraining fields at the Gauss quadrature points
  !=======================================================================
  Subroutine constraining_field()
    Implicit None
    Integer(ipr) :: ihli,lambda,icons
    Real(pr) :: z,rrr
    Real(pr), Dimension(0:8) :: Qval
    ! fields
    Do icons=1,numberCons
       lambda=multLambda(icons)
       ! Regular multipole moments
       If(lambda.Ge.1) Then
          Do ihli=1,nghl
             z=fh(ihli); rrr=fl(ihli)**2
             Call moments_valueMesh(z,rrr,Qval)
             qfield(ihli,lambda) = Qval(lambda)
          End do
       End If
       ! Gaussian neck operator
       If(lambda.Eq.0) Then
          Do ihli=1,nghl
             z=fh(ihli)
             qfield(ihli,lambda) = Exp(-((z-Z_NECK*bz)/AN_VAL)**2)
          End Do
       End If
    End Do
  End Subroutine constraining_field
  !=======================================================================
  !> Computes the mean field and pairing field of the HFB matrix in
  !> configuration space
  !=======================================================================
  Subroutine gamdel(WoodsSaxon,DoMixing)
    Implicit None
    logical, intent(in) :: WoodsSaxon,DoMixing
    Integer(ipr) :: i,ih,il,ib,ibx,nd,nd2,nza,nra,nla,nsa,nsb,nsab,icons,lambda
    Integer(ipr) :: ihil,laplus,im,JA,N1,N2,ndnd,n12,n21
    Integer(ipr) :: i1,i2,i3
    Real(pr) :: qla,yi,y,y2,qha,qhla,xmi,u2,un,up,xxx
    Real(pr) :: sml2,cnzaa,cnraa,SSU,SSD
    Real(pr) :: FITW1,FITW2,FITW3,FITW4
    Real(pr) :: fi1r,fi1z,fi2d,QHL1A,QH1LA
    Real(pr) :: vh,vdh,vsh,hbh,vsum
    Real(pr) :: SRFIh,SFIRh,SFIZh,SZFIh,SNABLARh,SNABLAZh
    Real(pr) :: xlam,xlam2,xlamy,xlamy2,xlap,xlap2,xlapy,xlapy2,XLAMPY
    Real(pr) :: FIUN1,FIDN1,FIURN1,FIDRN1,FIUZN1,FIDZN1,FIUD2N1,FIDD2N1
    Real(pr) :: FIUN2,FIDN2,FIURN2,FIDRN2,FIUZN2,FIDZN2,FIUD2N2,FIDD2N2
    Real(pr) :: FIUN12,FIDN12,FIURN12,FIDRN12,FIUZN12,FIDZN12
    Real(pr) :: vnhl,vrnhl,vznhl,vdnhl,vsnhl,vhbnhl,vSRFInhl,vSFIRnhl
    Real(pr) :: vSFIZnhl,vSZFInhl,vphl,vrphl,vzphl,vdphl,vsphl,vhbphl
    Real(pr) :: vSRFIphl,vSFIRphl,vSFIZphl,vSZFIphl,dvnhl,dvphl
    Integer(ipr) :: ibro
    Integer(ipr) :: ndxmax
    Parameter(ndxmax=(n00max+2)*(n00max+2)/4)
    Real(pr) :: OMPFIU(ndxmax),OMPFID(ndxmax),OMPFIUR(ndxmax),OMPFIDR(ndxmax),OMPFIUZ(ndxmax), &
                OMPFIDZ(ndxmax),OMPFIUD2N(ndxmax),OMPFIDD2N(ndxmax)
    integer(ipr) ::  t1,t2,countrate,countmax
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('gamdel',0)
    !
    If (neck_constraints) Call QNFIND()
    !
    !----------------------------------------------
    ! START BLOCKS
    !----------------------------------------------
    brout=zero; ibro=0
    If (.Not. Allocated(allibro)) Then
       Allocate(allibro(1:NB))
       allibro(1)=0
       Do ib=2,NB
          allibro(ib) = allibro(ib-1) + (ID(ib-1)*(ID(ib-1)+1)/2)
       End Do
    End If
!$OMP PARALLEL DO        &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(NB,ID,IA,NBX,NS,nghl, &
!$OMP&        NHHDIM2,NHHDIM3,NHHDIM4,allibro, &
!$OMP&        vSRFIn,vSFIRn,vSFIZn,vSZFIn, &
!$OMP&        vSRFIp,vSFIRp,vSFIZp,vSZFIp, &
!$OMP&        vn,vrn,vzn,vdn,vsn,vhbn,dvn, &
!$OMP&        vp,vrp,vzp,vdp,vsp,vhbp,dvp, &
!$OMP&        QHLA_opt,FI1R_opt, FI1Z_opt, FI2D_opt, y_opt, &
!$OMP&        nhhdim,kindhfb,ALA2,RK,brout,Parity) &
!$OMP& PRIVATE(I,ND,IB,IM,IBX,LAPLUS,XLAM,XLAP,XLAM2,IL,IH,IHIL,Y,Y2, &
!$OMP&         XLAMY,XLAMY2,XLAP2,XLAPY,XLAPY2,XLAMPY,N1,JA,NSA,SSU,SSD, &
!$OMP&         vnhl,vrnhl,vznhl,vdnhl,vsnhl,vhbnhl,dvnhl, &
!$OMP&         vphl,vrphl,vzphl,vdphl,vsphl,vhbphl,dvphl, &
!$OMP&         vSRFInhl,vSFIRnhl,vSFIZnhl,vSZFInhl,&
!$OMP&         vSRFIphl,vSFIRphl,vSFIZphl,vSZFIphl,&
!$OMP&         FI2D,i1,i2,i3,NSB,NSAB,SNABLARh, SNABLAZh,FI1R,FI1Z, &
!$OMP&         FIUD2N1,FIDD2N1,FIUD2N2,FIDD2N2,FITW3,FITW4,&
!$OMP&         OMPFIUD2N,OMPFIDD2N,OMPFIU,OMPFIUR,OMPFIUZ,OMPFID,OMPFIDR,OMPFIDZ, &
!$OMP&         FIUN1,FIDN1,FIURN1,FIDRN1,FIUZN1,FIDZN1,N2,FIUN2,FIDN2,FIURN2,  &
!$OMP&         FIDRN2,FIUZN2,FIDZN2,FIUN12,FIDN12,FIURN12,FIDRN12,FIUZN12,FIDZN12,VH,&
!$OMP&         HBH,VDH,VSH,SRFIH,SFIRH,SFIZH,SZFIH,UN,UP,N12,QHLA)
    Do ib=1,NB
       ND=ID(ib); IM=ia(ib); ibx=ib+nbx
       If(Parity) Then
          LAPLUS=(ib+1)/2 !Yesp
       Else
          LAPLUS=ib       !Nop
       End If
       XLAP=LAPLUS; XLAM=XLAP-ONE; xlap2=xlap*xlap; xlam2=xlam*xlam
       !----------------------------------------------
       ! SUM OVER GAUSS INTEGRATION POINTS
       !----------------------------------------------
       Do ihil=1,nghl
          y=y_opt(ihil); xlamy=xlam*y;     xlapy=xlap*y;   XLAMPY=XLAMY+XLAPY
          y2=y*y;        xlamy2=xlam2*y2;  xlapy2=xlap2*y2
          !
          vnhl=vn(ihil);         vrnhl=vrn(ihil);       vznhl=vzn(ihil);       vdnhl=vdn(ihil)
          vsnhl=vsn(ihil);       vhbnhl=vhbn(ihil);     vSRFInhl=vSRFIn(IHIL); vSFIRnhl=vSFIRn(IHIL)
          vSFIZnhl=vSFIZn(IHIL); vSZFInhl=vSZFIn(IHIL); vphl=vp(ihil);         vrphl=vrp(ihil)
          vzphl=vzp(ihil);       vdphl=vdp(ihil);       vsphl=vsp(ihil);       vhbphl=vhbp(ihil)
          vSRFIphl=vSRFIp(IHIL); vSFIRphl=vSFIRp(IHIL); vSFIZphl=vSFIZp(IHIL); vSZFIphl=vSZFIp(IHIL)
          dvnhl=dvn(ihil);       dvphl=dvp(ihil)
          !
          Do N1=1,ND
             JA=IM+N1;               NSA=NS(JA);             SSU=Max(NSA,0);         SSD=Max(-NSA,0)
             QHLA=QHLA_opt(JA,ihil); FI1R=FI1R_opt(JA,ihil); FI1Z=FI1Z_opt(JA,ihil); FI2D=FI2D_opt(JA,ihil)
             OMPFIU(N1)=QHLA*SSU;    OMPFIUR(N1)=fi1r*SSU
             OMPFIUZ(N1)=fi1z*SSU;   OMPFIUD2N(N1)=(FI2D-XLAMY2*QHLA)*SSU
             OMPFID(N1)=QHLA*SSD;    OMPFIDR(N1)=fi1r*SSD
             OMPFIDZ(N1)=fi1z*SSD;   OMPFIDD2N(N1)=(FI2D-XLAPY2*QHLA)*SSD
          End Do
          !
          i=allibro(ib)
          Do N1=1,ND
             JA=IM+N1;               NSA=NS(JA)
             FIUN1=OMPFIU(N1);       FIURN1=OMPFIUR(N1);
             FIUZN1=OMPFIUZ(N1);     FIUD2N1=OMPFIUD2N(N1)
             FIDN1=OMPFID(N1);       FIDRN1=OMPFIDR(N1);
             FIDZN1=OMPFIDZ(N1);     FIDD2N1=OMPFIDD2N(N1)
             Do N2=1,N1
                I=I+1; i1=i+nhhdim; i2=i+nhhdim2; i3=i+nhhdim3; NSB=NS(N2+IM); NSAB=NSA+NSB
                If (NSAB.Ne.0) Then
                   If (NSB.Gt.0) Then                                    !spin:UpUp
                      FIUN2    = OMPFIU(N2);    FIURN2 = OMPFIUR(N2)
                      FIUD2N2  = OMPFIUD2N(N2); FIUZN2 = OMPFIUZ(N2)
                      vh       = FIUN1*FIUN2
                      hbh      = vh*XLAMY2+FIURN1*FIURN2+FIUZN1*FIUZN2
                      vdh      = hbh+hbh+FIUN1*FIUD2N2+FIUN2*FIUD2N1
                      SNABLARh = FIURN1*FIUN2+FIURN2*FIUN1
                      SNABLAZh = FIUZN1*FIUN2+FIUZN2*FIUN1
                      vsh      = SNABLARh*XLAMY
                      SFIZh    = (vh+vh)*XLAMY ! =SFIZh (v103)
                   Else                                                  !spin:DoDo
                      FIDN2    = OMPFID(N2);  FIDRN2  = OMPFIDR(N2);
                      FIDZN2   = OMPFIDZ(N2); FIDD2N2 = OMPFIDD2N(N2)
                      vh       = FIDN1*FIDN2
                      hbh      = vh*XLAPY2+FIDRN1*FIDRN2+FIDZN1*FIDZN2
                      vdh      = hbh+hbh+FIDN1*FIDD2N2+FIDN2*FIDD2N1;
                      SNABLARh = FIDRN1*FIDN2+FIDRN2*FIDN1
                      SNABLAZh = FIDZN1*FIDN2+FIDZN2*FIDN1
                      vsh      =-SNABLARh*XLAPY
                      SFIZh    =-(vh+vh)*XLAPY ! =SFIZh (v103)
                   End If
                   brout(i )=brout(i )+vSFIZnhl*SFIZh+vh*vnhl+SNABLARh*vrnhl+SNABLAZh*vznhl+vdh*vdnhl+vsh*vsnhl+hbh*vhbnhl
                   brout(i1)=brout(i1)+vSFIZphl*SFIZh+vh*vphl+SNABLARh*vrphl+SNABLAZh*vzphl+vdh*vdphl+vsh*vsphl+hbh*vhbphl
                   brout(i2)=brout(i2)+vh*dvnhl
                   brout(i3)=brout(i3)+vh*dvphl
                Else
                   If (NSB.Gt.0) Then                                                                !spin:DoUp
                      !vh=ZERO; hbh=ZERO; vdh=ZERO; SNABLARh=ZERO; SNABLAZh=ZERO; SFIZh=ZERO
                      FIUN2   = OMPFIU(N2);    FIURN2 = OMPFIUR(N2);
                      FIUD2N2 = OMPFIUD2N(N2); FIUZN2 = OMPFIUZ(N2)
                      FITW3   =-FIDZN1*FIUN2; FITW4=FIUZN2*FIDN1
                      vsh     =-FIDRN1*FIUZN2+FIURN2*FIDZN1+FITW3*XLAMY-FITW4*XLAPY
                      SRFIh   =-FIDRN1*FIUN2+FIURN2*FIDN1
                      SFIRh   = FIDN1*FIUN2*XLAMPY
                      SZFIh   = FITW3+FITW4
                   Else                                                                             !spin:UpDo
                      !vh=ZERO; hbh=ZERO; vdh=ZERO; SNABLARh=ZERO; SNABLAZh=ZERO; SFIZh=ZERO
                      FIDN2  = OMPFID(N2);   FIDRN2  = OMPFIDR(N2);
                      FIDZN2 = OMPFIDZ(N2);  FIDD2N2 = OMPFIDD2N(N2)
                      FITW3  =-FIDZN2*FIUN1; FITW4=FIUZN1*FIDN2
                      vsh    = FIURN1*FIDZN2-FIDRN2*FIUZN1-FITW4*XLAPY+FITW3*XLAMY ! -vsh (v103)
                      SRFIh  = FIURN1*FIDN2-FIDRN2*FIUN1 !=SRFIh (v103)
                      SFIRh  = FIUN1*FIDN2*XLAMPY !=SFIRh(v103)
                      SZFIh  = FITW3+FITW4 !=SZFIh(v103)
                   End If
                   brout(i )=brout(i )+vsh*vsnhl+vSRFInhl*SRFIh+vSFIRnhl*SFIRh+vSZFInhl*SZFIh
                   brout(i1)=brout(i1)+vsh*vsphl+vSRFIphl*SRFIh+vSFIRphl*SFIRh+vSZFIphl*SZFIh
                End If
                !----------------------------------------------
                ! LN PH PART
                !----------------------------------------------
                If(kindhfb.Lt.0) Then
                   If(ihil.Eq.1) Then
                      un=zero; up=zero;
                      If(N1.Eq.N2) Then
                         un=-ala2(1); up=-ala2(2)
                      End If
                      n12=N1+(N2-1)*ND
                      brout(i )=brout(i )+two*(ala2(1)*rk(n12,ib )+un)
                      brout(i1)=brout(i1)+two*(ala2(2)*rk(n12,ibx)+up)
                   End If
                End If
             End Do !N2
          End Do !N1
       End Do !ihil
    End Do !IB
!$OMP End Parallel Do
    If (IDEBUG.Eq.1) Call get_CPU_time('gamdel',1)
    !
    if(.not.WoodsSaxon.and.(finite_range.or.coulomb_gaussian)) then
       call system_clock(t1,countrate,countmax)
       call gamdel_gogny
       call system_clock(t2,countrate,countmax)
       wct_gogny = (t2-t1)/real(countrate,kind=pr)
    else
       wct_gogny = 0._pr
    endif

    ! Lagrange parameters for the constraints
    Do lambda=1,lambdaMax
       brout(nhhdim4+lambda)=multLag(lambda)
    End Do
    If(neck_constraints) Then
       brout(nhhdim4+lambdaMax+1)=neckLag
    End If
    !----------------------------------------------
    ! BROYDEN/LINEAR MIXING
    !----------------------------------------------
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('broyden',0)
    !
    If(DoMixing) Then
       If(neck_constraints) Then
          Call broyden_min(nhhdim4+lambdaMax+1,brout,brin,alphamix,si,iiter,nbroyden,bbroyden)
       Else
          Call broyden_min(nhhdim4+lambdaMax,brout,brin,alphamix,si,iiter,nbroyden,bbroyden)
       End If
    Else
       brin=brout
    End If
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('broyden',1)
    !
    Do lambda=1,lambdaMax
       multLag(lambda)=brin(nhhdim4+lambda)
    End Do
    If(neck_constraints) Then
       neckLag=brin(nhhdim4+lambdaMax+1)
    End If
    !
  End Subroutine gamdel
  !=====================================================================
  !>  Calculates the HFB fields Gamma and Delta for a finite range
  !>  Gogny type functional using the axial and radial components of
  !>  two body potential matrix elements calculated with the subroutines
  !>  in the \ref hfbtho_gogny module.
  !>
  !>  If the preprocessor variable GOGNY_SYMMETRIES is set to 1 the
  !>  calculation takes advantage of the separability of the matrix
  !>  elements and first 'folds' the axial matrix elements with the
  !>  density and pairing matrices and creates the objects Zrho and
  !>  Zkappa. Subsequently this objects are 'folded' with the radial
  !>  matrix elements to obtain the Gamma and Delta fields. This
  !>  significantly reduces the number of operations and speeds up the
  !>  calculation. In this mode several other symmetries of the axial
  !>  matrix element are used to reduced the computation time.
  !>
  !>  For any other value of GOGNY_SYMMETRIES the loops to calculate
  !>  the HFB fields goes over every combination of quantum numbers.
  !>  This mode of course requires a larger computation time and is only
  !>  used when debugging the calculation of the fields.
  !======================================================================
  Subroutine gamdel_gogny
    Use UNEDF, Only: e2charg
    !--------------------------------------------------------------------
    ! ph- and pp- matrices in configurational space
    !--------------------------------------------------------------------
    Implicit None
    integer(ipr) :: ib,jb,ind,jnd,im,jm,ibx,jbx,i,ig
    integer(ipr) :: i1,i2,i3
    integer(ipr) :: n1,n2,n3,n4,j1,j2,j3,j4,n43,n12,n21,n34
    integer(ipr) :: nlb,nld,nzac
    integer(ipr) :: iz_abcd,iz_abdc,iz_acbd
    integer(ipr) :: iz_adcb,iz_adbc,iz_acdb
    integer(ipr) :: ir_abcd,ir_abdc,ir_acbd
    integer(ipr) :: ir_adcb,ir_adbc,ir_acdb
    integer(ipr) :: jr_abcd,jr_abdc,jr_acbd
    integer(ipr) :: jr_adcb,jr_adbc,jr_acdb
    real(pr) :: rho_n,rho_p,kap_n,kap_p
    real(pr) :: gamma_n_dir,gamma_n_exc,gamma_p_dir,gamma_p_exc
    real(pr) :: delta_n_dir,delta_n_exc,delta_p_dir,delta_p_exc
    real(pr) :: gamma_c_dir,gamma_c_exc,delta_c_dir
    real(pr) :: Vz_abcd, Vz_abdc
    real(pr) :: V_abcd,V_abdc,V_acbd,V_acdb
    real(pr) :: U_abcd,U_abdc,U_acbd,U_acdb
    real(pr) :: Vdir,Udir,Vexc,Uexc,Vdel,Udel
    real(pr) :: W,B,H,M
    integer(ipr) :: oa,ob,oc,od,la,lb,lc,ld,iba,ibb,ibc,ibd,fact
    real(pr) :: Zrho_n_dir(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    real(pr) :: Zrho_p_dir(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    real(pr) :: Zrho_n_exc(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    real(pr) :: Zrho_p_exc(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    real(pr) :: Zkap_n_dir(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    real(pr) :: Zkap_p_dir(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    real(pr) :: zrho_nd,zrho_pd,zrho_ne,zrho_pe,zkap_n,zkap_p
    integer(ipr) :: ka,pa,ra,za,lmina,lmaxa,sa,ja,na,ida,ima,ibax
    integer(ipr) :: kb,pb,rb,zb,lminb,lmaxb,sb,nb,idb,imb,ibbx
    integer(ipr) :: kc,pc,rc,zc,lminc,lmaxc,sc,jc,nc,idc,imc,ibcx
    integer(ipr) :: kd,pd,rd,zd,lmind,lmaxd,sd,jd,nd,idd,imd,ibdx
    integer(ipr) :: nac,ndb,sac,sdb,n_rlx
    integer(ipr) :: ita,itb,itc,itd,npa,npb,npc,npd,izb,izd
    if(allocated(gamma_g_dir)) deallocate(gamma_g_dir,gamma_g_exc,delta_g_dir,coulf_g_dir,coulf_g_exc,coulf_d_dir)
    allocate(gamma_g_dir(ndx**2,2*nbx),gamma_g_exc(ndx**2,2*nbx),delta_g_dir(ndx**2,2*nbx),coulf_g_dir(ndx**2,2*nbx),coulf_g_exc(ndx**2,2*nbx),coulf_d_dir(ndx**2,2*nbx))
    gamma_g_dir = zero
    gamma_g_exc = zero
    delta_g_dir = zero
    coulf_g_dir = zero
    coulf_g_exc = zero
    coulf_d_dir = zero
#if(GOGNY_SYMMETRIES==1)
    n_rlx = max(2*nrx,nlx)
    if(.not.force_is_dme) then
!$OMP Parallel Default(None)&
!$OMP& SHARED(nzx,nttx,nrr,nll,nss,noo,nzzx,nzz,ib_zrls,id,ia,nbx,&
!$OMP&        i_zrls,rk,ak,VzGogny,ntx,nz,nr,nl,ns,nox,nrx,nrlx,n_rlx,&
!$OMP&        W_g,B_g,H_g,M_g,VrGogny,allibro,nhhdim,nhhdim2,nhhdim3,&
!$OMP&        gamma_g_dir,gamma_g_exc,delta_g_dir,brout,n_g_all,n_g, &
!$OMP&        V_g_coul,coulf_g_dir,coulf_g_exc,coulf_d_dir,icou,e2charg) &
!$OMP& PRIVATE(nzac,Zrho_n_dir,Zrho_p_dir,Zrho_n_exc,Zrho_p_exc,&
!$OMP&         Zkap_n_dir,Zkap_p_dir,zc,za,itb,rb,lb,sb,npb,nlb,itd,npd,&
!$OMP&         rd,ld,sd,nld,izb,zb,ibb,idb,imb,ibbx,jb,nb,izd,zd,ibd,jd,&
!$OMP&         imd,nd,ndb,rho_n,rho_p,kap_n,kap_p,iz_abdc,iz_abcd,&
!$OMP&         iz_acbd,ig,V_abcd,V_abdc,V_acbd,ita,ra,la,sa,iba,ida,ima,&
!$OMP&         ibax,na,itc,rc,lc,sc,ibc,imc,nc,sac,gamma_n_dir,&
!$OMP&         gamma_p_dir,gamma_n_exc,gamma_p_exc,delta_n_dir,&
!$OMP&         delta_p_dir,ob,sdb,ir_abdc,ir_abcd,ir_acbd,jr_abdc,&
!$OMP&         jr_abcd,jr_acbd,W,B,H,M,zrho_nd,zrho_pd,zrho_ne,zrho_pe,&
!$OMP&         zkap_n,zkap_p,U_abcd,U_abdc,U_acbd,Vdir,Udir,Vexc,Uexc,&
!$OMP&         Vdel,Udel,i,i1,i2,i3,nac,&
!$OMP&         gamma_c_dir,gamma_c_exc,delta_c_dir)
!$OMP DO SCHEDULE(DYNAMIC)
    do nzac = 1,(nzx+1)**2
       Zrho_n_dir = 0
       Zrho_p_dir = 0
       Zrho_n_exc = 0
       Zrho_p_exc = 0
       Zkap_n_dir = 0
       Zkap_p_dir = 0
       zc = mod(nzac-1,nzx+1)
       za = (nzac-zc-1)/(nzx+1)
       do itb = 1,nttx
          rb = nrr(itb); lb = nll(itb); sb = nss(itb); npb = noo(itb)
          nlb = mod(npb-lb+1,2)
          do itd = 1,nttx
             npd = noo(itd)
             if(npb.ne.npd) cycle
             rd = nrr(itd); ld = nll(itd); sd = nss(itd)
             nld = mod(npd-ld+1,2)
             do izb = 1,nzzx(itb)
                zb = nzz(itb,izb)
                ibb = ib_zrls(zb,rb,lb,(sb+1)/2)
                idb = id(ibb)
                imb=ia(ibb); ibbx=ibb+nbx
                jb = i_zrls(zb,rb,lb,(sb+1)/2)
                nb = jb-imb
                do izd = 1,nzzx(itd)
                   zd = nzz(itd,izd)
                   if(mod(za+zb+zc+zd,2).ne.0) cycle
                   ibd = ib_zrls(zd,rd,ld,(sd+1)/2)
                   if(ibb.ne.ibd) cycle
                   jd = i_zrls(zd,rd,ld,(sd+1)/2)
                   imd=ia(ibd)
                   nd = jd - imd
                   ndb = nd + (nb-1)*idb
                   rho_n = rk(ndb,ibb )*0.5_pr !rho_db
                   rho_p = rk(ndb,ibbx)*0.5_pr
                   kap_n = ak(ndb,ibb )        !kappa_db
                   kap_p = ak(ndb,ibbx)
                   iz_abdc=zindex(za,zb,zd,zc)
                   if(zc.eq.zd) then
                      iz_abcd=iz_abdc
                   else
                      iz_abcd=zindex(za,zb,zc,zd)
                   endif
                   if(zd.eq.zb) then
                      iz_acbd=iz_abdc
                   elseif(zc.eq.zb) then
                      iz_acbd=iz_abcd
                   else
                      iz_acbd=zindex(za,zc,zb,zd)
                   endif
                   do ig = 1,n_g_all
                      V_abcd = VzGogny(ig,iz_abcd)
                      V_abdc = VzGogny(ig,iz_abdc)
                      V_acbd = VzGogny(ig,iz_acbd)
                      Zrho_n_dir(ig,rd,rb,nld,nlb,npb) = &
                           Zrho_n_dir(ig,rd,rb,nld,nlb,npb)+V_abcd*rho_n
                      Zrho_p_dir(ig,rd,rb,nld,nlb,npb) = &
                           Zrho_p_dir(ig,rd,rb,nld,nlb,npb)+V_abcd*rho_p
                      Zrho_n_exc(ig,rd,rb,nld,nlb,npb) = &
                           Zrho_n_exc(ig,rd,rb,nld,nlb,npb)+V_abdc*rho_n
                      Zrho_p_exc(ig,rd,rb,nld,nlb,npb) = &
                           Zrho_p_exc(ig,rd,rb,nld,nlb,npb)+V_abdc*rho_p
                      Zkap_n_dir(ig,rd,rb,nld,nlb,npb) = &
                           Zkap_n_dir(ig,rd,rb,nld,nlb,npb)+V_acbd*kap_n
                      Zkap_p_dir(ig,rd,rb,nld,nlb,npb) = &
                           Zkap_p_dir(ig,rd,rb,nld,nlb,npb)+V_acbd*kap_p
                   enddo !ig
                enddo ! izd
             enddo !izb
          enddo !itd
       enddo !itb
       do ita = 1,ntx
          if(za.ne.nz(ita)) cycle
          ra = nr(ita); la = nl(ita); sa = ns(ita)
          iba = ib_zrls(za,ra,la,(sa+1)/2)
          ida=id(iba); ima=ia(iba); ibax=iba+nbx
          na = ita-ima
          do itc = 1, ita
             if(zc.ne.nz(itc)) cycle
             rc = nr(itc); lc = nl(itc); sc = ns(itc)
             ibc = ib_zrls(zc,rc,lc,(sc+1)/2)
             if(ibc.ne.iba) cycle
             imc=ia(ibc)
             nc = itc-imc
             sac = sa + sc
             gamma_n_dir = zero
             gamma_p_dir = zero
             gamma_n_exc = zero
             gamma_p_exc = zero
             delta_n_dir = zero
             delta_p_dir = zero
             gamma_c_dir = zero
             gamma_c_exc = zero
             delta_c_dir = zero
             do ob = 0,nox
                ibb = ob + 1
                do nlb = 0,1
                   lb = ob+nlb
                   sb = -2*nlb+1
                   do nld = 0,1
                      ld = ob+nld
                      if(la+lb.ne.lc+ld.and.la-lb.ne.lc-ld) cycle
                      sd = -2*nld+1
                      sdb = sd + sb
                      do rb = 0,nrx
                         if(2*rb+lb.gt.nrlx) cycle
                         do rd = 0,nrx
                            if(2*rd+ld.gt.nrlx) cycle
                            if(sac.ne.0) then
                               ir_abdc=rindex(&
                                    ra,rb,rd,rc,la, lb, ld, lc,n_rlx)
                               if(rb.eq.rd.and.lb.eq.ld) then
                                  ir_acbd=ir_abdc
                               else
                                  ir_acbd=rindex(&
                                       ra,rc,rb,rd,la,-lc,lb,-ld,n_rlx)
                               endif
                               if(rc.eq.rd.and.lc.eq.ld) then
                                  ir_abcd = ir_abdc
                               elseif(rc.eq.rb.and.lc.eq.lb) then
                                  ir_abcd = ir_acbd
                               else
                                  ir_abcd=rindex(&
                                       ra,rb,rc,rd,la,lb,lc,ld,n_rlx)
                               endif
                               if(lb.eq.0.and.ld.eq.0) then
                                  jr_abdc = ir_abdc
                                  jr_acbd = ir_acbd
                                  jr_abcd = ir_abcd
                               else
                                  jr_abdc=rindex(&
                                       ra,rb,rd,rc,la,-lb,-ld,lc,n_rlx)
                                  if(rb.eq.rd .and.lb.eq.ld) then
                                     jr_acbd = jr_abdc
                                  else
                                     jr_acbd=rindex(&
                                          ra,rc,rb,rd,la,-lc,-lb,ld,n_rlx)
                                  endif
                                  if(rc.eq.rb.and.lb.eq.lc) then
                                     jr_abcd = ir_acbd
                                  elseif(rc.eq.rd.and.lc.eq.ld) then
                                     jr_abcd = ir_abdc
                                  else
                                     jr_abcd=rindex(&
                                          ra,rb,rc,rd,la,-lb,lc,-ld,n_rlx)
                                  endif
                               endif
                            else
                               if(sa.eq.sb) then
                                  jr_abdc=rindex(&
                                       ra,rb,rd,rc,la,-lb,-ld, lc,n_rlx)
                                  ir_acbd=rindex(&
                                       ra,rc,rb,rd,la,-lc, lb,-ld,n_rlx)
                               else
                                  ir_abdc=rindex(&
                                       ra,rb,rd,rc,la, lb, ld,lc,n_rlx)
                                  jr_acbd=rindex(&
                                       ra,rc,rb,rd,la,-lc,-lb,ld,n_rlx)
                               endif
                            endif
                            do ig = 1,n_g
                               W = W_g(ig)
                               B = B_g(ig)
                               H = H_g(ig)
                               M = M_g(ig)
                               zrho_nd = Zrho_n_dir(ig,rd,rb,nld,nlb,ibb)
                               zrho_pd = Zrho_p_dir(ig,rd,rb,nld,nlb,ibb)
                               zrho_ne = Zrho_n_exc(ig,rd,rb,nld,nlb,ibb)
                               zrho_pe = Zrho_p_exc(ig,rd,rb,nld,nlb,ibb)
                               zkap_n  = Zkap_n_dir(ig,rd,rb,nld,nlb,ibb)
                               zkap_p  = Zkap_p_dir(ig,rd,rb,nld,nlb,ibb)
                               if(sac.ne.0.and.sdb.ne.0) then
                                  V_abdc=VrGogny(ig,ir_abdc)
                                  U_abdc=VrGogny(ig,jr_abdc)
                                  V_abcd=VrGogny(ig,ir_abcd)
                                  U_abcd=VrGogny(ig,jr_abcd)
                                  V_acbd=VrGogny(ig,ir_acbd)
                                  U_acbd=VrGogny(ig,jr_acbd)
                                  if(sa.eq.sb) then
                                     !up up up up|or|down down down down
                                     Vdir = V_abcd
                                     Udir = U_abcd
                                     Vexc = V_abdc
                                     Uexc = U_abdc
                                     Vdel = V_acbd
                                     Udel = U_acbd
                                  else
                                     !up up down down|or|down down up up
                                     Vdir = U_abcd
                                     Udir = V_abcd
                                     Vexc = U_abdc
                                     Uexc = V_abdc
                                     Vdel = U_acbd
                                     Udel = V_acbd
                                  endif
                                  gamma_n_dir = gamma_n_dir &
                                       + Vdir*((W+B-M-H)*zrho_nd +&
                                       (W+B)*zrho_pd)&
                                       + Udir*(    (W-H)*zrho_nd +&
                                       W*zrho_pd)
                                  gamma_p_dir = gamma_p_dir &
                                       + Vdir*((W+B-M-H)*zrho_pd +&
                                       (W+B)*zrho_nd)&
                                       + Udir*(    (W-H)*zrho_pd +&
                                       W*zrho_nd)
                                  gamma_n_exc = gamma_n_exc &
                                       + Vexc*((M+H-W-B)*zrho_ne +&
                                       (M+H)*zrho_pe)&
                                       + Uexc*(    (M-B)*zrho_ne +&
                                       M*zrho_pe)
                                  gamma_p_exc = gamma_p_exc &
                                       + Vexc*((M+H-W-B)*zrho_pe +&
                                       (M+H)*zrho_ne)&
                                       + Uexc*(    (M-B)*zrho_pe +&
                                       M*zrho_ne)
                                  delta_n_dir = delta_n_dir &
                                       + (Vdel*(W-H)+Udel*(M-B))*zkap_n
                                  delta_p_dir = delta_p_dir &
                                       + (Vdel*(W-H)+Udel*(M-B))*zkap_p
                               elseif(sac.eq.0.and.sdb.eq.0) then
                                  if(sa.eq.sb) then
                                     !up down down up|or|down up up down
                                     Vexc = -VrGogny(ig,jr_abdc)
                                     Vdel =  VrGogny(ig,ir_acbd)
                                  else
                                     !up down up down|or|down up down up
                                     Vexc =  VrGogny(ig,ir_abdc)
                                     Vdel = -VrGogny(ig,jr_acbd)
                                  endif
                                  gamma_n_exc = gamma_n_exc &
                                       + Vexc*((H-W)*zrho_ne + H*zrho_pe)
                                  gamma_p_exc = gamma_p_exc &
                                       + Vexc*((H-W)*zrho_pe + H*zrho_ne)
                                  delta_n_dir = delta_n_dir &
                                       + Vdel*zkap_n*(W-H+B-M)
                                  delta_p_dir = delta_p_dir &
                                       + Vdel*zkap_p*(W-H+B-M)
                               endif
                            enddo !ig
                            do ig = n_g+1,n_g_all
                               W = V_g_coul(ig-n_g) * e2charg
                               zrho_pd = Zrho_p_dir(ig,rd,rb,nld,nlb,ibb)
                               if(icou.le.-4) then
                                  zrho_pe = 0._pr
                                  zkap_p  = 0._pr
                               elseif(icou.le.-3) then
                                  zrho_pe = Zrho_p_exc(ig,rd,rb,nld,nlb,ibb)
                                  zkap_p  = Zkap_p_dir(ig,rd,rb,nld,nlb,ibb)
                               elseif(icou.le.-2) then
                                  zrho_pe = Zrho_p_exc(ig,rd,rb,nld,nlb,ibb)
                                  zkap_p  = 0._pr
                               else
                                  zrho_pe = 0._pr
                                  zkap_p  = 0._pr
                               endif
                               if(sac.ne.0.and.sdb.ne.0) then
                                  if(sa.eq.sb) then
                                     !up up up up|or|down down down down
                                     Vdir = VrGogny(ig,ir_abcd)
                                     Udir = VrGogny(ig,jr_abcd)
                                     if(icou.le.-4) then
                                        Vexc = 0._pr
                                        Vdel = 0._pr
                                     elseif(icou.le.-3) then
                                        Vexc = VrGogny(ig,ir_abdc)
                                        Vdel = VrGogny(ig,ir_acbd)
                                     elseif(icou.le.-2) then
                                        Vexc = VrGogny(ig,ir_abdc)
                                        Vdel = 0._pr
                                     else
                                        Vexc = 0._pr
                                        Vdel = 0._pr
                                     endif
                                  else
                                     !up up down down|or|down down up up
                                     Vdir = VrGogny(ig,jr_abcd)
                                     Udir = VrGogny(ig,ir_abcd)
                                     if(icou.le.-4) then
                                        Vexc = 0._pr
                                        Vdel = 0._pr
                                     elseif(icou.le.-3) then
                                        Vexc = VrGogny(ig,jr_abdc)
                                        Vdel = VrGogny(ig,jr_acbd)
                                     elseif(icou.le.-2) then
                                        Vexc = VrGogny(ig,jr_abdc)
                                        Vdel = 0._pr
                                     else
                                        Vexc = 0._pr
                                        Vdel = 0._pr
                                     endif
                                  endif
                                  gamma_c_dir = gamma_c_dir &
                                       + (Vdir+Udir)*W*zrho_pd
                                  gamma_c_exc = gamma_c_exc &
                                       - Vexc*W*zrho_pe
                                  delta_c_dir = delta_c_dir &
                                       + Vdel*W*zkap_p
                               elseif(sac.eq.0.and.sdb.eq.0.and.icou.le.-2) then
                                  if(sa.eq.sb) then
                                     !up down down up|or|down up up down
                                     Vexc = -VrGogny(ig,jr_abdc)
                                     if(icou.eq.-3) then
                                        Vdel =  VrGogny(ig,ir_acbd)
                                     else
                                        Vdel = 0._pr
                                     endif
                                  else
                                     !up down up down|or|down up down up
                                     Vexc =  VrGogny(ig,ir_abdc)
                                     if(icou.eq.-3) then
                                        Vdel = -VrGogny(ig,jr_acbd)
                                     else
                                        Vdel = 0._pr
                                     endif
                                  endif
                                  gamma_c_exc = gamma_c_exc &
                                       - Vexc*W*zrho_pe
                                  delta_c_dir = delta_c_dir &
                                       + Vdel*zkap_p*W
                               endif
                            enddo !ig
                         enddo !rd
                      enddo !rb
                   enddo !nld
                enddo !nlb
             enddo !ob
             i=allibro(iba)+na*(na-1)/2+nc
             i1=i+nhhdim
             i2=i+nhhdim2
             i3=i+nhhdim3
             nac = na+(nc-1)*ida;
             gamma_g_dir(nac,iba )=gamma_n_dir
             gamma_g_dir(nac,ibax)=gamma_p_dir
             gamma_g_exc(nac,iba )=gamma_n_exc
             gamma_g_exc(nac,ibax)=gamma_p_exc
             delta_g_dir(nac,iba )=delta_n_dir
             delta_g_dir(nac,ibax)=delta_p_dir
             coulf_g_dir(nac,ibax)=gamma_c_dir
             coulf_g_exc(nac,ibax)=gamma_c_exc
             coulf_d_dir(nac,ibax)=delta_c_dir
             brout(i )=brout(i )+gamma_n_dir+gamma_n_exc
             brout(i1)=brout(i1)+gamma_p_dir+gamma_p_exc+gamma_c_dir+gamma_c_exc
             brout(i2)=brout(i2)+delta_n_dir
             brout(i3)=brout(i3)+delta_p_dir+delta_c_dir
          enddo !itc
       enddo !ita
    enddo !nzac
!$OMP End Do
!$OMP End Parallel
 else
!$OMP Parallel Default(None)&
!$OMP& SHARED(nzx,nttx,nrr,nll,nss,noo,nzzx,nzz,ib_zrls,id,ia,nbx,&
!$OMP&        i_zrls,rk,ak,VzGogny,ntx,nz,nr,nl,ns,nox,nrx,nrlx,n_rlx,&
!$OMP&        W_g,B_g,H_g,M_g,VrGogny,allibro,nhhdim,nhhdim2,nhhdim3,&
!$OMP&        gamma_g_dir,gamma_g_exc,delta_g_dir,brout,n_g,force_is_dme) &
!$OMP& PRIVATE(nzac,Zrho_n_dir,Zrho_p_dir,Zrho_n_exc,Zrho_p_exc,&
!$OMP&         Zkap_n_dir,Zkap_p_dir,zc,za,itb,rb,lb,sb,npb,nlb,itd,npd,&
!$OMP&         rd,ld,sd,nld,izb,zb,ibb,idb,imb,ibbx,jb,nb,izd,zd,ibd,jd,&
!$OMP&         imd,nd,ndb,rho_n,rho_p,kap_n,kap_p,iz_abdc,iz_abcd,&
!$OMP&         iz_acbd,ig,V_abcd,V_abdc,V_acbd,ita,ra,la,sa,iba,ida,ima,&
!$OMP&         ibax,na,itc,rc,lc,sc,ibc,imc,nc,sac,gamma_n_dir,&
!$OMP&         gamma_p_dir,gamma_n_exc,gamma_p_exc,delta_n_dir,&
!$OMP&         delta_p_dir,ob,sdb,ir_abdc,ir_abcd,ir_acbd,jr_abdc,&
!$OMP&         jr_abcd,jr_acbd,W,B,H,M,zrho_nd,zrho_pd,zrho_ne,zrho_pe,&
!$OMP&         zkap_n,zkap_p,U_abcd,U_abdc,U_acbd,Vdir,Udir,Vexc,Uexc,&
!$OMP&         Vdel,Udel,i,i1,i2,i3,nac)
!$OMP DO SCHEDULE(DYNAMIC)
    do nzac = 1,(nzx+1)**2
       Zrho_n_dir = 0
       Zrho_p_dir = 0
       zc = mod(nzac-1,nzx+1)
       za = (nzac-zc-1)/(nzx+1)
       do itb = 1,nttx
          rb = nrr(itb); lb = nll(itb); sb = nss(itb); npb = noo(itb)
          nlb = mod(npb-lb+1,2)
          do itd = 1,nttx
             npd = noo(itd)
             if(npb.ne.npd) cycle
             rd = nrr(itd); ld = nll(itd); sd = nss(itd)
             nld = mod(npd-ld+1,2)
             do izb = 1,nzzx(itb)
                zb = nzz(itb,izb)
                ibb = ib_zrls(zb,rb,lb,(sb+1)/2)
                idb = id(ibb)
                imb=ia(ibb); ibbx=ibb+nbx
                jb = i_zrls(zb,rb,lb,(sb+1)/2)
                nb = jb-imb
                do izd = 1,nzzx(itd)
                   zd = nzz(itd,izd)
                   if(mod(za+zb+zc+zd,2).ne.0) cycle
                   ibd = ib_zrls(zd,rd,ld,(sd+1)/2)
                   if(ibb.ne.ibd) cycle
                   jd = i_zrls(zd,rd,ld,(sd+1)/2)
                   imd=ia(ibd)
                   nd = jd - imd
                   ndb = nd + (nb-1)*idb
                   rho_n = rk(ndb,ibb )*0.5_pr !rho_db
                   rho_p = rk(ndb,ibbx)*0.5_pr
                   iz_abcd=zindex(za,zb,zc,zd)
                   do ig = 1,n_g
                      V_abcd = VzGogny(ig,iz_abcd)
                      Zrho_n_dir(ig,rd,rb,nld,nlb,npb) = &
                           Zrho_n_dir(ig,rd,rb,nld,nlb,npb)+V_abcd*rho_n
                      Zrho_p_dir(ig,rd,rb,nld,nlb,npb) = &
                           Zrho_p_dir(ig,rd,rb,nld,nlb,npb)+V_abcd*rho_p
                   enddo !ig
                enddo ! izd
             enddo !izb
          enddo !itd
       enddo !itb
       do ita = 1,ntx
          if(za.ne.nz(ita)) cycle
          ra = nr(ita); la = nl(ita); sa = ns(ita)
          iba = ib_zrls(za,ra,la,(sa+1)/2)
          ida=id(iba); ima=ia(iba); ibax=iba+nbx
          na = ita-ima
          do itc = 1, ita
             if(zc.ne.nz(itc)) cycle
             rc = nr(itc); lc = nl(itc); sc = ns(itc)
             ibc = ib_zrls(zc,rc,lc,(sc+1)/2)
             if(ibc.ne.iba) cycle
             imc=ia(ibc)
             nc = itc-imc
             sac = sa + sc
             if(sac.eq.0) cycle
             gamma_n_dir = zero
             gamma_p_dir = zero
             do ob = 0,nox
                ibb = ob + 1
                do nlb = 0,1
                   lb = ob+nlb
                   sb = -2*nlb+1
                   do nld = 0,1
                      ld = ob+nld
                      if(la+lb.ne.lc+ld.and.la-lb.ne.lc-ld) cycle
                      sd = -2*nld+1
                      sdb = sd + sb
                      if(sdb.eq.0) cycle
                      do rb = 0,nrx
                         if(2*rb+lb.gt.nrlx) cycle
                         do rd = 0,nrx
                            if(2*rd+ld.gt.nrlx) cycle
                            ir_abcd=rindex(&
                                 ra,rb,rc,rd,la,lb,lc,ld,n_rlx)
                            if(lb.eq.0.and.ld.eq.0) then
                               jr_abcd = ir_abcd
                            else
                               jr_abcd=rindex(&
                                    ra,rb,rc,rd,la,-lb,lc,-ld,n_rlx)
                            endif
                            do ig = 1,n_g
                               W = W_g(ig)
                               B = B_g(ig)
                               H = H_g(ig)
                               M = M_g(ig)
                               zrho_nd = Zrho_n_dir(ig,rd,rb,nld,nlb,ibb)
                               zrho_pd = Zrho_p_dir(ig,rd,rb,nld,nlb,ibb)
                               V_abcd=VrGogny(ig,ir_abcd)
                               U_abcd=VrGogny(ig,jr_abcd)
                               if(sa.eq.sb) then
                                  !up up up up|or|down down down down
                                  Vdir = V_abcd
                                  Udir = U_abcd
                               else
                                  !up up down down|or|down down up up
                                  Vdir = U_abcd
                                  Udir = V_abcd
                               endif
                               gamma_n_dir = gamma_n_dir &
                                    + Vdir*((W+B-M-H)*zrho_nd +&
                                    (W+B)*zrho_pd)&
                                    + Udir*(    (W-H)*zrho_nd +&
                                    W*zrho_pd)
                               gamma_p_dir = gamma_p_dir &
                                    + Vdir*((W+B-M-H)*zrho_pd +&
                                    (W+B)*zrho_nd)&
                                    + Udir*(    (W-H)*zrho_pd +&
                                    W*zrho_nd)
                            enddo !ig
                         enddo !rd
                      enddo !rb
                   enddo !nld
                enddo !nlb
             enddo !ob
             i=allibro(iba)+na*(na-1)/2+nc
             i1=i+nhhdim
             i2=i+nhhdim2
             i3=i+nhhdim3
             nac = na+(nc-1)*ida;
             gamma_g_dir(nac,iba )=gamma_n_dir
             gamma_g_dir(nac,ibax)=gamma_p_dir
             brout(i )=brout(i )+gamma_n_dir
             brout(i1)=brout(i1)+gamma_p_dir
          enddo !itc
       enddo !ita
    enddo !nzac
!$OMP End Do
!$OMP End Parallel
 endif
#else
    do ita = 1,ntx
       ra = nr(ita); za = nz(ita); la = nl(ita); sa = ns(ita)
       iba = ib_zrls(za,ra,la,(sa+1)/2)
       ida=id(iba); ima=ia(iba); ibax=iba+nbx
       na = ita-ima
       do itc = 1, ita
          rc = nr(itc); zc = nz(itc); lc = nl(itc); sc = ns(itc)
          ibc = ib_zrls(zc,rc,lc,(sc+1)/2)
          idc=id(ibc); imc=ia(ibc); ibcx=ibc+nbx
          nc = itc-imc
          sac = sa + sc
          if(ibc.ne.iba) cycle
          gamma_n_dir = zero
          gamma_p_dir = zero
          gamma_n_exc = zero
          gamma_p_exc = zero
          delta_n_dir = zero
          delta_p_dir = zero
          do itb = 1,ntx
             rb = nr(itb); zb = nz(itb); lb = nl(itb); sb = ns(itb)
             ibb = ib_zrls(zb,rb,lb,(sb+1)/2)
             idb=id(ibb); imb=ia(ibb); ibbx=ibb+nbx
             nb = itb-imb
             do itd = 1,ntx
                rd = nr(itd); zd = nz(itd); ld = nl(itd); sd = ns(itd)
                ibd = ib_zrls(zd,rd,ld,(sd+1)/2)
                idd=id(ibd); imd=ia(ibd); ibdx=ibd+nbx
                nd = itd-imd
                sdb = sd + sb
                if(ibb.ne.ibd) cycle
                ndb = nd + (nb-1)*idd
                rho_n = rk(ndb,ibb )*0.5_pr !rho_db
                rho_p = rk(ndb,ibbx)*0.5_pr
                kap_n = ak(ndb,ibb )        !kappa_db
                kap_p = ak(ndb,ibbx)
                do ig = 1,n_g
                   W = W_g(ig)
                   B = B_g(ig)
                   H = H_g(ig)
                   M = M_g(ig)
                   if(sac.ne.0.and.sdb.ne.0) then
                      V_abdc=Vr_Gogny(ig,ra,la,rb, lb,rd, ld,rc, lc)
                      U_abdc=Vr_Gogny(ig,ra,la,rb,-lb,rd,-ld,rc, lc)
                      V_abcd=Vr_Gogny(ig,ra,la,rb, lb,rc, lc,rd, ld)
                      U_abcd=Vr_Gogny(ig,ra,la,rb,-lb,rc, lc,rd,-ld)
                      V_acbd=Vr_Gogny(ig,ra,la,rc,-lc,rb, lb,rd,-ld)
                      U_acbd=Vr_Gogny(ig,ra,la,rc,-lc,rb,-lb,rd, ld)
                      if(sa.eq.sb) then
                         ! up up   up up   |or| down down down down
                         Vdir = V_abcd*Vz_Gogny(ig,za,zb,zc,zd)
                         Udir = U_abcd*Vz_Gogny(ig,za,zb,zc,zd)
                         Vexc = V_abdc*Vz_Gogny(ig,za,zb,zd,zc)
                         Uexc = U_abdc*Vz_Gogny(ig,za,zb,zd,zc)
                         Vdel = V_acbd*Vz_Gogny(ig,za,zc,zb,zd)
                         Udel = U_acbd*Vz_Gogny(ig,za,zc,zb,zd)
                      else
                         !up up  down down |or| down down up up
                         Vdir = U_abcd*Vz_Gogny(ig,za,zb,zc,zd)
                         Udir = V_abcd*Vz_Gogny(ig,za,zb,zc,zd)
                         Vexc = U_abdc*Vz_Gogny(ig,za,zb,zd,zc)
                         Uexc = V_abdc*Vz_Gogny(ig,za,zb,zd,zc)
                         Vdel = U_acbd*Vz_Gogny(ig,za,zc,zb,zd)
                         Udel = V_acbd*Vz_Gogny(ig,za,zc,zb,zd)
                      endif
                      gamma_n_dir = gamma_n_dir &
                           + Vdir*((W+B-M-H)*rho_n + (W+B)*rho_p)&
                           + Udir*(    (W-H)*rho_n +     W*rho_p)
                      gamma_p_dir = gamma_p_dir &
                           + Vdir*((W+B-M-H)*rho_p + (W+B)*rho_n)&
                           + Udir*(    (W-H)*rho_p +     W*rho_n)
                      gamma_n_exc = gamma_n_exc &
                           + Vexc*((M+H-W-B)*rho_n + (M+H)*rho_p)&
                           + Uexc*(    (M-B)*rho_n +     M*rho_p)
                      gamma_p_exc = gamma_p_exc &
                           + Vexc*((M+H-W-B)*rho_p + (M+H)*rho_n)&
                           + Uexc*(    (M-B)*rho_p +     M*rho_n)
                      delta_n_dir = delta_n_dir &
                           + (Vdel*(W-H)+Udel*(M-B))*kap_n
                      delta_p_dir = delta_p_dir &
                           + (Vdel*(W-H)+Udel*(M-B))*kap_p
                   elseif(sac.eq.0.and.sdb.eq.0) then
                      if(sa.eq.sb) then
                         !up down down up|or| down up up down
                         Vexc = -Vr_Gogny(ig,ra,la,rb,-lb,rd,-ld,rc,lc)&
                              *Vz_Gogny(ig,za,zb,zd,zc)
                         Vdel =  Vr_Gogny(ig,ra,la,rc,-lc,rb,lb,rd,-ld)&
                              *Vz_Gogny(ig,za,zc,zb,zd)
                      else
                         !up down up down |or| down up down up
                         Vexc = Vr_Gogny(ig,ra,la,rb,lb,rd,ld,rc,lc)&
                              *Vz_Gogny(ig,za,zb,zd,zc)
                         Vdel =-Vr_Gogny(ig,ra,la,rc,-lc,rb,-lb,rd,ld)&
                              *Vz_Gogny(ig,za,zc,zb,zd)
                      endif
                      gamma_n_exc = gamma_n_exc &
                           + Vexc*((H-W)*rho_n + H*rho_p)
                      gamma_p_exc = gamma_p_exc &
                           + Vexc*((H-W)*rho_p + H*rho_n)
                      delta_n_dir = delta_n_dir + Vdel*kap_n*(W-H+B-M)
                      delta_p_dir = delta_p_dir + Vdel*kap_p*(W-H+B-M)
                   endif
                enddo !ig
             enddo !itd
          enddo !itb
          if(force_is_dme) then
             !only the Hartree Part
             gamma_n_exc = zero
             gamma_p_exc = zero
             delta_n_dir = zero
             delta_p_dir = zero
          endif
          i=allibro(iba)+na*(na-1)/2+nc
          i1=i+nhhdim
          i2=i+nhhdim2
          i3=i+nhhdim3
          nac = na+(nc-1)*ida;
          gamma_g_dir(nac,iba )=gamma_n_dir
          gamma_g_dir(nac,ibax)=gamma_p_dir
          gamma_g_exc(nac,iba )=gamma_n_exc
          gamma_g_exc(nac,ibax)=gamma_p_exc
          delta_g_dir(nac,iba )=delta_n_dir
          delta_g_dir(nac,ibax)=delta_p_dir
          brout(i )=brout(i )+gamma_n_dir+gamma_n_exc
          brout(i1)=brout(i1)+gamma_p_dir+gamma_p_exc
          brout(i2)=brout(i2)+delta_n_dir
          brout(i3)=brout(i3)+delta_p_dir
       enddo !itc
    enddo !ita
#endif
  End Subroutine gamdel_gogny

  !=======================================================================
  !
  !======================================
  ! lib PnProjected specIfics Start >>>>>
  !======================================
  !=======================================================================
  !> Calculates expectation values (tz is the particle number) in the case
  !> of particle number projection
  !=======================================================================
  Subroutine expectpj(lpr)
    Use UNEDF
    Implicit None
    Logical :: lpr
    Integer(ipr) :: i,j,it,ihli,iw,iw1=901,iw2=902,icons,lambda,ign,igp,k
    Integer(ipr), Allocatable :: in_vec(:),ip_vec(:)
    Complex(pr) :: SZFIN,SFIZN,SRFIN,SFIRN,SZFIP,SFIZP,SRFIP,SFIRP
    Complex(pr) :: cekt(3),cdel(2),cept(3),cetot,etens,cq2pj(ilpj,ilpj)
    Complex(pr) :: cxn(2),crms(3),cq2(3),cq4(3),xnpj(2),rmspj(3),q2pj(3),q4pj(3)
    Complex(pr) :: evolpj,esurpj,ecdipj,ecexpj,ecoupj,ept1pj,ept2pj,epotpj,eki1pj, &
                   eki2pj,ekinpj,etotpj,espopj,epa1pj,epa2pj,epirpj,ede1pj,ede2pj, &
                   etenspj
    Complex(pr) :: eva,ev3,ev5,es5,eso,ecodi,ecoex,rn,rp,rnp1,rnp2,rt,rt2,tnt,tpt,tt, &
                   dn,dp,dt,akn,akp,akn2,akp2,adn,adp,evol,esurf,ecoul,pijk,row,cx,   &
                   dd1n,dd1p,ceptn,ceptp,cdeln,cdelp,cektn,cektp
    Complex(pr) :: rsa,rsa0
    Real(pr) :: whl,x,xn(3),q4(3),def(3),bet2(3),oct3(3),het4(3),r212,r222,rc,z,zz,rrr,p2,p3,p4 !MCedit oct3(3) for octupole beta3
    Real(pr) :: rdelta(2),repair(3),rekin(3),revolpj,resurpj,respopj,recdipj,recexpj,retenspj
    Real(pr) :: RpRMSsq,RnRMSsq,DarwinFoldy
    !
    Call densitpj ! calculates complex densities and the direct coulomb field
    !
    evolpj = zero; esurpj = zero; ecdipj = zero; ecexpj = zero; espopj  = zero;
    ept1pj = zero; ept2pj = zero; eki1pj = zero; eki2pj = zero; epj     = zero;
    etotpj = zero; epa1pj = zero; epa2pj = zero; ede1pj = zero; ede2pj  = zero;
    xnpj   = zero; rmspj  = zero; q2pj   = zero; q4pj   = zero; etenspj = zero; cq2pj = zero
    !
    Allocate(in_vec(1:ilpj*ilpj),ip_vec(1:ilpj*ilpj))
    i=0
    Do ign=1,ilpj
       Do igp=1,ilpj
          i=i+1
          in_vec(i) = ign
          ip_vec(i) = igp
       End Do
    End Do
    !
    ! Loop over both neutron and proton gauge angles. Could be multithreaded in the future
    Do k=1,ilpj*ilpj
       i = in_vec(k)
       j = ip_vec(k)

       pijk = pjk(i,1)*pjk(j,2)
       !
       cekt(:) = zero; cept(:) = zero; cdel(:) = zero
       cxn(:)  = zero; crms(:) = zero; cq2(:)  = zero; cq4(:) = zero
       !
       eva   = zero; ev3   = zero; ev5   = zero; es5   = zero
       eso   = zero; ecodi = zero; ecoex = zero; etens = zero
       ceptn = zero; ceptp = zero; cdeln = zero; cdelp = zero
       cektn = zero; cektp = zero
       !
       Do ihli = 1,nghl
          ! real
          whl = wdcor(ihli)
          z   = fh(ihli); zz = z*z; rrr = zz + fl(ihli)**2
          p2  = p32*zz   - half*rrr    !3/2 z*z-1/2 (z*z+r*r)=1/2(2 z*z-r*2)=1/2 Q
          p3  = p53*z*p2 - p23*rrr*z
          p4  = p74*z*p3 - p34*rrr*p2
          ! complex
          rn  = ropj(ihli,i,1); rp  = ropj(ihli,j,2); rnp2 = rn**2 + rp**2; rnp1=rn - rp
          ! ig - particle number, rms and deformations
          row = whl*rn; cxn(1)=cxn(1)+row; crms(1)=crms(1)+row*rrr; cq2(1)=cq2(1)+row*p2; cq4(1)=cq4(1)+row*p4
          row = whl*rp; cxn(2)=cxn(2)+row; crms(2)=crms(2)+row*rrr; cq2(2)=cq2(2)+row*p2; cq4(2)=cq4(2)+row*p4
          ! ig - energy contributions
          rt   = rn + rp;     rt2  = rt*rt
          tnt  = taupj(ihli,i,1); tpt  = taupj(ihli,j,2); tt = tnt + tpt
          dn   = dropj(ihli,i,1); dp   = dropj(ihli,j,2); dt = dn + dp
          akn  = akapj(ihli,i,1); akp  = akapj(ihli,j,2)
          akn2 = akn*akn;         akp2 = akp*akp
          adn  = akn*rn;          adp  = akp*rp
          ! ig-Pairing energy and delta
          rsa0=(rt/rho_c)
          dd1n=CpV0(0)*(ONE-rsa0*CpV1(0))*whl
          dd1p=CpV0(1)*(ONE-rsa0*CpV1(1))*whl
          !
          ceptn = ceptn + dd1n*akn2; ceptp = ceptp + dd1p*akp2
          cdeln = cdeln - dd1n*adn;  cdelp = cdelp - dd1p*adp
          !
          cektn = cektn + hb0n*whl*tnt; cektp = cektp + hb0p*whl*tpt   !kinetic (n and p)                              ! kinetic protons
          ev3     = ev3 + (tv1*rt2 - tv2*rnp2)*whl                     !volume
          eva     = eva + (tv3*rt2-tv4*rnp2)*rt**sigma*whl
          ev5     = ev5 + (tv5*rt*tt + tv6*(rn*tnt + rp*tpt))*whl
          es5     = es5 + (ts1*rt*dt + ts2*(rn*dn + rp*dp))*whl        !surface
          eso     = eso + (CrdJ(0)*rt+CrdJ(1)*rnp1)*djpj(ihli,i,1)*whl !spin-orbit n
          eso     = eso + (CrdJ(0)*rt-CrdJ(1)*rnp1)*djpj(ihli,j,2)*whl !spin-orbit p
          If(icou.Ge.1) ecodi = ecodi + half*coupj(ihli,j)*rp*whl      !Coulomb direct
          If(icou.Eq.2.Or.icou.Eq.-4) ecoex = ecoex - cex*rp**t4o3*whl !Coulomb exchange, Slater approximation
          If(use_j2terms) Then
             SFIZN=SFIZpj(IHLI,i,1); SFIRN=SFIRpj(IHLI,i,1); SZFIN=SZFIpj(IHLI,i,1); SRFIN=SRFIpj(IHLI,i,1)
             SFIZP=SFIZpj(IHLI,j,2); SFIRP=SFIRpj(IHLI,j,2); SZFIP=SZFIpj(IHLI,j,2); SRFIP=SRFIpj(IHLI,j,2)
             ETENS=ETENS+whl*(TA7*(SZFIN**2+SFIZN**2+SRFIN**2+SFIRN**2+SZFIP**2+SFIZP**2+SRFIP**2+SFIRP**2)&
                  +TA8*(SZFIN*SZFIP+SFIZN*SFIZP+SRFIN*SRFIP+SFIRN*SFIRP))
          End If
       End Do !ihli
       ! Intel compiler with -O2, -O3 does not like updating array element directly
       ! (like: cept(1)=cept(1)+...), hence the need for this construct
       cept(1) = ceptn; cept(2) = ceptp
       cdel(1) = cdeln; cdel(2) = cdelp
       cekt(1) = cektn; cekt(2) = cektp
       !
       evol     = ev3 + eva + ev5; esurf = es5; ecoul = ecodi  + ecoex*CExPar
       cekt(3)  = cekt(1) + cekt(2); cept(3)  = cept(1) + cept(2)
       cetot    = cekt(3) + evol + esurf + eso + ecoul + cept(3)+ ETENS
       cdel(1)  = cdel(1)/tz(1); cdel(2)  = cdel(2)/tz(2)
       !------------------------------------------------
       ! half-projected energies required for the matrix elements
       !------------------------------------------------
       epj(i,1) = epj(i,1) + cetot*pjk(j,2)
       epj(j,2) = epj(j,2) + cetot*pjk(i,1)
       !------------------------------------------------
       ! for constraint contributions to half-projected energies
       !------------------------------------------------
       If (icstr.Ne.0) cq2pj(i,j)=two*(cq2(1)+cq2(2))
       !------------------------------------------------
       ! projected energies
       !------------------------------------------------
       evolpj = evolpj + pijk*evol;    esurpj = esurpj + pijk*esurf;   espopj = espopj + pijk*eso
       epa1pj = epa1pj + pijk*cept(1); epa2pj = epa2pj + pijk*cept(2); epirpj = epa1pj + epa2pj
       ede1pj = ede1pj + pijk*cdel(1); ede2pj = ede2pj + pijk*cdel(2)
       ecdipj = ecdipj + pijk*ecodi;   ecexpj = ecexpj + pijk*ecoex;   ecoupj = ecdipj + ecexpj
       ept1pj = ept1pj + pijk*cept(1); ept2pj = ept2pj + pijk*cept(2); epotpj = ept1pj + ept2pj
       eki1pj = eki1pj + pijk*cekt(1); eki2pj = eki2pj + pijk*cekt(2); ekinpj = eki1pj + eki2pj
       !
       etotpj = etotpj + pijk*cetot
       etenspj= etenspj+ pijk*etens
       !------------------------------------------------
       ! unprojected hfb total energy and constraint
       !------------------------------------------------
       If(i.Eq.1.And.j.Eq.1) Then
          rehfbcan=Real(cetot,Kind=pr)
       End If
       !
       ! projected particle numbers, rms, deformations
       If(j.Eq.1) Then
          xnpj(1)  = xnpj(1)  + pjk(i,1)*cxn(1)
          rmspj(1) = rmspj(1) + pjk(i,1)*crms(1)
          q2pj(1)  = q2pj(1)  + pjk(i,1)*cq2(1)
          q4pj(1)  = q4pj(1)  + pjk(i,1)*cq4(1)
       End If
       If(i.Eq.1) Then
          xnpj(2)  = xnpj(2)  + pjk(j,2)*cxn(2)
          rmspj(2) = rmspj(2) + pjk(j,2)*crms(2)
          q2pj(2)  = q2pj(2)  + pjk(j,2)*cq2(2)
          q4pj(2)  = q4pj(2)  + pjk(j,2)*cq4(2)
       End If
       !
    End Do !k
    !
    ! Real quantities to the end
    !
    !------------------------------------------------
    ! Energies
    !------------------------------------------------
    rdelta(1) = Real(ede1pj,Kind=pr); rdelta(2) = Real(ede2pj,Kind=pr); retotpj   = Real(etotpj,Kind=pr);
    repair(1) = Real(epa1pj,Kind=pr); repair(2) = Real(epa2pj,Kind=pr); repair(3) = Real(epirpj,Kind=pr)
    rekin(1)  = Real(eki1pj,Kind=pr); rekin(2)  = Real(eki2pj,Kind=pr); rekin(3)  = Real(ekinpj,Kind=pr)
    revolpj   = Real(evolpj,Kind=pr); resurpj   = Real(esurpj,Kind=pr); respopj   = Real(espopj,Kind=pr);
    recdipj   = Real(ecdipj,Kind=pr); recexpj   = Real(ecexpj,Kind=pr); retenspj  = Real(etenspj,Kind=pr);
    depnp = retotpj - rehfbcan  !correlation energy due to projection
    !------------------------------------------------
    ! expectation values of multipole moments
    !------------------------------------------------
    Call moments_computeValue()
    !------------------------------------------------
    ! rms and deformations
    !------------------------------------------------
    Do it=itmin,itmax
       xn(it) = Real(xnpj(it),Kind=pr)
       rms(it)= Sqrt(Real(rmspj(it),Kind=pr)/xn(it))
       q2(it) = two*Real(q2pj(it),Kind=pr)    !Qnp=<2r^2P_2(teta)>=<2z^2-x^2-y^2>
       !MCcomment I believe the ffdef4 here should be replaced by "8.0_pr", the sqrt(117)/(4*pi) cannot be justified. 3/21/18
       !q4(it) = ffdef4*Real(q4pj(it),Kind=pr) !Hn=<8r^4P_4(teta)>=<8z^4-24z^2(x^2+y^2)+3(x^2+y^2)^2>
       q4(it) = eight*Real(q4pj(it),Kind=pr) !Hn=<8r^4P_4(teta)>=<8z^4-24z^2(x^2+y^2)+3(x^2+y^2)^2>
       def(it)= Sqrt(pi/5.0_pr)*q2(it)/(rms(it)**2*xn(it))
    End Do
    r212    = rms(1)**2; r222 = rms(2)**2
    rms(3)  = Sqrt((xn(1)*r212+xn(2)*r222)/amas)
    q2(3)   = q2(1) + q2(2)  ! quadrupole moment
    q4(3)   = q4(1) + q4(2)  ! hexadecapole moment
    def(3)  = Sqrt(pi/5.0_pr)*q2(3)/(rms(3)**2*amas) !deformation
        !------------------------------------------------
        ! other definitions of the same quantities
        !------------------------------------------------
    bet2(1) = ffdef6*q2(1)/(xn(1)*r02) !beta_n=Qn*Sqrt(5Pi)/(3N x^2)
    bet2(2) = ffdef6*q2(2)/(xn(2)*r02) !x=r0=1.2A^(1/3)
    bet2(3) = ffdef6*q2(3)/(amas*r02)
    het4(1) = ffdef7*q4(1)/(xn(1)*r04)
    het4(2) = ffdef7*q4(2)/(xn(2)*r04)
    het4(3) = ffdef7*q4(3)/(amas*r04)
    xn(3)   = xn(1) + xn(2)
    bet = def(3)
    !------------------------------------------------
    !  constraint constants and contributions to half-projected energies
    !------------------------------------------------
    If(icstr.Ne.0) Then
       cx=0.0_pr
       If (numberCons.Gt.0) Then
           Do icons=1,numberCons
              lambda=multLambda(icons)
              If(lambda.Ge.1) Then
                 cx = cx - multLag(lambda)*(qmoment(lambda,3)-multRequested(lambda))
              End If
              If(lambda.Eq.0) Then
                 cx = cx - neckLag*(neckValue-neckRequested)
              End If
           End Do
       End If
       !ty20=Sqrt(5.0_pr/pi)*hom/b0**2/two
       !cx=cqad*(cdef-bet)*ty20;
       Do i=1,ilpj
          Do j=1,ilpj
             epj(i,1) = epj(i,1) + cx*cq2pj(i,j)*pjk(j,2)
             epj(j,2) = epj(j,2) + cx*cq2pj(i,j)*pjk(i,1)
          End Do
       End Do
    End If
    !
    If (lpr) Then
       !rc=Sqrt(r222+0.640_pr)
       ! Charge radius, from Adv. Nucl. Phys. 8, 219 (1975)
       RpRMSsq=0.769_pr
       RnRMSsq=-0.1161_pr   ! J. Phys. G 33, 1 (2006)
       DarwinFoldy=0.033_pr ! Phys. Rev. A 56, 4579 (1997)
       rc = Sqrt(r222 + RpRMSsq + (xn(1)/xn(2))*RnRMSsq + DarwinFoldy)
       ! transitions to barn,barn^2,barn^4
       Do i=1,3
          q2(i)=q2(i)/100.0_pr; q4(i)=q4(i)/10000.0_pr
       End Do
       !
       ! STORE to projected buffer 'eresj'
       ! ieresj=50 from module definitions
       ! ' si ','JININ'
       eresj(1)=si; eresj(2)=inin;
       ! ' A','   N ','   Z '
       eresj(3)=npr(1)+npr(2); eresj(4)=npr(1); eresj(5)=npr(2);
       ! ' Jln ',' Jlp '
       eresj(6)=alast(1); eresj(7)=alast(2);
       ! ,'JEtot','Jbett','Jbetn','Jbetp',' JQt ',' JQn ',' JQp '  &
       eresj(8)=retotpj; eresj(9)=def(3); eresj(10)=def(1); eresj(11)=def(2);
       eresj(12)=q2(3); eresj(13)=q2(1); eresj(14)=q2(2);
       ! ' JpEn',' JpEp',' JpDn',' JpDp',' JAsn',' JAsp'  &
       eresj(15)=repair(1); eresj(16)=repair(2);
       eresj(17)=rdelta(1); eresj(18)=rdelta(2); eresj(19)=ass(1); eresj(20)=ass(2);
       ! ,' Jrt ',' Jrn ',' Jrp ',' Jrc ',' Jht ',' Jhn ',' Jhp '  &
       eresj(21)=rms(3); eresj(22)=rms(1); eresj(23)=rms(2); eresj(24)=rc;
       eresj(25)=het4(3); eresj(26)=het4(1); eresj(27)=het4(2);
       ! ,' Jqht',' Jqhn',' Jqhp'  &
       eresj(28)=q4(3); eresj(29)=q4(1); eresj(30)=q4(2);
       ! ,' JKINt',' JKINn','JKINp',' JSO ','JCDIR',' JCEX','JDisn','JDisp'  &
       eresj(31)=rekin(3); eresj(32)=rekin(1); eresj(33)=rekin(2); eresj(34)=respopj;
       eresj(35)=recdipj; eresj(36)=recexpj; eresj(37)=Dispersion(1); eresj(38)=Dispersion(2);
       ! ,'JV2Mn','JV2Mp','JILST','JKIND','  JL '  &
       eresj(39)=v2min(1); eresj(40)=v2min(2)
       eresj(41)=iLST; eresj(42)=kindhfb; eresj(43)=iLpj;
       !  ,'JECMPAV','JECMPAV','JECMPAV'
       eresj(44)=ECMPAV(3); eresj(45)=ECMPAV(1); eresj(46)=ECMPAV(2);
       ! 'JA','JN',JZ'
       eresj(47)=Nint(xn(3)); eresj(48)=Nint(xn(1)); eresj(49)=Nint(xn(2));
       ! 'iter'
       eresj(50)=iiter
       ! nucleus with wrong asymptotic
       If(iasswrong(3).Ne.0) eresj(21)=-eresj(21)
       !
       ! WRITE to screen 'lout' and tape akzout.dat 'lfile'
       Do iw=lout,lfile
          Write(iw,*)
          Write(iw,'(a,9x,a,/)')            '  NB! From expectpj (PNP PAV RESULTS)'
          Write(iw,*)
          If(iLST1.Ne.0)  &
               Write(iw,'(a,6f15.6)') '  hfb decay const. ass ',ass
          Write(iw,'(a,8f15.6)') '  pairing: CpV0,CpV1,pwi... ',CpV0,CpV1,pwi
          Write(iw,'(a,a,a,i3)') '  forces:   ',skyrme,',  Gauge points:',ilpj
          If(keyblo(1).Ne.0)  &
               Write(iw,'(a,i4,a,f10.3)')  '  Blocked neutron block    ',  &
               bloblo(keyblo(1),1)
          If(keyblo(2).Ne.0)  &
               Write(iw,'(a,i4,a,f10.3)')  '  Blocked proton  block    ',  &
               bloblo(keyblo(2),2)
          Write(iw,*)
          Write(iw,'(/,28x,a,8x,a,9x,a)') ' neutrons ','protons','total'
          Write(iw,'(a,6f15.6)') '  Requested part.numbs.',tz,Sum(tz)
          Write(iw,'(a,6f15.6)') '  Projected part.numbs.',xn
          Write(iw,'(a,3f15.6)') '  Dispersion dN2 ......',Dispersion
          Write(iw,'(a,6f15.6)') '  b0, bz, bp ..........',b0,bz,bp
          Write(iw,*)
          Write(iw,'(a,6f15.6)') '  lambda (ala) ........',ala
          Write(iw,'(a,6f15.6)') '  Lambda (alast) ......',alast
          Write(iw,'(a,6f15.6)') '  delta(n,p) ..........',rdelta
          Write(iw,'(a,6f15.6)') '  pairing energy ......',repair
          Write(iw,*)
          Write(iw,'(a,6f15.6)') '  rms-radius ..........',rms
          Write(iw,'(a,15x,2f15.6)') '  charge-radius, r0 ...',rc,r00
          Write(iw,'(a,6f15.6)') '  deformation beta2 ...',def
          Write(iw,'(a,6f15.6)') '  quadrupole moment[b] ',q2
          Write(iw,'(a,6f15.6)') '  hexadecapole moment .',q4
          Write(iw,*)
          Write(iw,'(a,6f15.6)')     '  kinetic energy ......',rekin
          Write(iw,'(a,6f15.6)')     '  cmc-diagonal part ...',rekin(1)/hb0n*hbzeron-rekin(1),&
               rekin(2)/hb0p*hbzerop-rekin(2),rekin(1)/hb0n*hbzeron-rekin(1)+rekin(2)/hb0p*hbzerop-rekin(2)
          Write(iw,'(a,6f15.6)')     '  cmc-PAV .............',ECMPAV
          Write(iw,*)
          Write(iw,'(a,30x,6f15.6)') '  volume energy .......',revolpj
          Write(iw,'(a,30x,6f15.6)') '  surface energy ......',resurpj
          Write(iw,'(a,30x,6f15.6)') '  spin-orbit energy ...',respopj
          Write(iw,'(a,30x,6f15.6)') '  coulomb direct ......',recdipj
          Write(iw,'(a,30x,6f15.6)') '  coulomb exchange ....',recexpj
          Write(iw,'(a,30x,6f15.6)') '  tensor energy .......',retenspj
          Write(iw,*)
          Write(iw,'(a,30x,f15.6)')  '  Energy: ehfb(qp) ....',ehfb
          Write(iw,'(a,30x,f15.6)')  '  Energy: ehfb(can,pj).',rehfbcan
          Write(iw,'(a,30x,f15.6)')  '  ehfb(qp)-ehfb(can,pj)',ehfb-rehfbcan
          Write(iw,'(a,30x,f15.6)')  '  Epj-ehfb(can,pj) ....',depnp
          Write(iw,'(a,30x,6f15.6)') '  Energy: Epj=E(PAV) ..',retotpj
          Write(iw,*)
       End Do
       !
       ! APPEND the results to file 'thodef.dat'
       ! ieres=ieresu+ieresl+ieresj+ierebl from module definitions
       If(iappend.Ne.0) Then
          ierest=0
          ! charge buffers
          Do i=1,ieresj       !charge projected buffer
             ierest=ierest+1
             eres(ierest)=eresj(i)
          End Do
          Do i=1,ieresu       !charge unprojected buffer
             ierest=ierest+1
             eres(ierest)=eresu(i)
          End Do
          Do i=1,ieresl       !charge LN  buffer
             ierest=ierest+1
             eres(ierest)=eresl(i)
          End Do
          Do i=1,ieresbl      !charge Blocking buffer
             ierest=ierest+1
             eres(ierest)=eresbl(i)
          End Do
          If(ierest.Ne.ieres) Then
             ierror_flag=ierror_flag+1
             ierror_info(ierror_flag)='STOP: In expectpj: ierest wrong'
             Return
          End If
          If(Print_Screen) Then
             ! recording results
100          Continue                        ! complications are due to eagle_ornl
             If(iLST1.Le.0) Then
#if(DO_MASSTABLE==1 || DRIP_LINES==1 || DO_PES==1)
                Open (unit=iw2,file='hodef'//row_string//'.dat',err=100,iostat=i,position='append')
#else
                Open (unit=iw2,file='hodef.dat',err=100,iostat=i,position='append')
#endif
                Write(iw2,'(3(1x,a,1x),160(1x,f14.6))') nucname,ereslbl,eres(1:ierest)
                Close(iw2)
             Else
                If (iasswrong(3).Eq.0) Then
#if(DO_MASSTABLE==1 || DRIP_LINES==1 || DO_PES==1)
                   Open (unit=iw1,file='thodef'//row_string//'.dat',err=100,iostat=i,position='append')
#else
                   Open (unit=iw1,file='thodef.dat',err=100,iostat=i,position='append')
#endif
                   Write(iw1,'(3(1x,a,1x),160(1x,f14.6))') nucname,ereslbl,eres(1:ierest)
                   Close(iw1)
                End If
             End If
          End If
       End If
    End If
    !
  End Subroutine expectpj
  !=======================================================================
  !> Calculate gauge-dependent densities
  !=======================================================================
  Subroutine densitpj
    Use UNEDF
    Implicit None
    !
    Complex(pr) :: tpfiu1,tpfid1,v2ig,dig,sumsum
    Complex(pr), Allocatable :: ank1(:,:),pfiun1(:,:),pfidn1(:,:)
    Complex(pr), Allocatable :: pakapj(:),propj(:), pdjpj(:), ptaupj(:),pdropj(:)
    Complex(pr), Allocatable :: pszfipj(:),psfizpj(:),psrfipj(:),psfirpj(:)
    Complex(pr), Pointer:: ppjk(:),pcpj(:,:),prpj(:,:),pypj(:,:)
    Real(pr) :: f,s,sd,su,sud,y,y2,sml2,cnzaa,cnraa,u,v2,tauin,xxx,yyy
    Real(pr) :: aav,anik,anik2,fi1r,fi1z,fi2d,qhla
    !,qhla,qh1la,qhl1a,qla,qha,fi1r,fi1z,fi2d
    Real(pr) :: xlam,xlam2,xlamy,xlamy2,xlap,xlap2,xlapy,xlapy2,xlampy
    Real(pr) :: tfiu,tfid,tfiur,tfidr,tfiuz,tfidz,tfiud2,tfidd2,tpfiu2,tpfid2,TW_T
    Real(pr), Allocatable :: an2(:),ank2(:),pfiun2(:),pfidn2(:)
    Integer(ipr) :: iw,nsa,nza,nra,nla,k,i,nd,il,ih,ihil,laplus,kkymu,n12
    Integer(ipr) :: imen,ib,m,im,ig,it,j,jj,ja,jn,ILIHLI,k1,k2,kkk,kky,mu,kkkmu
    Integer(ipr) :: k0(2),ky(2),kk(nqx),kyk(nqx)
    !
    Allocate(ank1(nqx,ilpj),pfiun1(ndx,ilpj),pfidn1(ndx,ilpj))
    Allocate(pfiun2(ndx),pfidn2(ndx),an2(nqx),ank2(nqx))
    Allocate(pakapj(ilnghl),propj(ilnghl), pdjpj(ilnghl), ptaupj(ilnghl),pdropj(ilnghl),  &
             pSZFIpj(ilnghl),pSFIZpj(ilnghl),pSRFIpj(ilnghl),pSFIRpj(ilnghl))
    !
    ! Projection grid points
    ! keypj=max(1,keypj); ilpj=keypj;  ilpj2=ilpj**2 !all
    ! when a value two*pi is used the results are precisely the same
    ! but the accuracy for even L is slow with increasing L.
    ! when 'pi' is used it gives regular and better convergence
    ! with respect to both, odd and even, L.
    ! Write(*,'(2x,a,i2,a,f12.8,a,f12.8)') 'point ig= ',i,' phi= ',yyy,' pi/2= ',pi/two
    xxx = pi/Real(ilpj,Kind=pr) ! equivalent to xxx = two*pi/Real(ilpj)
    Do i=1,ilpj
       yyy          = Real(i-1,Kind=pr)*xxx
       phypj(i)     = onei*yyy
       sinphy(i)    = onei*Sin(yyy)
       exp1iphy(i)  = Exp(onei*yyy)
       exp1iphym(i) = Exp(-onei*yyy)
       exp2iphy(i)  = Exp(two*onei*yyy)
       exp2iphym(i) = Exp(-two*onei*yyy)
    End Do
    !
    ! initialize parameters
    varmas = zero
    !
    Do it=itmin,itmax
       !
       ! zero for densities
       Do J=1,ilnghl
          pakapj(J)=zero; propj(J)=zero; pdjpj(J)=zero; ptaupj(J)=zero; pdropj(J)=zero;
       End Do
       Do J=1,ilnghl
          pszfipj(J)=zero; psfizpj(J)=zero; psrfipj(J)=zero; psfirpj(J)=zero;
       End Do
       !
       ! it-pointers
       prpj => rpj(:,:,it);  pcpj => cpj(:,:,it);
       pypj => ypj(:,:,it);  ppjk => pjk(:,it);
       !
       ! null for all pointers
       pypj=zero; prpj=zero; pcpj=zero; ppjk=one;
       !
       ! particle-init (kkk-even: 2 x number of pairs)
       kkk=npr(it); If(kkk.Ne.2*(kkk/2)) kkk=npr(it)-1
       ppjk(1:ilpj)=exp1iphym(1:ilpj)**kkk
       !
       ! start blocks
       k0(it)=0; ky(it)=0
       Do ib=1,nb
          nd=id(ib); im=ia(ib)
          If(Parity) Then
             LAPLUS=(ib+1)/2 !Yesp
          Else
             LAPLUS=ib       !Nop
          End If
          xlap=laplus; xlap2=xlap*xlap; xlam=xlap-one; xlam2=xlam*xlam
          !
          ! charge block can quantities
          m=ib+(it-1)*nbx; k1=ka(ib,it)+1; k2=ka(ib,it)+kd(ib,it); imen=0
          If(k1.Le.k2) Then
             ! below the pwi cut-off
             imen = nd
             !lcanon(ib,it)=lc
             Do k = 1,nd
                k0(it) = k0(it) + 1; kk(k)  = k0(it); kkk = k0(it)
                ky(it) = ky(it) + 1; kyk(k) = ky(it); kky = ky(it)
                aav    = vk(kkk,it)                                         ! v^2
                Do ig=1,ilpj
                   v2ig = exp2iphy(ig)*aav                                  ! gauged v^2
                   dig  = one - aav + v2ig                                  ! denominator
                   If(kkk.Ne.blocanon(it)) Then
                      ppjk(ig)  = ppjk(ig)*dig                              ! y(ig,it) <<<<<
                   End If
                   prpj(kkk,ig) = v2ig/dig                                  ! rho(mu,ig,it)
                   pcpj(kky,ig) = exp2iphy(ig)/dig                          ! c(mu,ig,it)
                   pypj(kky,ig) = exp1iphy(ig)/dig*onei*Sin(phypj(ig)/onei) ! sinphy(ig) !Y(mu,ig,it)
                End Do
             End Do
             ! At this point density (and related) are strictly equivalent in qp- and can-representation
             ! (up to 10^-14). Pairing density is not so strict (up to 10^-5) due to uv from v^2 but
             ! pairing density is taken directly in qp representation so both representations
             ! qp and can are strictly exact (up to 10^-14).
             j=0
             Do jj = 1,nd
                Do k = 1,nd
                   j=j+1; n12 = jj+(k-1)*nd;
                   an2(j)  = ddc(jj,kk(k),it)
                   ank2(j) = ak(n12,m)                 ! half \tilde{\rho} in q.p. basis
                   Do ig=1,ilpj
                      ank1(j,ig) = zero
                   End Do
                   Do mu=1,nd                          ! for half e^(-i\phy)*C(\phy)*\tilde{\rho} in q.p. basis
                      kkkmu = kk(mu); kkymu=kyk(mu)
                      Do ig=1,ilpj                     ! e^(-i\phy)*C in q.p. basis
                         ank1(j,ig) = ank1(j,ig) + ddc(jj,kkkmu,it)*ddc(k,kkkmu,it)*pcpj(kkymu,ig)*exp1iphym(ig)
                      End Do
                   End Do
                End Do
             End Do
          Else
             ! above the pwi cut-off (NB! Attention)
             ! here imem=0 and the contribution does
             ! not enter the densities but the Hamiltonian matrix
             ! used only in VAP regime
             ky(it)=ky(it)+1; kky = ky(it)
             Do ig=1,ilpj
                pcpj(kky,ig) = exp2iphy(ig)
                pypj(kky,ig) = exp1iphy(ig)*sinphy(ig)
             End Do
          End If
          !
          ! calculate the densities only below the PWI cutoff
          If (imen.Gt.0) Then
             ! gauss integration points
             Do il=1,ngl
                v2 = half/xl(il)
                Do ih=1,ngh
                   ihil = ih + (il-1)*ngh; ilihli=(ihil-1)*ilpj
                   !u = xh(ih); y = fli(ihil); y2=y*y
                   u = xh(ih); y = y_opt(ihil); y2=y*y
                   xlamy=xlam*y; xlamy2=xlam2*y2;
                   xlapy=xlap*y; xlapy2=xlap2*y2;
                   xlampy=xlamy+xlapy
                   !
                   ! initialize spin up/down funct
                   Do k=1,nd
                      fiu(k)=zero; fiuz(k)=zero; fiur(k)=zero; fiud2n(k)=zero; pfiun2(k)=zero;
                      fid(k)=zero; fidz(k)=zero; fidr(k)=zero; fidd2n(k)=zero; pfidn2(k)=zero;
                      Do ig=1,ilpj
                         pfiun1(k,ig)=zero; pfidn1(k,ig)=zero
                      End Do
                   End Do
                   !
                   ! scan over basis states
                   jn=0
                   Do i=1,nd
                      ja = i+im; nla = nl(ja); nra = nr(ja); nza = nz(ja); nsa = ns(ja);
                      sml2  = nla*nla; cnzaa = nza+nza+1; cnraa = nra+nra+nla+1
                      QHLA=QHLA_opt(JA,ihil); FI2D=FI2D_opt(JA,ihil)
                      FI1Z=FI1Z_opt(JA,ihil); FI1R=FI1R_opt(JA,ihil)

                      !qha   = qh(nza,ih); qla = ql(nra,nla,il); qhla = qha*qla
                      !qhl1a = qha*ql1(nra,nla,il)*v2; qh1la = qh1(nza,ih)*qla
                      !fi1z  = fp1(ihil)*qhla+fp2(ihil)*qh1la+fp3(ihil)*qhl1a
                      !fi1r  = fp4(ihil)*qhla+fp5(ihil)*qh1la+fp6(ihil)*qhl1a
                      !fi2d  = (fs1(ihil)*qh1la**2 + four*fs4(ihil)*qh1la*qhl1a        &
                      !      +  fs2(ihil)*qhl1a**2 + two*(fs5(ihil)*qh1la              &
                      !      +  fs6(ihil)*qhl1a)*qhla + ((u*u - cnzaa)*fs1(ihil)       &
                      !      +  (p14-cnraa*v2+sml2*v2*v2)*fs2(ihil)+fs3(ihil))*qhla**2 &
                      !      -  two*(fi1r**2+fi1z**2))/(two*qhla)
                      !
                      ! wave function(spin:up,down; grad:r,z,d2)
                      If (nsa.Gt.0) Then
                         Do k=1,nd
                            jn = jn+1; anik = an2(jn); anik2 = ank2(jn)
                            Do ig=1,ilpj
                               pfiun1(k,ig) = pfiun1(k,ig) + ank1(jn,ig)*qhla
                            End Do
                            pfiun2(k) = pfiun2(k) + anik2*qhla
                            fiu(k)    = fiu(k)    + anik*qhla
                            fiur(k)   = fiur(k)   + anik*fi1r
                            fiuz(k)   = fiuz(k)   + anik*fi1z
                            fiud2n(k) = fiud2n(k) + anik*fi2d
                            !
                         End Do
                      Else
                         Do k=1,nd
                            jn = jn+1; anik = an2(jn); anik2 = ank2(jn)
                            Do ig=1,ilpj
                               pfidn1(k,ig) = pfidn1(k,ig) + ank1(jn,ig)*qhla
                            End Do
                            pfidn2(k) = pfidn2(k) + anik2*qhla
                            fid(k)    = fid(k)    + anik*qhla
                            fidr(k)   = fidr(k)   + anik*fi1r
                            fidz(k)   = fidz(k)   + anik*fi1z
                            fidd2n(k) = fidd2n(k) + anik*fi2d
                            !
                         End Do
                      End If
                   End Do ! i
                   !
                   ! calculate densities
                   Do k=1,nd
                      kkk =kk(k)
                      tfiu=fiu(k); tfiuz=fiuz(k); tfiur=fiur(k); tfiud2=fiud2n(k); tpfiu2=pfiun2(k);
                      tfid=fid(k); tfidz=fidz(k); tfidr=fidr(k); tfidd2=fidd2n(k); tpfid2=pfidn2(k);
                      Do ig=1,ilpj
                         I=ig+ilihli; v2ig=prpj(kkk,ig); tpfiu1=pfiun1(k,ig); tpfid1=pfidn1(k,ig)
                         !
                         pakapj(I)  = pakapj(I)  +  (tpfiu1*tpfiu2+tpfid1*tpfid2)
                         propj(I)   = propj(I)   +  (tfiu**2+tfid**2)*v2ig
                         pdjpj(I)   = pdjpj(I)   +  (tfiur*tfidz-tfidr*tfiuz+xlamy*tfiu*(tfiur-tfidz) &
                                                 -   xlapy*tfid*(tfidr+tfiuz))*v2ig
                         TW_T=(tfiur**2+tfidr**2+tfiuz**2+tfidz**2)
                         tauin      = (xlamy2*tfiu**2+xlapy2*tfid**2+TW_T)
                         ptaupj(I)  = ptaupj(I)  +   tauin*v2ig
                         pdropj(I)  = pdropj(I)  +  (TW_T + tfiu*tfiud2 + tfid*tfidd2)*v2ig
                         psrfipj(I) = psrfipj(I) + (tfiur*tfid - tfidr*tfiu)*v2ig
                         psfirpj(I) = psfirpj(I) + (tfiu*tfid*xlampy)*v2ig
                         psfizpj(I) = psfizpj(I) + (xlamy*tfiu**2 - xlapy*tfid**2)*v2ig
                         pszfipj(I) = pszfipj(I) + (tfiuz*tfid - tfidz*tfiu)*v2ig
                         !
                      End Do !ig
                   End Do !k
                End Do !ih
             End Do !il
          End If
       End Do !ib
       !
       ! normalized pjk
       sumsum = Sum(ppjk(1:ilpj)); ppjk(1:ilpj) = ppjk(1:ilpj)/sumsum
       !
       ! Y minus second term of Y
       Do k=1,ky(it)
          sumsum = Sum(ppjk(1:ilpj)*pypj(k,1:ilpj))
          pypj(k,1:ilpj) = pypj(k,1:ilpj) - sumsum
       End Do
       !
       ! norm of the projected/unprojected density
       s = zero; sd = zero; su = zero; sud = zero;
       Do ihil=1,nghl
          ilihli=(ihil-1)*ilpj
          Do ig=1,ilpj
             I=ig+ilihli
             s=s+Real(two*propj(I)*ppjk(ig)); sd=sd+Real(four*pdropj(I)*ppjk(ig))
          End Do
          I=1+ilihli
          su=su+Real(two*propj(I)); sud=sud+Real(four*pdropj(I))
       End Do
       !
       ! print unprojected normalization
       Do iw=lout,lfile
          Write(iw,'(2(a,2(2x,D15.8)),(a,D15.8),a,i3)') &
               '   pj/unpj  s= ',s,su,'   pj/unpj sd= ',sd,sud,' ala1= ',ala1(it),' inner= ',inner(it)
       End Do
       varmas = varmas + su
       varmasNZ(it) = su; pjmassNZ(it) = s
       !
       s = Real(npr(it),Kind=pr)/s; dnfactor(it) = s; drhoi(it) = sd
       !
       Do ihil = 1,nghl
          ilihli=(ihil-1)*ilpj
          ! wdcor moves out the int.weight and multiply by the jacobian
          f = two*wdcori(ihil)
          Do ig=1,ilpj
             I=ig+ilihli
             ropj (ihil,ig,it)  = f*propj (I)
             taupj(ihil,ig,it)  = f*ptaupj(I)
             dropj(ihil,ig,it)  = f*pdropj(I)*two
             djpj (ihil,ig,it)  = f*pdjpj (I)*two
             akapj(ihil,ig,it)  = f*pakapj(I)*half
             SRFIpj(ihil,ig,it) = f*psrfipj(I)
             SFIRpj(ihil,ig,it) = f*psfirpj(I)
             SFIZpj(ihil,ig,it) = f*psfizpj(I)
             SZFIpj(ihil,ig,it) = f*pszfipj(I)
          End Do !ig
       End Do !ihil
       !
    End Do !it
    !
    dnfactor(3)=dnfactor(1)+dnfactor(2)
    !
    Deallocate(ank1,pfiun1,pfidn1)
    Deallocate(pakapj,propj, pdjpj, ptaupj,pdropj,pSZFIpj,pSFIZpj,pSRFIpj,pSFIRpj)
    !
    Call coulompj !complex coulomb fields
    !
  End Subroutine densitpj
  !=======================================================================
  !> Computes the Coulomb field (direct part) with Gauge-dependent densities
  !=======================================================================
  Subroutine coulompj
    Implicit None
    Integer(ipr) :: i,j,k
    Real(pr) :: zd2,y1,y2,xx1,s1,vik,f,r,r1,rr2,z,z1,zd1,t
    Real(pr) :: bb,r2,r12,rrr,rz1,rz2,rrz1,rrz2,xx
    Real(pr) :: bb1=3.5156229_pr,g1=0.39894228_pr,g7=0.02635537_pr,  &
         bb2=3.0899424_pr,g2=0.01328592_pr,g8=0.01647633_pr,  &
         bb3=1.2067492_pr,g3=0.00225319_pr,g9=0.00392377_pr,  &
         bb4=0.2659732_pr,g4=0.00157565_pr,bbxx=3.750_pr,     &
         bb5=0.0360768_pr,g5=0.00916281_pr,                   &
         bb6=0.0045813_pr,g6=0.02057706_pr
    If(icacoupj.Eq.0) Then
       icacoupj = 1; bb = Max(bp,bz)**4; f = chargee2/Sqrt(pi);  ! f=e^2/Sqrt(pi)
       Do i = 1,nghl
          r = fl(i); z = fh(i); r2 = r*r
          Do k = 1,i
             r1   = fl(k);       z1   = fh(k);      r12 = r1*r1
             rrr  = two*r*r1;    rr2  = (r - r1)**2
             zd1  = (z - z1)**2; zd2  = (z + z1)**2
             rz1  = r2+r12+zd1;  rz2  = r2+r12+zd2
             rrz1 = rr2+zd1;     rrz2 = rr2+zd2
             !
             xx1=zero
             Do j=1,nleg
                xx=Sqrt(one-xleg(j)**2); y1=(xleg(j)/(bb*xx))**2; s1=y1*rrr
                If(s1.Le.bbxx) Then
                   t=(s1/bbxx)**2; y2=one+t*(bb1+t*(bb2+t*(bb3+t*(bb4+t*(bb5+t*bb6)))))
                   y2=y2*(Exp(-rz1*y1)+Exp(-rz2*y1))
                Else
                   t=(bbxx/s1); y2=g1+t*(g2+t*(g3+t*(-g4+t*(g5+t*(-g6+t*(g7+t*(-g8+t*g9)))))))
                   y2=y2/Sqrt(s1)*(Exp(-rrz1*y1)+Exp(-rrz2*y1))
                End If
                xx1 = xx1 + wleg(j)*y2/(bb*xx**3)
             End Do
             vik=f*xx1; vc(i,k)=vik*wdcor(k); vc(k,i)=vik*wdcor(i)  !wdcor=pi*wh*wl*bz*bp*bp/fd
          End Do  !k
       End Do  !i
    End If
    ! calculation of the coulomb field
    coupj = zero
    Do i = 1,nghl
       Do k=1,ilpj
          coupj(:,k) = coupj(:,k) + vc(:,i)*ropj(i,k,2)
       End Do
    End Do
  End Subroutine coulompj
  !=======================================================================
  !> Modified Broyden's method: D.D.Johnson, PRB 38, 12807 (1988)
  !> Adopted from: (C) 2001 PWSCF group
  !=======================================================================
  Subroutine broyden_min(N,vout,vin,alpha,si,iter,M,bbroyden)
    !---------------------------------------------------------------------
    ! Input :
    !  N      dimension of arrays vin,vout
    !  vin    outpu at previous iteration
    !  vout   output at current iteration
    !  alpha  mixing factor (0 < alpha <= 1)
    !  iter   current iteration number
    !  M      number of iterations in Broyden history
    !  M=0    Linear mixing
    ! Output:
    !  si     MaxVal(|vout-vin|)
    !  vin    Broyden/Linear mixing result
    !  vout   vout-vin
    !  bbroyden='B' Broyden mixing, curvature>0
    !  bbroyden='L' Linear mixing,  curvature<0
    !---------------------------------------------------------------------
    Use HFBTHO_utilities, Only: pr,ipr
    Use HFBTHO, Only: ierror_flag,ierror_info
    Implicit None
    Integer(ipr),     Intent(In)    :: N,iter,M
    Real(pr),        Intent(In)     :: alpha
    Real(pr),        Intent(Out)    :: si
    Character(1),      Intent(Out)  :: bbroyden
    Real(pr),        Intent(InOut)  :: vout(N),vin(N)
    Integer(ipr)                    :: i,j,iter_used,ipos,inext,info
    Integer(ipr), Allocatable, Save :: iwork(:)
    Real(pr),    Allocatable, Save  :: beta(:,:),work(:)
    Real(pr),    Allocatable, Save  :: df(:,:),dv(:,:),curv(:)
    Real(pr),                 Save  :: w0
    Real(pr)                        :: DDOT,DNRM2,normi,gamma,curvature,sf
    !
    sf=-1.0_pr; Call DAXPY(N,sf,vin,1,vout,1)
    si=Maxval(Abs(vout))
    ! Linear mixing
    If(M.Eq.0.Or.iter.Eq.0) Then
       bbroyden='L'; Call DAXPY(N,alpha,vout,1,vin,1)
       !If(iter.Eq.0) Write(6,*) '  Linear mixing (alpha) : ',alpha
       Return
    End If
    ! Broyden mixing
    iter_used=Min(iter-1,M)
    ipos=iter-1-((iter-2)/M)*M
    inext=iter-((iter-1)/M)*M
    If (iter.Eq.1) Then
       w0=0.010_pr
       If(Allocated(df)) Deallocate(curv,df,dv,beta,work,iwork)
       Allocate(curv(N),df(N,M),dv(N,M),beta(M,M),work(M),iwork(M))
    Else
       df(:,ipos)=vout(:)-df(:,ipos); dv(:,ipos)=vin(:)-dv(:,ipos)
       Normi=1.0_pr/Sqrt((DNRM2(N,df(1,ipos),1))**2)
       Call DSCAL(N,Normi,df(1,ipos),1)
       Call DSCAL(N,Normi,dv(1,ipos),1)
    End If
    Do i=1,iter_used
       Do j=i+1,iter_used
          beta(i,j)=DDOT(N,df(1, j),1,df(1,i),1)
       End Do
       beta(i,i)=1.0_pr+w0*w0
    End Do
#if(SWITCH_ESSL==0)
    Call DSYTRF('U', iter_used, beta, M, iwork, work, M, info)
#else
    Call DPOTRF('U', iter_used, beta, M, info)
#endif
    If(info.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='STOP: In Broyden: info at DSYTRF '
       Return
    End If
#if(SWITCH_ESSL==0)
    Call DSYTRI('U', iter_used, beta, M, iwork, work, info)
#else
    Call DPOTRI('U', iter_used, beta, M, info)
#endif
    If(info.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='STOP: In Broyden: info at DSYTRI '
       Return
    End If
    Do i=1,iter_used
       Do j=i+1,iter_used
          beta(j,i)=beta(i,j)
       End Do
       work(i)=DDOT(N,df(1,i),1,vout,1)
    End Do
    curv=alpha*vout
    Do i=1,iter_used
       gamma=0.0_pr
       Do j=1,iter_used
          gamma=gamma+beta(j,i)*work(j)
       End Do
       curv=curv-gamma*(dv(:,i)+alpha*df(:,i))
    End Do
    Call DCOPY(N,vout,1,df(1,inext),1)
    Call DCOPY(N,vin,1,dv(1,inext),1)
    curvature=DDOT(N,vout,1,curv,1)
    If(curvature.Gt.-1.0_pr) Then
       bbroyden='B'; sf=+1.0_pr; Call DAXPY(N,sf,curv,1,vin,1)
    Else
       bbroyden='L'; sf=alpha*0.50_pr; Call DAXPY(N,sf,vout,1,vin,1)
    End If
  End Subroutine broyden_min
  !=======================================================================
  !> Calculates expectation values (xn is the particle number)
  !> for lpr=.true. also calculates PAV corrections
  !=======================================================================
  Subroutine expect(lpr)
    Use UNEDF, Only: FunctionalName   !EOedit: this subroutine needs to know the name of the functional
    Implicit None
    Logical :: lpr
    Integer(ipr) :: i,it,ihli,iw,LAMACT,I_TYPE
    Real(pr) :: ekt(3),xn(3),q4(3),def(3),bet2(3),oct3(3),het4(3),econst
    Real(pr) :: q3(3)   !EOedit for LN corrected octupole deformations
    Real(pr) :: z,zz,rrr,p2,p3,p4,row,r212,r222,rc
    Real(pr) :: eso,ecodi,ecoex,rn,rp,rnp1,rnp2,rt,whl,tnt,tpt,tt
    Real(pr) :: dn,dp,dt,akn,akp,akn2,akp2,adn,adp,evol,esurf,ecoul
    Real(pr) :: etens,dd1n,dd1p,rt1,tt1,dt1,djn,djp,djt,djt1
    Real(pr) :: RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,DJ_0,DJ_1,J2_0,J2_1,JabJba_0,JabJba_1
    Real(pr) :: SZFIN,SFIZN,SRFIN,SFIRN,SZFIP,SFIZP,SRFIP,SFIRP
    Real(pr) :: SZFI_0,SFIZ_0,SRFI_0,SFIR_0,SZFI_1,SFIZ_1,SRFI_1,SFIR_1
    Real(pr) :: SNABLARN,SNABLAZN,SNABLARP,SNABLAZP
    Real(pr) :: SNABLAR_0,SNABLAZ_0,SNABLAR_1,SNABLAZ_1
    Real(pr) :: xn1,xn2,rms1,rms2,q21,q22,q41,q42,EKIN_N,EKIN_P,ept1,ept2,del1,del2,rsa0
    Real(pr) :: ESURF_rho_DELTA_rho,ESURF_NABLA_rho_NABLA_rho,ESO_rho_NABLA_J,ESO_NABLA_rho_J
    Real(pr) :: EVOL_rho_tau,EVOL_rho_rho,EExtra,E_HARTREE_DIR,tempE_Crho0,tempREARR
    Real(pr) :: E_EXT_FIELD
    Real(pr) :: ezp     ! EOedit: approx cm correction for SV-min
    Real(pr) :: QLMLEF,QLMRIG
    Real(pr) :: RpRMSsq,RnRMSsq,DarwinFoldy
    Real(pr) :: e_gamma_fr_dir,e_gamma_fr_exc
    Real(pr) :: delrho1,delrho2,delkap1,delkap2
    Real(pr) :: cdelrho1,cdelrho2,cdelkap1,cdelkap2
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
    !integer ::  t1,t2,countrate,countmax
    !------------------------------------------------
    ! Part called during iterations (lpr=F)
    !------------------------------------------------
      !
    Call DENSIT
    If(ierror_flag.Ne.0) Return
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('expect',0)
    !
    !------------------------------------------------
    ! zero energy variables
    !------------------------------------------------
    EKIN_N=zero;              EKIN_P=zero;
    EVOL_rho_tau=zero;        EVOL_rho_rho=zero;
    ESURF_rho_DELTA_rho=zero; ESURF_NABLA_rho_NABLA_rho=zero;
    ESO_rho_NABLA_J=zero;     ESO_NABLA_rho_J=zero; E_HARTREE_DIR=zero
    ept1=zero; ept2=zero;     del1=zero; del2=zero;
    ecodi=zero; ecoex=zero; etens=zero
    EExtra=zero ; E_EXT_FIELD = zero ;
    xn1=zero; xn2=zero; rms1=zero; rms2=zero
    q21=zero; q22=zero; q41=zero; q42=zero
    tempE_Crho0=zero; tempREARR=zero
    DEROT=zero; SQUJ=zero; CRAN=zero; ERIGHFB=zero
    !------------------------------------------------
    ! zero optimization variables
    !------------------------------------------------
    If(DO_FITT) Then
       efit_0=zero; efitV0=zero; dfitV0=zero
       efit_rhorho=zero; efit_rhorhoD=zero;
       efit_rhotau=zero; efit_rhoDrho=zero;
       efit_rhonablaJ=zero; efit_JJ=zero;
    End If
    !------------------------------------------------
    ! Integration in coordinate space
    !------------------------------------------------
    Do ihli=1,nghl
       whl=wdcor(ihli)
       !------------------------------------------------
       ! np-representation
       !------------------------------------------------
       rn=ro(ihli,1);      rp=ro(ihli,2); rnp2=rn**2+rp**2; rnp1=rn-rp
       tnt=tau(ihli,1);    tpt=tau(ihli,2);
       dn=dro(ihli,1);     dp=dro(ihli,2);
       djn=dj(ihli,1);     djp=dj(ihli,2);
       akn=aka(ihli,1);    akp=aka(ihli,2)
       akn2=akn*akn;       akp2=akp*akp
       adn=akn*rn;         adp=akp*rp
       SFIZN=SFIZ(IHLI,1); SFIZP=SFIZ(IHLI,2);
       SFIRN=SFIR(IHLI,1); SFIRP=SFIR(IHLI,2);
       SZFIN=SZFI(IHLI,1); SZFIP=SZFI(IHLI,2);
       SRFIN=SRFI(IHLI,1); SRFIP=SRFI(IHLI,2);
       SNABLARN=NABLAR(IHLI,1); SNABLARP=NABLAR(IHLI,2);
       SNABLAZN=NABLAZ(IHLI,1); SNABLAZP=NABLAZ(IHLI,2);
       !------------------------------------------------
       ! t-representation
       !------------------------------------------------
       RHO_0=rn+rp;        RHO_1=rn-rp;
       TAU_0=tnt+tpt;      TAU_1=tnt-tpt;
       DRHO_0=dn+dp;       DRHO_1=dn-dp;
       DJ_0=djn+djp;       DJ_1=djn-djp;
       SFIZ_0=SFIZN+SFIZP; SFIZ_1=SFIZN-SFIZP;
       SFIR_0=SFIRN+SFIRP; SFIR_1=SFIRN-SFIRP;
       SZFI_0=SZFIN+SZFIP; SZFI_1=SZFIN-SZFIP;
       SRFI_0=SRFIN+SRFIP; SRFI_1=SRFIN-SRFIP;
       SNABLAR_0=SNABLARN+SNABLARP; SNABLAR_1=SNABLARN-SNABLARP;
       SNABLAZ_0=SNABLAZN+SNABLAZP; SNABLAZ_1=SNABLAZN-SNABLAZP;
       J2_0=SFIZ_0**2+SFIR_0**2+SZFI_0**2+SRFI_0**2
       J2_1=SFIZ_1**2+SFIR_1**2+SZFI_1**2+SRFI_1**2
       JabJba_0=2*(SFIZ_0*SZFI_0+SFIR_0*SRFI_0)
       JabJba_1=2*(SFIZ_1*SZFI_1+SFIR_1*SRFI_1)
       !
       Call calculate_U_parameters(RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,  &
            (SNABLAR_0**2+SNABLAZ_0**2),(SNABLAR_1**2+SNABLAZ_1**2))
       !------------------------------------------------
       ! rms and deformations
       !------------------------------------------------
       z=fh(ihli); zz=z*z; rrr=zz+fl(ihli)**2
       p2=p32*zz-half*rrr; p3=p53*z*p2-p23*rrr*z; p4=p74*z*p3-p34*rrr*p2
       row=whl*rn; xn1=xn1+row; rms1=rms1+row*rrr; q21=q21+row*p2; q41=q41+row*p4
       row=whl*rp; xn2=xn2+row; rms2=rms2+row*rrr; q22=q22+row*p2; q42=q42+row*p4
       !------------------------------------------------
       ! PH energies
       !------------------------------------------------
       EKIN_N=EKIN_N+hb0n*(TAU_0+TAU_1)*HALF*whl*facECM                         ! kinetic, n
       EKIN_P=EKIN_P+hb0p*(TAU_0-TAU_1)*HALF*whl*facECM                         ! kinetic, p
       EVOL_rho_tau=EVOL_rho_tau+(Urhotau(0,0)*RHO_0*TAU_0  &                   ! volume rho tau
            +Urhotau(1,0)*RHO_1*TAU_1+Urhotau(2,0)*RHO_0*TAU_1  &
            +Urhotau(3,0)*RHO_1*TAU_0 )*whl
       EVOL_rho_rho=EVOL_rho_rho+(Urhorho(0,0)*RHO_0**2  &                      ! volume density dependent
            +Urhorho(1,0)*RHO_1**2+(Urhorho(3,0)+Urhorho(2,0))*RHO_0*RHO_1)*whl
       ESURF_rho_DELTA_rho =ESURF_rho_DELTA_rho+(UrhoDrho(0,0)*RHO_0*DRHO_0  &  ! surface: rho delta rho
            +UrhoDrho(1,0)*RHO_1*DRHO_1+UrhoDrho(2,0)*RHO_0*DRHO_1  &
            +UrhoDrho(3,0)*RHO_1*DRHO_0 )*whl
       ESURF_NABLA_rho_NABLA_rho=ESURF_NABLA_rho_NABLA_rho  &                   ! surface: (nabla rho)**2
            +(Unablarho(0,0)*(SNABLAR_0*SNABLAR_0+SNABLAZ_0*SNABLAZ_0)  &
            +Unablarho(1,0)*(SNABLAR_1*SNABLAR_1+SNABLAZ_1*SNABLAZ_1)  &
            +(Unablarho(3,0)+Unablarho(2,0))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1) )*whl
       ESO_rho_NABLA_J=ESO_rho_NABLA_J+(UrhonablaJ(0,0)*RHO_0*DJ_0  &           ! spin-orbit rho Nabla . J
            +UrhonablaJ(1,0)*RHO_1*DJ_1+UrhonablaJ(2,0)*RHO_0*DJ_1  &
            +UrhonablaJ(3,0)*RHO_1*DJ_0 )*whl
       ESO_NABLA_rho_J=ESO_NABLA_rho_J  &
            +(UJnablarho(0,0)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))  &  ! spin-orbit J . Nabla rho
            + UJnablarho(1,0)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))  &
            + UJnablarho(2,0)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))  &
            + UJnablarho(3,0)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1)) )*whl
       ETENS=ETENS+(UJJ(0,0)*J2_0+UJJ(1,0)*J2_1  &                              ! tensor J^2
            +(UJJ(3,0)+UJJ(2,0))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)&
            +UJabJba(0,0)*JabJba_0 + UJabJba(1,0)*JabJba_1 &
            +(UJabJba(3,0)+UJabJba(2,0))*(SFIZ_0*SZFI_1+SFIR_0*SRFI_1+SZFI_0*SFIZ_1+SRFI_0*SFIR_1))*whl
       EExtra=EExtra+(UEnonstdr(0)+UEnonstdr(1))*whl                            ! extra field if needed
       E_EXT_FIELD=E_EXT_FIELD + ( Vexternal(0,zero,fl(ihli),z)*RHO_0 &         ! external field
            +Vexternal(1,zero,fl(ihli),z)*RHO_1 )*whl
       !------------------------------------------------
       ! Coulomb & Hartree
       !------------------------------------------------
       If (icou.Ge.1) ecodi=ecodi+half*cou(ihli)*rp*whl
       If (icou.Eq.2.Or.icou.Eq.-4) ecoex=ecoex-CExPar*cex*rp**p43*whl ! Slater approximation
       E_HARTREE_DIR=E_HARTREE_DIR +half*vDHartree(ihli,1)*RHO_0*whl+half*vDHartree(ihli,2)*RHO_1*whl
       ! just for printing
       tempE_Crho0=tempE_Crho0+RHO_0**2*whl
       tempREARR=tempREARR+(Cdrho(0)*RHO_0**2+Cdrho(1)*RHO_1**2)*RHO_0**sigma*whl
       !------------------------------------------------
       ! pairing energy and delta
       !------------------------------------------------
       rsa0=(RHO_0/rho_c)
       If(pairing_regularization .and. All(geff_inv .Ne. 0.0_pr)) Then
          ept1=ept1+akn2*whl/geff_inv(ihli,1); del1=del1-adn*whl/geff_inv(ihli,1)
          ept2=ept2+akp2*whl/geff_inv(ihli,2); del2=del2-adp*whl/geff_inv(ihli,2)
       Else
          dd1n=CpV0(0)*(ONE-rsa0*CpV1(0))*whl
          dd1p=CpV0(1)*(ONE-rsa0*CpV1(1))*whl
          ept1=ept1+dd1n*akn2; del1=del1-dd1n*adn
          ept2=ept2+dd1p*akp2; del2=del2-dd1p*adp
       End If
       !------------------------------------------------
       ! optimization quantities
       !------------------------------------------------
       If(DO_FITT) Then
          efitV0(0)=efitV0(0)+(ONE-rsa0*CpV1(0))*akn2*whl
          efitV0(1)=efitV0(1)+(ONE-rsa0*CpV1(1))*akp2*whl
          dfitV0(0)=dfitV0(0)-(ONE-rsa0*CpV1(0))*adn*whl
          dfitV0(1)=dfitV0(1)-(ONE-rsa0*CpV1(1))*adp*whl
          !
          efit_rhotau(0)=efit_rhotau(0)+RHO_0*TAU_0*whl              ! rho tau
          efit_rhotau(1)=efit_rhotau(1)+RHO_1*TAU_1*whl              ! rho tau
          efit_rhorho(0)=efit_rhorho(0)+RHO_0**2*whl                 ! rho^2
          efit_rhorho(1)=efit_rhorho(1)+RHO_1**2*whl                 ! rho^2
          efit_rhorhoD(0)=efit_rhorhoD(0)+RHO_0**sigma*RHO_0**2*whl  ! rho^2
          efit_rhorhoD(1)=efit_rhorhoD(1)+RHO_0**sigma*RHO_1**2*whl  ! rho^2
          efit_rhoDrho(0)=efit_rhoDrho(0)+RHO_0*DRHO_0*whl           ! rho Delta rho
          efit_rhoDrho(1)=efit_rhoDrho(1)+RHO_1*DRHO_1*whl           ! rho Delta rho
          efit_rhonablaJ(0)=efit_rhonablaJ(0)+RHO_0*DJ_0*whl         ! rho nabla J J
          efit_rhonablaJ(1)=efit_rhonablaJ(1)+RHO_1*DJ_1*whl         ! rho nabla J J
          efit_JJ(0)=efit_JJ(0)+J2_0*whl                             ! J.J
          efit_JJ(1)=efit_JJ(1)+J2_1*whl                             ! J.J
       End If
    End Do !ihli
    !------------------------------------------------
    ! Finite range traces
    !------------------------------------------------
    if(finite_range) then
       e_gamma_fr_dir = trace_product(gamma_g_dir,rk)*0.5_pr
       e_gamma_fr_exc = trace_product(gamma_g_exc,rk)*0.5_pr
       call trace_product_2(delta_g_dir,rk,delrho1,delrho2)
       call trace_product_2(delta_g_dir,ak,delkap1,delkap2)
       delrho1 = -delrho1; delrho2 = -delrho2
    else
       e_gamma_fr_dir = zero; e_gamma_fr_exc = zero
       delrho1 = zero; delrho2 = zero
       delkap1 = zero; delkap2 = zero
    endif
    !------------------------------------------------
    ! Coulomb Gaussians traces
    !------------------------------------------------
    if(coulomb_gaussian) then
       if(icou.le.-1) ecodi = trace_product(coulf_g_dir,rk)*0.5_pr
       if(icou.eq.-2.or.icou.eq.-3) ecoex = trace_product(coulf_g_exc,rk)*0.5_pr
       call trace_product_2(coulf_d_dir,rk,cdelrho1,cdelrho2)
       call trace_product_2(coulf_d_dir,ak,cdelkap1,cdelkap2)
       cdelrho1 = -cdelrho1; cdelrho2 = -cdelrho2
    else
       cdelrho1 = zero; cdelrho2 = zero
       cdelkap1 = zero; cdelkap2 = zero
    endif
    !------------------------------------------------
    ! after the integration
    !------------------------------------------------
    xn(1)=xn1;              xn(2)=xn2;         xn(3)=xn1+xn2;
    rms(1)=rms1;            rms(2)=rms2
    q2(1)=q21;              q2(2)=q22;
    q4(1)=q41;              q4(2)=q42
    ekt(1)=EKIN_N;          ekt(2)=EKIN_P;     ekt(3)=ekt(1)+ekt(2)
    ept(1)=ept1;            ept(2)=ept2;       ept(3)=ept(1)+ept(2)
    frept(1)=delkap1;       frept(2)=delkap2;  frept(3)=frept(1)+frept(2)
    coept(1)=cdelkap1;      coept(2)=cdelkap2; coept(3)=coept(1)+coept(2)
    del(1)=del1/xn(1);      del(2)=del2/xn(2);
    frdel(1)=delrho1/xn(1); frdel(2)=delrho2/xn(2);
    codel(1)=cdelrho1/xn(1);codel(2)=cdelrho2/xn(2);
    !
    ept = ept + frept + coept
    del = del + frdel + codel
    !
    EVOL=EVOL_rho_tau+EVOL_rho_rho+E_HARTREE_DIR
    esurf=ESURF_rho_DELTA_rho+ESURF_NABLA_rho_NABLA_rho
    ESO=ESO_rho_NABLA_J+ESO_NABLA_rho_J
    ecoul=ecodi+ecoex
    !--------------------------------------------
    !EOedit: Center of mass correction for SV-min
    !--------------------------------------------
    ezp = 17.3_pr/amas**0.2_pr
    !
    If(FunctionalName=='SV-min') Then
      etot=ekt(3)+evol+esurf+eso+ecoul+ept(3)+ETENS+EExtra+E_EXT_FIELD&
         + e_gamma_fr_dir+e_gamma_fr_exc-ezp
    Else
    etot=ekt(3)+evol+esurf+eso+ecoul+ept(3)+ETENS+EExtra+E_EXT_FIELD&
         + e_gamma_fr_dir+e_gamma_fr_exc
    End If
    !etot=ekt(3)+evol+esurf+eso+ecoul+ept(3)+ETENS+EExtra+E_EXT_FIELD&
    !     + e_gamma_fr_dir+e_gamma_fr_exc
    ehfb=etot
    entropy(3)=entropy(1)+entropy(2)
    !------------------------------------------------
    ! rms and deformations
    !------------------------------------------------
    Do it=itmin,itmax
       rms(it)=Sqrt(rms(it)/xn(it))
       q2(it)=two*q2(it)       !Qnp=<2r^2P_2(teta)>=<2z^2-x^2-y^2>
       !MCcomment I believe the ffdef4 here should be replaced by "8.0_pr", the sqrt(117)/(4*pi) cannot be justified. 3/21/18
       !q4(it)=ffdef4*q4(it)    !Hn=8r^4P_4(teta)=8z^4-24z^2(x^2+y^2)+3(x^2+y^2)^2
       q4(it)=eight*q4(it)    !Hn=8r^4P_4(teta)=8z^4-24z^2(x^2+y^2)+3(x^2+y^2)^2
       def(it)=Sqrt(pi/5.0_pr)*q2(it)/(rms(it)**2*xn(it))
    End Do
    r212=rms(1)**2; r222=rms(2)**2
    rms(3)=Sqrt((xn(1)*r212+xn(2)*r222)/amas)
    q2(3)=q2(1)+q2(2)          ! quadrupole moment
    q4(3)=q4(1)+q4(2)          ! hexadecapole moment
    def(3)=Sqrt(pi/5.0_pr)*q2(3)/(rms(3)**2*amas) !deformation
    bet=def(3)
    !bet=ffdef6*q2(3)/(amas*r02)  ! bet=Q2*Sqrt(5 Pi)/(3A x^2);  x=r0 A^(1/3)
    !------------------------------------------------
    ! Lipkin-Nogami energy
    !------------------------------------------------
    If(kindhfb.Lt.0) Then
       Call tracesln
       If(ierror_flag.Ne.0) Return
       etot=etot+etr(3)
    End If
    !------------------------------------------------
    ! optimization quantities
    !------------------------------------------------
    If(DO_FITT) Then
       efV_0=0.0_pr
       If(kindhfb.Lt.0) Then
          efV_0(0)=ala2(1)
          efV_0(1)=ala2(2)
       End If
       dfitV0(0)=dfitV0(0)/xn(1)
       dfitV0(1)=dfitV0(1)/xn(2)
       efit_0=etot-efitV0(0)*CpV0(0)-efitV0(1)*CpV0(1)  &
            -efit_rhotau(0)*Ctau(0)-efit_rhotau(1)*Ctau(1)  &
            -efit_rhorho(0)*Crho(0)-efit_rhorho(1)*Crho(1)  &
            -efit_rhorhoD(0)*Cdrho(0)-efit_rhorhoD(1)*Cdrho(1)  &
            -efit_rhoDrho(0)*CrDr(0)-efit_rhoDrho(1)*CrDr(1)  &
            -efit_rhonablaJ(0)*CrdJ(0)-efit_rhonablaJ(1)*CrdJ(1)  &
            -efit_JJ(0)*CJ(0)-efit_JJ(1)*CJ(1)
    End If
    !------------------------------------------------
    ! expectation values of multipole moments
    !------------------------------------------------
    Call moments_computeValue()
    !
    ! ATDHFB and GCM+GOA collective mass at the perturbative cranking approximation
    If (collective_inertia .And. lpr) Then
        Call calculate_collective_mass()
    End If
    !
    ! Fission fragment characteristics at convergence
    If (fission_fragments .And. lpr) Then
        ! Number of particles in the neck
        Call QNFIND()
        ! Position of the fragment centers of mass
        Call center_of_mass(Z_NECK,CENLEF,CENRIG)
        ! Mass multipole moments in the fragment intrinsic frame
        If(.Not.Allocated(QLMTOT)) Allocate(QLMTOT(0:lambdaMax,0:1),QLMPRO(0:lambdaMax,0:1))
        I_TYPE=1
        Do LAMACT=0,lambdaMax
           Call QLMFRA(Z_NECK,LAMACT,QLMLEF,QLMRIG,CENLEF,CENRIG,I_TYPE)
           QLMTOT(LAMACT,0) = QLMLEF
           QLMTOT(LAMACT,1) = QLMRIG
        End Do
        ! Charge of the fission fragments
        I_TYPE=2
        Do LAMACT=0,lambdaMax
           Call QLMFRA(Z_NECK,LAMACT,QLMLEF,QLMRIG,CENLEF,CENRIG,I_TYPE)
           QLMPRO(LAMACT,0) = QLMLEF
           QLMPRO(LAMACT,1) = QLMRIG
        End Do
    End If
    !------------------------------------------------
    ! debug
    !------------------------------------------------
    If(Print_Screen.And.IDEBUG.Gt.10) Then
       Write(*,'(4(a12,g13.6))')  &
            ' Tn=     ',ekt(1),           ' Tp=     ',ekt(2), &
            ' EPn=    ',ept(1),           ' EPp=    ',ept(2),  &
            ' EVOL=   ',EVOL,             ' Esurf=  ',esurf,  &
            ' NrNr=   ',ESURF_NABLA_rho_NABLA_rho,' rDr=    ',ESURF_rho_DELTA_rho,  &
            ' Etens=  ',ETENS,            ' Eso=   ',eso,  &
            ' rNJ=    ',ESO_rho_NABLA_J,  ' NrJ=   ',ESO_NABLA_rho_J,  &
            ' ECd=    ',ecodi,            ' ECex=  ',ecoex, &
            ' EHd=    ',E_HARTREE_DIR,    ' Ir0^2= ',tempE_Crho0, &
            ' Eextra= ',EExtra,           ' Ext.Fl= ',E_EXT_FIELD, &
            ' Etot=  ',etot
       If(DO_FITT) Then
          Write(*,'(4(a12,g13.6))')
          Write(*,'(4(a12,g13.6))')  &
               ' efrr0= ',efit_rhorho(0),     ' efrr1= ',efit_rhorho(1), &
               ' efrrD0=   ',efit_rhorhoD(0),          ' efrr1D=  ',efit_rhorhoD(1),  &
               ' efrt0= ',efit_rhotau(0),     ' efrt1= ',efit_rhotau(1), &
               ' efrDr0=   ',efit_rhoDrho(0),          ' efrDr1=  ',efit_rhoDrho(1),  &
               ' efrDj0=',efit_rhonablaJ(0),  ' efrDj1=',efit_rhonablaJ(1), &
               ' efjj0=    ',efit_JJ(0),               ' efjj1=   ',efit_JJ(1),  &
               ' efV0_0=',efitV0(0),          ' efV0_1=',efitV0(1), &
               ' dfV0_0=   ',dfitV0(0),                ' dfV0_1=  ',dfitV0(1),  &
               ' efV0=  ',efV_0(0),           ' efV_1= ',efV_0(1), &
               ' ef0=      ',efit_0,                   ' etot=    ',etot
       End If
    End If
    !------------------------------------------------
    ! Part called at the very end only (lpr=T)
    !------------------------------------------------
    If (lpr) Then
       !------------------------------------------------
       ! other definitions of deformations  (ffdef6=Sqrt(5.0_pr*pi)/3.0_pr)
       !------------------------------------------------
       bet2(1)=ffdef6*q2(1)/(xn(1)*r02) ! beta_n=Qn*Sqrt(5 Pi)/(3N x^2)
       bet2(2)=ffdef6*q2(2)/(xn(2)*r02) ! x=r0 A^(1/3)
       bet2(3)=ffdef6*q2(3)/(amas*r02)
       het4(1)=ffdef7*q4(1)/(xn(1)*r04)
       het4(2)=ffdef7*q4(2)/(xn(2)*r04)
       het4(3)=ffdef7*q4(3)/(amas*r04)
       ! Charge radius, from Adv. Nucl. Phys. 8, 219 (1975)
       RpRMSsq=0.769_pr
       RnRMSsq=-0.1161_pr   ! J. Phys. G 33, 1 (2006)
       DarwinFoldy=0.033_pr ! Phys. Rev. A 56, 4579 (1997)
       rc = Sqrt(r222 + RpRMSsq + (xn(1)/xn(2))*RnRMSsq + DarwinFoldy)
       ! transitions to barn,barn^2,barn^4
       Do i=1,3
          q2(i)=q2(i)/100.0_pr; q4(i)=q4(i)/10000.0_pr
       End Do
       !------------------------------------------------
       ! STORE to unprojected buffer 'eresu'
       !------------------------------------------------
       ! ieresu=50 from module definitions
       ! ,'UEtot','Ubett','Ubetn','Ubetp',' UQt ',' UQn ',' UQp '  &
       eresu(1)=etot; eresu(2)=def(3); eresu(3)=def(1); eresu(4)=def(2);
       eresu(5)=q2(3); eresu(6)=q2(1); eresu(7)=q2(2);
       ! ,' Uln ',' Ulp ',' UpEn',' UpEp',' UpDn',' UpDp',' UAsn',' UAsp'  &
       eresu(8)=alast(1); eresu(9)=alast(2); eresu(10)=ept(1); eresu(11)=ept(2);
       eresu(12)=del(1); eresu(13)=del(2); eresu(14)=ass(1); eresu(15)=ass(2);
       ! ,' Urt ',' Urn ',' Urp ',' Urc ',' Uht ',' Uhn ',' Uhp '  &
       eresu(16)=rms(3); eresu(17)=rms(1); eresu(18)=rms(2); eresu(19)=rc;
       eresu(20)=het4(3); eresu(21)=het4(1); eresu(22)=het4(2);
       ! ,' Uqht',' Uqhn',' Uqhp'  &
       eresu(23)=q4(3); eresu(24)=q4(1); eresu(25)=q4(2);
       ! ,'UKINT','UKINN','UKINP',' USO ','UCDIR',' UCEX','UDisn','UDisp'  &
       eresu(26)=ekt(3); eresu(27)=ekt(1); eresu(28)=ekt(2); eresu(29)=eso;
       eresu(30)=ecodi; eresu(31)=ecoex; eresu(32)=Dispersion(1); eresu(33)=Dispersion(2);
       ! ,'UV2Mn','UV2Mp'
       eresu(34)=v2min(1); eresu(35)=v2min(2);
       !  ,'UECMT','UECMN','UECMP'
       eresu(36)=ECMHFB(3); eresu(37)=ECMHFB(1); eresu(38)=ECMHFB(2);
       !  ,'UROTT','UROTN','UROTP'
       eresu(39)=DEROT(3); eresu(40)=DEROT(1); eresu(41)=DEROT(2);
       !  ,'USQUJT','USQUJN','USQUJP'
       eresu(42)=SQUJ(3); eresu(43)=SQUJ(1); eresu(44)=SQUJ(2);
       !  ,'UCRANT','UCRANN','UCRANP'
       eresu(45)=CRAN(3); eresu(46)=CRAN(1); eresu(47)=CRAN(2);
       !  ,'UERIGT','UERIGN','UERIGP'
       eresu(48)=ERIGHFB(3); eresu(49)=ERIGHFB(1); eresu(50)=ERIGHFB(2);
       !
       ! nucleus with wrong assymptotic
       If(iasswrong(3).Ne.0) eresu(16)=-eresu(16)
       !------------------------------------------------
       ! WRITE UNPROJECTED OUTPUT
       !------------------------------------------------
       Do iw=lout,lfile
          Write(iw,*)
          Write(iw,'(a,9x,a)')            '  NB! From expect (UNPROJECTED RESULTS)'
          Write(iw,*)
          If(iLST1.Ne.0)  &
               Write(iw,'(a,3f15.6)') '  hfb decay const. ass ',ass
          Write(iw,'(a,5f15.6)') '  pairing: CpV0,CpV1,...    ',CpV0,CpV1
          Write(iw,'(a,a)')      '  forces:   ',skyrme
          If(keyblo(1).Ne.0)  &
               Write(iw,'(a,i4,a,f10.3)')  '  Blocked neutron block    ',  &
               bloblo(keyblo(1),1)
          If(keyblo(2).Ne.0)  &
               Write(iw,'(a,i4,a,f10.3)')  '  Blocked proton  block    ',  &
               bloblo(keyblo(2),2)
          Write(iw,*)
          Write(iw,'(/,28x,a,8x,a,9x,a)') ' neutrons ','protons','total'
          Write(iw,'(a,6f15.6)')          '  Requested part.numbs.',tz,Sum(tz)
          Write(iw,'(a,6f15.6)')          '  UnPj(av) part.numbs .',xn
          Write(iw,'(a,3f15.6)')          '  b0, bz, bp ..........',b0,bz,bp
          Write(iw,*)
          Write(iw,'(a,3f15.6)') '  lambda (ala) ........',ala
          Write(iw,'(a,3f15.6)') '  Lambda (alast) ......',alast
          Write(iw,'(a,3f15.6)') '  delta(n,p), pwi .....',del,pwi
          If(finite_range) then
             Write(iw,'(a,3f15.6)') '  pairing energy ......',frept(3)
          Else
             Write(iw,'(a,3f15.6)') '  pairing energy ......',ept
          End If
          If(kindhfb.Lt.0) Then
             Write(iw,'(a,3f15.6)') '  LN lambda_2 ... ala2 ',ala2
             Write(iw,'(a,3f15.6)') '  LN energies .........',etr
             Write(iw,'(a,3f15.6)') '  delta(n,p)+ala2 .....',del+ala2
             Write(iw,'(a,3f15.6)') '  Geff(n,p) ...........',Geff
          End If
          Write(iw,*)
          Write(iw,'(a,3f15.6)') '  rms-radius ..........',rms
          Write(iw,'(a,15x,2f15.6)') '  charge-radius, r0 ...',rc,r00
          Write(iw,'(a,3f15.6)') '  deformation beta2....',def
          Write(iw,'(a,3f15.6)') '  dipole moment[fm] ...',(qmoment(1,it),it=1,3)
          Write(iw,'(a,3f15.6)') '  quadrupole moment[b] ',(qmoment(2,it),it=1,3)
          Write(iw,'(a,3f15.6)') '  octupole moment .....',(qmoment(3,it),it=1,3)
          Write(iw,'(a,3f15.6)') '  hexadecapole moment .',(qmoment(4,it),it=1,3)
          Write(iw,'(a,3f15.6)') '  q5 ..................',(qmoment(5,it),it=1,3)
          Write(iw,'(a,3f15.6)') '  q6 ..................',(qmoment(6,it),it=1,3)
          Write(iw,'(a,3f15.6)') '  q7 ..................',(qmoment(7,it),it=1,3)
          Write(iw,'(a,3f15.6)') '  q8 ..................',(qmoment(8,it),it=1,3)
          Write(iw,*)
          Write(iw,'(a,3f15.6)')    '  kinetic energy ......',ekt
          Write(iw,'(a,30x,f15.6)') '  volume energy .......',evol
          Write(iw,'(a,30x,f15.6)') '        rho_tau .......',EVOL_rho_tau
          Write(iw,'(a,30x,f15.6)') '        rho_rho .......',EVOL_rho_rho
          Write(iw,'(a,30x,f15.6)') '  surface energy ......',esurf
          Write(iw,'(a,30x,f15.6)') '   rho_DELTA_rho ......',ESURF_rho_DELTA_rho
          Write(iw,'(a,30x,f15.6)') '   (NABLA_rho)^2 ......',ESURF_NABLA_rho_NABLA_rho
          Write(iw,'(a,30x,f15.6)') '  spin-orbit energy ...',eso
          Write(iw,'(a,30x,f15.6)') '        rho_NABLA_J ...',ESO_rho_NABLA_J
          Write(iw,'(a,30x,f15.6)') '        NABLA_rho_J ...',ESO_NABLA_rho_J
          if(finite_range) then
             Write(iw,'(a,30x,f15.6)') ' finite range direct...',e_gamma_fr_dir
             Write(iw,'(a,30x,f15.6)') ' finite range exchange.',e_gamma_fr_exc
             Write(iw,'(a,30x,f15.6)') ' finite range total....',e_gamma_fr_dir+e_gamma_fr_exc
             Write(iw,'(a,30x,f15.6)') ' finite range pairing..',frept(3)
          endif
          Write(iw,'(a,30x,f15.6)') '  coulomb energy ......',ecodi+ecoex+coept(3)
          Write(iw,'(a,30x,f15.6)') '          direct ......',ecodi
          Write(iw,'(a,30x,f15.6)') '          exchange ....',ecoex
          Write(iw,'(a,30x,f15.6)') '          pairing .....',coept(3)
          Write(iw,'(a,30x,f15.6)') '  tensor energy .......',etens
          Write(iw,'(a,30x,f15.6)') '  direct Hartree E  ...',E_HARTREE_DIR
          Write(iw,'(a,30x,f15.6)') '  Extra E .............',EEXTRA
          Write(iw,'(a,30x,f15.6)') '  External field E ....',E_EXT_FIELD
          Write(iw,'(a,3f15.6)')    '  Entropy .............',entropy
          Write(iw,*)
          Write(iw,'(a,30x,f15.6)')    '  tEnergy: ehfb (qp)...',ehfb
          If(kindhfb.Lt.0) Then
             Write(iw,'(a,30x,f15.6)') '  tEnergy: ehfb(qp)+LN ',etot
          End If
          Write(iw,*)
          Write(iw,'(a,6f15.6)')    '  Calculated but not added corrections'
          Write(iw,'(a,6f15.6)')    '  ===================================='
          Write(iw,'(a,6f15.6)')    '  cmc-diagonal part ...',ekt(1)/hb0n*hbzeron-ekt(1),&
               ekt(2)/hb0p*hbzerop-ekt(2),ekt(1)/hb0n*hbzeron-ekt(1)+ekt(2)/hb0p*hbzerop-ekt(2)
          Write(iw,'(a,6f15.6)')    '  cmc-hfb .............',ECMHFB
          Write(iw,'(a,6f15.6)')    '  cranking rot corr ...',DEROT
          Write(iw,*)
          Write(iw,'(a,6f15.6)')    '  SQUJ ................',SQUJ
          Write(iw,'(a,6f15.6)')    '  CRAN x 4 ............',4.0_pr*CRAN
          Write(iw,'(a,6f15.6)')    '  Rigit Body ..........',ERIGHFB
          Write(iw,*)
          ! Print collective inertia mass tensor
          If(collective_inertia) Then
             Write(iw,'(a,6f15.6)')    '  ZPE (ATDHFB) ........',E0_ATD
             Write(iw,'(a,6f15.6)')    '  ZPE (GCM) ...........',E0_GCM
             Write(iw,*)
             Call print_collective(iw)
          End If
          ! Printing fision fragment characteristics
          If(fission_fragments) Call print_moments(iw)
       End Do
       !------------------------------------------------
       ! START corrected Lipkin-Nogami characteristics
       !------------------------------------------------
       If(kindhfb.Lt.0) Then
          Call densitln       !density LN corrections
          If(ierror_flag.Ne.0) Return
          Do it=itmin,itmax
             xn(it)=zero
             rms(it)=zero; q2(it)=zero; q4(it)=zero
             q3(it)=zero  !EOedit for LN corrected octupole deformations
          End Do
          !
          Do ihli=1,nghl
             whl=wdcor(ihli)
             rn=ro(ihli,1); rp=ro(ihli,2); rnp2=rn**2+rp**2
             ! rms and deformations
             z=fh(ihli); zz=z*z; rrr=zz+fl(ihli)**2
             !MCcomment, the coefficients for p2~4 are matched with legendre polynomials, thus later on the q2=two*q2 etc. is to get the numerator portion of spherical harmonics in spherical coordinates, which are the physical moments. rrr = x^2+y^2+z^2
             p2=p32*zz   -half*rrr
             p3=p53*z*p2 -p23*rrr*z
             p4=p74*z*p3 -p34*rrr*p2
             row=whl*rn
             xn(1)=xn(1)+row
             rms(1)=rms(1)+row*rrr
             q2(1)=q2(1)+row*p2
             q3(1)=q3(1)+row*p3  !EOedit for LN corrected octupole deformations
             q4(1)=q4(1)+row*p4
             row=whl*rp
             xn(2)=xn(2)+row
             rms(2)=rms(2)+row*rrr
             q2(2)=q2(2)+row*p2
             q3(2)=q3(2)+row*p3  !EOedit for LN corrected octupole deformations
             q4(2)=q4(2)+row*p4
          End Do !ihli
          !------------------------------------------------
          ! rms and deformations
          !------------------------------------------------
          Do it=itmin,itmax
             rms(it)=Sqrt(rms(it)/xn(it))
             !MCedit 3/21/18, q2~q4 should be the numerator portion of spherical harmonics Y_20~Y_40, in spherical coord.
             q2(it)=two*q2(it)     !Qnp=<2r^2P_2(teta)>=<2z^2-x^2-y^2>
             q3(it)=two*q3(it)     !On=<(2z^2-3x^2-3y^2)*z>
             !MCcomment I believe the ffdef4 here should be replaced by "8.0_pr" or "eight", the sqrt(117)/(4*pi) cannot be justified. 3/21/18
             !q4(it)=ffdef4*q4(it)  !Hn=<8r^4P_4(teta)>=<8z^4-24z^2(x^2+y^2)+3(x^2+y^2)^2>
             q4(it)=eight*q4(it)  !Hn=<8r^4P_4(teta)>=<8z^4-24z^2(x^2+y^2)+3(x^2+y^2)^2>
             def(it)=Sqrt(pi/5.0_pr)*q2(it)/(rms(it)**2*xn(it))
          End Do
          r212=rms(1)**2; r222=rms(2)**2
          r03=r00*r02
          rms(3)=Sqrt((xn(1)*r212+xn(2)*r222)/amas)
          q2(3)=q2(1)+q2(2)  ! quadrupole moment
          q3(3)=q3(1)+q3(2)  !EOedit for LN corrected octupole deformations
          q4(3)=q4(1)+q4(2)  ! hexadecapole moment
          def(3)=Sqrt(pi/5.0_pr)*q2(3)/(rms(3)**2*amas) !deformation
          ! other definitions of the same quantitsies, beta2~4
          ffdef65=Sqrt(seven*pi)/three !MCedit ffdef65 between ffdef6, ffdef7 for oct3 use
          !MCcomment, for future record, Q_n0 = C_n * 3/(4*pi) * A * r0^n * beta_n, A is neutron/proton/mass number, C_n is
          !the overall numerical portion of spherical harmonics
          !for example: C_2= 1/4 * sqrt(5/pi), C_3 = 1/4 * sqrt(7/pi), C_4 = 3/16 * sqrt(1/pi)
          !https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Spherical_harmonics
          bet2(1)=ffdef6*q2(1)/(xn(1)*r02) !beta_2=Q2*Sqrt(5Pi)/(3N x^2)
          bet2(2)=ffdef6*q2(2)/(xn(2)*r02) !x=r0=1.2A^(1/3)
          bet2(3)=ffdef6*q2(3)/(amas*r02)
          oct3(1)=ffdef65*q3(1)/(xn(1)*r03)    !MCedit beta3 for octupole deformation
          oct3(2)=ffdef65*q3(2)/(xn(2)*r03)
          oct3(3)=ffdef65*q3(3)/(amas*r03)
          het4(1)=ffdef7*q4(1)/(xn(1)*r04)
          het4(2)=ffdef7*q4(2)/(xn(2)*r04)
          het4(3)=ffdef7*q4(3)/(amas*r04)
          xn(3)=xn(1)+xn(2)
          bet=def(3)
          RpRMSsq=0.769_pr
          RnRMSsq=-0.1161_pr   ! J. Phys. G 33, 1 (2006)
          DarwinFoldy=0.033_pr ! Phys. Rev. A 56, 4579 (1997)
          rc = Sqrt(r222 + RpRMSsq + (xn(1)/xn(2))*RnRMSsq + DarwinFoldy)
          ! transitions to barn,barn^2,barn^4
          Do i=1,3
             q2(i)=q2(i)/100.0_pr; q4(i)=q4(i)/10000.0_pr
             q3(i)=q3(i)/1000.0_pr   !EOedit for LN corrected octupole deformations
          End Do
          !------------------------------------------------
          ! STORE to unprojected LN buffer 'eresl'
          !------------------------------------------------
          ! ieresl=20 from module definitions
          ! ,' EHFBLN',' EHFB',' LNEt','LNbet','LNben','LNbep',' LNQt',' LNQn',' LNQp'  &
          eresl(1)=etot; eresl(2)=etot-etr(3);
          eresl(3)=def(3); eresl(4)=def(1); eresl(5)=def(2)
          eresl(6)=q2(3); eresl(7)=q2(1); eresl(8)=q2(2);
          ! ,'LNpEn','LNpEp','LNpDn','LNpDp',' LNrt',' LNrn',' LNrC'  &
          eresl(9)=ept(1); eresl(10)=ept(2); eresl(11)=del(1)+ala2(1); eresl(12)=del(2)+ala2(2);
          eresl(13)=rms(3); eresl(14)=rms(1); eresl(15)=rms(2); eresl(16)=rc;
          ! ,' LNam2n',' LNam2p',' LNe2n',' LNe2p'
          eresl(17)=ala2(1); eresl(18)=ala2(2); eresl(19)=etr(1); eresl(20)=etr(2)
          !------------------------------------------------
          ! WRITE UNPROJECTED LN OUTPUT
          !------------------------------------------------
          Do iw=lout,lfile
             Write(iw,'(a,3f15.6)')
             Write(iw,'(a,3f15.6)') '  With Lipkin-Nogami Corrections'
             Write(iw,'(a,3f15.6)') '================================'
             Write(iw,'(a,3f15.6)') '  rms-radius ..........',rms
             Write(iw,'(a,15x,2f15.6)') '  charge-radius, r0 ...',rc,r00
             Write(iw,'(a,3f15.6)') '  deformation beta ....',def
             Write(iw,'(a,3f15.6)') '  q2 deformation,  beta',bet2
             Write(iw,'(a,3f15.6)') '  q3 deformation,  octa',oct3
             Write(iw,'(a,3f15.6)') '  q4 deformation,  heta',het4
             Write(iw,'(a,3f15.6)') '  quadrupole moment[b] ', q2
             Write(iw,'(a,3f15.6)') '  octupole moment[fm^3] ', q3
             Write(iw,'(a,3f15.6)') '  hexadecapole moment .', q4
             Write(iw,'(a,3f15.6)') '================================'
             Write(iw,'(a,3f15.6)')
          End Do
       End If
       !------------------------------------------------
       ! WRITE all blocking candidates
       !------------------------------------------------
       If(keyblo(3).Eq.0) Then
          Do iw=lout,lfile
             Write(iw,*)
             Do it=itmin,itmax
                If(it.Eq.1) Then
                   EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN
                Else
                   EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP
                End If
                !
                Write(iw,*) ' ',' Blocking candidates are:'
                Write(iw,*) '  ',protn(it),' eqpmin=',eqpmin(it),' pwiblo=',pwiblo
                Do i=1,blomax(it)
                   Write(iw,'(a,i4,a,i4,a,i4,2x,i4,3(a,1x,f12.8,1x),a)') '    num=',i,  &
                        ' block=',bloblo(i,it),  &
                        ' state=',blo123(i,it),blok1k2(i,it),  &
                        ' Eqp=',EqpPo(KqpPo(blok1k2(i,it))),  &
                        ' (1-2N)E=',(one-two*uk(blok1k2(i,it),it))*EqpPo(KqpPo(blok1k2(i,it))),  &
                        ' Ovlp=',vkmax(blok1k2(i,it),it),  &
                        tb(numax(blok1k2(i,it),it))
                End Do
                Write(iw,*)
             End Do
             Write(iw,*)
          End Do
       End If
       !
       !------------------------------------------------
       ! PAV
       !------------------------------------------------
       ! Projecting on different nucleus
       If(iproj.Ne.0) Then
          npr(1)=Int(npr1pj); npr(2)=Int(npr2pj)
          !tz(1)=Real(npr(1),Kind=pr); tz(2)=Real(npr(2),Kind=pr)
          Call expectpj(.True.)
       End If
    End If
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('expect',1)
    !
  End Subroutine expect
  !=======================================================================
  !> Determines if calculations include constraints on multipole moments
  !=======================================================================
  Subroutine Constraint_or_not(inin_INI0,inin0,icstr0)
    Implicit None
    Integer(ipr), Intent(in) :: inin_INI0
    Integer(ipr), Intent(inout) :: inin0,icstr0
    Integer(ipr) :: icount,l
    icount=0
    Do l=1,lambdaMax
       If(lambda_active(l).Gt.0) icount=icount+1
    End Do
    If(icount.Gt.0) Then
       icstr0=1; inin0=inin_INI0
    Else
       icstr0=0; inin0=inin_INI0
    End If
  End Subroutine Constraint_or_not
  !=======================================================================
  !> This routine updates the Lagrange multipliers of the multi-dimensional
  !> linear constraints based on the variation of the generalized density
  !> matrix and the QRPA matrix at the cranking approximation; see
  !> \cite berger1980 \cite younes2009-a \cite schunck2012
  !=======================================================================
  Subroutine getLagrange(ite)
    Implicit None
    Character(Len=1) :: trans
    Integer(ipr) :: ite,icons,lambda,icount,it,i,j,l,ierror
    Integer(ipr) :: ib,nd,nd2,nhfb,i0,m,k1,k2,n1,n2,nd1,k,kk,ll
    Integer(ipr) :: i_uvN,i_uvP,incx,incy
    Integer(ipr), allocatable :: ipivot(:),iftN(:),iftP(:)
    Real(pr) :: minu,hla,t_term,temp_k,temp_l,result,brakev,epsilo
    Real(pr), allocatable :: EqpN(:),EqpP(:)
    Real(pr), allocatable :: vecold(:),qmultt(:),veclam(:),veccns(:)
    Real(pr), allocatable :: cnsorg(:,:),cnsmat(:,:),cnsvec(:)
    Real(pr), allocatable :: fn12pl(:,:,:),fp12pl(:,:,:)
    Real(pr), allocatable :: fn11pl(:,:,:),fp11pl(:,:,:),fn11mi(:,:,:),fp11mi(:,:,:)
    Real(pr), allocatable :: doubln(:,:),doublp(:,:),dsum_n(:,:),dsum_p(:,:)
    Real(pr), allocatable :: workcn(:),dblmul(:,:),Umatr(:,:),Vmatr(:,:)
    !
    minu=-one
    epsilo=1.E-14
    !
    ! initializing the multipole moment template array
    Allocate(qmultt(0:lambdaMax));qmultt=zero
    Do lambda=0,lambdaMax
       qmultt(lambda)=qmoment(lambda,3)
    End Do
    !
    ! constructing the vector of the deviations of the current constraint
    ! from the requested values
    !
    Allocate(vecold(1:numberCons));vecold=zero
    Allocate(cnsvec(1:numberCons));cnsvec=zero
    Allocate(veclam(1:numberCons));veclam=zero
    !
    Do icons=1,numberCons
       lambda=multLambda(icons)
       ! regular multipole
       If(lambda.Ge.1) Then
          cnsvec(icons)=multRequested(lambda)-qmultt(lambda)
          If (nbroyden.lt.1) Then
              vecold(icons)=multLag(lambda)
          Else
              vecold(icons)=brin(nhhdim4+lambda)
          End If
       End If
       ! Gaussian neck
       If(lambda.Eq.0) Then
          cnsvec(icons)=neckRequested-neckValue
          If (nbroyden.lt.1) Then
              vecold(icons)=neckLag
          Else
              vecold(icons)=brin(nhhdim4+lambdaMax+1)
          End If
       End If
       veclam(icons)=vecold(icons)
    End Do
    !
    ! proceeding to determine the matrix of the constraint operators
    ! in the q.p. basis
    ! loop over the K blocks
    Allocate(cnsmat(numberCons,numberCons));cnsmat=zero
    Allocate(cnsorg(numberCons,numberCons));cnsorg=zero
    !
    i_uvN=0 ! new index referring to all q.p. vectors
    i_uvP=0 ! new index referring to all q.p. vectors
    !
    Do ib = 1,nb
       !
       !------------------------------------------------------
       ! matrix of the constraint in q.p. basis
       !------------------------------------------------------
       !
       !------------------------------------------------------
       ! neutron sector
       !------------------------------------------------------
       !
       it=1
       !
       nd=id(ib); nd2=nd*nd; nhfb=nd+nd; i0=ia(ib); m=ib+(it-1)*nbx
       !
       If(kd(ib,it).Gt.0) Then
          Allocate(doubln(nd,kd(ib,it))); doubln=zero
          Allocate(fn12pl(kd(ib,it),kd(ib,it),numberCons)); fn12pl=zero
          Allocate(Umatr(nd,kd(ib,it))); Umatr=zero
          Allocate(Vmatr(nd,kd(ib,it))); Vmatr=zero
          Allocate(EqpN(kd(ib,it))); EqpN=zero
          Allocate(ifTN(kd(ib,it))); ifTN=1
          ! temperature
          If(switch_on_temperature) Then
             Allocate(fn11pl(kd(ib,it),kd(ib,it),numberCons)); fn11pl=zero
          End If
          !
          ! U and V for this block (v. 101)
          Do k=1,kd(ib,it)
             ifTN(k)=ka(ib,it)+k; kk=KqpN(ka(ib,it)+k); EqpN(k)=REqpN(kk)
             Do n1=1,nd
                i_uvN=i_uvN+1
                Vmatr(n1,k)=RVqpN(i_uvN)
                Umatr(n1,k)=RUqpN(i_uvN)
             End Do
          End Do
          !
          Do icons=1,numberCons
             !
             lambda=multLambda(icons)
             If(lambda.Ge.1) Then
                Allocate(multMatElems(1:nd2)); multMatElems=zero
                Call moments_computeField(lambda,ib)
                ! matrix of the constraints in HO basis (size nd x nd)
                Allocate(dblmul(nd,nd));dblmul=zero
                j=0
                Do n1=1,nd
                   Do n2=1,n1
                      j=j+1;hla=multMatElems(j)
                      dblmul(n1,n2)=hla;dblmul(n2,n1)=hla
                   End Do
                End Do
             End If
             If(lambda.Eq.0) Then
                Allocate(gaussian_neck(1:nd2)); gaussian_neck=zero
                Call neck_computeField(ib)
                ! matrix of the constraints in HO basis (size nd x nd)
                Allocate(dblmul(nd,nd));dblmul=zero
                j=0
                Do n1=1,nd
                   Do n2=1,n1
                      j=j+1;hla=gaussian_neck(j)
                      dblmul(n1,n2)=hla;dblmul(n2,n1)=hla
                   End Do
                End Do
             End If
             !
             ! matrix of the constraint operator in the qp basis. due to
             ! the q.p. cut-off the actual size of the q.p. basis is not
             ! the same as the s.p. (ho) basis, and it is  not the  same
             ! for protons and neutrons.  the formulas implemented below
             ! differ from the 2 references for 3 reasons:
             !  - different  phase convention for the bogoliubov  matrix
             !  - block structure of the bogoliubov matrix in hfodd
             !  - storage in a() and b() arrays correspond to complex
             !    conjugate of the actual matrices
             !
             ! second term: v^{+} f^{*} u^{*} = v^{T} f u
             Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Umatr,nd,zero,doubln,nd)
             Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Vmatr,nd,doubln,nd,zero,fn12pl(1,1,icons),kd(ib,it))
             !
             ! first term:  u^{+} f v^{*} = u^{T} f v
             Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Vmatr,nd,zero,doubln,nd)
             Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Umatr,nd,doubln,nd,minu,fn12pl(1,1,icons),kd(ib,it))
             !
             ! temperature - computing \tilde{f}^{11}
             If(switch_on_temperature) Then
               !
               ! second term: v^{+} f^{*} v = v^{T} f v
               Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Vmatr,nd,zero,doubln,nd)
               Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Vmatr,nd,doubln,nd,zero,fn11pl(1,1,icons),kd(ib,it))
               ! first term:  u^{+} f u = u^{T} f u
               Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Umatr,nd,zero,doubln,nd)
               Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Umatr,nd,doubln,nd,minu,fn11pl(1,1,icons),kd(ib,it))
               !
             End If
             !
             If(Allocated(multMatElems)) Deallocate(multMatElems)
             If(Allocated(gaussian_neck)) Deallocate(gaussian_neck)
             Deallocate(dblmul)
             !
          End Do ! end icons (neutrons)
          !
          Deallocate(doubln,Umatr,Vmatr)
       End If
       !
       !------------------------------------------------------
       ! Proton sector
       !------------------------------------------------------
       !
       it=2
       !
       nd=id(ib); nd2=nd*nd; nhfb=nd+nd; i0=ia(ib); m=ib+(it-1)*nbx
       !
       If(kd(ib,it).Gt.0) Then
          Allocate(doublp(nd,kd(ib,it))); doublp=zero
          Allocate(fp12pl(kd(ib,it),kd(ib,it),numberCons)); fp12pl=zero
          Allocate(Umatr(nd,kd(ib,it))); Umatr=zero
          Allocate(Vmatr(nd,kd(ib,it))); Vmatr=zero
          Allocate(EqpP(kd(ib,it))); EqpP=zero
          Allocate(ifTP(kd(ib,it))); ifTP=1
          ! temperature
          If(switch_on_temperature) Then
             Allocate(fp11pl(kd(ib,it),kd(ib,it),numberCons)); fp11pl=zero
          End If
          !
          ! U and V for this block
          Do k=1,kd(ib,it)
             ifTP(k)=ka(ib,it)+k; kk=KqpP(ka(ib,it)+k); EqpP(k)=REqpP(kk)
             Do n1=1,nd
                i_uvP=i_uvP+1
                Vmatr(n1,k)=RVqpP(i_uvP)
                Umatr(n1,k)=RUqpP(i_uvP)
             End Do
          End Do
          !
          Do icons=1,numberCons
             !
             lambda=multLambda(icons)
             If(lambda.Ge.1) Then
                Allocate(multMatElems(1:nd2)); multMatElems=zero
                Call moments_computeField(lambda,ib)
                Allocate(dblmul(nd,nd));dblmul=zero
                j=0
                Do n1=1,nd
                   Do n2=1,n1
                      j=j+1;hla=multMatElems(j)
                      dblmul(n1,n2)=hla;dblmul(n2,n1)=hla
                   End Do
                End Do
             End If
             If(lambda.Eq.0) Then
                Allocate(gaussian_neck(1:nd2)); gaussian_neck=zero
                Call neck_computeField(ib)
                Allocate(dblmul(nd,nd));dblmul=zero
                j=0
                Do n1=1,nd
                   Do n2=1,n1
                      j=j+1;hla=gaussian_neck(j)
                      dblmul(n1,n2)=hla;dblmul(n2,n1)=hla
                   End Do
                End Do
             End If
             !
             ! matrix of the constraint operator in the qp basis. due to
             ! the q.p. cut-off the actual size of the q.p. basis is not
             ! the same as the s.p. (ho) basis, and it is  not the  same
             ! for protons and neutrons.  the formulas implemented below
             ! differ from the 2 references for 3 reasons:
             !  - different  phase convention for the bogoliubov  matrix
             !  - block structure of the bogoliubov matrix in hfodd
             !  - storage in a() and b() arrays correspond to complex
             !    conjugate of the actual matrices
             !
             ! second term: v^{+} f^{*} u^{*} = v^{t} f u
             Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Umatr,nd,zero,doublp,nd)
             Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Vmatr,nd,doublp,nd,zero,fp12pl(1,1,icons),kd(ib,it))
             !
             ! first term:  u^{+} f v^{*} = u^{t} f v
             Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Vmatr,nd,zero,doublp,nd)
             Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Umatr,nd,doublp,nd,minu,fp12pl(1,1,icons),kd(ib,it))
             !
             ! temperature - computing \tilde{f}^{11}
             If(switch_on_temperature) Then
               !
               ! second term: v f^{*} v = v^{T} f v
               Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Vmatr,nd,zero,doublp,nd)
               Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Vmatr,nd,doublp,nd,zero,fp11pl(1,1,icons),kd(ib,it))
               ! first term:  u f u = u^{T} f u
               Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Umatr,nd,zero,doublp,nd)
               Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Umatr,nd,doublp,nd,minu,fp11pl(1,1,icons),kd(ib,it))
               !
             End If
             !
             If(Allocated(multMatElems)) Deallocate(multMatElems)
             If(Allocated(gaussian_neck)) Deallocate(gaussian_neck)
             Deallocate(dblmul)
             !
          End Do ! end icons (protons)
          !
          Deallocate(doublp,Umatr,Vmatr)
       End If
       !
       !------------------------------------------------------
       ! constraint correlation matrix
       !------------------------------------------------------
       !
       !------------------------------------------------------
       ! neutron sector
       !------------------------------------------------------
       !
       it=1
       !
       If(kd(ib,it).Gt.0) Then
          Allocate(doubln(kd(ib,it),kd(ib,it))); doubln=zero
          Allocate(dsum_n(kd(ib,it),kd(ib,it))); dsum_n=zero
          !
          Do i=1,numberCons
             Do j=1,numberCons
                !
                ! temperature
                If((.Not.switch_on_temperature)) Then
                   !
                   Do l=1,kd(ib,it)
                      Do k=1,kd(ib,it)
                         If(Abs(EqpN(k)+EqpN(l)).Gt.Epsilo) Then
                            doubln(k,l)=fn12pl(k,l,i)/(EqpN(k)+EqpN(l))
                         Else
                            doubln(k,l)=zero
                         End If
                      End Do
                   End do
                   !
                   Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fn12pl(1,1,j),kd(ib,it),&
                                                         doubln,kd(ib,it),zero,dsum_n,kd(ib,it))
                Else
                   !
                   ! term corresponding to f^12
                   Do l=1,kd(ib,it)
                      Do k=1,kd(ib,it)
                         kk=iftN(k);ll=iftN(l)
                         temp_k=fn_T(kk)
                         temp_l=fn_T(ll)
                         If(Abs(EqpN(k)+EqpN(l)).Gt.Epsilo) Then
                            doubln(k,l)=fn12pl(k,l,i)*(one+temp_k+temp_l)/(EqpN(k)+EqpN(l))
                         Else
                            doubln(k,l)=zero
                         End If
                      End Do
                   End Do
                   !
                   Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fn12pl(1,1,j),kd(ib,it),&
                                                         doubln,kd(ib,it),zero,dsum_n,kd(ib,it))
                   !
                   ! first term: positive simplex
                   Do l=1,kd(ib,it)
                      Do k=1,kd(ib,it)
                         kk=iftN(k);ll=iftN(l)
                         temp_k=fn_T(kk)
                         temp_l=fn_T(ll)
                         If(k.ne.l.And.(Abs(EqpN(k)-EqpN(l)).Gt.Epsilo)) Then
                            t_term=-(temp_k-temp_l)/(EqpN(k)-EqpN(l))
                         Else
                            t_term=-temp_k*(temp_k-one)/temper
                         End If
                         doubln(k,l)=half*t_term*fn11pl(k,l,i)
                      End Do
                   End Do
                   !
                   Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fn11pl(1,1,j),kd(ib,it),&
                                                          doubln,kd(ib,it),one,dsum_n,kd(ib,it))
                   !
                End If
                !
                ! taking the trace of the resulting matrix
                !
                result=zero
                Do l=1,kd(ib,it)
                   result=result+dsum_n(l,l)
                End Do
                !
                cnsmat(i,j)=cnsmat(i,j)+result
                !
             End Do ! end of loop over j constraint
          End Do ! end of loop over i constraint
          !
          Deallocate(doubln,dsum_n,fn12pl,EqpN,ifTN)
          If(switch_on_temperature) Deallocate(fn11pl)
       End If
       !
       !------------------------------------------------------
       ! proton sector
       !------------------------------------------------------
       !
       it=2
       !
       If(kd(ib,it).Gt.0) Then
          Allocate(doublp(kd(ib,it),kd(ib,it))); doublp=zero
          Allocate(dsum_p(kd(ib,it),kd(ib,it))); dsum_p=zero
          !
          Do i=1,numberCons
             Do j=1,numberCons
                !
                ! temperature
                If((.Not.switch_on_temperature)) Then
                   !
                   Do l=1,kd(ib,it)
                      Do k=1,kd(ib,it)
                         If(Abs(EqpP(k)+EqpP(l)).Gt.Epsilo) Then
                            doublp(k,l)=fp12pl(k,l,i)/(EqpP(k)+EqpP(l))
                         Else
                            doublp(k,l)=zero
                         End If
                      End Do
                   End do
                   !
                   Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fp12pl(1,1,j),kd(ib,it),&
                                                         doublp,kd(ib,it),zero,dsum_p,kd(ib,it))
                Else
                   !
                   ! term corresponding to f^12
                   Do l=1,kd(ib,it)
                      Do k=1,kd(ib,it)
                         kk=iftP(k);ll=iftP(l)
                         temp_k=fp_T(kk)
                         temp_l=fp_T(ll)
                         If(Abs(EqpP(k)+EqpP(l)).Gt.Epsilo) Then
                            doublp(k,l)=fp12pl(k,l,i)*(one+temp_k+temp_l)/(EqpP(k)+EqpP(l))
                         Else
                            doublp(k,l)=zero
                         End If
                      End Do
                   End Do
                   !
                   Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fp12pl(1,1,j),kd(ib,it),&
                                                         doublp,kd(ib,it),zero,dsum_p,kd(ib,it))
                   !
                   ! first term: positive simplex
                   Do l=1,kd(ib,it)
                      Do k=1,kd(ib,it)
                         kk=iftP(k);ll=iftP(l)
                         temp_k=fp_T(kk)
                         temp_l=fp_T(ll)
                         If(k.ne.l.And.(Abs(EqpP(k)-EqpP(l)).Gt.Epsilo)) Then
                            t_term=-(temp_k-temp_l)/(EqpP(k)-EqpP(l))
                         Else
                            t_term=-temp_k*(temp_k-one)/temper
                         End If
                         doublp(k,l)=half*t_term*fp11pl(k,l,i)
                      End Do
                   End Do
                   !
                   Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fp11pl(1,1,j),kd(ib,it),&
                                                          doublp,kd(ib,it),one,dsum_p,kd(ib,it))
                   !
                End If
                !
                ! taking the trace of the resulting matrix
                !
                result=zero
                Do l=1,kd(ib,it)
                   result=result+dsum_p(l,l)
                End Do
                !
                cnsmat(i,j)=cnsmat(i,j)+result
                !
             End Do ! end of loop over j constraint
          End Do ! end of loop over i constraint
          !
          Deallocate(doublp,dsum_p,fp12pl,EqpP,ifTP)
          If(switch_on_temperature) Deallocate(fp11pl)
       End If
       !
    End Do ! end of loop over blocks ib
    !
    ! computing the inverse of the correlation matrix
    cnsorg=cnsmat
    !
    ierror=0
    Allocate(ipivot(numberCons))
    Call dgetrf(numberCons,numberCons,cnsmat,numberCons,ipivot,ierror)
    !
    ierror=0
    Allocate(workcn(numberCons))
    Call dgetri(numberCons,cnsmat,numberCons,ipivot,workcn,numberCons,ierror)
    Deallocate(ipivot)
    !
    ! constructing the vector of variations of the linear constraints
    trans='N'; incx=1; incy=1
    Call dgemv(trans,numberCons,numberCons,one,cnsmat,numberCons,cnsvec,incx,zero,workcn,incy)
    !
    ! updating the linear constraint vector (mixing has to be done simultaneously).
    If (ite.Eq.0) Then
        brakev=zero
    Else
        brakev=xmix
    End If
    !
    Allocate(veccns(numberCons))
    Do i=1,numberCons
       veccns(i)=veclam(i)+workcn(i)
       lambda=multLambda(i)
       ! regular multipole
       If(lambda.Ge.1) Then
          If(nbroyden.lt.1) Then
             multLag(lambda)=brakev*vecold(i)+(1.0-brakev)*veccns(i)
          Else
             multLag(lambda)=veccns(i)
             brout(nhhdim4+lambda)=multLag(lambda)
          End If
       End If
       If(lambda.Eq.0) Then
          If(nbroyden.lt.1) Then
             neckLag=brakev*vecold(i)+(1.0-brakev)*veccns(i)
          Else
             neckLag=veccns(i)
             brout(nhhdim4+lambdaMax+1)=neckLag
          End If
       End If
    End Do
    !
    Deallocate(veccns,vecold,workcn)
    Deallocate(cnsmat,cnsorg)
    Deallocate(qmultt,cnsvec,veclam)
    !
    Return
  End Subroutine getLagrange
  !=======================================================================
  !> Search for the requested state to block
  !=======================================================================
  Subroutine requested_blocked_level(ib,it)
    Implicit None
    Integer(ipr), Intent(in)  :: ib,it
    Integer(ipr) :: nd,im,k,ndk,na2,nad2,iqn,k0,LAPLUS,OMEGA,n1,n2,n3
    Real(pr) :: s1,s2,UUk,VVk
    k0=0
    If(nkblo(it,2).Eq.0)      Return
    If(Parity) Then
       LAPLUS=(ib+1)/2 !Yesp
    Else
       LAPLUS=ib       !Nop
    End If
    OMEGA=2*LAPLUS-1
    If(nkblo(it,1).Ne.OMEGA) Return
    nd=ID(ib); im=ia(ib);
    Do k=1,nd
       ndk=k+nd; s1=zero
       Do na2=1,nd
          nad2=na2+nd
          UUk=allhfb(ib)%arr(na2,ndk)
          VVk=allhfb(ib)%arr(nad2,ndk)
          s2=Max(s1,Abs(UUk),Abs(VVk))
          If(s2.Gt.s1) Then
             s1=s2
             iqn=na2+im  ! the position in [123] numbering
          End If
       End Do
       ! quantum numbers: Omega,P[n1,n2,n3]=>OMEGA,tpar(npar(iqn))[nz(iqn)+2*nr(iqn)+nl(iqn),nz(iqn),nl(iqn)]
       If(nkblo(it,2).Ne.tpar(npar(iqn))) Cycle
       n3=nl(iqn);          If(nkblo(it,5).Ne.n3) Cycle
       n2=nz(iqn);          If(nkblo(it,4).Ne.n2) Cycle
       n1=n2+2*nr(iqn)+n3;  If(nkblo(it,3).Ne.n1) Cycle
       k0=iqn
       keyblo(it)=1
       bloblo(keyblo(it),it)=ib
       blo123(keyblo(it),it)=k
       Exit
    End Do
  End Subroutine requested_blocked_level
  !=======================================================================
  !
  !=======================================================================
End Module HFBTHO_solver

