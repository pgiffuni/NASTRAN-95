SUBROUTINE sma2
! ******
! THIS ROUTINE IS A DRIVER AND INITIALIZATION PROGRAM FOR MODULE
! 2.4.2 OF THE NASTRAN SYSTEM.  IT GENERATES THE MASS MATRIX, MGG, AND
! THE DAMPING MATRIX, BGG.
! ******
 
 
 
 
 LOGICAL :: heat
 
 
 DOUBLE PRECISION :: dz                 ,zzzzzz
 
 
 
 INTEGER :: iz(1)              ,eor  &
     ,                  clsrw              ,clsnrw ,                  frowic  &
     ,                  tnrows             ,outrw ,                  bggind
 
 
 
 DIMENSION nmsma2(2)
 COMMON /BLANK/      wtmass             ,nomgg ,                  nobgg
 
 
 
 COMMON   /system/  isys,isew1(53),iprec,itherm
 
! SMA2 I/O PARAMETERS
 
 COMMON   /sma2io/ ifcstm             ,ifmpt  &
     ,                  ifdit              ,idum1  &
     ,                  ifecpt             ,igecpt  &
     ,                  ifgpct             ,iggpct  &
     ,                  idum2              ,idum3  &
     ,                  ifmgg              ,igmgg  &
     ,                  ifbgg              ,igbgg  &
     ,                  idum4              ,idum5  &
     ,                  inrw               ,outrw  &
     ,                  clsnrw             ,clsrw  &
     ,                  neor               ,eor  &
     ,                  mcbmgg(7)          ,mcbbgg(7)
 
! SMA2 VARIABLE CORE
 
 COMMON   /zzzzzz /  z(1)
 
! SMA2 VARIABLE CORE BOOKKEEPING PARAMETERS.
 
 COMMON   /sma2bk/ icstm              ,ncstm  &
     ,                  igpct              ,ngpct  &
     ,                  ipoint             ,npoint  &
     ,                  i6x6m              ,n6x6m  &
     ,                  i6x6b              ,n6x6b
 
! SMA2 PROGRAM CONTROL PARAMETERS
 
 COMMON   /sma2cl/ ioptb              ,bggind  &
     ,                  npvt               ,left  &
     ,                  frowic             ,lrowic  &
     ,                  nrowsc             ,tnrows  &
     ,                  jmax               ,nlinks  &
     ,                  link(10)           ,nogo
 
! ELEMENT DATA
 
 COMMON   /gpta1/ nelems, last, incr, NE(1)
 
! ECPT COMMON BLOCK
 
 COMMON   /sma2et/ ecpt(100)
 
! SCRATCH BLOCK FOR ELEMENT ROUTINES
 
 COMMON   /sma2dp/ zzzzzz(300)
 
 COMMON   /sma2ht/  heat
 
 COMMON   /hmatdd/  ihmat, nhmat, mptmpt, idit
 
 
 EQUIVALENCE (z(1),iz(1),dz)
 
 
 
 DATA nmsma2(1) /4HSMA2/ ,nmsma2(2) /4H    /
 
!*****
!  SET HEAT FLAG
!*****
 heat = .false.
 IF( itherm /= 0 )   heat = .true.
 
 
 CALL delset
 izmax = korsz (z)
 
! SET PURGE FLAGS FOR BGG AND NO PURGE FLAG FOR MGG.
 
 bggind = -1
 nobgg  = -1
 nomgg = 1
 
! ATTEMPT TO OPEN THE OUTPUT FILE FOR THE MASS MATRIX.  IF IT IS NOT
! IN THE OSCAR, EXECUTION WILL BE TERMINATED SINCE WE DO NOT ALLOW
! THE USER TO GENERATE ONLY A BGG. (EXCEPT IN A HEAT TRANSER PROBLEM)
 
 igmgg = izmax - isys
 IF( heat ) GO TO 5
 CALL OPEN(*100,ifmgg,z(igmgg),outrw)
 
! WRITE A TWO WORD BCD HEADER AND CLOSE THE MGG FILE WITHOUT REWIND.
 
 CALL fname (ifmgg,z(1))
 CALL WRITE (ifmgg,z(1),2,eor)
 CALL CLOSE (ifmgg,clsnrw)
 
! ATTEMPT TO OPEN THE BGG FILE.
 
 5 igbgg = igmgg
 ioptb = 0
 CALL OPEN(*10,ifbgg,z(igbgg),outrw)
 ioptb = 1
 igbgg = igbgg - isys
 CALL fname (ifbgg,z(1))
 CALL WRITE (ifbgg,z(1),2,eor)
 CALL CLOSE(ifbgg,clsnrw)
 
! SET UP POINTERS TO GINO BUFFERS AND SET UP MATRIX CONTROL BLOCKS.
 
 10 igecpt = igbgg  - isys
 iggpct = igecpt - isys
 mcbmgg(1) = ifmgg
 mcbmgg(2) = 0
 mcbmgg(3) = 0
 mcbmgg(4) = 6
 mcbmgg(5) = iprec
 mcbmgg(6) = 0
 mcbmgg(7) = 0
 IF (ioptb == 0) GO TO 30
 mcbbgg(1) = ifbgg
 DO  i = 2,7
   mcbbgg(i) = mcbmgg(i)
 END DO
 
! ATTEMPT TO READ THE CSTM INTO CORE.
 
 30 ncstm = 0
 icstm = 0
 left = iggpct - 1
 CALL OPEN(*50,ifcstm,z(igmgg),inrw)
 CALL fwdrec(*9020,ifcstm)
 CALL READ(*9030,*40,ifcstm,z(1),left,eor,ncstm)
 
! IF CORE WAS FILLED WITHOUT HITTING AN EOR CALL MESAGE
 
 CALL mesage (-8,ifcstm,ifcstm)
 40 left = left - ncstm
 
! PRETRD SETS UP FUTURE CALLS TO TRANSD.
 
 CALL pretrd (z(icstm+1),ncstm)
 CALL pretrs(z(icstm+1),ncstm)
 CALL CLOSE (ifcstm,clsrw)
 50 imat1 = ncstm
 nmat1 = 0
 nmat2 = 0
 nmat3 = 0
 nmat4 = 0
 imat11 = imat1 + 1
!*****
!  IF -HEAT- PROBLEM THEN HMAT IS USED FOR MAT4 AND MAT5 CARDS.
!*****
 IF( .NOT. heat ) GO TO 56
 ihmat = imat11 + 1
 nhmat = imat11 + left - 2
 mptmpt = ifmpt
 idit = ifdit
 CALL hmat( 0 )
 left = left - nhmat + ihmat
 igpct = nhmat +1
 GO TO 58
 56 CALL premat(iz(imat11),z(imat11),z(igmgg),left,matcr,ifmpt,ifdit)
 left = left - matcr
 igpct = ncstm + matcr
 
! OPEN THE ECPT INPUT FILE AND THE GPCT INPUT FILE.
 
 58 CALL OPEN(*9070,ifecpt,z(igecpt),inrw)
 CALL fwdrec(*9080,ifecpt)
 CALL OPEN(*9090,ifgpct,z(iggpct),inrw)
 CALL fwdrec(*9100,ifgpct)
 
! REOPEN THE MGG OUTPUT FILE WITHOUT REWIND, AND THE BGG, IF CALLED FOR.
 
 IF(.NOT.heat)CALL OPEN(*9110,ifmgg,z(igmgg),3)
 IF(ioptb /= 0)CALL OPEN(*9120,ifbgg,z(igbgg),3)
 
! CALL SUBROUTINE SMA2A WHICH WILL PERFORM ALL THE COMPUTATIONS.
 
 CALL sma2a
 
! CLOSE FILES AND WRITE TRAILERS.
 
 CALL CLOSE (ifgpct,clsrw)
 CALL CLOSE (ifecpt,clsrw)
 CALL CLOSE (ifmgg ,clsrw)
 mcbmgg(3) = mcbmgg(2)
 IF (mcbmgg(6) /= 0) GO TO 70
 DO  i = 2,7
   mcbmgg(i) = 0
 END DO
 nomgg = -1
 70 IF( .NOT. heat ) CALL wrttrl( mcbmgg )
 IF (ioptb == 0) GO TO 100
 CALL CLOSE (ifbgg ,clsrw)
 IF (mcbbgg(6) == 0) GO TO 80
 mcbbgg(3) = mcbbgg(2)
 CALL wrttrl (mcbbgg(1))
 nobgg = 1
 GO TO 100
 80 DO  i = 2,7
   mcbbgg(i) = 0
 END DO
 nobgg = -1
 100 RETURN
 
! SUBROUTINE SMA2 ERROR EXITS.
 
 9020 ifile = ifcstm
 GO TO 10002
 9030 ifile = - ifcstm
 GO TO 10002
 9070 ifile = ifecpt
 GO TO 10001
 9080 ifile = ifecpt
 GO TO 10002
 9090 ifile = ifgpct
 GO TO 10001
 9100 ifile = ifgpct
 GO TO 10002
 9110 ifile = ifmgg
 GO TO 10002
 9120 ifile = ifbgg
 GO TO 10002
 10001 iparm = -1
 GO TO 10010
 10002 iparm = -2
 10010 CALL mesage (iparm,ifile,nmsma2(1))
 RETURN
END SUBROUTINE sma2
