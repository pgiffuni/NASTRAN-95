SUBROUTINE sma1
!*****
! THIS ROUTINE IS A DRIVER AND INITIALIZATION PROGRAM FOR MODULE
! 2.4.1 OF THE NASTRAN SYSTEM.  IT GENERATES THE STIFFNESS MATRIX, KGG,
! THE STRUCTURAL DAMPING MATRIX, K4GG, AND THE GRID POINT SINGULARITY
! TABLE, GPST.
!*****
 DOUBLE PRECISION :: dz                 ,dpdum
 
 INTEGER :: iz(1)              ,eor  &
     ,                  clsrw              ,clsnrw ,                  frowic  &
     ,                  tnrows             ,outrw ,                  option
 
 LOGICAL :: anytab        ,linear
 LOGICAL :: dodet              ,heat
 
 DIMENSION nmsma1(2)
 DIMENSION ibuf(7)
 
 COMMON /BLANK/  nogenl             ,nok4gg   ,option(2)
 COMMON   /system/  isys,skip(53),iprec,itherm
 
! SMA1 I/O PARAMETERS
 
 COMMON   /sma1io/ ifcstm             ,ifmpt  &
     ,                  ifdit              ,idum1  &
     ,                  ifecpt             ,igecpt  &
     ,                  ifgpct             ,iggpct  &
     ,                  ifgei              ,iggei  &
     ,                  ifkgg              ,igkgg  &
     ,                  if4gg              ,ig4gg  &
     ,                  ifgpst             ,iggpst  &
     ,                  inrw               ,outrw  &
     ,                  clsnrw             ,clsrw  &
     ,                  neor               ,eor  &
     ,                  mcbkgg(7)          ,mcb4gg(7)
 
! SMA1 VARIABLE CORE
 
 COMMON   /zzzzzz /  z(1)
 
! SMA1 VARIABLE CORE BOOKKEEPING PARAMETERS
 
 COMMON   /sma1bk/ icstm              ,ncstm  &
     ,                  igpct              ,ngpct  &
     ,                  ipoint             ,npoint  &
     ,                  i6x6k              ,n6x6k  &
     ,                  i6x64              ,n6x64
 
! SMA1 PROGRAM CONTROL PARAMETERS
 
 COMMON   /sma1cl/ iopt4              ,k4ggsw  &
     ,                  npvt               ,left  &
     ,                  frowic             ,lrowic  &
     ,                  nrowsc             ,tnrows  &
     ,                  jmax               ,nlinks  &
     ,                  link(10)           ,idetck  &
     ,                  dodet              ,nogo
 
! ELEMENT DATA
 
 COMMON /gpta1/ nelems, last, incr, NE(1)
 
! ECPT COMMON BLOCK
 
 COMMON   /sma1et/ ecpt(100)
 
! SCRATCH COMMON BLOCK USED BY ELEMENT ROUTINES.
 
 COMMON   /sma1dp/ dpdum(300)
 
! COMMON INTERFACE FOR HMAT -HEAT- MATERIAL ROUTINE.
 
 COMMON /hmatdd/ ihmat,nhmat,mptmpt,idit,linear,anytab
 
 COMMON   /sma1ht/  heat
 
 EQUIVALENCE (z(1),iz(1),dz)
 
 DATA nmsma1(1) /4HSMA1/ ,nmsma1(2) /4H    /
!*****
!  SET THE LOGICAL HEAT FLAG IF THIS IS A -HEAT- FORMULATION
!*****
 CALL delset
 linear =.true.
 option(1) = -1
 heat = .false.
 IF( itherm /= 0 )   heat = .true.
 
 izmax = korsz(z)
 
! IF NOGENL .GT. 0, GENERAL ELEMENTS EXIST AND HENCE THE GPST IS NOT
! CREATED AND SO DETCK WILL NOT BE CALLED.
 
 dodet = .true.
 IF (nogenl > 0) dodet = .false.
 ibuf(1) = ifecpt
 CALL rdtrl(ibuf(1))
 IF (ibuf(3) == 1) dodet = .false.
 
! SET K4GG PURGE FLAGS
 
 nok4gg = -1
 k4ggsw = -1
 
! ATTEMPT TO OPEN THE OUTPUT FILE FOR THE KGG  MATRIX.  IF IT IS NOT
! IN THE OSCAR, EXECUTION WILL BE TERMINATED SINCE WE DO NOT ALLOW
! THE USER TO GENERATE ONLY A K4GG.
 
 igkgg = izmax - isys
 CALL OPEN(*100,ifkgg,z(igkgg),outrw)
 
! WRITE A TWO WORD BCD HEADER AND CLOSE THE KGG FILE WITHOUT REWIND.
 
 CALL fname (ifkgg,z(1))
 CALL WRITE (ifkgg,z(1),2,eor)
 CALL CLOSE (ifkgg,clsnrw)
 
! ATTEMPT TO OPEN THE K4GG FILE.
 
 ig4gg = igkgg
 iopt4 = 0
 CALL OPEN(*10,if4gg,z(ig4gg),outrw)
 iopt4 = 1
 ig4gg = ig4gg - isys
 CALL fname (if4gg,z(1))
 CALL WRITE (if4gg,z(1),2,eor)
 CALL CLOSE(if4gg,clsnrw)
 
! SET UP POINTERS TO GINO BUFFERS AND SET UP MATRIX CONTROL BLOCKS.
 
 10 igecpt = ig4gg - isys
 iggpct = igecpt - isys
 iggpst = iggpct - isys
 IF (.NOT. dodet) iggpst = iggpst + isys
 mcbkgg(1) = ifkgg
 mcbkgg(2) = 0
 mcbkgg(3) = 0
 mcbkgg(4) = 6
 mcbkgg(5) = iprec
 mcbkgg(6) = 0
 mcbkgg(7) = 0
 IF (iopt4 == 0) GO TO 30
 mcb4gg(1) = if4gg
 DO  i = 2,7
   mcb4gg(i) = mcbkgg(i)
 END DO
 
! ATTEMPT TO READ THE CSTM INTO CORE.
 
 30 ncstm = 0
 icstm = 0
 left = iggpst - 1
 CALL OPEN(*50,ifcstm,z(igkgg),inrw)
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
 
! CALL PREMAT TO READ MPT AND THE DIT INTO CORE
 
 imat11 = imat1 + 1
!*****
!  IF THIS IS A -HEAT- PROBLEM THE HMAT ROUTINE IS USED TO READ MAT4 AND
!  MAT5 CARDS INTO CORE.
!*****
 IF( .NOT. heat ) GO TO 56
 ihmat = imat11 + 1
 nhmat = imat11 + left - 2
 mptmpt = ifmpt
 idit = ifdit
 CALL hmat( 0 )
 left = left - nhmat + ihmat
 igpct = nhmat + 1
 GO TO 58
!*****
!  NORMAL PREMAT PROCESSING.
!*****
 56 CALL premat (iz(imat11),z(imat11),z(igkgg),left,matcr,ifmpt,ifdit)
 left = left - matcr
 igpct = ncstm + matcr
 
! OPEN THE ECPT AND GPCT INPUT FILES AND THE GPST OUTPUT FILE.
 
 58 CALL OPEN(*9070,ifecpt,z(igecpt),inrw)
 CALL fwdrec(*9080,ifecpt)
 CALL OPEN(*9090,ifgpct,z(iggpct),inrw)
 CALL fwdrec(*9100,ifgpct)
 IF (.NOT. dodet) GO TO 60
 CALL OPEN(*9110,ifgpst,z(iggpst),outrw)
 CALL fname(ifgpst,ecpt(1))
 CALL WRITE(ifgpst,ecpt(1),2,eor)
 
! REOPEN THE KGG OUTPUT FILE WITHOUT REWIND, AND THE K4GG, IF CALLED FOR
 
 60 CALL OPEN(*9120,ifkgg,z(igkgg),3)
 IF(iopt4 /= 0)CALL OPEN(*9130,if4gg,z(ig4gg),3)
 
! CALL SUBROUTINE SMA1A WHICH WILL PERFORM ALL THE COMPUTATIONS.
 
 CALL sma1a
 IF(.NOT. linear) option(1)= 1
 
! CLOSE FILES AND WRITE TRAILERS.
 
 CALL CLOSE(ifecpt,clsrw)
 CALL CLOSE(ifgpct,clsrw)
 IF (.NOT. dodet) GO TO 70
 CALL CLOSE (ifgpst,clsrw)
 CALL wrttrl (ifgpst)
 70 CALL CLOSE (ifkgg,clsrw)
 mcbkgg(3) = mcbkgg(2)
 CALL wrttrl (mcbkgg(1))
 IF (iopt4 == 0) GO TO 100
 CALL CLOSE(if4gg,clsrw)
 IF (mcb4gg(6) == 0) GO TO 80
 mcb4gg(3) = mcb4gg(2)
 CALL wrttrl (mcb4gg(1))
 nok4gg = 1
 GO TO 100
 80 DO  i = 2,7
   mcb4gg(i) = 0
 END DO
 nok4gg = -1
 100 RETURN
 
! SUBROUTINE SMA1 ERROR EXITS.
 
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
 9110 ifile = ifgpst
 GO TO 10001
 9120 ifile = ifkgg
 GO TO 10001
 9130 ifile = if4gg
 10001 iparm = -1
 GO TO 10010
 10002 iparm = -2
 10010 CALL mesage (iparm,ifile,nmsma1(1))
 RETURN
END SUBROUTINE sma1
