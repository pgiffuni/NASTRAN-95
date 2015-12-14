SUBROUTINE emg
     
!     ELEMENT-MATRIX-GENERATOR MAIN DRIVING ROUTINE.
 
!     DMAP SEQUENCE
 
!     EMG, EST,CSTM,MPT,DIT,GEOM2, /KMAT,KDICT, MMAT,MDICT, BMAT,BDICT/
!          V,N,NOKGG/V,N,NOMGG/V,N,NOBGG/V,N,NOK4GG/V,N,NOKDGG/
!          C,Y,COUPMASS/C,Y,CPBAR/C,Y,CPROD/C,Y,CPQUAD1/C,Y,CPQUAD2/
!          C,Y,CPTRIA1/C,Y,CPTRIA2/C,Y,CPTUBE/C,Y,CQDPLT/C,Y,CPTRPLT/
!          C,Y,CPTRBSC/V,Y,VOLUME/V,Y,SURFACE $
 
 LOGICAL :: error, anycon, nogo, heat, linear
 INTEGER :: z, est, cstm, dit, geom2, dictn
 INTEGER :: precis, cmass, flags, NAME(2)
 DIMENSION       ibuf(7),mcb(7)
 COMMON /BLANK / nok, nom, nob, nok4gg, nokdgg, cmass
 COMMON /emgprm/ icore, jcore, ncore, icstm, ncstm, imat, nmat,  &
     ihmat, nhmat, idit, ndit, icong, ncong, lcong,  &
     anycon, flags(3), precis, error, heat, icmbar,  &
     lcstm, lmat, lhmat, kflags(3), l38
 COMMON /zzzzzz/ z(1)
 COMMON /emgfil/ est, cstm, mpt, dit, geom2, mats(3), dictn(3)
 COMMON /hmatdd/ skp(4), linear
 COMMON /system/ ksystm(65)
 COMMON /machin/ mach
 EQUIVALENCE     (ksystm(3),nogo), (ksystm(55),ipreci),  &
     (ksystm(2),nout), (ksystm(56),noheat)
 DATA    NAME  / 4HEMG ,4H     /
 
!     SET EMG PRECISION FLAG TO SYSTEM PRECISION FLAG
 
 precis = ipreci
 
!     IF .NOT.1 .AND. .NOT.2 DEFAULT EMG PRECISION TO SINGLE
 
 IF (precis < 1 .OR. precis > 2) precis = 1
 
!     HEAT  FORMULATION
 
 heat   = .false.
 IF (noheat <= 0) GO TO 2
 heat   = .true.
 linear = .true.
 nokdgg = -1
 
!     TEST FOR NO SIMPLE ELEMENTS
 
 2 nogo   = .false.
 mcb(1) = 101
 CALL rdtrl (mcb)
 IF (mcb(1) < 0) GO TO 3
 IF (mcb(2) /= 0 .OR. mcb(5) /= 0 .OR. mcb(6) /= 0 .OR. mcb(7) /= 0) GO TO 5
 3 nok = -1
 nom = -1
 nob = -1
 nok4gg = -1
 RETURN
 
!     SET OPEN CORE
 
 5 ncore = korsz(z(1))
 icore = 3
 IF (mach == 3 .OR. mach == 4) CALL emgsoc (icore,ncore,heat)
 ncore = ncore - 1
 jcore = icore
 
!     SET WORKING CORE TO ALL ZEROS
 
 DO  i = icore,ncore
   z(i) = 0
 END DO
 
!     THIS MODULE WILL SET NOK4GG = -1 . IF DURING EXECUTION A NON-ZERO
!     DAMPING CONSTANT IS DETECTED IN A DICTIONARY BY EMGOUT, NOK4GG
!     WILL BE SET TO 1
 
!     A DMAP DETERMINATION CAN THEN BE MADE WHETHER OR NOT TO HAVE EMA
!     FORM THE K4GG MATRIX
 
 nok4gg = -1
 
!     SET GINO FILE NUMBERS
 
 est    = 101
 cstm   = 102
 mpt    = 103
 dit    = 104
 geom2  = 105
 DO  i = 1,3
   mats(i) = 199 + 2*i
   dictn(i) = mats(i) + 1
 END DO
 error = .false.
 
!     IF DIAG 38 IS ON, PRINT TOTAL TIME (IN SECONDS) USED BY EMGPRO
!     AND MESSAGES 3113 AND 3107  WHILE PRPCESSING ELEMENTS
 
 CALL sswtch (38,l38)
 
!     READ AND SETUP INTO CORE MISC. TABLES.
!     E.G. MPT, CSTM, DIT, ETC.
 
 CALL emgtab
 
!     PROCESS ANY CONGRUENT DATA CARDS AND BUILD TABLE IN OPEN CORE.
 
 CALL emgcng
 
!     SETUP BALANCE OF CORE WITH REQUIRED BUFFERS AND OPEN
!     REQUIRED DATA BLOCKS.
 
 CALL emgcor (ibuf)
 
!     PASS THE EST AND WRITE THE OUTPUT DATA BLOCKS.
 
 IF (l38 == 1) CALL klock (i)
 CALL emgpro (ibuf)
 IF (l38 == 0) GO TO 40
 CALL klock (j)
 j = j - i
 WRITE  (nout,30) j
 30 FORMAT (///,34H *** emg element processing time =,i10,8H seconds)
 
!     WRAP-UP OPERATIONS.
 
 40 CALL emgfin
 IF (nogo .OR. error) CALL mesage (-37,0,NAME)
 RETURN
END SUBROUTINE emg
