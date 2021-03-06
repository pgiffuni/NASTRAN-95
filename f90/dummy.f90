SUBROUTINE dummy
     
!     NOTE:
!     THIS DUMMY.MIS ROUTINE CONTAINS 4 MACHINE VERSIONS (IBM,CDC,VAX,
!     AND UNIVAC). MOVE THIS SUBROUTINE TO THE MDS GROUP AND
!     REPLACE ALL THE 'C+' BY 2 SPACES IF MACHINE IS IBM, OR
!     REPLACE ALL THE 'C-' BY 2 SPACES IF MACHINE IS CDC, OR
!     REPLACE ALL THE 'C=' BY 2 SPACES IF MACHINE IS VAX, AND UNIX, OR
!     REPLACE ALL THE 'C*' BY 2 SPACES IF MACHINE IS UNIVAC
!     REPLACE ALL THE 'C.' BY 2 SPACES IF MACHINE TYPE IS 1, AND 11-20
 
! ****
!     IBM VERSION
 
!     THIS SUBROUTINE PROVIDES ENTRIES FOR THE DUMMY ROUTINES
!     USED BY OTHER COMPUTER MACHINES, AND ARE REFERENCED IN
!     VARIOUS NASTRAN LINKS
 
!     THIS SUBROUTINE INCLUDES ALSO SOME DUMMY ROUTINES NOT YET
!     WRITTEN
 
!     THIS ROUTINE SHOULD BE MOVED TO NASTRAN MACHINE-DEPENDENT
!     SECTION (MDS)
! ****
 
!+    DIMENSION       N(1)
!+    CHARACTER*8     NAME
 
!+    COMMON /MACHIN/ MACH
!+    COMMON /SYSTEM/ ISYSBF, NOUT
 
!+    IF (MACH .EQ. 2) GO TO 250
!+    WRITE  (NOUT,20) MACH
!+ 20 FORMAT (/,' MACH =',I7)
!+    NAME = 'DUMMY'
!+    GO TO 100
 
! ****
!     ROUTINES USED ONLY IN UNIVAC MACHINE
! ****
 
 
!+    ENTRY NTRAN (I,J,K)
!+    NAME = 'NTRAN'
!+    GO TO 100
 
!+    ENTRY CONTIN
!+    NAME = 'CONTIN'
!+    GO TO 100
 
!+    ENTRY FACIL (I,J)
!+    NAME = 'FACIL'
!+    GO TO 100
 
!+    ENTRY FACSF (I)
!+    NAME = 'FACSF'
!+    GO TO 100
 
!+    ENTRY UNVOPN (I)
!+    NAME = 'UNVOPN'
!+    GO TO 100
 
!+    ENTRY UNVCLS (I)
!+    NAME = 'UNVCLS'
!+    GO TO 100
 
!+    ENTRY ADDCRD (I,J)
!+    NAME = 'ADDCRD'
!+    GO TO 100
 
! ****
!     ROUTINES USED BY UNIVAC AND IBM
! ****
 
!     ENTRY RETURN
!     GO TO 250
 
!+    ENTRY MSGUNI
!+    IF (MACH .EQ. 2) GO TO 250
!+    NAME = 'MSGUNI'
!+    GO TO 100
 
!+    ENTRY XEOT (I,J,K,L)
!+    IF (MACH .EQ. 2) GO TO 250
!+    NAME = 'XEOT'
!+    GO TO 100
 
!     ENTRY TPSWIT (I,J,K,L)
!     NAME = 'TPSWIT'
!     GO TO 100
 
! ****
!     ROUTINES USED ONLY IN IBM MACHINE
! ****
 
!     ENTRY UMFTRN (I)
!     NAME = 'UMFTRN'
!     GO TO 100
 
!     ENTRY TAPSWI (I,J,K,L)
!     NAME = 'TAPSWI'
!     GO TO 100
 
!     ENTRY SOFIOI
!     NAME = 'SOFIOI'
!     GO TO 100
 
!     ENTRY SEARCH (I)
!     NAME = 'SEARCH'
!     GO TO 100
 
! ... NEXT THREE ARE SYSTEM ROUTINES THAT OPEN FILE DYNAMICALLY WITHOUT
!     THE USE OF JCL. THESE ROUTINES ARE COMMONLY 'LOCAL INSTALLED'.
 
!     IQADDN CHECKS WHETHER A FILE EXISTS OR NOT
!     QQDCBF BUILDS AN ATTRIBUTE LIST BY DDNAME
!     QQGETF ALLOCATES FILE IN TSO OR BATCH
 
!     ENTRY IQZDDN (I)
!     NAME = 'IQZDDN'
!     GO TO 100
 
!     ENTRY QQDCBF (I,J,K,L,M,N)
!     NAME = 'QQDCBF'
!     GO TO 100
 
!     ENTRY QQGETF (I,J,K)
!     NAME = 'QQGETF'
!     GO TO 100
 
! ****
!     ROUTINE USED ONLY BY IBM AND VAX
! ****
 
!     ENTRY SOFIOF
!     NAME = 'SOFIOF'
!     GO TO 100
 
!     THE FOLLOWING THREE ARE FUNCTIONS FOR QUAD WORD OPERATIONS
!                                           (REAL*16)
!     ENTRY QABS (I)
!     NAME = 'QABS'
!     GO TO 100
 
!     ENTRY SNGLQ (I)
!     NAME = 'SNGLQ'
!     GO TO 100
 
!     ENTRY DBLEQ (I)
!     NAME = 'DBLEQ'
!     GO TO 100
 
!     ENTRY QSQRT (I)
!     NAME = 'QSQRT'
!     GO TO 100
 
!     ENTRY QLOG (I)
!     NAME = 'QLOG'
!     GO TO 100
 
!     ENTRY QEXTD (I)
!     NAME = 'QEXTD'
!     GO TO 100
 
! ****
!     ROUTINE USED BY UNIVAC AND VAX
! ****
 
!+    ENTRY DEFCOR
!+    NAME = 'DEFCOR'
!+    GO TO 100
 
! ****
!     ROUTINES USED BY ALL MACHINES, EXCEPT VAX
! ****
 
!     ENTRY GPERR
!     NAME = 'GPERR'
!     GO TO 100
 
!     ENTRY PDUMP
!     GO TO 250
 
!     ENTRY MPY1
!     NAME = 'MPY1'
!     GO TO 100
 
!     ENTRY MPY2NT
!     NAME = 'MPY2NT'
!     GO TO 100
 
!     ENTRY MPY2T
!     NAME = 'MPY2T'
!     GO TO 100
 
! ****
!     ROUTINES USED ONLY IN CDC MACHINE
! ****
 
!+    ENTRY LINK (I,J,K)
!+    NAME = 'LINK'
!+    GO TO 100
 
!+    ENTRY REMARK (I)
!+    NAME = 'REMARK'
!+    GO TO 100
 
!+    ENTRY CDCBUG (I,J,K,L)
!+    NAME = 'CDCBUG'
!+    GO TO 100
 
!+    ENTRY CDCOPN (I)
!+    NAME = 'CDCOPN'
!+    GO TO 100
 
!+    ENTRY CDCCLS (I)
!+    NAME = 'CDCCLS'
!+    GO TO 100
 
!+    ENTRY CDCKSZ (I)
!+    NAME = 'CDCKSZ'
!+    GO TO 100
 
!+    ENTRY PF (I,J,K)
!+    NAME = 'PF'
!+    GO TO 100
 
!+    ENTRY ISWAP (I)
!+    NAME = 'ISWAP'
!+    GO TO 100
 
! ****
!     ROUTINES USED ONLY IN VAX MACHINE
! ****
 
!+    ENTRY VAXEND
!+    NAME = 'VAXEND'
!+    GO TO 100
 
!+    ENTRY VAXERR (L)
!+    WRITE (NOUT,50) L
!+ 50 FORMAT (/,' *** GINO ERROR AT LOC',I5)
!+    GO TO 220
 
!+    ENTRY VAXSCH
!+    NAME = 'VAXSCH'
!+    GO TO 100
 
!+    ENTRY VAXBRK
!+    NAME = 'VAXBRK'
!+    GO TO 100
 
!+    ENTRY MPY1V (I,J,K)
!+    NAME = 'MPY1V'
!+    GO TO 100
 
!+    ENTRY MPY2NV (I,J,K)
!+    NAME = 'MPY2NV'
!+    GO TO 100
 
!+    ENTRY MPY2TV (I,J,K)
!+    NAME = 'MPY2TV'
!+    GO TO 100
 
! ****
!     ROUTINES THAT PERFORM NO PARTICULAR FUNCTIONS, BUT THEY
!     ARE STILL CALLED BY NASTRAN
! ****
 
!+    ENTRY UNLOAD (I)
!     CALLED BY INPTT1
!+    GO TO 250
 
! ****
!     THE FOLLOWING ROUTINES SEEM TO BE NO LONGER USED IN NASTRAN
! ****
 
!+    ENTRY JIDINT (I)
!+    NAME = 'JIDINT'
!+    GO TO 100
 
!+    ENTRY OPMESG
!+    NAME = 'OPMESG'
!+    GO TO 100
 
!     ENTRY PDUM1,PDUM2,...,PDUM9 HAD BEEN REPLACED BY PDUMI
!     ENTRY QDMM3, SQDM31, AND SQDM32 ARE NOW OBSOLETE
 
!+    ENTRY SEMTRN
!+    NAME = 'SEMTRN'
!+    GO TO 100
 
! ****
!     DUMMY ROUTINES REFERENCED ONLY IN LINK 2, ALL MACHINES
! ****
 
!+    ENTRY PDUMI (*,*,*,I,J,K,L,M,N,O)
!+    NAME = 'PDUMI'
!+    GO TO 100
 
! ****
!     DUMMY ROUTINES REFERENCED ONLY IN LINK 5, ALL MACHINES
! ****
 
!+    ENTRY PLBAR1 (I,J)
!+    NAME = 'PLBAR1'
!+    GO TO 100
 
!+    ENTRY PLOADX
!+    NAME = 'PLOADX'
!+    GO TO 100
 
!+    ENTRY ERRTRC (NAM)
!     ==================
!     ERROR TRACEBACK
 
!+    GO TO 220
 
!+100 WRITE  (NOUT,150) NAME
!+150 FORMAT ('0*** SYSTEM FATAL ERROR  ---  JOB TERMINATED',
!+   1        ' DUE TO CALL TO DUMMY SUBROUTINE.  ENTRY NAME IS ', A8)
!+    GO TO 220
 
! ****
!     TO FORCE A SYSTEM FATAL ERROR FOR TRACEBACK
! ****
 
!+220 WRITE  (NOUT,230)
!+230 FORMAT ('0*** ERROR TRACEBACK IN SYSTEM LOG FILE')
!+    I = 987654321
!+    N(I) = 1
!+250 RETURN
 
 
!     SUBROUTINE DUMMY
 
! ****
!     CDC VERSION
 
!     THIS SUBROUTINE PROVIDES ENTRIES FOR THE DUMMY ROUTINES
!     USED BY OTHER COMPUTER MACHINES, AND ARE REFERENCED IN
!     VARIOUS NASTRAN LINKS
 
!     THIS SUBROUTINE INCLUDES ALSO SOME DUMMY ROUTINES NOT YET
!     WRITTEN
 
!     THIS ROUTINE SHOULD BE MOVED TO NASTRAN MACHINE-DEPENDENT
!     SECTION (MDS)
! ****
 
!-    CHARACTER*8     NAME
 
!-    COMMON /MACHIN/ MACH
!-    COMMON /SYSTEM/ ISYSBF, NOUT
 
!-    IF (MACH .EQ. 4) GO TO 250
!-    WRITE  (NOUT,20) MACH
!- 20 FORMAT (/,' MACH =',I7)
!-    NAME = 'DUMMY'
!-    GO TO 100
 
! ****
!     ROUTINES USED ONLY IN UNIVAC MACHINE
! ****
 
!-    ENTRY NTRAN (I,J,K)
!-    NAME = 'NTRAN'
!-    GO TO 100
 
!-    ENTRY CONTIN
!-    NAME = 'CONTIN'
!-    GO TO 100
 
!-    ENTRY FACIL (I,J)
!-    NAME = 'FACIL'
!-    GO TO 100
 
!-    ENTRY FACSF (I)
!-    NAME = 'FACSF'
!-    GO TO 100
 
!-    ENTRY UNVOPN (I)
!-    NAME = 'UNVOPN'
!-    GO TO 100
 
!-    ENTRY UNVCLS (I)
!-    NAME = 'UNVCLS'
!-    GO TO 100
 
!-    ENTRY ADDCRD (I,J)
!-    NAME = 'ADDCRD'
!-    GO TO 100
 
! ****
!     ROUTINES USED BY UNIVAC AND IBM
! ****
 
!-    ENTRY RETURN
!-    GO TO 250
 
!-    ENTRY MSGUNI
!-    IF (MACH .EQ. 2) GO TO 250
!-    NAME = 'MSGUNI'
!-    GO TO 100
 
!-    ENTRY XEOT (I,J,K,L)
!-    IF (MACH .EQ. 2) GO TO 250
!-    NAME = 'XEOT'
!-    GO TO 100
 
!-    ENTRY TPSWIT (I,J,K,L)
!-    NAME = 'TPSWIT'
!-    GO TO 100
 
! ****
!     ROUTINES USED ONLY IN IBM MACHINE
! ****
 
!-    ENTRY UMFTRN (I)
!-    NAME = 'UMFTRN'
!-    GO TO 100
 
!-    ENTRY TAPWSI (I,J,K,L)
!-    NAME = 'TAPSWI'
!-    GO TO 100
 
!-    ENTRY SEARCH (I)
!-    NAME = 'SEARCH'
!-    GO TO 100
 
!-    ENTRY SOFIOI
!-    NAME = 'SOFIOI'
!-    GO TO 100
 
!-    ENTRY IQZDDN (I)
!-    NAME = 'IQZDDN'
!-    GO TO 100
 
!-    ENTRY QQDCBF (I,J,K,L,M,N)
!-    NAME = 'QQDCBF'
!-    GO TO 100
 
!-    ENTRY QQGETF (I,J,K)
!-    NAME = 'QQGETF'
!-    GO TO 100
 
! ****
!     ROUTINE USED ONLY BY IBM AND VAX
! ****
 
!-    ENTRY SOFIOF
!-    NAME = 'SOFIOF'
!-    GO TO 100
 
!     THE FOLLOWING THREE ARE FUNCTIONS FOR QUAD WORD OPERATIONS
!                                           (REAL*16)
!-    ENTRY QABS (I)
!-    NAME = 'QABS'
!-    GO TO 100
 
!-    ENTRY SNGLQ (I)
!-    NAME = 'SNGLQ'
!-    GO TO 100
 
!-    ENTRY DBLEQ (I)
!-    NAME = 'DBLEQ'
!-    GO TO 100
 
!-    ENTRY QSQRT (I)
!-    NAME = 'QSQRT'
!-    GO TO 100
 
!-    ENTRY QLOG (I)
!-    NAME = 'QLOG'
!-    GO TO 100
 
!-    ENTRY QEXTD (I)
!-    NAME = 'QEXTD'
!-    GO TO 100
 
! ****
!     ROUTINE USED BY UNIVAC AND VAX
! ****
 
!-    ENTRY DEFCOR
!-    NAME = 'DEFCOR'
!-    GO TO 100
 
! ****
!     ROUTINES USEDS BY ALL MACHINES, EXCEPT VAX
! ****
 
!     ENTRY GPERR (I,J)
!     NAME = 'GPERR'
!     GO TO 100
 
!     ENTRY PDUMP
!     NAME = 'PDUMP'
!     GO TO 250
 
!     ENTRY MPY1
!     NAME = 'MPY1'
!     GO TO 100
 
!     ENTRY MPY2NT
!     NAME = 'MPY2NT'
!     GO TO 100
 
!     ENTRY MPY2T
!     NAME = 'MPY2T'
!     GO TO 100
 
! ****
!     ROUTINES USED ONLY IN CDC MACHINE
! ****
 
!     ENTRY LINK (I,J,K)
!     NAME = 'LINK'
!     GO TO 100
 
!     ENTRY REMARK (I)
!     NAME = 'REMARK'
!     GO TO 100
 
!     ENTRY CDCBUG (I,J,K,L)
!     NAME = 'CDCBUG'
!     GO TO 100
 
!     ENTRY CDCOPN (I)
!     NAME = 'CDCOPN'
!     GO TO 100
 
!     ENTRY CDCCLS (I)
!     NAME = 'CDCCLS'
!     GO TO 100
 
!     ENTRY PF (I,J,K)
!     NAME = 'PF'
!     GO TO 100
 
!     ENTRY ISWAP (I)
!     NAME = 'ISWAP'
!     GO TO 100
 
!-    ENTRY CDCKSZ (I)
!-    ENCODE (20,30,A) I
!- 30 FORMAT ('OPEN CORE =',I7,2X)
!-    CALL REMARK (A)
!-    GO TO 250
 
! ****
!     ROUTINES USED ONLY IN VAX MACHINE
! ****
 
!-    ENTRY VAXEND
!-    NAME = 'VAXEND'
!-    GO TO 100
 
!-    ENTRY VAXERR (L)
!-    WRITE  (NOUT,50) L
!- 50 FORMAT (/,' *** GINO ERROR AT LOC',I5)
!-    GO TO 220
 
!-    ENTRY VAXSCH
!-    NAME = 'VAXSCH'
!-    GO TO 100
 
!-    ENTRY VAXBRK
!-    NAME = 'VAXBRK'
!-    GO TO 100
 
!-    ENTRY MPY1V (I,J,K)
!-    NAME = 'MPY1V'
!-    GO TO 100
 
!-    ENTRY MPY2NV (I,J,K)
!-    NAME = 'MPY2NV'
!-    GO TO 100
 
!-    ENTRY MPY2TV (I,J,K)
!-    NAME = 'MPY2TV'
!-    GO TO 100
 
! ****
!     ROUTINES THAT PERFORM NO PARTICULAR FUNCTIONS, BUT THEY
!     ARE STILL CALLED BY NASTRAN
! ****
 
!-    ENTRY UNLOAD (I)
!     CALLED BY INPTT1
!-    GO TO 250
 
! ****
!     THE FOLLOWING ROUTINES SEEM TO BE NO LONGER USED IN NASTRAN
! ****
 
!-    ENTRY JIDINT (I)
!-    NAME = 'JIDINT'
!-    GO TO 100
 
!-    ENTRY OPMESG
!-    NAME = 'OPMESG'
!-    GO TO 100
 
!     ENTRY PDUM1,PDUM2,...,PDUM9 HAD BEEN REPLACED BY PDUMI
!     ENTRY QDMM3, SQDM31, AND SQDM32 ARE NOW OBSOLETE
 
!-    ENTRY SEMTRN
!-    NAME = 'SEMTRN'
!-    GO TO 100
 
! ****
!     DUMMY ROUTINES REFERENCED ONLY IN LINK 2, ALL MACHINES
! ****
 
!-    ENTRY PDUMI (*,*,*,I,J,K,L,M,N,O)
!-    NAME = 'PDUMI'
!-    GO TO 100
 
! ****
!     DUMMY ROUTINES REFERENCED ONLY IN LINK 5, ALL MACHINES
! ****
 
!-    ENTRY PLBAR1 (I,J)
!-    NAME = 'PLBAR1'
!-    GO TO 100
 
!-    ENTRY PLOADX
!-    NAME = 'PLOADX'
!-    GO TO 100
 
!-    ENTRY ERRTRC (NAM)
!     ==================
!     ERROR TRACEBACK
 
!-    GO TO 220
 
!-100 WRITE  (NOUT,150) NAME
!-150 FORMAT ('0*** SYSTEM FATAL ERROR  ---  JOB TERMINATED',
!-   1        ' DUE TO CALL TO DUMMY SUBROUTINE.  ENTRY NAME IS ', A8)
!-    GO TO 220
 
! ****
!     TO FORCE A SYSTEM FATAL ERROR FOR TRACEBACK
! ****
 
!-220 WRITE  (NOUT,230)
!-230 FORMAT ('0*** ERROR TRACEBACK IN SYSTEM LOG FILE')
!-    I =-3
!-    READ (I) J,K,M,N,O
!-250 RETURN
 
 
!     SUBROUTINE DUMMY
 
! ****
!     VAX VERSION  (MODIFIED FOR DEC/ULTRIX)
 
!     THIS SUBROUTINE PROVIDES ENTRIES FOR THE DUMMY ROUTINES
!     USED BY OTHER COMPUTER MACHINES, AND ARE REFERENCED IN
!     VARIOUS NASTRAN LINKS
 
!     THIS SUBROUTINE INCLUDES ALSO SOME DUMMY ROUTINES NOT YET
!     WRITTEN
 
!     THIS ROUTINE SHOULD BE MOVED TO NASTRAN MACHINE-DEPENDENT
!     SECTION (MDS)
! ****
 
 DIMENSION       n(1)
 CHARACTER (LEN=8) :: NAME
 
 COMMON /machin/ mach
 COMMON /system/ isysbf, nout
 
 IF (mach == 6) GO TO 250
 WRITE  (nout,20) mach
 20 FORMAT (/,' MACH =',i7)
 NAME = 'DUMMY'
 GO TO 100
 
! ****
!     ROUTINES USED ONLY IN UNIVAC MACHINE
! ****
 
 ENTRY zcorsz (i)
 NAME = 'ZCORSZ'
 GO TO 100
 
 ENTRY mvbits (i1,i2,i3,i4,i5)
 NAME = 'MVBITS'
 GO TO 100
 
 ENTRY codkey (code,key)
 NAME = 'CODKEY'
 GO TO 100
 
 ENTRY kconeq
 NAME = 'KCONEQ'
 GO TO 100
 
!hgs      ENTRY FNXTVQ (V1,V2,V3,V4,V5,ZB,I)
!hgs      NAME = 'FNXTVQ'
!hgs      GO TO 100
 
 ENTRY ntran (i,j,k)
 NAME = 'NTRAN'
 GO TO 100
 
 ENTRY contin
 NAME = 'CONTIN'
 GO TO 100
 
 ENTRY facil (i,j)
 NAME = 'FACIL'
 GO TO 100
 
 ENTRY facsf (i)
 NAME = 'FACSF'
 GO TO 100
 
 ENTRY unvopn (i)
 NAME = 'UNVOPN'
 GO TO 100
 
 ENTRY unvcls (i)
 NAME = 'UNVCLS'
 GO TO 100
 
 ENTRY addcrd (i,j)
 NAME = 'ADDCRD'
 GO TO 100
 
! ****
!     ROUTINES USED BY UNIVAC AND IBM
! ****
 
 ENTRY RETURN
 GO TO 250
 
 ENTRY msguni
 IF (mach == 2) GO TO 250
 NAME = 'MSGUNI'
 GO TO 100
 
 ENTRY xeot (i,j,k,l)
 IF (mach == 2) GO TO 250
 NAME = 'XEOT'
 GO TO 100
 
 ENTRY tpswit (i,j,k,l)
 NAME = 'TPSWIT'
 GO TO 100
 
! ****
!     ROUTINES USED ONLY IN IBM MACHINE
! ****
 
 ENTRY umftrn (i)
 NAME = 'UMFTRN'
 GO TO 100
 
 ENTRY tapswi (i,j,k,l)
 NAME = 'TAPSWI'
 GO TO 100
 
 ENTRY search (i)
 NAME = 'SEARCH'
 GO TO 100
 
 ENTRY sofioi
 NAME = 'SOFIOI'
 GO TO 100
 
 ENTRY iqzddn (i)
 NAME = 'IQZDDN'
 GO TO 100
 
 ENTRY qqdcbf (i,j,k,l,m,n)
 NAME = 'QQDCBF'
 GO TO 100
 
 ENTRY qqgetf (i,j,k)
 NAME = 'QQGETF'
 GO TO 100
 
! ****
!     ROUTINE USED ONLY BY IBM AND VAX
! ****
 
!     ENTRY SOFIOF
!     NAME = 'SOFIOF'
!     GO TO 100
 
!     THE FOLLOWING THREE ARE FUNCTIONS FOR QUAD WORD OPERATIONS
!                                           (REAL*16)
 ENTRY qabs (i)
 NAME = 'QABS'
 GO TO 100
 
 ENTRY snglq (i)
 NAME = 'SNGLQ'
 GO TO 100
 
 ENTRY dbleq (i)
 NAME = 'DBLEQ'
 GO TO 100
 
 ENTRY qsqrt (i)
 NAME = 'QSQRT'
 GO TO 100
 
 ENTRY qlog (i)
 NAME = 'QLOG'
 GO TO 100
 
 ENTRY qextd (i)
 NAME = 'QEXTD'
 GO TO 100
 
! ****
!     ROUTINE USED BY UNIVAC AND VAX
! ****
 
!     ENTRY DEFCOR
!     NAME = 'DEFCOR'
!     GO TO 100
 
! ****
!     ROUTINES USED BY ALL MACHINES, EXCEPT VAX
! ****
 
 ENTRY gperr (i,j)
 NAME = 'GPERR'
 GO TO 100
 
 ENTRY pdump
 GO TO 250
 
 ENTRY mpy1
 NAME = 'MPY1'
 GO TO 100
 
 ENTRY mpy2nt
 NAME = 'MPY2NT'
 GO TO 100
 
 ENTRY mpy2t
 NAME = 'MPY2T'
 GO TO 100
 
! ****
!     ROUTINES USED ONLY IN CDC MACHINE
! ****
 
 ENTRY link (i,j,k)
 NAME = 'LINK'
 GO TO 100
 
 ENTRY remark (i)
 NAME = 'REMARK'
 GO TO 100
 
 ENTRY cdcbug (i,j,k,l)
 NAME = 'CDCBUG'
 GO TO 100
 
 ENTRY cdcopn (i)
 NAME = 'CDCOPN'
 GO TO 100
 
 ENTRY cdccls (i)
 NAME = 'CDCCLS'
 GO TO 100
 
 ENTRY cdcksz (i)
 NAME = 'CDCKSZ'
 GO TO 100
 
 ENTRY pf (i,j,k)
 NAME = 'PF'
 GO TO 100
 
 ENTRY iswap (i)
 NAME = 'ISWAP'
 GO TO 100
 
! ****
!     ROUTINES USED ONLY IN VAX MACHINE
! ****
 
 ENTRY vaxend
 NAME = 'VAXEND'
 GO TO 100
 
 ENTRY vaxerr (l)
 WRITE  (nout,50) l
 50 FORMAT (/,' *** GINO ERROR AT LOC',i5)
 GO TO 220
 
!     ENTRY VAXSCH
!     NAME = 'VAXSCH'
!     GO TO 100
 
 ENTRY vaxbrk
 NAME = 'VAXBRK'
 GO TO 100
 
!     ENTRY MPY1V (I,J,K)
!     NAME = 'MPY1V'
!     GO TO 100
 
!     ENTRY MPY2NV (I,J,K)
!     NAME = 'MPY2NV'
!     GO TO 100
 
!     ENTRY MPY2TV (I,J,K)
!     NAME = 'MPY2TV'
!     GO TO 100
 
! ****
!     ROUTINES THAT PERFORM NO PARTICULAR FUNCTIONS, BUT THEY
!     ARE STILL CALLED BY NASTRAN
! ****
 
 ENTRY unload (i)
!     CALLED BY INPTT1
 GO TO 250
 
! ****
!     THE FOLLOWING ROUTINES SEEM TO BE NO LONGER USED IN NASTRAN
! ****
 
 ENTRY jidint (i)
 NAME = 'JIDINT'
 GO TO 100
 
 ENTRY opmesg
 NAME = 'OPMESG'
 GO TO 100
 
!     ENTRY PDUM1,PDUM2,...,PDUM9 HAD BEEN REPLACED BY PDUMI
!     ENTRY QDMM3, SQDM31, AND SQDM32 ARE NOW OBSOLETE
 
 ENTRY semtrn
 NAME = 'SEMTRN'
 GO TO 100
 
! ****
!     DUMMY ROUTINES REFERENCED ONLY IN LINK 2, ALL MACHINES
! ****
 
 ENTRY pdumi (*,*,*,i,j,k,l,m,n,o)
 NAME = 'PDUMI'
 GO TO 100
 
! ****
!     DUMMY ROUTINES REFERENCED ONLY IN LINK 5, ALL MACHINES
! ****
 
!HGS      ENTRY PLBAR1 (I,J)
!HGS      NAME = 'PLBAR1'
!HGS      GO TO 100
 
!HGS      ENTRY PLOADX
!HGS      NAME = 'PLOADX'
!HGS      GO TO 100
 
!WKBD ENTRY ERRTRC (NAM)
!     ==================
!     ERROR TRACEBACK
 
!WKBD GO TO 220
 
 100 WRITE  (nout,150) NAME
 150 FORMAT ('0*** SYSTEM FATAL ERROR  ---  JOB TERMINATED',  &
     ' DUE TO CALL TO DUMMY SUBROUTINE.  ENTRY NAME IS ', a8)
 GO TO 220
 
! ****
!     TO FORCE A SYSTEM FATAL ERROR FOR TRACEBACK (VAX ONLY, NOT UNIX)
! ****
 
 220 IF (mach /= 5) GO TO 240
 WRITE  (nout,230)
 230 FORMAT ('0*** ERROR TRACEBACK IN SYSTEM LOG FILE')
 i = 987654321
 n(i) = 0
 240 STOP
 250 RETURN
 
 
!     SUBROUTINE DUMMY
 
! ****
!     UNIVAC  VERSION
 
!     THIS SUBROUTINE PROVIDES ENTRIES FOR THE DUMMY ROUTINES
!     USED BY OTHER COMPUTER MACHINES, AND ARE REFERENCED IN
!     VARIOUS NASTRAN LINKS
 
!     THIS SUBROUTINE INCLUDES ALSO SOME DUMMY ROUTINES NOT YET
!     WRITTEN
 
!     THIS ROUTINE SHOULD BE MOVED TO NASTRAN MACHINE-DEPENDENT
!     SECTION (MDS)
! ****
 
!*    CHARACTER*8     NAME
 
!*    COMMON /MACHIN/ MACH
!*    COMMON /SYSTEM/ ISYSBF, NOUT
 
!*    IF (MACH .EQ. 3) GO TO 250
!*    WRITE  (NOUT,20) MACH
!* 20 FORMAT (/,' MACH =',I7)
!*    NAME = 'DUMMY'
!*    GO TO 100
 
! ****
!     ROUTINES USED ONLY IN UNIVAC MACHINE
! ****
 
!     ENTRY NTRAN (I,J,K)
!     NAME = 'NTRAN'
!     GO TO 100
 
!     ENTRY CONTIN
!     NAME = 'CONTIN'
!     GO TO 100
 
!     ENTRY FACIL (I,J)
!     NAME = 'FACIL'
!     GO TO 100
 
!     ENTRY FACSF (I)
!     NAME = 'FACSF'
!     GO TO 100
 
!     ENTRY UNVOPN (I)
!     NAME = 'UNVOPN'
!     GO TO 100
 
!     ENTRY UNVCLS (I)
!     NAME = 'UNVCLS'
!     GO TO 100
 
!     ENTRY ADDCRD (I,J)
!     NAME = 'ADDCRD'
!     GO TO 100
 
! ****
!     ROUTINES USED BY UNIVAC AND IBM
! ****
 
!*    ENTRY RETURN
!*    GO TO 250
 
!     ENTRY MSGUNI
!     IF (MACH .EQ. 2) GO TO 250
!     NAME = 'MSGUNI'
!     GO TO 100
 
!     ENTRY XEOT (I,J,K,L)
!     IF (MACH .EQ. 2) GO TO 250
!     NAME = 'XEOT'
!     GO TO 100
 
!     ENTRY TPSWIT (I,J,K,L)
!     NAME = 'TPSWIT'
!     GO TO 100
 
! ****
!     ROUTINES USED ONLY IN IBM MACHINE
! ****
 
!*    ENTRY UMFTRN (I)
!*    NAME = 'UMFTRN'
!*    GO TO 100
 
!*    ENTRY TAPSWI (I,J,K,L)
!*    NAME = 'TAPWSI'
!*    GO TO 100
 
!*    ENTRY SEARCH (I)
!*    NAME = 'SEARCH'
!*    GO TO 100
 
!*    ENTRY SOFIOI
!*    NAME = 'SOFIOI'
!*    GO TO 100
 
!*    ENTRY IQZDDN (I)
!*    NAME = 'IQZDDN'
!*    GO TO 100
 
!*    ENTRY QQDCBF (I,J,K,L,M,N)
!*    NAME = 'QQDCBF'
!*    GO TO 100
 
!*    ENTRY QQGETF (I,J,K)
!*    NAME = 'QQGETF'
!*    GO TO 100
 
! ****
!     ROUTINE USED ONLY BY IBM AND VAX
! ****
 
!*    ENTRY SOFIOF
!*    NAME = 'SOFIOF'
!*    GO TO 100
 
!     THE FOLLOWING THREE ARE FUNCTIONS FOR QUAD WORD OPERATIONS
!                                           (REAL*16)
!*    ENTRY QABS (I)
!*    NAME = 'QABS'
!*    GO TO 100
 
!*    ENTRY SNGLQ (I)
!*    NAME = 'SNGLQ'
!*    GO TO 100
 
!*    ENTRY DBLEQ (I)
!*    NAME = 'DBLEQ'
!*    GO TO 100
 
!*    ENTRY QSQRT (I)
!*    NAME = 'QSQRT'
!*    GO TO 100
 
!*    ENTRY QLOG (I)
!*    NAME = 'QLOG'
!*    GO TO 100
 
!*    ENTRY QEXTD (I)
!*    NAME = 'QEXTD'
!*    GO TO 100
 
! ****
!     ROUTINE USED BY UNIVAC AND VAX
! ****
 
!     ENTRY DEFCOR
!     NAME = 'DEFCOR'
!     GO TO 100
 
! ****
!     ROUTINES USED BY ALL MACHINES, EXCEPT VAX
! ****
 
!     ENTRY GPERR (I,J)
!     NAME = 'GPERR'
!     GO TO 100
 
!     ENTRY PDUMP
!     GO TO 250
 
!     ENTRY MPY1
!     NAME = 'MPY1'
!     GO TO 100
 
!     ENTRY MPY2NT
!     NAME = 'MPY2NT'
!     GO TO 100
 
!     ENTRY MPY2T
!     NAME = 'MPY2T'
!     GO TO 100
 
! ****
!     ROUTINES USED ONLY IN CDC MACHINE
! ****
 
!*    ENTRY LINK (I,J,K)
!*    NAME = 'LINK'
!*    GO TO 100
 
!*    ENTRY REMARK (I)
!*    NAME = 'REMARK'
!*    GO TO 100
 
!*    ENTRY CDCBUG (I,J,K,L)
!*    NAME = 'CDCBUG'
!*    GO TO 100
 
!*    ENTRY CDCOPN (I)
!*    NAME = 'CDCOPN'
!*    GO TO 100
 
!*    ENTRY CDCCLS (I)
!*    NAME = 'CDCCLS'
!*    GO TO 100
 
!*    ENTRY CDCKSZ (I)
!*    NAME = 'CDCKSZ'
!*    GO TO 100
 
!*    ENTRY PF (I,J,K)
!*    NAME = 'PF'
!*    GO TO 100
 
!*    ENTRY ISWAP (I)
!*    NAME = 'ISWAP'
!*    GO TO 100
 
! ****
!     ROUTINES USED ONLY IN VAX MACHINE
! ****
 
!*    ENTRY VAXEND
!*    NAME = 'VAXEND'
!*    GO TO 100
 
!*    ENTRY VAXERR (L)
!*    WRITE  (NOUT,50) L
!* 50 FORMAT (/,' *** GINO ERROR AT LOC',I5)
!*    GO TO 220
 
!*    ENTRY VAXSCH
!*    NAME = 'VAXSCH'
!*    GO TO 100
 
!*    ENTRY VAXBRK
!*    NAME = 'VAXBRK'
!*    GO TO 100
 
!*    ENTRY MPY1V (I,J,K)
!*    NAME = 'MPY1V'
!*    GO TO 100
 
!*    ENTRY MPY2NV (I,J,K)
!*    NAME = 'MPY2NV'
!*    GO TO 100
 
!*    ENTRY MPY2TV (I,J,K)
!*    NAME = 'MPY2TV'
!*    GO TO 100
 
! ****
!     ROUTINES THAT PERFORM NO PARTICULAR FUNCTIONS, BUT THEY
!     ARE STILL CALLED BY NASTRAN
! ****
 
!*    ENTRY UNLOAD (I)
!     CALLED BY INPTT1
!*    GO TO 250
 
! ****
!     THE FOLLOWING ROUTINES SEEM TO BE NO LONGER USED IN NASTRAN
! ****
 
!*    ENTRY JIDINT (I)
!*    NAME = 'JIDINT'
!*    GO TO 100
 
!*    ENTRY OPMESG
!*    NAME = 'OPMESG'
!*    GO TO 100
 
!     ENTRY PDUM1,PDUM2,...,PDUM9 HAD BEEN REPLACED BY PDUMI
!     ENTRY QDMM3, SQDM31, AND SQDM32 ARE NOW OBSOLETE
 
!*    ENTRY SEMTRN
!*    NAME = 'SEMTRN'
!*    GO TO 100
 
! ****
!     DUMMY ROUTINES REFERENCED ONLY IN LINK 2, ALL MACHINES
! ****
 
!*    ENTRY PDUMI (*,*,*,I,J,K,L,M,N,O)
!*    NAME = 'PDUMI'
!*    GO TO 100
 
! ****
!     DUMMY ROUTINES REFERENCED ONLY IN LINK 5, ALL MACHINES
! ****
 
!*    ENTRY PLBAR1 (I,J)
!*    NAME = 'PLBAR1'
!*    GO TO 100
 
!*    ENTRY PLOADX
!*    NAME = 'PLOADX'
!*    GO TO 100
 
!*    ENTRY ERRTRC (NAM)
!     ==================
!     ERROR TRACEBACK
 
!*    GO TO 220
 
!*100 WRITE  (NOUT,150) NAME
!*150 FORMAT ('0*** SYSTEM FATAL ERROR  ---  JOB TERMINATED',
!*   1        ' DUE TO CALL TO DUMMY SUBROUTINE.  ENTRY NAME IS ', A8)
!*    GO TO 220
 
! ****
!     TO FORCE A SYSTEM FATAL ERROR FOR TRACEBACK
! ****
 
!*220 WRITE  (NOUT,230)
!*230 FORMAT ('0*** ERROR TRACEBACK IN SYSTEM LOG FILE')
!*    X =-1.0
!*    X = SQRT(X)
!*250 RETURN
 
 
!     SUBROUTINE DUMMY
 
! ****
!     MACHINES 1, AND 6 THRU 20 VERSION
 
!     THIS SUBROUTINE PROVIDES ENTRIES FOR THE DUMMY ROUTINES
!     USED BY OTHER COMPUTER MACHINES, AND ARE REFERENCED IN
!     VARIOUS NASTRAN LINKS
 
!     THIS SUBROUTINE INCLUDES ALSO SOME DUMMY ROUTINES NOT YET
!     WRITTEN
 
!     THIS ROUTINE SHOULD BE MOVED TO NASTRAN MACHINE-DEPENDENT
!     SECTION (MDS)
! ****
 
!.    DIMENSION       N(1)
!.    CHARACTER*8     NAME
 
!.    COMMON /MACHIN/ MACH
!.    COMMON /SYSTEM/ ISYSBF, NOUT
 
!.    IF (MACH.EQ.1 .AND. MACH.GE.6) GO TO RETURN
!.    WRITE  (NOUT,150) NAME,MACH
!.150 FORMAT ('0*** SYSTEM FATAL ERROR  ---  JOB TERMINATED', /5X,
!.   1       'SUBROUTINE DUMMY FOR MACHINE TYPE',I4,' IS NOT AVAILABLE')
!.    I = 987654321
!.    N(I) = 0
!.    STOP
 
END SUBROUTINE dummy
