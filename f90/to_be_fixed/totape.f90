SUBROUTINE totape (caller,z)
     
!     THIS ROUTINE IS CALLED ONLY BY DPLTST (CALLER=1), DPLOT (CALLER=2)
!     AND/OR OFP (CALLER=3) TO COPY NUMBERS OF GINO INPUT FILES TO A
!     SAVE FILE, INP9, FOR NEXT INTERACTIVE NASTRAN RUN (INTRA.LT.0)
!     THE SAVE FILE CAN BE A TAPE OR DISC.
 
!     WRITTEN BY G.CHAN/SPERRY     NOV. 1985
 
!     FILE STRUCTURE IN SAVE TAPE
 
!     RECORD NO.   CONTENT
!     ----------   -----------------------------------------------
!        1           6-WORDS (3 CALLER ID WORDS AND 3 DATE WORDS)
!                   96-WORD HEADING
!                  100-SYSTEM WORDS
!        2         MARK
!        3         CALLER ID, NO. OF FILES, NO. OF PARAMETERS
!        4         7-WORD TRAILER OF FIRST GINO INPUT FILE
!      5 TO N      FIRST GINO INPUT FILE (IF FILE IS NOT PURGED)
!       N+1        MARK
!       N+2        7-WORD TRAILER OF SECOND GINO INPUT FILE
!    N+2 TO M      SECOND GINO INPUT FILE (IF FILE IS NOT PURGED)
!       M+1        MARK
!    M+2 TO ..R    REPEAT FOR ADDITION FILES, TRAILER, AND MARK
!       R+1        PARAMETERS IN /BLANK/ OF CURRENT CALLER
!       R+2        MARK
!       R+3        NASTRAN EOF MARK
!    R+4 TO LAST   REPEAT 3 TO R+3 AS MANY TIMES AS NEEDED FROM THE
!                  SAME OR A DIFFERENT CALLER AT DIFFERENT TIME
!     LAST+1       SYSTEM EOF MARK
 
!     THE INTERACTIVE FLAG, INTRA, IN /SYSTEM/ WAS SET BY XCSA TO
!         1 FOR PLOT ONLY,
!         2 FOR OUTPUT PRINT ONLY
!      OR 3 FOR BOTH
 
 
 INTEGER, INTENT(IN)                      :: caller
 INTEGER, INTENT(OUT)                     :: z(3)
 IMPLICIT INTEGER (a-z)
 LOGICAL :: disc,     tapbit
 DIMENSION  tab(3,3), mark(3),  sub(2),   fn(2), date(3),  who(2)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,      uwm,      uim
 COMMON /BLANK / param(1)
 COMMON /system/ ksystm(100)
 COMMON /output/ head(96)
 COMMON /names / rd,       rdrew,    wrt,      wrtrew,   rew, norew
 EQUIVALENCE     (ksystm( 1),ibuf), (ksystm(15),date(1)),  &
     (ksystm( 2),nout), (ksystm(86),intra  ), (tab(2,3) ,BLANK)
 DATA    tab   / 4HPLTS,   4HET  ,   2, 4HPLOT,   4H    ,   5,  &
     4HOFP ,   4H    ,   3/
 DATA    FILE,   nfile,    mark  /   4HINP9,23,  2*65536,11111  /
 DATA    sub /   4HTOTA,   4HPE  /
 
 IF (intra >= 0 .OR. caller < 1 .OR. caller > 3) RETURN
 IF (caller <= 2 .AND. intra == -2) RETURN
 IF (caller == 3 .AND. intra == -1) RETURN
 who(1) = tab(1,caller)
 who(2) = tab(2,caller)
 nparam = tab(3,caller)
 kore   = korsz(z(1))
 ibuf1  = kore  - ibuf
 ibuf2  = ibuf1 - ibuf
 kore   = ibuf2 - 1
 fn(1)  = FILE
 fn(2)  = BLANK
 disc   = .true.
 IF (tapbit(fn(1))) disc = .false.
 IF (.NOT.disc .OR. intra > 0) GO TO 30
 
 CALL OPEN (*120,FILE,z(ibuf2),rdrew)
 10   CALL READ (*20,*20,FILE,z(1),2,0,m)
 CALL skpfil (FILE,1)
 GO TO 10
 20   CALL CLOSE (FILE,norew)
 30   CALL OPEN (*120,FILE,z(ibuf2),wrt)
 IF (intra < 0) GO TO 40
 DO  i = 1,2
   IF (intra /= i .AND. intra /= 3) CYCLE
   FILE = tab(3,i+1)
   CALL WRITE (FILE,tab(1,caller),3,0)
   CALL WRITE (FILE,  date(1),  3,0)
   CALL WRITE (FILE,  head(1), 96,0)
   CALL WRITE (FILE,ksystm(1),100,1)
   CALL WRITE (FILE,  mark(1),  3,1)
 END DO
 intra = -intra
 FILE  = tab(3,caller)
 40   z(1)  = caller
 z(2)  = nfile
 z(3)  = nparam
 CALL WRITE (FILE,z(1),3,1)
 WRITE  (nout,50) uim,who,FILE
 50   FORMAT (a29,', THE FOLLOWING FILES WERE COPIED FROM DMAP ',a4,a2,  &
     4H TO ,a4,5H FILE,/)
 DO  i = 1,nfile
   infil = 100 + i
   CALL OPEN (*100,infil,z(ibuf1),rdrew)
   z(1)  = infil
   CALL rdtrl (z(1))
   CALL WRITE (FILE,z(1),7,1)
   IF (z(1) <= 0) GO TO 80
   60   CALL READ (*80,*70,infil,z(1),kore,1,m)
   CALL mesage (-8,0,sub)
   70   CALL WRITE (FILE,z(1),m,1)
   GO TO 60
   80   CALL CLOSE (infil,rew)
   CALL fname (infil,fn)
   WRITE  (nout,90) fn
   90   FORMAT (5X,2A4)
   100  CALL WRITE (FILE,mark(1),3,1)
 END DO
 CALL WRITE (FILE,param(1),nparam,1)
 CALL WRITE (FILE,mark(1),3,1)
 IF (.NOT.disc) CALL CLOSE (FILE,norew)
 IF (     disc) CALL CLOSE (FILE,  rew)
 RETURN
 
 120  CALL mesage (-1,FILE,sub)
 RETURN
END SUBROUTINE totape
