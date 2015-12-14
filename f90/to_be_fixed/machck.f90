SUBROUTINE machck (*)
     
!     NEW MACHINE COMPATIBILITY CHECK
!     THIS ROUTINE IS CALLED ONCE ONLY BY XSEM01 IF DEBUG FLAG IS ON
 
!     FOR LINK1 DEBUG PURPOSE, PRINT OUT GOES TO UNIT 6, NOT NOUT
 
!     WRITTEN BY G.CHAN/UNISYS  5/1991
 
!     NEXT LINE IS NEEDED FOR HP WORKSTATION. THE $ STARTS ON COLUMN 1
 
!  $MIXED_FORMATS
 
 
 , INTENT(OUT)                            :: *
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,rshift
 INTEGER :: ar(5)
 REAL :: xx
 COMPLEX :: e,d,f
 CHARACTER (LEN=8) :: lrv
 CHARACTER (LEN=8) :: 2
 CHARACTER (LEN=29) :: uimx
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim
 COMMON /sem   / mask
 COMMON /system/ sysbuf,nout,nogo,dumm1(11),date(3),dumm2(13),  &
     hicore,timew,dumm3(62),sperlk
 COMMON /machin/ machx,ijhalf(2),lqro,mchnam
 COMMON /lhpwx / lowpw,highpw
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (s,z(1))
 EQUIVALENCE     (l21(1),l2),(xx,ix)
 DATA    i,j,k / 4HABCD, 4H1234, 4HA3CD  /,  &
     d,e   / (1.0,-2.0), (-3.0,4.0)  /, l1,bt / 4HWORD, 4HBYTE          /,  &
     ia,ir / 4HA   , 4HR             /, l2,rv / 'NATURAL ', 'REVERSED'  /,  &
     as,ds / 3H as,  3HDES           /
 DATA    uimx  / '0*** USER INFORMATION MESSAGE' /
 
 nogo = 0
 IF (uimx == uim) GO TO 20
 nogo = 1
 WRITE  (6,10)
 10 FORMAT (/,' -LINK1 DEBUG- SEMDBD DATA BLOCK NOT LOADED CORRECTLY')
 
!     CALL BTSTRP TO INITIALIZE MACHINE CONSTANTS
 
 20 WRITE  (6,30)
 30 FORMAT (/,' -LINK1 DEBUG- MACHCK CALLING BTSTRP NEXT')
 CALL btstrp
 
 WRITE  (6,40)
 40 FORMAT ( ' -LINK1 DEBUG-  CHARACTER, SHIFT AND COMPLEX CHECKS')
 l  = khrfn1(i,2,j,3)
 IF (l == k) GO TO 60
 nogo = 2
 WRITE  (6,50) i,j,k,l
 50 FORMAT (' * KHRFN1 ERROR   I,J,K,L =',4(1X,a4))
 60 i  = 128
 j  = lshift(i,2)
 k  = rshift(i,2)
 IF (j == 512 .AND. k == 32) GO TO 80
 nogo = 3
 WRITE  (6,70)
 70 FORMAT (' * LSHIFT AND/OR RSHIFT ERROR')
 
!     JUMP TO 100 IF MACHINE DOES NOT HAVE ISHFT FUNCTION
 
 80 j  = ishft(i,+2)
 k  = ishft(i,-2)
 IF (j == 512 .AND. k == 32) GO TO 110
 nogo = 4
 IF (j /= 512) WRITE (6,90)
 IF (k /=  32) WRITE (6,100)
 90 FORMAT (' * ISHFT(+) NOT SAME AS LSHIFT')
 100 FORMAT (' * ISHFT(-) NOT SAME AS RSHIFT')
 
!     CHECK ISHFT IS ZERO-FILL
 
 110 i  = -1
 j  = ishft(i,-1)
 k  = ishft(i,+1)
 IF (j > 0 .AND. MOD(k,2) == 0) GO TO 130
 nogo = 5
 WRITE  (6,120)
 120 FORMAT (' * SYSTEM ISHFT IS NOT ZERO-FILL')
 
!     CHECK K2B SUBROUTINE
 
 130 CALL k2b (l21,ar,5)
 IF (ar(2) == ia .AND. ar(5) == ir) GO TO 150
 nogo = 6
 WRITE  (6,140) ar(2),ar(5)
 140 FORMAT (' * K2B ERROR   A,R ==',2A4)
 
!     COMPLEX NUMBER CHECK
 
 150 f  = d*e
 dr = REAL (d)
 di = AIMAG(d)
 er = REAL (e)
 ei = AIMAG(e)
 a  = dr*er - di*ei
 b  = dr*ei + di*er
 IF (ABS(a-REAL(f)) <= .01 .AND. ABS(b-AIMAG(f)) <= .01) GO TO 170
 nogo = 7
 WRITE  (6,160)
 160 FORMAT (' * COMPLEX ERROR')
 170 IF (mask == 65535) GO TO 190
 nogo = 8
 WRITE  (6,180)
 180 FORMAT (' * LABEL COMMON /SEM/ ERROR')
 190 IF (sperlk == 1 .OR. sperlk == 0) GO TO 210
 nogo = 9
 WRITE  (6,200)
 200 FORMAT (' * LABEL COMMON /SYSTEM/ ERROR')
 210 IF (nogo == 0) WRITE (6,220)
 220 FORMAT ('  OK')
 
!     LOGICAL 'AND' AND 'OR' CHECK.
!     SYSTEM MAY NAME THESE FUNCTIONS - 'IAND', 'IOR', OR 'AND', 'OR'
!     IF UNSATISFIED EXTERNALS OCCUR, FIX THEM HERE AND IN MAPFNS.MDS
 
 WRITE  (6,230)
 230 FORMAT ('  LOGICAL "AND" AND "OR" CHECK.  IF ERROR OCCURS, ',  &
     'SEE MACHCK')
 k  = IAND(i,j)
 k  =  ior(i,j)
 WRITE (6,220)
 
!     CHECK DATE AND TIME, ALREADY SAVED IN /SYSTEM/ BY NASTRN OR NAST01
 
!     TIME IS SYSTEM CPU TIME, COMMONLY IN 1/60 SECONDS ACCURCY
!     IF UNSATISFIED EXTERNALS OCCUR, FIX THEM IN TDATE, KLOCK, WALTIM,
!     CPUTIM, AND/OR SECNDS(IN MAPFNS) SUBROUTINES
 
 WRITE  (6,240)
 240 FORMAT ('  DATE AND TIME CHECKS.  IF ERROR OCCURS, SEE MACHCK')
 i  = timew/3600
 j  = (timew-i*3600)/60
 k  = timew - i*3600 - j*60
 WRITE  (6,250) date,i,j,k
 250 FORMAT (' -MONTH/DAY/YEAR = ',i2,1H/,i2,1H/,i2, 3X,  &
     ' -HOUR:MIN:SEC = ',i2,':',i2,':',i2)
 IF (date(1) > 12 .OR. date(3) > 1000) WRITE (6,260)
 260 FORMAT (' * SYSTEM DATE SHOULD BE IN mm,dd,yy', ' ORDER  <===')
 
 r  = MOD(lqro,100)/10
 IF (r > 1) l1 = bt
 IF (MOD(lqro,10) == 1) l2 = rv
 WRITE  (6,270) machx,mchnam,l1,l2,sysbuf,lqro
 270 FORMAT (/,' -MACHINE =',i3,2H, ,a4,', RECL BY ',a4,  &
     'S, BCD WORD IN ',a8,' ORDER,', /3X,'SYSBUF =',i7, ' WORDS,  LQRO =',i7)
 
!     OPEN A DIRECT FILE, FORTRAN UNIT 41, AND TEST FOR RECORD LENGTH
 
 i  = sysbuf - 3
 IF (machx == 3 .OR. machx >= 5) i = sysbuf - 4
 i  = i*r
 OPEN  (UNIT=41,ACCESS='DIRECT',RECL=i,STATUS='SCRATCH',ERR=310)
 WRITE (41,REC=1,ERR=280) (z(j),j=1,i)
 GO TO 300
 280 IF (r > 1) GO TO 300
 nogo = 10
 WRITE  (6,290) r
 290 FORMAT (' * FORTRAN I/O RECORD LENGTH IN BTSTRP MAY BE IN ERROR.',  &
     5X,'R =',i4)
 300 CLOSE  (UNIT=41)
 
!     CHECK OPEN CORE IN MEMORY                          ** NEW, NEXT 28
 
 310 j  = 11
 i  = locfx(z(j))
 j  = locfx(z(1))
 k  = as
 IF (i < j) k = ds
 WRITE  (6,320) k,j,i
 320 FORMAT (' * SYSTEM MEMORY IN ',a3,'CENDING ORDER',i15,'==>',i12)
 
!     CHECK WHETHER NUMTYP.MIS IS SET UP FOR THIS CURRENT MACHINE
 
 k  = 123
 IF (numtyp(k) /= 1) nogo = 11
 
!     CHECK /SOFPTR/ LOCATION WITH RESPECT TO /ZZZZZZ/ LOCATION HERE IF
!     AND ONLY IF CURRENT NASTRAN VERSION STILL USES /SOFPTR/, AND
!     SET K = 1
!                  K J          I           I          J K
!               ---+-+----------+   OR   ---+----------+-+
!                    ASCENDING                 DECENDING
 
 k  = 0
 IF (k /= 1) GO TO 340
 
 k  = locfx(s)
 IF (i > j .AND. k > j) WRITE (6,330)
 IF (i < j .AND. k < j) WRITE (6,330)
 330 FORMAT (' * COMMONS /SOFPTR/ AND /ZZZZZZ/ POSITIONS SHOULD BE ',  &
     'REVERSED IN OPNCOR.MDS')
 
!     CHECK S.P. NUMERIC RANGE
 
 340 IF (10.0**(lowpw+1) >= 0.0 .AND. 10.0**(highpw-1) > 10.0**36) GO TO 360
 nogo = 12
 WRITE  (6,350) lowpw,highpw
 350 FORMAT (' * MACHINE NUMERIC RANGE, 10.**',i3,' THRU 10.**',i2,  &
     ' SET BY BTSTRP, EXCEEDS MACHINE LIMIT.')
 
!     CHECK FORTRAN MIXED FORMAT WRITE USED IN SUBOURINTE OFPPNT OF THE
!     OFP MODULE.
!     DEC/ULTRIX FORTRAN 3.0 (1992) FAILS ON THIS TEST.
 
 360 ix = 123456
 WRITE  (6,370,ERR=380) xx
 370 FORMAT (/,i10)
 380 WRITE  (6,390)
 390 FORMAT (' IF 123456 IS NOT PRINTED ON ABOVE LINE, MIXED FORMAT ',  &
     'PRINT OUT IS NOT', /1X,  &
     'ALLOWED, AND NASTRAN OFP MODULE MAY NOT WORK PROPERLY')
 
!     CHECK OPEN CORE
 
 j  = 5000
 z(j) = 1
 z(hicore) = 2
 
 IF (nogo /= 0) GO TO 410
 WRITE  (6,400) uim
 400 FORMAT (a29,', MACHINE COMPATIBILITY CHECK ROUTINE FINDS NO ',  &
     /5X,'SIGNIFICANT SYSTEM ERROR')
 RETURN 1
 
 410 WRITE  (6,420) uim,nogo
 420 FORMAT (a29,' * ERROR IN MACHCK.  NOGO =',i3)
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE machck
