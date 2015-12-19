SUBROUTINE mplprt
     
!     PRINTS MPL TABLE FOR DOCUMENTATION PURPOSES
!     AND CHECKS VALIDITY OF MANY ITEMS.
 
 DOUBLE PRECISION :: xx
 REAL :: x(2,1)
 INTEGER :: kp(6),flag,flagb,flags,tot,flgtot,add(2),  &
     t1,t2,t3,h1,h2,h3,h1x(32),h2x(32),h3x(32)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm,swm
 COMMON /system/ nb,no,junk1(6),nlpp,junk2(2),line
 COMMON /xfist / nfist
 COMMON /xpfist/ npfist
 COMMON /output/ t1(32),t2(32),t3(32),h1(32),h2(32),h3(32)
 COMMON /xgpi2 / lmpl,mplpnt,mpl(1)
 COMMON /xgpi2x/ xx(1)
 
 EQUIVALENCE    (xx(1),x(1,1))
 
 DATA kp    / 1,1,2,2,2,4 / ,  add   /4HADD ,4H    /
 DATA flagb / 1H  /  ,  flags /4H ***/
 DATA h1x/4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    &
     ,4H   m,4H o d,4H u l,4H e  ,4H p r,4H o p,4H e r,4H t i    &
     ,4H e s,4H   l,4H i s,4H t  ,4H    ,4H    ,4H    ,4H        &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    /
 DATA h2x/4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H   -,4H - -,4H - -    &
     ,4H - -,4H p a,4H r a,4H m e,4H t e,4H r s,4H - -,4H - -    &
     ,4H - -,4H - -,4H    ,4H    ,4H    ,4H    ,4H    ,4H    /
 DATA h3x/4H  mp,4HLID ,4HNWDS,4H  wd,4H1  m,4HOD-n,4HAME ,4HTYP &
     ,4H in ,4H out,4H  sc,4HR  t,4HOT  ,4H   i,4HD ty,4HP       &
 ,4HP   ,4H   d,4HEFAU,4HLT (,4HIF a,4HNY) ,4H    ,4H   w        &
     ,4H1-w2,4H flg,4H    ,4H    ,4H    ,4H    ,4H    ,4H    /
 
 21 FORMAT (7H0      ,3I5,2X ,2A4 ,i3,4I5  ,10X,a3)
 22 FORMAT (7H0      ,3I5,2X ,2A4 ,i3,4I5  &
     ,10X,50H----- n o   p a r a m e t e r s   e x i s t ----- ,10X,a3)
 23 FORMAT (7H0      ,3I5,2X,8H (NONE) )
 24 FORMAT (7H0      ,3I5,2X,2A4,i3    )
 31 FORMAT (59X,i2,5H. INT,i5,7X,16H-- no default --,6X ,i2          )
 32 FORMAT (59X,i2,5H. rsp,i5,7X,16H-- no default --,6X ,i2          )
 33 FORMAT (59X,i2,5H. bcd,i5,7X,16H-- no default --,6X ,i2,1H-,i2   )
 34 FORMAT (59X,i2,5H. rdp,i5,7X,16H-- no default --,6X ,i2,1H-,i2,a4)
 35 FORMAT (59X,i2,5H. csp,i5,7X,16H-- no default --,6X ,i2,1H-,i2,a4)
 36 FORMAT (59X,i2,5H. cdp,i5,7X,16H-- no default --,6X ,i2,1H-,i2,a4)
 41 FORMAT (59X,i2,5H. INT,i5, i15,14X                  ,i2          )
 42 FORMAT (59X,i2,5H. rsp,i5, 1P,e20.4,9X              ,i2          )
 43 FORMAT (59X,i2,5H. bcd,i5, 11X,2A4,10X              ,i2,1H-,i2   )
 44 FORMAT (59X,i2,5H. rdp,i5, 1P,d20.4,9X              ,i2,1H-,i2,a4)
 45 FORMAT (59X,i2,5H. csp,i5, 3H  (,1P,e11.4,1H, ,1P,e11.4,3H   ,  &
     i2,1H-,i2,a4)
 46 FORMAT (59X,i2,5H. cdp,i5, 3H  (,1P,d11.4,1H, ,1P,d11.4,3H)  ,  &
     i2,1H-,i2,a4)
 47 FORMAT (10X,'NOTE - THE ABOVE PARAMETER DEFAULTS WILL BE CHANGED',  &
 ' TO ALL ZEROS BY THE ADD MODULE.  HOWEVER, IF ALL 4 PARAMETERS',  &
     ' ARE NOT', /10X,'SPECIFIED, THEY WILL BE CHANGED TO 2*(1.,0.),',  &
     ' 2*(0.D0,0.D0), OR 2*(0.,0.), 2*(1.D0,0.D0) DEPENDING ON ',  &
     'MATRICES INVOLVED')
 
!     INITIALIZATION
 
 CALL page
 mplid = 0
 npad  = 0
 i2    = 0
 DO  i = 1,32
   h1(i) = h1x(i)
   h2(i) = h2x(i)
   h3(i) = h3x(i)
 END DO
 CALL page
 
!     PROCESS NEXT ENTRY
 
 100 CONTINUE
 IF (i2-lmpl < 0) THEN
   GO TO   110
 ELSE IF (i2-lmpl == 0) THEN
   GO TO   900
 ELSE
   GO TO  9901
 END IF
 110 i0 = i2
 i1 = i2 + 1
 i2 = i0 + mpl(i1)
 mplid = mplid + 1
 
!     TEST FOR MODULE TYPE
 
 IF (mpl(i1)-1 < 0) THEN
   GO TO  9904
 ELSE IF (mpl(i1)-1 == 0) THEN
   GO TO   120
 END IF
 112 IF (mpl(i1+1) == 0) GO TO 120
 IF (mpl(i0+4) < 3) GO TO 130
 
!     EXECUTIVE MODULE
 
 CALL page2 (-2)
 l1 = i0 + 2
 l2 = l1 + 2
 WRITE (no,24) mplid,mpl(i1),i1,(mpl(l),l=l1,l2)
 GO TO 100
 
!     PAD SPACE
 
 120 CALL page2 (-2)
 WRITE (no,23) mplid,mpl(i1),i1
 npad = npad + 1
 GO TO 100
 
!     FUNCTIONAL MODULE
 
 130 IF (mpl(i1) > 7) GO TO 140
 
!     NO PARAMETERS EXIST FOR THIS FUNCTIONAL MODULE
 
 CALL page2 (-2)
 l1  = i0 + 5
 l2  = l1 + 2
 tot = 0
 DO  l = l1,l2
   tot = tot + mpl(l)
 END DO
 flgtot = flagb
 IF (tot > nfist-npfist) flgtot = flags
 l1  = i0 + 2
 WRITE (no,22) mplid,mpl(i1),i1,(mpl(l),l=l1,l2),tot,flgtot
 GO TO 100
 
!     PARAMETERS EXIST FOR THIS FUNCTIONAL MODULE
 
 140 CONTINUE
 
!     DETERMINE THE NUMBER OF PARAMETERS FOR FUNCTIONAL MODULE
 
 np = 0
 i  = i0 + 8
 150 CONTINUE
 IF ((i-1)-(i2) < 0) THEN
   GO TO   151
 ELSE IF ((i-1)-(i2) == 0) THEN
   GO TO   160
 ELSE
   GO TO  9903
 END IF
 151 ip = IABS(mpl(i))
 IF (ip > 6) GO TO 9902
 IF (mpl(i) < 0) THEN
   GO TO   152
 ELSE IF (mpl(i) == 0) THEN
   GO TO  9902
 ELSE
   GO TO   154
 END IF
 152 np = np + 1
 i  = i  + 1
 GO TO 150
 154 np = np + 1
 i  = i + 1 + kp(ip)
 GO TO 150
 160 IF (np <= 0) GO TO 9903
 
 CALL page2 (-2-np)
 l1 = i0 + 5
 l2 = l1 + 2
 tot= 0
 DO  l = l1,l2
   tot = tot + mpl(l)
 END DO
 flgtot = flagb
 IF (tot > nfist-npfist) flgtot = flags
 l1 = i0 + 2
 WRITE (no,21) mplid,mpl(i1),i1,(mpl(l),l=l1,l2),tot,flgtot
 
!     PRINT PARAMETERS
 
 np = 0
 i  = i0 + 8
 j2 = 0
 170 CONTINUE
 np = np + 1
 j1 = j2 + 1
 IF ((i-1)-(i2) < 0) THEN
   GO TO   175
 ELSE IF ((i-1)-(i2) == 0) THEN
   GO TO   200
 ELSE
   GO TO  9903
 END IF
 175 ip = IABS(mpl(i))
 IF (ip > 6) GO TO 9902
 IF (mpl(i) < 0) THEN
   GO TO   180
 ELSE IF (mpl(i) == 0) THEN
   GO TO  9902
 ELSE
   GO TO   190
 END IF
 
!     PARAMETER HAS NO DEFAULT VALUE
 
 180 CONTINUE
 j2 = j1
 SELECT CASE ( ip )
   CASE (    1)
     GO TO 181
   CASE (    2)
     GO TO 182
   CASE (    3)
     GO TO 183
   CASE (    4)
     GO TO 184
   CASE (    5)
     GO TO 185
   CASE (    6)
     GO TO 186
 END SELECT
 
!     INTEGER
 
 181 WRITE (no,31) np,i,j1
 GO TO 188
 
!     REAL SINGLE-PRECISION
 
 182 WRITE (no,32) np,i,j1
 GO TO 188
 
!     ALPHANUMERIC (BCD)
 
 183 j2 = j2 + 1
 WRITE (no,33) np,i,j1,j2
 GO TO 188
 
!     REAL DOUBLE-PRECISION
 
 184 j2 = j2 + 1
 flag = flagb
 IF (MOD(j1,2) == 0) flag = flags
 WRITE (no,34) np,i,j1,j2,flag
 GO TO 188
 
!     COMPLEX SINGLE-PRECISION
 
 185 j2   = j2 + 1
 flag = flagb
 IF (MOD(j1,2) == 0) flag = flags
 WRITE (no,35) np,i,j1,j2,flag
 GO TO 188
 
!     COMPLEX DOUBLE-PRECISION
 
 186 j2   = j2 + 3
 flag = flagb
 IF (MOD(j1,2) == 0) flag = flags
 WRITE (no,36) np,i,j1,j2,flag
 GO TO 188
 188 CONTINUE
 i = i + 1
 GO TO 170
 
!     PARAMETER HAS A DEFAULT VALUE
 
 190 CONTINUE
 SELECT CASE ( ip )
   CASE (    1)
     GO TO 191
   CASE (    2)
     GO TO 192
   CASE (    3)
     GO TO 193
   CASE (    4)
     GO TO 194
   CASE (    5)
     GO TO 195
   CASE (    6)
     GO TO 196
 END SELECT
 
!     INTEGER
 
 191 j2 = j1
 WRITE (no,41) np,i,mpl(i+1),j1
 i = i + 2
 GO TO 198
 
!     REAL SINGLE-PRECISION
 
 192 j2 = j1
 m  = mpl(i+1)
 WRITE (no,42) np,i,x(1,m),j1
 i  = i + 2
 GO TO 198
 
!     ALPHANUMERIC (BCD)
 
 193 j2 = j1 + 1
 WRITE (no,43) np,i,mpl(i+1),mpl(i+2),j1,j2
 i = i + 3
 GO TO 198
 
!     REAL DOUBLE-PRECISION
 
 194 j2 = j1 + 1
 m  = mpl(i+1)
 flag = flagb
 IF (MOD(j1,2) == 0) flag = flags
 WRITE (no,44) np,i,xx(m),j1,j2,flag
 i = i + 3
 GO TO 198
 
!     COMPLEX SINGLE-PRECISION
 
 195 j2 = j1 + 1
 m  = mpl(i+1)
 flag = flagb
 IF (MOD(j1,2) == 0) flag = flags
 WRITE (no,45) np,i,x(1,m),x(2,m),j1,j2,flag
 i  = i + 3
 GO TO 198
 
!     COMPLEX DOUBLE-PRECISION
 
 196 j2 = j1 + 3
 m1 = mpl(i+1)
 m2 = mpl(i+3)
 flag = flagb
 IF (MOD(j1,2) == 0) flag = flags
 WRITE (no,46) np,i,xx(m1),xx(m2),j1,j2,flag
 i = i + 5
 GO TO 198
 198 CONTINUE
 GO TO 170
 
 200 CONTINUE
 IF (mpl(l1) /= add(1) .OR. mpl(l1+1) /= add(2)) GO TO 100
 CALL page2 (-2)
 WRITE (no,47)
 GO TO 100
 
!     TERMINATION
 
 900 CONTINUE
 CALL page2 (-4)
 WRITE  (no,901)
901 FORMAT ('0*** END OF MPL PRINTOUT')
WRITE  (no,902) mplid,npad
902 FORMAT ('0*** THE MPL CONTAINS ',i3,' ENTRYS.  OF THESE, ',i3,  &
    ' ARE PAD ENTRYS.')

RETURN

!     ERROR MESSAGES

9901 WRITE  (no,9951) swm,i2,lmpl
9951 FORMAT (a27,' 65, POINTER I2 =',i10,' DOES NOT AGREE WITH LMPL =', i11)
GO TO 9995

9902 WRITE  (no,9952) swm
9952 FORMAT (a27,' 66, ILLEGAL PARAMETER TYPE CODE.')
GO TO 9995

9903 WRITE  (no,9953) swm
9953 FORMAT (a27,' 67, ERROR IN PARAMETER SEQUENCE.')
GO TO 9995

9904 WRITE  (no,9954) swm
9954 FORMAT (a27,' 68, ILLEGAL WORD COUNT.')


9995 CALL page2 (4)
WRITE  (no,9996)
9996 FORMAT (5X,'MPL TABLE LISTING CANCELLED.')

RETURN
END SUBROUTINE mplprt
