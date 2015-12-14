SUBROUTINE ifpmdc
     
!     IFPMDC MODIFIES BULK DATA CARDS GIVEN THE INFORMATION ON IFIL
 
 EXTERNAL        lshift,rshift,andf,orf
 LOGICAL :: abort,cf,diag
 INTEGER :: ret,andf,orf,rshift,t1,cnt,dum,x,exi,test,apprch,  &
     ick(6),ivc(2),inc(2),xi(2),con(38)
 DIMENSION       rm(1),rm1(1),cd(6)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /ifpx1 / ncds,t1(2,1)
 COMMON /machin/ mach
 COMMON /system/ ibuf,nout,abort,dum(17),apprch,dum1(17),nbits,x,  &
     ncpw,dum2(41),jrun
 COMMON /ifpdta/ id(2),kn,d1(52),m(50),mf(50),m1(35),m1f(35),  &
     d2(3),nopen,d3(6),knt,d4(18)
 COMMON /two   / itwo(32)
 COMMON /ifpx0 / lbd,lcc,ibits(18),iparpt
 COMMON /zzzzzz/ kor(1)
 COMMON /xsrtcm/ im1(6),im2(5),im3(4),im4(4),im5(2),im6,im7(8),im8,  &
     isft(4),im9(7),isfim,im10(3),mis
 EQUIVALENCE (rm(1),m(1)),(rm1(1),m1(1)),(ick(1),cd(1)),(k,ick(1))
 DATA con/4H    ,4H   0,4H   1,4H   2,4H   3,4H   4,4H   5,4H   6,  &
     4H   7,4H   8,4H   9,4H   a,4H   b,4H   c,4H   d,4H   e,4H   f,  &
     4H   g,4H   h,4H   i,4H   j,4H   k,4H   l,4H   m,4H   n,4H   o,  &
     4H   p,4H   q,4H   r,4H   s,4H   t,4H   u,4H   v,4H   w,4H   x,  &
     4H   y,4H   z,4H    /
 DATA    ifil, ieof, icycl,   iefm,   iend,   diag  /  &
     213,    0,     0, -32767, 4HZZZZ, .false. /
 
 IF (ieof == -1) GO TO 190
 cnt = 0
 IF (ieof ==  1) GO TO 10
 
!     FIRST CALL INITIALIZE OPEN FILE ADJUST CORE
 
 ibuf1 = nopen + 2*ibuf
 nopen = nopen - ibuf
 DO  i = 1,38
   con(i) = andf(con(i),im3(4))
 END DO
 IF (nopen <= 0) GO TO 1001
 cf   = .false.
 iod  = 0
 isc  = 0
 nf   = 0
 ionf = 0
 ilst = 0
 ieof = 1
 CALL OPEN (*1002,ifil,kor(ibuf1+1),0)
 5 CALL READ (*180,*180,ifil,ick,6,0,nw)
 
!     CHECK INCOMING  CALL FOR VARY MATCH SORT, UNSORT AND/OR CONT
 
 10 IF (k == kn) GO TO 20
 
!     NOT CARD WE ARE WORKING ON CHECK ALPH POSITION
 
 IF (cf .OR. iod == kn) GO TO 190
 iod = kn
 isc = 0
 ASSIGN 15 TO exi
 xi(1) = t1(1,k)
 xi(2) = t1(2,k)
 GO TO 100
 15 ivc(1) = xi(1)
 ivc(2) = xi(2)
 ASSIGN 16 TO exi
 xi(1) = t1(1,kn)
 xi(i) = t1(2,kn)
 GO TO 100
 16 inc(1) = xi(1)
 inc(2) = xi(2)
 IF (mach == 2) GO TO 18
 inc(1) = rshift(inc(1),1)
 ivc(1) = rshift(ivc(1),1)
 18 IF (inc(1) < ivc(1)) GO TO 190
 IF (inc(1) > ivc(1)) GO TO 1004
 
!     SHIFT IN CASE OF STAR
 
 inc(2) = rshift(inc(2),nbits)
 ivc(2) = rshift(ivc(2),nbits)
 IF (inc(2)-ivc(2) < 0) THEN
   GO TO   190
 ELSE
   GO TO  1004
 END IF
 
!     CARD TYPE FOUND TRY ID
 
 20 IF (ick(2) < 0) GO TO 70
 IF (cf .AND. nf /= 0 .AND. ilst == ick(2) .AND. cnt == 1) GO TO 31
 IF (cf .AND. nf /= 0 .AND. ilst == ick(2)) GO TO 25
 IF (cf) GO TO 190
 nf   = 0
 ionf = 0
 ASSIGN 5 TO ret
 IF (m(1) < ick(2)) GO TO 190
 IF (m(1) > ick(2)) GO TO 1004
 ilst = ick(2)
 
!     FIND FIELD FORMAT DOES NOT COUNT FOR FIELD  1 OR 10 K1=COUNT
 
 25 DO  i = 1,50
   IF (mf(i) == iefm) GO TO 30
 END DO
 GO TO 1002
 30 nf  = nf + i - 1
 cnt = 1
 31 k1  = ick(3)
 
!     FIND NUMBER OF FIELDS TO PITCH
 
 i  = k1/10
 j  = (k1-1)/10
 k1 = k1 - i - j - 1
 
!     CHECK TO SEE IF WE HAVE IT NOW
 
 IF (k1 > nf) GO TO 60
 
!     CHECK FORMAT FIELD FOR TYPE
 
 k1 = k1 - ionf
 IF (mf(k1) /= 2 .AND. mf(k1) /= 0) GO TO 1003
 j = 0
 DO  i = 1,k1
   j = j + 1
   IF (mf(i) > 2) j = j +1
 END DO
 
!     PERFORM VARY
 
 IF (cd(6) == 0.0) GO TO 38
 rm(j) = rm(j)*(1.0 + cd(4)*cd(5))**cd(6)
 IF (diag) WRITE (nout,1000) uim,t1(1,k),t1(2,k),knt,ick(2),ick(3), rm(j)
 GO TO 40
 38 rm(j ) = rm(j) + cd(4)*cd(5)
 mf(k1) = 2
 IF (diag) WRITE (nout,1000) uim,t1(1,k),t1(2,k),knt,ick(2),ick(3), rm(j)
 GO TO 40
 
!     SET RESTART BITS
 
 40 IF (apprch >= 0) GO TO 50
 
!     CHECK FOR PARAM CARDS (82)
 
 IF (kn /= 82) GO TO 45
 DO  i = iparpt,ncds
   IF (m(1) == t1(1,i) .AND. m(2) == t1(2,i)) GO TO 42
 END DO
 GO TO 50
 42 j = i - 1
 GO TO 46
 45 j = kn - 1
 46 karl = 1
 IF (icycl == 0) ibits(karl) = orf(ibits(karl),rshift(1,(x-1)))
 icycl = (j/31) + karl
 ipos  = MOD(j,31) + 2
 ibits(icycl) = orf(ibits(icycl),itwo(ipos))
 50 GO TO ret, (5,90)
 60 IF (m1(1) /= 0 .AND. m1(2) /= 0) GO TO 1004
 GO TO 190
 
!     SORTED TYPE OF IDS NEED TO COUNT PARENTS IN THE GROUP
 
 70 CONTINUE
 IF (cf .AND. nf /= 0 .AND. isc == ick(2) .AND. cnt == 1) GO TO 31
 IF (cf .AND. nf /= 0 .AND. isc == ick(2)) GO TO 25
 IF (cf) GO TO 190
 IF (cnt == 1) GO TO 80
 cnt = 1
 nf  = 0
 ionf= 0
 ASSIGN 90 TO ret
 isc = isc - 1
 80 CONTINUE
 IF (isc > ick(2)) GO TO 190
 IF (isc-ick(2) < 0) THEN
   GO TO  1004
 ELSE
   GO TO    25
 END IF
 
!     FOUND ID FIND FIELD
 
 90 CALL READ (*180,*180,ifil,ick,6,0,nw)
 IF (k == kn .AND. nf /= 0 .AND. isc == ick(2)) GO TO 31
 GO TO 10
 
!     CHANGE EXTERNAL BCD TO INTERNAL BCD FOR SORT TEST
 
 100 DO  i = 1,2
   itm  = xi(i)
   DO  j = 1,4
     ji   = 5 - j
     ists = isft(ji)
     test = rshift(andf(itm,im3(j)),ists)
     DO  l = 1,37
       IF (test == con(l)) GO TO 120
     END DO
     l = 1
     EXIT
     120 itm = orf(andf(itm,im4(j)),lshift(l,ists +isfim))
     IF (l == 1) EXIT
   END DO
   140 xi(i) = itm
   IF (l == 1) EXIT
 END DO
 160 GO TO exi, (15,16)
 
!     IFP IS DONE BUT VARY IS NOT   MESSAGES FOR ANY LEFT
 
 170 WRITE (nout,1014) ufm,t1(1,k),t1(2,k),ick(2),ick(3)
 CALL READ (*180,*180,ifil,ick,6,0,nw)
 GO TO 170
 
!     END OF IFIL
 
 180 CALL CLOSE (ifil,1)
 ieof  = -1
 ncore = ncore + ibuf
 
 190 cf = .false.
 ionf = nf
 IF (m1(1) == 0 .AND. m1(2) == 0) cf = .true.
 IF (m1(1) /= iend) GO TO 200
 
!     LAST TIME ENTERED MAKE SURE FILE IS USED UP
 
 IF (ieof >= 0) GO TO 170
 200 RETURN
 
!     ERROR MESSAGES
 
 1000 FORMAT (a29,' 3310, CARD TYPE ',2A4,' SORTED',i9,' ID',i9,  &
     ' FIELD',i9,' CHANGED TO ',e16.8)
 1001 WRITE  (nout,1011) ufm
 1011 FORMAT (a23,' 303, NO OPEN CORE IFP')
 GO TO  1111
 1002 WRITE  (nout,1012) sfm
 1012 FORMAT (a25,' 3037, ERROR IN IFPMDC')
 GO TO  1111
 1003 WRITE  (nout,1013) ufm,t1(1,k),t1(2,k),knt,ick(2),ick(3)
 1013 FORMAT (a23,' 0301, FIELD TO VARY IS NOT A REAL NUMBER. CARD ',  &
     2A4,'SORTED',i9,' ID',i9,' FIELD',i9)
 abort = .true.
 GO TO ret, (5,90)
 1004 WRITE  (nout,1014) ufm,t1(1,k),t1(2,k),ick(2),ick(3)
 1014 FORMAT (a23,' 520, CARD TO VARY NOT FOUND. CARD ',2A4,' ID',i9,  &
     ' FIELD',i9)
 GO TO ret, (5,90)
 1111 abort = .true.
 nopen = nopen + ibuf
 ieof  = -1
 GO TO 190
END SUBROUTINE ifpmdc
