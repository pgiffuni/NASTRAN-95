SUBROUTINE sma3
     
!     THIS ROUTINE, FOR EACH GENERAL ELEMENT, READS THE GENERAL ELEMENT
!     INPUT FILE, GEI, CALLS SMA3A OR SMA3B, DEPENDING UPON WHETHER OR
!     NOT THE ORDERS OF THE K OR Z AND S MATRICES WILL ALLOW THE IN CORE
!     MATRIX ROUTINES (CALLED BY SMA3A) TO BE USED, AND THEN CALLS THE
!     MATRIX ADD ROUTINE TO ADD THE KGGX MATRIX TO THE GENERAL ELEMENT
!     MATRIX.
 
 LOGICAL :: even,onlyge
 INTEGER :: iq(1),eor,outrw,clsrw,clsnrw
 DOUBLE PRECISION :: dq(1)
 DIMENSION        ibuff3(3),NAME(2),mcbid(7),BLOCK(11),iblock(11)
 COMMON /BLANK /  luset,ngenel,noecpt
 COMMON /system/  ksystm(65)
 COMMON /zzzzzz/  q(1)
 COMMON /genely/  ifgei,ifkggx,ifout,ifa,ifb,ifc,ifd,ife,iff,inrw,  &
     outrw,clsrw,clsnrw,eor,neor,mcba(7),mcbb(7),  &
     mcbc(7),mcbd(7),mcbe(7),mcbf(7),mcbkgg(7),  &
     iui,iud,izi,is,izis,istzis,ibuff3,left
 EQUIVALENCE      (ksystm(1),isys),(ksystm(55),iprec),  &
     (iq(1),dq(1),q(1)),(ibuff3(2),m),(ibuff3(3),n),  &
     (mcbid(1),mcbc(1)),(BLOCK(1),iblock(1))
 DATA    NAME  /  4HSMA3,4H    /
 
!     GENERAL INITIALIZATION
 
 ifgei  = 101
 ifkggx = 102
 if201  = 201
 if301  = 301
 if302  = 302
 if303  = 303
 if304  = 304
 if305  = 305
 if306  = 306
 ifout  = if201
 ifa    = if301
 ifb    = if302
 ifc    = if303
 ifd    = if304
 ife    = if305
 iff    = if306
 ifg    = 307
 inrw   = 0
 outrw  = 1
 clsrw  = 1
 clsnrw = 2
 eor    = 1
 neor  = 0
 
!     DETERMINE THE SIZE OF VARIABLE CORE AVAILABLE AND SET IUI TO THE
!     ZEROTH LOCATION OF VARIABLE CORE.
 
 iqmax = korsz (q)
 iui   = 0
 
!     OPEN THE GENERAL ELEMENT INPUT FILE AND SKIP OVER THE HEADER
!     RECORD.
 
 iggei = iqmax - isys + 1
 CALL gopen (ifgei,q(iggei),0)
 iga   = iggei - isys
 
!     DETERMINE IF THE NUMBER OF GENERAL ELEMENTS IS EVEN OR ODD.
 
 even = .true.
 IF ((ngenel/2)*2 /= ngenel) even = .false.
 ipass = 0
 
!     COMPUTE LENGTH OF OPEN CORE
 
 left = iga - 1
 nz   = left
 
!     READ THE TRAILER FOR KGGX TO SEE IF IT EXISTS.
 
 onlyge = .false.
 mcbkgg(1) = ifkggx
 CALL rdtrl (mcbkgg(1))
 IF (mcbkgg(1) < 0) GO TO 12
 ifb = mcbkgg(1)
 DO  i = 1,7
   mcbb(i) = mcbkgg(i)
   mcbkgg(i) = 0
 END DO
 GO TO 14
 12 onlyge = .true.
 
!     INITIALIZATION PRIOR TO LOOP
 
 14 IF (onlyge) GO TO 21
 ifout = if201
 IF (even) ifout = if302
 GO TO 30
 21 ifa = ifout
 IF (even) ifa = if302
 
!     BEGIN MAIN LOOP OF THE PROGRAM
 
 30 ipass = ipass + 1
 
!     READ THE ELEMENT ID, THE LENGTH OF THE UI SET, M, AND THE LENGTH
!     OF THE UD SET, N
 
 CALL READ (*200,*210,ifgei,ibuff3(1),3,neor,idummy)
 needed = 2*(m+n+m**2 + n**2 + 2*m*n)
 itemp1 = 2*(m+n+ m**2) + 3*m
 IF (itemp1 > needed) needed = itemp1
 
!     DETERMINE IF THERE IS ENOUGH CORE STORAGE AVAILABLE TO USE THE IN
!     CORE MATRIX ROUTINES.
 
 IF (needed > left) GO TO 140
 
 
!     **********  IN CORE VERSION  ****************
 
!     USE THE IN CORE MATRIX ROUTINES.  CALL SMA3A.
 
 CALL makmcb (mcba,ifa,0,6,iprec)
 
!     OPEN THE FILE ON WHICH THE CURRENT GENERAL ELEMENT WILL BE OUTPUT.
 
 CALL gopen (ifa,q(iga),1)
 CALL sma3a (mcba)
 
!     STORE THE CORRECT NUMBER OF ROWS IN THE 3RD WORD OF THE MATRIX
!     CONTROL BLOCK AND CLOSE THE FILE WITH REWIND.
 
 mcba(3) = mcba(2)
 CALL wrttrl (mcba)
 CALL CLOSE (ifa,clsrw)
 
!     SUMATION
 
!     JUMP TO 100 ONLY IF THIS IS THE FIRST PASS AND KGGX DOES NOT EXIST
 
 60 IF (ipass == 1 .AND. onlyge) GO TO 100
 CALL makmcb (mcbkgg,ifout,0,6,iprec)
 iblock(1) = 1
 BLOCK (2) = 1.0
 BLOCK (3) = 0.0
 BLOCK (4) = 0.0
 BLOCK (5) = 0.0
 BLOCK (6) = 0.0
 iblock(7) = 1
 BLOCK (8) = 1.0
 BLOCK (9) = 0.0
 BLOCK(10) = 0.0
 BLOCK(11) = 0.0
 
!     CLOSE GEI WITH NO REWIND SO SUBROUTINE ADD CAN HAVE THE BUFFER
 
 CALL CLOSE (ifgei,2)
 
!     CALL SSG2C TO PERFORM SUMMATION - OUTPUT ON IFOUT
 
 CALL ssg2c (ifa,ifb,ifout,0,BLOCK)
 IF (ipass == ngenel) GO TO 160
 CALL rdtrl (mcbkgg)
 
!     RESTORE GEI AFTER SUMATION
 
 CALL gopen (ifgei,q(iggei),2)
 IF (ipass  > 1) GO TO 130
 100 IF (ngenel == 1) GO TO 160
 ifa   = if301
 ifb   = if302
 ifout = if201
 IF (.NOT.even) GO TO 130
 ifb   = if201
 ifout = if302
 
!     SWITCH FILES IFB AND IFOUT FOR NEXT GENEL PROCESSING
 
 130 DO  i = 1,7
   ii = mcbkgg(i)
   mcbkgg(i) = mcbb(i)
   mcbb(i) = ii
 END DO
 ii = ifout
 ifout = ifb
 ifb = ii
 
!     RETURN TO BEGIN LOOP
 
 GO TO 30
 
!     ***********  OUT OF CORE VERSION  *************
 
!     IFOUT MUST CONTAIN THE RESULTS OF THE LAST GENEL PROCESSED
!     SWITCH FILES IFB AND IFOUT FOR OUT OF CORE VERSION
 
 140 IF (ipass == 1 .AND. onlyge .AND. .NOT.even) GO TO 142
 DO  i = 1,7
   ii = mcbkgg(i)
   mcbkgg(i) = mcbb(i)
   mcbb(i) = ii
 END DO
 ii = ifout
 ifout = ifb
 ifb = ii
 
!     THE IN CORE MATRIX ROUTINES CANNOT BE USED.SUBROUTINE SMA3B BUILDS
!     THE ZE IF Z IS INPUT OR THE ZINYS IF K IS INPUT AND IF PRESENT THE
!     SE MATRICES. IF THE SE MATRIX IS PRESENT ISE IS POSITIVE.
!     NOTE - SE(T) IS ON THE SE FILE.
 
 142 CALL sma3b (ise,izk)
 IF (izk == 2) GO TO 145
 
!     FACTOR DECOMPOSES THE ZE MATRIX INTO ITS UPPER AND LOWER
!     TRIANGULAR FACTORS.  TWO SCRATCH FILES ARE NEEDED.
 
 CALL factor (ifa,ife,iff,ifd,ifc,ifg)
 
!     CONVERT IFB INTO THE IDENTITY MATRIX.  (MCBID HAS BEEN SET UP BY
!     SMA3B)
 
 CALL wrttrl (mcbid)
 
!     COMPUTE Z INVERSE
 
 CALL ssg3a (ifa,ife,ifc,ifd,0,0,-1,0)
 145 CONTINUE
 
!     GO TO 150 IF NO SE MATRIX IS PRESENT.
 
 IF (ise < 0) GO TO 150
 
!               T        T  -1
!     COMPUTE -S XK OR -S XZ  AND STORE ON IFF
!               E  E     E  E
 
 CALL ssg2b (ifb,ifd,0,iff,0,iprec,1,ifc)
 
!     TRANSPOSE THE SE FILE ONTO IFA.  HENCE IFA CONTAINS THE -SE MATRIX
 
 CALL tranp1 (ifb,ifa,1,ifc,0,0,0,0,0,0,0)
 
!                       -1
!     COMPUTE K X-S OR Z  X-S AND STORE ON IFE
!              E   E    E    E
 
 CALL ssg2b (ifd,ifa,0,ife,0,iprec,1,ifc)
 
!              T          T  -1
!     COMPUTE S XK XS OR S XZ  XS AND STORE ON IFC
!              E  E  E    E  E   E
 
 CALL ssg2b (ifb,ife,0,ifc,0,iprec,1,ifa)
 
!     SMA3C BUILDS THE FINAL MATRIX OF G (LUSET) SIZE.
 
 mcba(1) = ifa
 150 CALL sma3c (ise,mcba)
 
!     RETURN FILES IFB AND IFOUT TO ORIGIONAL FILES AFTER OUT OF CORE
 
 IF (ipass == 1 .AND. onlyge .AND. .NOT.even) GO TO 60
 DO  i = 1,7
   ii = mcbkgg(i)
   mcbkgg(i) = mcbb(i)
   mcbb(i) = ii
 END DO
 ii = ifout
 ifout = ifb
 ifb = ii
 
!     RETURN TO SUMATION
 
 GO TO 60
 
!     WRAP-UP
 
 160 CALL CLOSE (ifgei, clsrw)
 IF (ifout /= if201) CALL mesage (-30,28,5)
 RETURN
 
!     FATAL ERROR MESSAGES
 
 200 CALL mesage (-2,ifgei,NAME)
 210 CALL mesage (-3,ifgei,NAME)
 RETURN
END SUBROUTINE sma3
