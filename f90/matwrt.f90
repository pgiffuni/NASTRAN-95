SUBROUTINE matwrt (ifile,xname,xitem,lcore)
     
 
 INTEGER, INTENT(IN)                      :: ifile
 REAL, INTENT(IN)                         :: xname(2)
 REAL, INTENT(IN)                         :: xitem
 INTEGER, INTENT(IN)                      :: lcore
 INTEGER :: otpe,sysbuf
 DOUBLE PRECISION :: dcol
 DIMENSION       ia(7),TYPE(10),FORM(18),dcol(1)
 COMMON /zzzzzz/ col(1)
 COMMON /unpakx/ it,k,l,incr
 COMMON /system/ sysbuf,otpe,inx(6),nlpp,inx1(2),line
 COMMON /output/ head1(96),head2(96)
 EQUIVALENCE     (col(1),dcol(1))
 DATA    TYPE  / 4HREAL,4H    ,4HDB  ,4HPREC,4HCOMP,4HLEX ,4HCMP ,  &
     4HD.p.,4HILL ,4HDEFN/
 DATA    FORM  / 4HSQUA,4HRE  ,4HRECT,4HANG ,4HDIAG,4HONAL,4HLOW ,  &
     4HTRI ,4HUPP ,4HTRI ,4HSYME,4HTRIC,4HVECT,4HOR  ,  &
     4HIDEN,4HITY ,4HILL ,4HDEFN/
 DATA    BLANK , su    ,bstr  ,uctu  ,re    ,xit   ,em    ,cont  /  &
     4H    , 4H  su,4HBSTR,4HUCTU,4HRE  ,4H  it,4HEM  ,4HCONT/
 DATA    xinue , dx    / 4HINUE, 4HD   /
 
 
!     TRANSFER MATRIX FORM SOF TO GINO
 
 CALL mtrxi (ifile,xname,xitem,0,itest)
 IF (itest /= 1) RETURN
 ia(1) = ifile
 CALL rdtrl (ia(1))
 
 DO  i = 1,96
   head2(i)  = BLANK
 END DO
 head2( 1) = su
 head2( 2) = bstr
 head2( 3) = uctu
 head2( 4) = re
 head2( 5) = xname(1)
 head2( 6) = xname(2)
 head2( 7) = xit
 head2( 8) = em
 head2( 9) = xitem
 head2(11) = cont
 head2(12) = xinue
 head2(13) = dx
 namea = ifile
 lcol  = lcore - sysbuf
 incr  = 1
 CALL gopen (namea,col(lcol+1),0)
 it = ia(5)
 IF (it <= 0 .OR. it > 4) it = 5
 IF = ia(4)
 IF (IF <= 0 .OR. IF > 8) IF = 9
 ncol = ia(2)
 nrow = ia(3)
 IF (IF == 7) ncol = ia(3)
 CALL page1
 WRITE (otpe,20) xname,xitem,TYPE(2*it-1),TYPE(2*it),ncol,nrow,  &
     FORM(2*IF-1),FORM(2*IF)
 20 FORMAT (1H0,6X,13HSUBSTRUCTURE ,2A4,6H item ,a4,6H is a ,2A4,  &
     1X,i6,10H column x ,i6,5H row ,2A4,8H matrix. )
 IF (it == 5 .OR. IF == 9 .OR. ncol == 0 .OR. nrow == 0) GO TO 320
 IF (IF-8) 30,300,320
 30 IF (IF /= 3 .AND. IF /= 7) GO TO 40
 ncol = 1
 nrow = ia(3)
 40 inull= 0
 it1  = 5
 IF (it == 1 .OR. it == 3) it1 = 9
 ASSIGN 60 TO ihop
 jj = 1
 50 k  = 0
 l  = 0
 CALL unpack (*190,namea,col)
 IF (inull == 1) GO TO 330
 60 nrow = l - k + 1
 SELECT CASE ( IF )
   CASE (    1)
     GO TO 80
   CASE (    2)
     GO TO 80
   CASE (    3)
     GO TO 220
   CASE (    4)
     GO TO 80
   CASE (    5)
     GO TO 80
   CASE (    6)
     GO TO 80
   CASE (    7)
     GO TO 240
 END SELECT
 80 WRITE (otpe,90) jj,k,l
 line = line + 3
 IF (line >= nlpp) CALL page
 90 FORMAT (8H0COLUMN ,i6,5X,6H rows ,i6,6H thru ,i6,5X,50(1H-),/1H )
 IF (it > 2) nrow = 2*nrow
 91 k = 0
 100 j = k + 1
 IF (j > nrow) GO TO 200
 k = j + it1
 IF (k > nrow) k = nrow
 SELECT CASE ( it )
   CASE (    1)
     GO TO 110
   CASE (    2)
     GO TO 130
   CASE (    3)
     GO TO 150
   CASE (    4)
     GO TO 170
 END SELECT
 
!     REAL SINGLE PRECISION
 
 110 WRITE  (otpe,120) (col(i),i=j,k)
 120 FORMAT (1X,1P,10E13.5)
 121 line = line + 1
 IF (line >= nlpp) CALL page
 GO TO 100
 
!     REAL DOUBLE PRECISION
 
 130 WRITE  (otpe,140) (dcol(i),i=j,k)
 140 FORMAT (1P,6D22.14)
 GO TO 121
 
!     COMPLEX SINGLE
 
 150 WRITE  (otpe,160) (col(i),i=j,k)
 160 FORMAT (5(1P,e12.4,1H+,1P,e12.4,1HI))
 GO TO 121
 
!     COMPLEX DOUBLE
 
 170 WRITE  (otpe,180) (dcol(i),i=j,k)
 180 FORMAT (3(1P,d20.12,1H+,1P,d20.12,2HI ))
 GO TO 121
 190 IF (inull == 1) GO TO 200
 ibegn = jj
 inull = 1
 200 jj =  jj + 1
 IF (jj  <= ncol) GO TO 50
 ASSIGN 210 TO ihop
 IF (inull == 1) GO TO 330
 210 CALL CLOSE (namea,1)
 GO TO 270
 220 WRITE (otpe,230)k,l
 line = line + 2
 230 FORMAT (30H0DIAGONAL elements for columns,i6,3H TO,i7,4H are,/1H0)
 GO TO 91
 240 WRITE (otpe,250) k,l
 line = line + 2
 250 FORMAT (25H0ROW elements for columns,i6,4H TO ,i6,4H are ,/1H0 )
 GO TO 91
 270 WRITE  (otpe,280) ia(6)
 280 FORMAT (53H0THE NUMBER of non-zero words in the longest record =, i8 )
 ia7a = ia(7)/100
 ia7c = ia(7) - 100*ia7a
 ia7b = ia7c/10
 ia7c = ia7c - 10*ia7b
 WRITE  (otpe,285) ia7a,ia7b,ia7c
 285 FORMAT (31H0THE density of this matrix is ,i3,1H.,i1,i1, 9H percent.)
 290 RETURN
 
 300 WRITE  (otpe,310)
 310 FORMAT (16H0IDENTITY matrix)
 320 CALL CLOSE (namea,1)
 
!     FUNNY MATRIX -- TABLE PRINT IT
 
 CALL tabprt (namea)
 GO TO 290
 330 ifin = jj - 1
 WRITE (otpe,340) ibegn,ifin
 inull = 0
 line  = line + 2
 IF (line >= nlpp) CALL page
 340 FORMAT (9H0COLUMNS ,i7,6H thru ,i7,10H are null.)
 GO TO ihop, (60,210)
END SUBROUTINE matwrt
