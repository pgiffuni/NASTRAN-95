SUBROUTINE stqme2( ntype )
     
!     PHASE TWO STRESS DATA RECOVERY TRIANGULAR MEMBRANE
 
!     NTYPE = 1 TRI-MEMBRANE
!     NTYPE = 2 QUAD-MEMBRANE
 
!     PH1OUT CONTAINS THE FOLLOWING
!     *** NTYPE = 1 ***
!     ELEMENT ID
!     3 SILS
!     1 DUMMY
!     T SUB 0
!     S SUB T 3X1
!     3 S ARRAYS EACH 3X3
 
!     *** NTYPE = 2 ***
!     ELEMENT ID
!     4 SILS
!     T SUB 0
!     S SUB T 3X1
!     4 S ARRAYS EACH 3X3
 
 
 INTEGER, INTENT(IN)                      :: ntype
 REAL :: frlast(2)
 INTEGER :: eject    ,ishd(7)  ,istyp(2) ,typ(3)
 DIMENSION nsil(4), si(36), ph1out(45), st(3)
 COMMON /zzzzzz/ z(1)
 COMMON /sdr2x4/ dummy(35),ivec,ivecn,ldtemp,deform
 COMMON /sdr2x7/ est(100),stres(100),forvec(25)
 COMMON /sdr2x8/ stress(3),vec(3),tem,temp,npoint,delta,nsize, cvec(3),cstr(4)
 COMMON /sdr2x9/ nchk,isub,ild,frtmei(2),twotop,fnchk
 COMMON   /system/  ibfsz    ,nout     ,idm(9)   ,line
 
 EQUIVALENCE (ph1out(1),est(1))  &
     ,(nsil(1),ph1out(2)) ,(tsub0,ph1out(6))  &
     ,(st(1),ph1out(7)) ,(si(1),ph1out(10))  &
     ,(ftemp,ldtemp) , (ishd(1),lsub) , (ishd(2),lld) , (ishd(6),frlast(1) )
 
 DATA lsub,lld,frlast / 2*-1 , -1.0E30, -1.0E30 /
 DATA typ / 2HTR, 2HQD, 3HMEM /
!     ******************************************************************
!     ZERO OUT THE STRESS VECTOR
 DO  i = 1,3
   stress(i) = 0.0E0
   cstr(i+1) = 0.0E0
 END DO
 
!                          I=NSIZE
!        STRESS VECTOR =(SUMMATION (S )(U )) - (S )(LDTEMP - T SUB 0)
!                          I=1       I   I       T
 nsize = ntype + 2
 DO  i = 1,nsize
   
!     POINTER TO DISPLACEMENT VECTOR IN VARIABLE CORE
   
   npoint = ivec + nsil(i) - 1
   
   CALL smmats (si(9*i-8),3,3,0, z(npoint),3,1,0, vec,cvec)
   DO  j = 1,3
     cstr(j+1) = cstr(j+1) + cvec(j)
     stress(j) = stress(j) + vec(j)
   END DO
   
 END DO
 
 stres(1) = ph1out(1)
 stres(2) = stress(1)
 stres(3) = stress(2)
 stres(4) = stress(3)
 
!     ADD IN TEMPERATURE EFFECTS
 
 IF( ldtemp == (-1) ) GO TO 200
 tem = ftemp - t_sub_0
 DO  i = 2,4
   stres(i) = stres(i) - st(i-1) * tem
 END DO
!     STRESS VECTOR COMPLETE AND CONTAINS SIGMA X ,  SIGMA Y ,  SIGMA XY
 
!     ******************************************************************
 
!     PRINCIPAL STRESSES AND ANGLE OF ACTION PHI
 200 temp = stres(2) - stres(3)
 stres(8) = SQRT( (temp/2.0E0)**2 + stres(4)**2 )
 delta = (stres(2) + stres(3))/2.0E0
 stres(6) = delta + stres(8)
 stres(7) = delta - stres(8)
 delta = 2.0E0 * stres(4)
 IF( ABS(delta) < 1.0E-15 .AND. ABS(temp) < 1.0E-15)GO TO 101
 stres(5) = ATAN2( delta,temp ) * 28.6478898 e00
 RETURN
 101 stres(5) = 0.0E0
 IF (nchk <= 0) GO TO 250
 
!  . CHECK PRECISION...
 
 cstr(1) = ph1out(1)
 k = 0
 
 CALL sdrchk (stres(2),cstr(2),3,k)
 IF (k == 0) GO TO 250
 
!  . LIMITS EXCEEDED...
 istyp(1) = typ(ntype)
 istyp(2) = typ(3)
 j = 0
 IF (lsub == isub .AND. frlast(1) == frtmei(1) .AND.  &
     lld == ild  .AND. frlast(2) == frtmei(2) ) GO TO 220
 lsub = isub
 lld = ild
 frlast(1) = frtmei(1)
 frlast(2) = frtmei(2)
 j = 2
 CALL page1
 205 CALL sd2rhd (ishd,j)
 line = line + 1
 WRITE(nout,210)
 210 FORMAT (7X,4HTYPE,5X,3HEID,5X,2HSX,5X,2HSY,4X,3HSXY)
 GO TO 230
 220 IF (eject(2) /= 0) GO TO 205
 230 WRITE(nout,240) istyp,cstr
 240 FORMAT (1H0,6X,a2,a3,i7,3F7.1)
 
 250 CONTINUE
 RETURN
END SUBROUTINE stqme2
