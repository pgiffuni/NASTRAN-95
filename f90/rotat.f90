SUBROUTINE rotat (ect2,b1,gplst,x)
     
 
 INTEGER, INTENT(IN OUT)                  :: ect2
 INTEGER, INTENT(IN OUT)                  :: b1
 INTEGER, INTENT(IN OUT)                  :: gplst(1)
 REAL, INTENT(IN)                         :: x(3,1)
 INTEGER :: types(13),esym,elid,gpts(12), bar,offset
 REAL :: magtud(2),normal(3)
 DIMENSION       rec1(146),rec2(17),a(3,3),t(2,2),v(2,3),cross(3),  &
     shear(3),isym(13)
 COMMON /BLANK / skip(23),oes1,scr1,scr2,newoes
 COMMON /xxparm/ skppar(211),icase,flag,DATA,skparm
 EQUIVALENCE     (rec1(3),itype),(rec1(4),isub),(rec1( 5),time),  &
     (rec1(6),eigen),(rec1(10),nwds)
 DATA    types / 6,7,8,9,15,16,17,18,19,62,63,64,83 /
 DATA    isym  / 2HT1,2HTB,2HTP,2HTM,2HQP,2HQM,2HT2,2HQ2,2HQ1,2HM1,  &
     2HM2,2HQ4,2HT3/,    esym/2H  /,    bar /2HBR/
 
 twopi = 8.0*ATAN(1.0)
 irdect= 0
 sum   = 0.0
 CALL OPEN (*110,newoes,gplst(b1),1)
 irec  = 0
 elid  = 0
 10 CALL READ (*110,*110,oes1,rec1,146,1,m)
 IF (isub /= icase) GO TO 13
 IF (flag ==   0.0) GO TO 14
 IF (flag == 1.0 .AND. time == DATA) GO TO 14
 figen = SQRT(ABS(eigen))/twopi
 IF (flag == 2.0 .AND. ABS(figen-DATA) > 1.0E-5) GO TO 14
 13 IF (irec == 0) GO TO 18
 GO TO 100
 
!     CHECK ELEMENT TYPE
 
 14 irec = irec + 2
 DO  it = 1,13
   IF (itype == types(it)) GO TO 20
 END DO
 
!     SKIP SUBCASE
 
 18 CALL fwdrec (*110,oes1)
 GO TO 10
 
!     CHECK ELEMENT TYPE
 
 20 IF (elid /= 0) GO TO 21
 CALL READ (*22,*22,ect2,esym,1,0,n)
 CALL fread (ect2,ngppe,1,0)
 irdect = 1
 offset = 0
 IF (esym == bar) offset = 6
 IF (esym == isym(12) .OR. esym == isym(13)) offset = 1
 IF (esym == isym(it)) GO TO 23
 21 CALL fread (ect2,elid,1,0)
 IF (elid == 0) GO TO 20
 j = 1 + ngppe + offset
 CALL fread (ect2,0,-j,0)
 GO TO 21
 22 CALL bckrec (ect2)
 irdect = 0
 GO TO 18
 
!     PROCESS SUBCASE
 
 23 CALL WRITE (newoes,rec1,146,1)
 nwds = nwds - 1
 25 CALL READ  (*100,*56,oes1,ielmt,1,0,m)
 CALL fread (oes1,rec2,nwds,0)
 30 CALL fread (ect2,elid,1,0)
 IF (elid == 0) GO TO 55
 CALL fread (ect2,0,-1,0)
 CALL fread (ect2,gpts,ngppe,0)
 IF (offset /= 0) CALL fread (ect2,0,-offset,0)
 29 IF (elid == ielmt/10) GO TO 31
 IF (elid > ielmt/10) GO TO 60
 GO TO 30
 31 ig1 = gpts(1)
 ig2 = gpts(2)
 ig1 = IABS(gplst(ig1))
 ig2 = IABS(gplst(ig2))
 ig3 = gpts(3)
 ig3 = IABS(gplst(ig3))
 DO  i = 1,3
   v(1,i) = x(i,ig1) - x(i,ig2)
   v(2,i) = x(i,ig1) - x(i,ig3)
 END DO
 magtud(1) = SQRT(v(1,1)**2 + v(1,2)**2 + v(1,3)**2)
 magtud(2) = SQRT(v(2,1)**2 + v(2,2)**2 + v(2,3)**2)
 DO  i = 1,3
   v(1,i) = v(1,i)/magtud(1)
   v(2,i) = v(2,i)/magtud(2)
   a(1,i) = v(1,i)
 END DO
 a(2,1) = a(1,2)
 a(3,1) = a(1,3)
 a(3,3) =   v(1,1)*v(2,2) - v(2,1)*v(1,2)
 cross(1) = v(1,2)*v(2,3) - v(2,2)*v(1,3)
 cross(2) = v(2,1)*v(1,3) - v(1,1)*v(2,3)
 cross(3) = a(3,3)
 a(2,2) = cross(1)*v(1,3) - v(1,1)*cross(3)
 a(2,3) = v(1,1)*cross(2) - cross(1)*v(1,2)
 a(3,2) = a(2,3)
 iel    = 0
 DO  more = 1,2
   IF (itype == 9 .OR. itype == 16) GO TO 34
   norm   = iel + 2
   ishear = iel + 4
   GO TO 35
   34 norm   = iel + 1
   ishear = iel + 3
   35 t(1,1) = rec2(norm  )
   t(2,2) = rec2(norm+1)
   t(1,2) = rec2(ishear)
   t(2,1) = t(1,2)
   DO  i = 1,3
     sum = 0.0
     DO  j = 1,2
       DO  k = 1,2
         sum = sum + a(i,j)*a(i,k)*t(j,k)
       END DO
     END DO
     normal(i) = sum
   END DO
   shear(1) = a(2,1)*a(1,1)*t(1,1) + a(2,1)*a(1,2)*t(1,2)  &
       + a(2,2)*a(1,2)*t(2,1) + a(2,2)*a(1,2)*t(2,2)
   shear(2) = a(3,1)*a(1,1)*t(1,1) + a(3,1)*a(1,2)*t(1,2)  &
       + a(3,2)*a(1,2)*t(2,1) + a(3,2)*a(1,2)*t(2,2)
   shear(3) = a(3,1)*a(2,1)*t(1,1) + a(3,1)*a(2,2)*t(1,2)  &
       + a(3,2)*a(2,1)*t(2,1) + a(3,2)*a(2,2)*t(2,2)
   DO  i = 1,3
     ishear = ishear + 1
     rec2(norm  ) = normal(i)
     rec2(ishear) = shear(i)
     norm = norm + 1
   END DO
   iel  = iel + 8
   IF (itype == 9 .OR. itype == 16) EXIT
 END DO
 50 CALL WRITE (newoes,ielmt,1,0)
 CALL WRITE (newoes,rec2,nwds,0)
 GO TO 25
 
!     CLOSE RECORD
 
 55 CALL fread (oes1,0,0,1)
 56 CALL WRITE (newoes,0,0,1)
 GO TO 10
 
!     SKIP ELEMENT
 
 60 CALL READ  (*100,*56,oes1,ielmt,1,0,m)
 CALL fread (oes1,rec2,nwds,0)
 GO TO 29
 100 CONTINUE
 110 IF (irdect > 0) CALL bckrec (ect2)
 CALL bckrec (oes1)
 CALL CLOSE  (newoes,1)
 RETURN
END SUBROUTINE rotat
