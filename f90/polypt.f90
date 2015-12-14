SUBROUTINE polypt( loctof,stedge,tr, ngridf,fledge,fl,locfos, eps,  &
        npoly,p)
     
!     POLYPT DETERMINES PERIMETER POINTS OF AREA COMMON TO STRUCTURAL
!        TRIANGLE BOUNDED BY TR POINTS AND FLUID ELEMENT BOUNDED BY
!        (3 OR 4) FL POINTS
 
 
 INTEGER, INTENT(IN)                      :: loctof(3)
 INTEGER, INTENT(IN)                      :: stedge(2,3)
 DOUBLE PRECISION, INTENT(IN)             :: tr(3,3)
 INTEGER, INTENT(IN)                      :: ngridf
 INTEGER, INTENT(IN)                      :: fledge(2,4)
 DOUBLE PRECISION, INTENT(IN)             :: fl(3,4)
 INTEGER, INTENT(IN OUT)                  :: locfos(4)
 DOUBLE PRECISION, INTENT(IN OUT)         :: eps(2)
 INTEGER, INTENT(OUT)                     :: npoly
 DOUBLE PRECISION, INTENT(OUT)            :: p(2,7)
 
 DOUBLE PRECISION :: ss(2), p1(2)
 INTEGER :: kedge(2,5), jedge(2,7)
 
 ip= 0
 npoly= 0
 
 DO  i=1,2
   DO  j=1,7
     p(i,j)= 0.d0
   END DO
 END DO
 
 DO  k=1,3
   IF ( loctof(k) < 0)  GO TO 40
 END DO
 
!     STRUCTURAL TRIANGLE IS COMMON AREA WHEN NO STR PTS LIE OUTSIDE
!        FLUID ELEMENT BOUNDRY
 ip= 3
 DO  k=1,3
   DO  i=1,2
     p(i,k)= tr(i,k)
   END DO
 END DO
 GO TO 9000
 
 40 CONTINUE
 
 k= ngridf -1
 DO  i=1,2
   DO  j=1,k
     jedge(i,j)= fledge(i,j)
     jedge(i,j+ngridf)= fledge(i,j)
   END DO
   jedge(i,ngridf)= fledge(i,ngridf)
 END DO
 
 DO  i=1,2
   DO  j=1,2
     kedge(i,j)= stedge(i,j)
     kedge(i,j+3)= stedge(i,j)
   END DO
   kedge(i,3)= stedge(i,3)
 END DO
 
 
 DO  k=1,3
   k1= kedge(1,k)
   k2= kedge(2,k)
   DO  j=1,ngridf
     j1= jedge(1,j)
     j2= jedge(2,j)
     CALL ptintr( tr(1,k1),tr(1,k2), fl(1,j1),fl(1,j2), ss, inter, eps)
     IF (inter == 1)  GO TO 200
   END DO
 END DO
 
! - - AREAS ARE DISJOINT
 GO TO 9000
 
 
 200 jlast= j
 jj1= j
 jj2= j +ngridf -1
 klast= k
 kk1= k +1
 kk2= k+2
 
 IF (loctof(k1) == 1)  GO TO 1800
!     1ST TRI POINT IS OUTSIDE FLUID BOUNDRY
 p1(2)= ss(2)
 p1(1)= ss(1)
 ap1= (p1(1)-tr(1,k1))**2 +(p1(2)-tr(2,k1))**2
 jp1= jlast
 jj1= jlast+1
 
 DO  j=jj1,jj2
   j1= jedge(1,j)
   j2= jedge(2,j)
   CALL ptintr( tr(1,k1),tr(1,k2), fl(1,j1),fl(1,j2), ss, inter, eps)
   IF (inter == 1)  GO TO 400
 END DO
 
 ip= ip+1
 p(1,ip)= p1(1)
 p(2,ip)= p1(2)
 GO TO 1000
 
 400 ap2= (ss(1)-tr(1,k1))**2 + (ss(2)-tr(2,k1))**2
 IF (ap1 < ap2)  GO TO 500
 
 p(1,ip+1)= ss(1)
 p(2,ip+1)= ss(2)
 p(1,ip+2)= p1(1)
 p(2,ip+2)= p1(2)
 ip= ip+2
 jlast= jp1
 GO TO 600
 
 500 p(1,ip+1)= p1(1)
 p(2,ip+1)= p1(2)
 p(1,ip+2)= ss(1)
 p(2,ip+2)= ss(2)
 ip= ip+2
 jlast= j
 
 600 CONTINUE
 IF ( jlast > ngridf)  jlast= jlast -ngridf
 jj1= jlast
 jj2= jj1 +ngridf -1
 j2= jedge(2,jlast)
 GO TO 2000
 
!     SEARCH ALONG LAST STRUCTURAL TRIANGLE EDGE FOR NEXT PTINTR
 
 1000 IF ( loctof(k2) < 0)  GO TO 1100
 IF (tr(1,k2) == p(1,1)  .AND. tr(2,k2) == p(2,1))  GO TO 9000
 ip= ip+1
 p(1,ip)= tr(1,k2)
 p(2,ip)= tr(2,k2)
 klast= klast +1
 IF ( klast == kk2)   GO TO 9000
 k2= kedge(2,klast)
 GO TO 1000
 
 1100 CONTINUE
 jj1= jlast
 IF (jj1 > jj2)  GO TO 9000
 DO  j= jj1,jj2
   j1= jedge(1,j)
   j2= jedge(2,j)
   CALL ptintr( p(1,ip),tr(1,k2), fl(1,j1),fl(1,j2), ss, inter, eps)
   IF ( inter == 1)  GO TO 1200
 END DO
 
 GO TO 9000
 
 1200 IF (ss(1) == p(1,1)  .AND.  ss(2) == p(2,1))  GO TO 9000
 ip= ip +1
 p(1,ip)= ss(1)
 p(2,ip)= ss(2)
 jlast= j
 GO TO 2000
 
 1800 p(1,ip+1)= tr(1,k1)
 p(2,ip+1)= tr(2,k1)
 p(1,ip+2)= ss(1)
 p(2,ip+2)= ss(2)
 ip= ip+2
 
!     SEARCH ALONG LAST FLUID EDGE FOR NEXT PTINTR
 
 2000 IF ( locfos(j2) < 0)  GO TO 2100
 IF (fl(1,j2) == p(1,1)  .AND.  fl(2,j2) == p(2,1))  GO TO 9000
 ip= ip+1
 p(1,ip)= fl(1,j2)
 p(2,ip)= fl(2,j2)
 jlast= jlast +1
 IF ( jlast > jj2)   GO TO 9000
 j2= jedge(2,jlast)
 GO TO 2000
 
 2100 CONTINUE
 kk1= klast
 IF (kk1 > kk2)  GO TO 9000
 DO  k=kk1,kk2
   k1= kedge(1,k)
   k2= kedge(2,k)
   CALL ptintr( p(1,ip),fl(1,j2), tr(1,k1),tr(1,k2), ss, inter, eps)
   IF ( inter == 1)  GO TO 2200
 END DO
 
 GO TO 9000
 
 2200 IF (ss(1) == p(1,1)  .AND.  ss(2) == p(2,1))  GO TO 9000
 ip= ip +1
 p(1,ip)= ss(1)
 p(2,ip)= ss(2)
 klast= k
 GO TO 1000
 
 
 9000 CONTINUE
 npoly= ip
 RETURN
END SUBROUTINE polypt
