SUBROUTINE stpaic(bloc,dy,nsize,gap,bm,gm,pm,ns,cla,ajjl)
     
 REAL, INTENT(IN)                         :: bloc(1)
 REAL, INTENT(IN)                         :: dy(1)
 INTEGER, INTENT(IN)                      :: nsize(1)
 REAL, INTENT(IN)                         :: gap(1)
 REAL, INTENT(IN)                         :: bm(4,4,ns)
 REAL, INTENT(IN)                         :: gm(4,3,ns)
 REAL, INTENT(IN OUT)                     :: pm(37,ns)
 INTEGER, INTENT(IN)                      :: ns
 REAL, INTENT(IN)                         :: cla(1)
 REAL, INTENT(IN OUT)                     :: ajjl
 COMPLEX :: ch,cdum,ekm
 
 
 
 DIMENSION ch(3,3),cdum(4,4)
 COMMON /stripc/nns,bref,clam,fm,ncirc,nncirc,ekr(1),  &
     dum,       bb(4),beta(4),ekm(4,4)
 COMMON /amgmn / mcb(7),nrow,nd,NE,refc,fmach,rfk,tskj(7),isk,nsk
 COMMON /packx / iti,ito,ii,nn,incr
 
 k = 1
 ii = nrow +1
 nn = nrow
 IF(ekr(1) <= .00001) ekr(1) = 0.0
 nsted=0
 IF(ekr(k) == 0.0) nsted=1
 DO  n=1,ns
   bob=bloc(n)/bref
   ekl=ekr(k)*bob
   const=  cla(n)*dy(n)*clam
   cr = fm
   IF ( ncirc /= 0 ) cr = bb(1)
   ci = 0.
   nopen = 0
   IF(nsize(n) == 3.AND.gap(n) == 0.0) nopen = 1
   tsr= 0.5*gap(n)/bloc(n)
   im=nsize(n)
   IF(im-3 < 0) THEN
     GO TO    31
   ELSE
     GO TO    32
   END IF
   31 jm=2
   j1=2
   GO TO 33
   32 jm=4
   j1=3
   33 CONTINUE
   CALL stpk(ekl,n,nsize(n),nopen,nsted,tsr,pm(1,n),cr,ci,im,j1)
   DO  i=1,im
     DO  j=1,jm
       cdum(i,j)=CMPLX(0.0,0.0)
       DO  m=1,jm
         cdum(i,j) = cdum(i,j) + bm(i,m,n)*ekm(m,j)
       END DO
     END DO
   END DO
   DO  i=1,im
     DO  j=1,j1
       ch(i,j)  =CMPLX(0.0,0.0)
       DO  m=1,jm
         ch(i,j)   = ch(i,j)   + cdum(i,m)*gm(m,j,n)
       END DO
       ch(i,j)   = const * ch(i,j)
     END DO
   END DO
   nn = nn + im
   DO  i=1,im
     CALL pack(ch(1,i),ajjl,mcb)
   END DO
   ii = ii + im
 END DO
 RETURN
END SUBROUTINE stpaic
