SUBROUTINE tiger (ig,list,inv,ii3,norig,kg,jg)
     
!     THIS ROUTINE MAKES ADDITIONS TO THE CONNECTION TABLE IG TO REFLECT
!     THE PRESENCE OF MPC'S AND STORES THE DEPENDENT POINTS IN LIST.
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     NEQ =NUMBER OF MPC EQUATIONS.
!     NEQR=NUMBER OF MPC EQUATIONS COMING FROM RIGID ELEMENTS
 
 
 INTEGER, INTENT(IN OUT)                  :: ig(1)
 INTEGER, INTENT(OUT)                     :: list(1)
 INTEGER, INTENT(IN OUT)                  :: inv(ii3,1)
 INTEGER, INTENT(IN OUT)                  :: ii3
 INTEGER, INTENT(IN OUT)                  :: norig(1)
 INTEGER, INTENT(IN)                      :: kg(1)
 INTEGER, INTENT(IN)                      :: jg(1)
 INTEGER :: scr1,     bunpk,    rdrew,    rd,       rew
 DIMENSION  sub(2),
 COMMON /banda /  ibuf1,    nompc,    nodep
 COMMON /bandb /  dum6b(6), kdim
 COMMON /bandd /  dum(7),   neq,      neqr
 COMMON /bands /  nn,       mm,       dum2s(2), maxgrd,   maxdeg,  &
     dum3s(3), nedge
 COMMON /geomx /  gdum(3),  scr1
 COMMON /system/  ibuf,     nout
 COMMON /names /  rd,       rdrew,    ndum(2),  rew
 COMMON /zzzzzz/  iz(1)
 DATA             sub /     4HTIGE, 4HR   /
 
 IF (neq+neqr == 0) GO TO 170
 kdim4=kdim*4
 CALL OPEN (*200,scr1,iz(ibuf1),rdrew)
 
!     GENERATE NEW CONNECTIONS.
!     TWO PASSES.   FIRST PASS FOR MPC CARDS, AND SECOND FOR RIGID ELEM.
 
 DO  jj=1,2
   IF (jj == 1) nq=neq
   IF (jj == 2) nq=neqr
   IF (nq == 0) CYCLE
   
!     READ MPC EQUATIONS AND RIGID ELEMENT GRIDS
!     AND CONVERT ORIGINAL GRID NOS. TO INTERNAL LABELS.
   
   DO  ii=1,nq
     CALL READ (*210,*210,scr1,nterm,1,0,m)
     kk=1
     j2=2
     IF (jj == 1) GO TO 10
     k=MOD(nterm,1000)
     nterm=nterm/1000
     kk=nterm-k
     j2=nterm
     10 IF (nterm > kdim4) GO TO 70
     CALL READ (*210,*210,scr1,kg,nterm,1,m)
     CALL scat (kg,nterm,inv,ii3,norig)
     
     loop40:  DO  k=1,kk
       igrid=kg(k)
       IF (nodep == +1) list(igrid)=igrid
       
!     IGRID=DEPENDENT GRID POINT IN AN MPC EQUATION.
       
       CALL bunpak(ig,igrid,maxdeg,jg)
       DO  i=1,maxdeg
         l=jg(i)
         IF (l <= 0) CYCLE loop40
         
!     L= A GRID POINT THAT IGRID IS CONNECTED TO BEFORE THE MPC IS APPLI
         
         IF (nterm < 2) CYCLE
         DO  j=j2,nterm
           CALL setig (l,kg(j),ig,norig)
         END DO
       END DO
     END DO loop40
   END DO
 END DO
 GO TO 90
 
 70 WRITE (nout,80)
 80 FORMAT (72H0*** mpc cards NOT processed in bandit due TO insuffici  &
     ent SCRATCH SPACE,//)
 neq =0
 neqr=0
 90 CALL CLOSE (scr1,rew)
 
!     QUIT HERE IF MPC DEPENDENT POINTS ARE NOT TO BE DELETED FROM THE
!     CONNECTION TABLE IG.
 
 IF (nodep /= +1) GO TO 170
 
!     COMPRESS OUT ZEROS FORM LIST
 
 n=0
 DO  i=1,nn
   IF (list(i) == 0) CYCLE
   n=n+1
   list(n)=list(i)
 END DO
 
!     DELETES ALL REFERENCE IN THE CONNECTION TABLE IG TO THOSE POINTS
!     IN LIST
 
 IF (n <= 0) GO TO 170
 mm1=mm-1
 loop160:  DO  ii=1,n
   i=list(ii)
   CALL bunpak (ig,i,mm,jg)
   DO  j=1,mm
     l=jg(j)
     IF (l == 0) CYCLE loop160
     nedge=nedge-1
     k=0
     120 k=k+1
     m=bunpk(ig,l,k)
     IF (m /= i) GO TO 120
     IF (k >= mm) GO TO 140
     DO  np=k,mm1
       is=bunpk(ig,l,np+1)
       CALL bpack (ig,l,np,is)
     END DO
     140 CALL bpack (ig,l,mm1+1,0)
     CALL bpack (ig,i,j,0)
   END DO
 END DO loop160
 170 RETURN
 
!     SCR1 FILE ERROR
 
 200 k=-1
 GO TO 220
 210 k=-2
 220 CALL mesage (k,scr1,sub)
 RETURN
END SUBROUTINE tiger
