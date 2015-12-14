SUBROUTINE rod
     
!     ELEMENT TEMPERATURE AND DEFORMATION LOADING FOR THE ROD, CONROD,
!     TUBE
 
 INTEGER :: eltype,eid,gpida,gpidb
 REAL :: arry(3),gpida1(1),gpidb1(1)
 COMMON /condas/ pi,twopi,radeg,degra,s4pisq
 COMMON /zzzzzz/ core(1)
 COMMON /trimex/ eid, gpida, gpidb, iarry(97)
 COMMON /matin / matid,inflag,temp,stress,sinth,costh
 COMMON /matout/ e1,g,nu,rho,alpha,to1,GE,sigmat,sigmac,sigmas, SPACE(10)
 COMMON /ssgwrk/ ti(16),vect(3),force(3),bgpdt(9),vmag,in,l,tbar, delta,xl
 COMMON /ssgett/ eltype,oldel,eorflg,endid,bufflg,itemp,ideft,idefm
 EQUIVALENCE     (iarry(1),arry(1)), (icstma,bgpdt(1)),(icstmb,bgpdt(5)),  &
     (gpida1(1),bgpdt(2)),(gpidb1(1),bgpdt(6))
 
 nept = 5
 IF (eltype == 3) nept = 4
 a = arry(2)
 
!     RECOMPUTE AREA IF ELEMENT IS TUBE
 
 IF (nept == 4) a = pi*(a-arry(3))*arry(3)
 
 DO  i = 1,9
   nept = nept + 1
   bgpdt(i) = arry(nept)
 END DO
 
!     OBTAIN THE MATERIAL DATA
 
 inflag = 1
 matid  = iarry(1)
 temp   = bgpdt(9)
 CALL mat (eid)
 IF (itemp == 0) THEN
   GO TO   250
 END IF
 240 CALL ssgetd (eid,ti,0)
 tbar = ti(1) - to1
 GO TO 260
 250 tbar = 0.0
 260 IF (ideft == 0) THEN
   GO TO   280
 END IF
 270 CALL fedt (eid,delta,idefm)
 GO TO 290
 280 delta = 0.0
 290 DO  i = 1,3
   vect(i) = gpida1(i) - gpidb1(i)
 END DO
 CALL norm (vect(1),xl)
 vmag = e1*a*(delta + alpha*xl*tbar)/xl
 DO  i = 1,3
   vect(i) = -vect(i)*vmag
   force (i) = -vect(i)
 END DO
 IF (icstmb == 0) THEN
   GO TO   340
 END IF
 330 CALL basglb (vect(1),vect(1),gpidb1,icstmb)
 340 in = gpidb - 1
 DO  i = 1,3
   l  = in + i
   core(l) = core(l) + vect(i)
 END DO
 IF (icstma == 0) THEN
   GO TO   380
 END IF
 370 CALL basglb (force(1),force(1),gpida1,icstma)
 380 in = gpida - 1
 DO  i = 1,3
   l  = in + i
   core(l) = core(l) + force(i)
 END DO
 RETURN
END SUBROUTINE rod
