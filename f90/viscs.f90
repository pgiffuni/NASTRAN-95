SUBROUTINE viscs
     
!     THIS SUBROUTINE COMPUTES THE 12X12 MATRIX BGG FOR A VISCOUS
!     (DASHPOT) ELEMENT
 
!     SINGLE PRECISION VERSION
 
!     THE ECPT ENTRIES FOR THE VISC ELEMENT ARE
 
!         ECPT
!     ECPT( 1)   ELEMENT ID
!     ECPT( 2)   SIL NUMBER FOR GRID POINT A
!     ECPT( 3)   SIL NUMBER FOR GRID POINT B
!     ECPT( 4)   EXTENSIONAL DAMPING CONSTANT  - C1
!     ECPT( 5)   TORSIONAL DAMPING COEFFICIENT - C2
!     ECPT( 6)   COORD. SYSTEM ID FOR POINT A
!     ECPT( 7)   X1
!     ECPT( 8)   Y1
!     ECPT( 9)   Z1
!     ECPT(10)   COORD. SYSTEM ID FOR POINT B
!     ECPT(11)   X2
!     ECPT(12)   Y2
!     ECPT(13)   Z2
!     ECPT(14)   ELEMENT TEMPERATURE (NOT USED)
 
 
 LOGICAL :: nogo,idbug
 INTEGER :: iecpt(14),elid,estid,dict(7),indx(4),kx(4),kbx(4)
 REAL :: vec(3),d(64),b(144),ta(9),tb(9)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ skp,ioutpt,ksystm(53),iheat
 COMMON /emgest/ ecpt(14)
 COMMON /emgprm/ ixtra,jcore,ncore,dum(12),istif,imass,idamp,  &
     iprec,nogo,heat,icmbar,lcstm,lmat,lhmat
 COMMON /zzzzzz/ xx(1)
 COMMON /emgdic/ idm,ldict,ngrids,elid,estid
 EQUIVALENCE     (ecpt(1),iecpt(1),ielid), (dict(5),dict5),  &
     (indx(1),ia), (indx(2),iab), (indx(3),iba), (indx(4),ib)
 DATA    kx    / 1 ,7 ,73 ,79 /
 DATA    kbx   / 40,46,112,118/
 
!     INITIALIZE EMGOUT PARAMETERS
 
 idbug   = .true.
 ngrids  = 2
 ldict   = 5 + ngrids
 dict(1) = estid
 dict(2) = 1
 dict(3) = 12
 dict(4) = 63
 dict5   = 0.
 ifile   = 3
 ip      = iprec
 
!     NOW COMPUTE THE LENGTH OF THE ROD AND NORMALIZE
 
 fl = 0.
 DO  i = 1,3
   vec(i) = ecpt(i+6) - ecpt(i+10)
   fl = fl + vec(i)**2
 END DO
 fl = SQRT(fl)
 
 IF (fl <= 0) GO TO 7770
 DO  i = 1,3
   vec(i) = vec(i)/fl
 END DO
 
!     SET UP THE N MATRIX
 
 DO  i = 1,3
   DO  j = 1,3
     ix = (i-1)*3 + j
     d(ix) = vec(i)*vec(j)
   END DO
 END DO
 
!     INITIALIZE THE B MATRIX
 
 DO  i = 1,144
   b(i) = 0.
 END DO
 
!     SWAP INDICES A AND B IF NECESSARY SO MATRIX WILL BE ORDERED
!     BY INCREASING SIL VALUE
 
 ipa = 6
 ipb = 10
 IF (iecpt(2) < iecpt(3)) GO TO 60
 ix  = ipa
 ipa = ipb
 ipb = ipa
 
!     CONVERT GRID POINTS TO BASIC COORDINATES IF NECESSARY
 
 60 ia  = 1
 iab = 1
 IF (iecpt(ipa) == 0) GO TO 70
 ia  = 19
 iab = 10
 CALL transs (ecpt(ipa),ta(1))
 CALL gmmats (ta(1), 3,3,1, d(1), 3,3,0, d(10))
 CALL gmmats (d(10), 3,3,0, ta(1),3,3,0, d(19))
 
 70 ib  = 1
 iba = 1
 IF (iecpt(ipb) == 0) GO TO 80
 ib  = 28
 iba = 37
 CALL transs (ecpt(ipb), tb(1))
 CALL gmmats (tb(1),3,3,1, d(1), 3,3,0, d(37))
 CALL gmmats (d(37),3,3,0, tb(1),3,3,0, d(28))
 
 CALL gmmats (d(iab),3,3,0, tb(1), 3,3,0, d(46))
 iab = 46
 
 80 IF (iecpt(ipa) == 0) GO TO 90
 CALL gmmats (d(iba),3,3,0, ta(1),3,3,0, d(55))
 iba = 55
 
!     CALCULATE THE DAMPING MATRIX B
 
!                       ****                    ****
!                       *      /     /      /      *
!                       * C D  /   0 /-C D  /  0   *
!                       *  1 AA/     /  1 AB/      *
!                       *--------------------------*
!                       *  0   /C D  /   0  /-C D  *
!                       *      / 2 AA/      /  2 AB*
!         B    =        *--------------------------*
!                       *-C D  /   0 / C D  /  0   *
!                       *  1 BA/     /  1 BB/      *
!                       *------------/-------------*
!                       *  0   /-C D /   0  / C D  *
!                       *      /  2 BA      /  2 BB*
!                       *      /     /      /      *
!                       ****                    ****
 
 90 c1 = ecpt (4)
 c2 = ecpt (5)
 
 DO  jtj = 1,4
   kb  = kx(jtj)
   kbb = kbx(jtj)
   j   = 0
   i1  = indx(jtj)
   i2  = i1 + 8
   IF (MOD(jtj,2) /= 0) GO TO 100
   c1  = -c1
   c2  = -c2
   
   100 DO  i = i1,i2
     b(kb)  = c1*d(i)
     b(kbb) = c2*d(i)
     IF (MOD(i,3) == 0) j = 9
     kb  = kb  + 1 + j
     kbb = kbb + 1 + j
     j = 0
   END DO
   
 END DO
 
!     OUTPUT THE MATRIX
 
 CALL emgout (b,b,144,1,dict,ifile,ip)
 RETURN
 
!     ERROR EXITS
 
 7770 WRITE  (ioutpt,7775) ufm,ielid
 7775 FORMAT (a23,' 31XX, ILLEGAL GEOMETRY OR CONNECTIONS FOR VISC ',  &
     'ELEMENT',i10)
 nogo = .true.
 RETURN
END SUBROUTINE viscs
