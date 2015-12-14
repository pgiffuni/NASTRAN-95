SUBROUTINE sptchk
     
!     THIS ROUTINE IS CALLED ONLY BY BANDIT TO CHECK THE PRESENCE OF ANY
!     UNDEFINED SPOINT. RESET NGRID AND RETURN FOR ONE MORE COMPUTATION
!     IF THAT IS THE CASE
 
 INTEGER :: geom1,    geom2,    rd,       rdrew,    rew,  &
     z,        spoint(2),NAME(2),  kg(200)
 COMMON /system/ ibuf,     nout
 COMMON /banda / ibuf1,    dum6(6),  npt(2)
 COMMON /bandb / dum3(3),  ngrid,    dum4(4),  irept
 COMMON /bandd / ndd(9)
 COMMON /bands / skip(4),  maxgrd
 COMMON /geomx / geom1,    geom2
 COMMON /names / rd,       rdrew,    dum2(2),  rew
 COMMON /gpta1 / NE,       last,     incr,     ke(1)
 COMMON /zzzzzz/ z(1)
 DATA            spoint,   NAME /    5551,49,  4HBMIS, 4HS      /
 
!     LIST ALL SPOINTS IN Z(1) THRU Z(NS)
 
 IF (irept == 3) GO TO 160
 ns=1
 CALL preloc (*160,z(ibuf1),geom2)
 CALL locate (*40,z(ibuf1),spoint,k)
 30   CALL READ (*180,*40,geom2,z(ns),1,0,k)
 ns=ns+1
 GO TO 30
 40   ns=ns-1
 CALL REWIND (geom2)
 
!     CHECK THE PRESENCE OF ELAST, DAMP AND MASS CARDS (ELEMENT TYPES
!     201 THRU 1301).  THEY MAY SPECIFY SCALAR POINTS WITHOUT USING
!     SPOINT CARDS.
 
 nss=ns
 DO  ielem=26,350,incr
   CALL locate (*100,z(ibuf1),ke(ielem+3),j)
   nwds =ke(ielem+5)
   ngpt1=ke(ielem+12)
   ngpts=ke(ielem+9)+ngpt1-1
   50   CALL READ (*180,*100,geom2,kg(1),nwds,0,j)
   DO  i=ngpt1,ngpts
     IF (ns == 0) GO TO 70
     CALL bisloc (*70,kg(i),z(1),1,ns,k)
     CYCLE
     70   nss=nss+1
     IF (nss >= ibuf1) EXIT
     z(nss)=kg(i)
   END DO
   GO TO 50
 END DO
 110  CALL CLOSE (geom2,rew)
 k=nss-ns-1
 IF (k < 0) THEN
   GO TO   160
 ELSE IF (k == 0) THEN
   GO TO   140
 END IF
 
!     SOME SCALAR POINTS ARE USED, BUT NOT SPECIFIED BY SPOINT CARDS.
!     SORT THEM, AND THROW OUT DUPLICATES
 
 120  ns1=ns+1
 CALL sort (0,0,1,1,z(ns1),nss-ns)
 k  =nss
 nss=ns1
 j  =ns+2
 DO  i=j,k
   IF (z(i) == z(i-1)) CYCLE
   nss=nss+1
   z(nss)=z(i)
 END DO
 
!     RE-COMPUTE THE TOTAL NO. OF GRID POINTS, NGRID, AND RETURN FOR
!     ONE MORE BANDIT COMPUTATION
 
 140  npt(2)=nss-ns
 ngrid =npt(1)+npt(2)
 DO  i=1,9
   ndd(i)=0
 END DO
 irept =2
 RETURN
 
 160  WRITE (nout,170) maxgrd
 170  FORMAT (120H1*** user fatal error 2007,  this structure model uses  &
     more grid points than the total no. of grid cards in bulk DATA (= ,i6,1H),/)
 ngrid=0
 GO TO 190
 
 180  CALL mesage (-3,geom2,NAME)
 190  RETURN
END SUBROUTINE sptchk
