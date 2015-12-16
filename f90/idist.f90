FUNCTION idist (ns,ml,maxlev,ig,ic,ideg,idis,iw,icc,jg)
     
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     THIS FUNCTION HAS AS ITS VALUE THE MAXIMUM DISTANCE OF ANY NODE
!     IN COMPONENT IC(NS) FROM THE NODE NS.
!     THE DISTANCE OF EACH NODE IN THIS COMPONENT IS STORED IN THE ARRAY
!     IDIS.
!     THE MAXIMUM NUMBER OF NODES AT THE SAME DISTANCE FROM NS IS
!     STORED IN ML.
 
!     INTEGER          BUNPK
 
 INTEGER, INTENT(IN)                      :: ns
 INTEGER, INTENT(OUT)                     :: ml
 INTEGER, INTENT(IN)                      :: maxlev
 INTEGER, INTENT(IN OUT)                  :: ig(1)
 INTEGER, INTENT(IN)                      :: ic(1)
 INTEGER, INTENT(IN)                      :: ideg(1)
 INTEGER, INTENT(OUT)                     :: idis(1)
 INTEGER, INTENT(OUT)                     :: iw(1)
 INTEGER, INTENT(IN)                      :: icc(1)
 INTEGER, INTENT(IN)                      :: jg(1)
 
 COMMON /bands /  nn
 
 icn = ic(ns)
 nnc = icc(icn+1) - icc(icn)
 DO  i = 1,nn
   IF (ic(i)-ic(ns) == 0) THEN
     GO TO    40
   ELSE
     GO TO    50
   END IF
   40 idis(i) = 0
   50 CONTINUE
 END DO
 ll = 1
 l  = 0
 ki = 0
 ko = 1
 ml = 0
 iw(1) = ns
 idis(ns) = -1
 130 ki = ki + 1
 IF (ki-ll == 0) THEN
   GO TO   132
 ELSE
   GO TO   135
 END IF
 132 l  = l + 1
 ll = ko + 1
 k  = ko - ki + 1
 IF (k-ml > 0) THEN
   GO TO   133
 ELSE
   GO TO   135
 END IF
 133 ml = k
 IF (ml-maxlev > 0) THEN
   GO TO   220
 END IF
 135 ii = iw(ki)
 n  = ideg(ii)
 IF (n == 0) THEN
   GO TO   215
 END IF
 140 CALL bunpak (ig,ii,n,jg)
 DO  i = 1,n
   ia = jg(i)
   IF (idis(ia) == 0) THEN
     GO TO   150
   ELSE
     GO TO   200
   END IF
   150 idis(ia) = l
   ko = ko + 1
   iw(ko) = ia
   200 CONTINUE
 END DO
 IF (ko-nnc < 0) THEN
   GO TO   130
 END IF
 205 idist = l
 idis(ns) = 0
 k = ko - ll + 1
 IF (k-ml > 0) THEN
   GO TO   207
 ELSE
   GO TO   206
 END IF
 207 ml = k
 206 CONTINUE
 RETURN
 
 215 l = 0
 GO TO 205
 220 idist = 1
 
 RETURN
END FUNCTION idist
