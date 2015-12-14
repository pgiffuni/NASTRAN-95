SUBROUTINE relabl (ns,nodes,ig,ic,ideg,idis,iw,NEW,icc,ild,iaj,  &
        jg,idim)
     
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     GENERATE A RELABELING SCHEME STARTING WITH NS NODES FOR WHICH
!     LABELS HAVE BEEN STORED IN ARRAY NODES.
!     SET UP ILD AND NEW.
!     ILD(OLD) = NEW
!     NEW(NEW) = OLD, THE INVERSE OF ILD
!     IAJ IS DIMENSIONED TO IDIM
 
 
 INTEGER, INTENT(IN)                      :: ns
 INTEGER, INTENT(IN)                      :: nodes(1)
 INTEGER, INTENT(IN OUT)                  :: ig(1)
 INTEGER, INTENT(IN)                      :: ic(1)
 INTEGER, INTENT(IN)                      :: ideg(1)
 INTEGER, INTENT(OUT)                     :: idis(1)
 INTEGER, INTENT(OUT)                     :: iw(1)
 INTEGER, INTENT(OUT)                     :: NEW(1)
 INTEGER, INTENT(IN)                      :: icc(1)
 INTEGER, INTENT(OUT)                     :: ild(1)
 INTEGER, INTENT(OUT)                     :: iaj(1)
 INTEGER, INTENT(IN)                      :: jg(1)
 INTEGER, INTENT(IN OUT)                  :: idim
 INTEGER :: x
 
 COMMON /bands / nn,       dums(3),  maxgrd
 COMMON /system/ ibuf,     nout
 
 i   = nodes(1)
 icn = ic(i)
 nt  = icc(icn) - 1
 DO  i = 1,nn
   IF (ic(i)-icn == 0) THEN
     GO TO    80
   ELSE
     GO TO    90
   END IF
   80   idis(i) = 0
 END DO
 DO  j = 1,ns
   jj = nodes(j)
   idis(jj) =-1
   jt = j + nt
   NEW(jt) = jj
   ild(jj) = jt
 END DO
 ki = nt
 ko = ns + nt
 ll = ko
 l  = 1
 j  = ko
 nnc= icc(icn+1) - 1
 110  ki = ki + 1
 IF (ki-ll == 0) THEN
   GO TO   120
 ELSE
   GO TO   130
 END IF
 120  l  = l  + 1
 ll = ko + 1
 130  ii = NEW(ki)
 n  = ideg(ii)
 IF (n == 0) THEN
   GO TO   270
 END IF
 140  ij = 0
 CALL bunpak (ig,ii,n,jg)
 DO  i = 1,n
   ia = jg(i)
   IF (idis(ia) == 0) THEN
     GO TO   150
   ELSE
     GO TO   170
   END IF
   150  ij = ij + 1
   IF (ij <= idim) GO TO 160
   
!     DIMENSION EXCEEDED.  STOP JOB.
   
   ngrid = -2
   RETURN
   
   160  idis(ia) = l
   ko       = ko + 1
   iaj(ij)  = ia
   iw(ij)   = ideg(ia)
 END DO
 IF (ij-1 < 0) THEN
   GO TO   260
 ELSE IF (ij-1 == 0) THEN
   GO TO   180
 ELSE
   GO TO   190
 END IF
 180  j  = ko
 iz = iaj(1)
 NEW(ko) = iz
 ild(iz) = ko
 GO TO 260
 190  x = 0
 DO  i = 2,ij
   IF (iw(i)-iw(i-1) < 0) THEN
     GO TO   210
   ELSE
     GO TO   230
   END IF
   210  CONTINUE
   x = iw(i)
   iw(i  ) = iw(i-1)
   iw(i-1) = x
   x = iaj(i)
   iaj(i  ) = iaj(i-1)
   iaj(i-1) = x
 END DO
 IF (x > 0.0) THEN
   GO TO   190
 END IF
 240  DO  i = 1,ij
   j  = j + 1
   iz = iaj(i)
   NEW(j ) = iz
   ild(iz) = j
 END DO
 260  IF (ko-nnc < 0) THEN
   GO TO   110
 END IF
 270  CONTINUE
 
!     REVERSE SEQUENCE FOR THIS COMPONENT (ICN).
 
!     ICC IS AN ARRAY USED FOR IDENTIFYING COMPONENTS IN THE NEW ARRAY.
!     ICC(N1) CONTAINS THE INDEX FOR THE NEW ARRAY AT WHICH COMPONENT
!         N1 STARTS.
 
 n1 = icc(icn) - 1
 n2 = nn - icc(icn+1) + 1
 IF (n2 > nn) n2 = 0
 
!     REVERSE THE NODAL CM SEQUENCE, OMITTING THE FIRST N1 AND THE LAST
!     N2 POINTS.
 
!     NEW(N1) = OLD LABEL FOR NODE NOW LABELLED N1.
!     ILD(N1) = NEW LABEL FOR NODE ORIGINALLY LABELED N1.
!     N1      = NUMBER OF POINTS AT BEGINNING OF SEQUENCE TO OMIT FROM
!               REVERSAL.
!     N2      = NUMBER OF POINTS AT END OF SEQUENCE TO OMIT FROM
!               REVERSAL.
!     NN      = NUMBER OF NODES.
!     J       = NUMBER OF INTERCHANGES TO MAKE.
 
 j  = (nn-n1-n2)/2
 IF (j <= 0) RETURN
 ll = nn - n2 + 1
 
!     MAKE INTERCHANGES IN NEW ARRAY.
 
 DO  i = 1,j
   l = ll - i
   k = NEW(l)
   m = n1 + i
   NEW(l) = NEW(m)
   NEW(m) = k
 END DO
 
!     CORRECT ILD, THE INVERSE OF NEW.
 
 l = 1  + n1
 m = nn - n2
 DO  i = l,m
   k = NEW(i)
   ild(k) = i
 END DO
 
 RETURN
END SUBROUTINE relabl
