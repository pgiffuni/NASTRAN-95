SUBROUTINE genvec (*,ibuf,filea,nx,ix,ncol,b,bbar,c,cbar,r,ientry)
     
!     GENVEC WILL PICK THE OPTIMUM VALUE OF B AND BBAR FOR A GIVEN
!     MATRIX
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN OUT)                  :: ibuf(1)
 INTEGER, INTENT(IN)                      :: filea(1)
 INTEGER, INTENT(IN)                      :: nx
 INTEGER, INTENT(OUT)                     :: ix(2)
 INTEGER, INTENT(IN)                      :: ncol
 INTEGER, INTENT(OUT)                     :: b
 INTEGER, INTENT(OUT)                     :: bbar
 INTEGER, INTENT(OUT)                     :: c
 INTEGER, INTENT(OUT)                     :: cbar
 INTEGER, INTENT(OUT)                     :: r
 INTEGER, INTENT(IN OUT)                  :: ientry
 INTEGER :: NAME(2)  ,bmax     ,cmax     ,  rsp      ,eol      ,sysbuf   ,  &
       bb       ,cc       ,bbr      ,  &
     ccr       ,rrr      ,bbr1     ,ccr1     ,  &
     bbr2      ,ccr2     ,rr1      ,rr2      ,  &
     p         ,dbname(2),findc    ,namin(2,2)
 DIMENSION  xmb(2)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm       ,uwm      ,uim
 COMMON /names / rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw   ,rsp
 COMMON /zntpkx/ ia(4)     ,ii       ,eol
!     COMMON /DESCRP/ LENGTH    ,MAJOR(1)
 COMMON /ntime / lntime    ,tcons(15)
 COMMON /system/ istv(65)
 COMMON /dcompx/ dum(35)   ,isym
 EQUIVALENCE     (istv( 1),sysbuf)   ,(istv( 2),nout  )  ,  &
     (istv(55),p     )   ,(tcons(8),xmb(1))
 DATA    NAME  / 4HGENV,4HEC  /  ,cmax / 200 /,  &
     namin / 4H rea,1HL   ,4HCOMP,3HLEX  /
 
 
 CALL fname (filea,dbname)
 CALL sswtch (11,l11)
 IF (l11 /= 0) WRITE (nout,6) filea
 6 FORMAT ('O*** DIAG 11 OUTPUT FROM GENVEC (UNSYMMETRIC DECOMP) FOR'  &
     ,       ' FILE',i6 , /9X,1HB,6X,4HBBAR,9X,1HC,6X,4HCBAR,9X,1HR,3X, 4HTIME )
 
 bmax = MIN0(IFIX(1.0E+05/SQRT(FLOAT(ncol)*xmb(p))),ncol)
 ifile= filea(1)
 CALL OPEN (*280,filea(1),ibuf,rdrew)
 i1   = ncol
 i4   = 4*ncol + 2*cmax
 icrq = i4 - nx + sysbuf
 IF (i4 > nx-sysbuf) GO TO 300
 DO  i = 1,i4
   ix(i) = 0
 END DO
 nmax  = 0
 mmax  = 0
 CALL fwdrec (*290,filea(1))
 
!     GENERATE THE ROW AND COLUMN VECTORS
 
 DO  i = 1,ncol
   CALL intpk (*320,filea(1),0,rsp,0)
   CALL zntpki
   in1  = i1 + i
   ix(in1) = ii
   nmax = MAX0(nmax,i-ii+1)
   20 IF (ix(ii) == 0) THEN
     GO TO    30
   ELSE
     GO TO    40
   END IF
   30 ix(ii) = i
   mmax = MAX0(mmax,ii-i+1)
   40 IF (eol == 0.0) THEN
     GO TO    50
   ELSE
     GO TO    60
   END IF
   50 CALL zntpki
   GO TO 20
 END DO
 CALL CLOSE (filea(1),rew)
 i2   = i1 + ncol + 1
 i3   = i2 + 2*ncol
 nmax = MIN0(nmax,bmax)
 mmax = MIN0(mmax,bmax)
 mmax = MAX0(mmax,2)
 
!     SET UP ACTIVE COLUMN BANDWIDTH VECTOR
 
 DO  i = 2,ncol
   j = ncol - i + 1
   icount = 0
   DO  k = 1,j
     l = i2 - k
     IF (ix(l)-i < 0) THEN
       GO TO    70
     ELSE
       GO TO    80
     END IF
     70 icount = icount + 1
     80 l = i2 + (j-k)*2
     ix(l) = MAX0(ix(l),icount)
   END DO
 END DO
 
!     REDUCE LIST TO UNIQUE PAIRS
 
 i = i2
 j = i2 + 2
 k = 2
 110 IF (ix(j) == 0) GO TO 140
 IF (ix(j) -ix(i) == 0) THEN
   GO TO   130
 END IF
 120 i = i + 2
 ix(i  ) = ix(j)
 ix(i+1) = k
 130 j = j + 2
 k = k + 1
 GO TO 110
 140 CONTINUE
 i = i + 2
 ix(i  ) = 0
 ix(i+1) = k
 ilast   = 0
 
!     BEGIN SEARCH FOR B,BBAR
 
 time = 1000000.
 b    = 0
 bbar = 0
 c    = 0
 cbar = 0
 150 bb = ix(i+1)
 IF (bb <= bmax) GO TO 155
 i  = i - 2
 GO TO 150
 155 CONTINUE
 
!    MAKE PRELIMINARY SEARCH
 
 tt1 = 1000000.
 156 CONTINUE
 bb  = ix(i+1)
 cc  = ix(i) + 1
 IF (cc == 1) cc = 0
 bbr = bb
 ccr = cc
 CALL rcore (bb,bbr,cc,ccr,ncol,ientry,nx,rrr)
 rrr = MIN0(rrr,bb+bbr-1,ncol-1)
 IF (rrr < 2) GO TO 157
 CALL timeeq (FLOAT(bb),FLOAT(bbr),FLOAT(cc),FLOAT(ccr),FLOAT(rrr),  &
     ientry,ncol,tt)
 IF (ilast == 0) ilast = i
 IF (l11   == 0) GO TO 1500
 WRITE  (nout,151) bb,bbr,cc,ccr,rrr,tt
 151 FORMAT (5I10,f10.2)
 1500 CONTINUE
 IF (tt > tt1) GO TO 157
 tt1  = tt
 bbr1 = bbr
 ccr1 = ccr
 rr1  = rrr
 157 i    = i - 2
 IF (bb <   3) GO TO 158
 IF (i >= i2+2) GO TO 156
 158 CONTINUE
 i  = i + 2
 IF (tt1 == 1000000.)GO TO 300
 bb = bbr1
 cc = ccr1
 tt1= 1000000.
 
!     SEARCH ON INCREASING BBAR
 
 159 bbr = bb
 incrxx = MAX1(.02*FLOAT(bb),1.)
 160 ccr = findc(bb,bbr,ncol,ix(1),ix(i3))
 CALL rcore (bb,bbr,cc,ccr,ncol,ientry,nx,rrr)
 rrr = MIN0(rrr,bb+bbr-1)
 rrr = MIN0(rrr,ncol-1)
 IF (rrr < 2) GO TO 170
 CALL timeeq (FLOAT(bb),FLOAT(bbr),FLOAT(cc),FLOAT(ccr),FLOAT(rrr),  &
     ientry,ncol,tt)
 IF (l11 == 0) GO TO 1600
 WRITE (nout,151) bb,bbr,cc,ccr,rrr,tt
 1600 CONTINUE
 IF (tt1 == 1000000.) tt1 = tt
 IF (tt  > tt1) GO TO 170
 tt1  = tt
 bbr1 = bbr
 ccr1 = ccr
 rr1  = rrr
 170 CONTINUE
 bbr = bbr + incrxx
 IF (tt > 1.2*tt1) GO TO 180
 IF (ccr  ==     0) GO TO 180
 IF (bbr  <  bmax) GO TO 160
 
!     BEGIN SEARCH ON DECREASING BBAR
 
 180 tt2 = 1000000.
 bbr = bb - incrxx
 190 IF (bbr <= 2 ) GO TO 210
 ccr = findc(bb,bbr,ncol,ix(1),ix(i3))
 CALL rcore (bb,bbr,cc,ccr,ncol,ientry,nx,rrr)
 rrr = MIN0(rrr,bb+bbr-1)
 rrr = MIN0(rrr,ncol-1)
 IF (rrr < 2) GO TO 200
 CALL timeeq (FLOAT(bb),FLOAT(bbr),FLOAT(cc),FLOAT(ccr),FLOAT(rrr),  &
     ientry,ncol,tt)
 IF (l11 == 0) GO TO 195
 WRITE (nout,151) bb,bbr,cc,ccr,rrr,tt
 195 CONTINUE
 IF (tt2 == 1000000.) tt2 = tt
 IF (tt  > tt2) GO TO 200
 tt2  = tt
 bbr2 = bbr
 ccr2 = ccr
 rr2  = rrr
 200 CONTINUE
 bbr  = bbr - incrxx
 IF (tt > 1.20*tt2) GO TO 210
 GO TO 190
 210 CONTINUE
 IF (tt1 >= time) GO TO 220
 time = tt1
 b    = bb
 c    = cc
 bbar = bbr1
 cbar = ccr1
 r    = rr1
 220 IF (tt2 >= time) GO TO 230
 time = tt2
 b    = bb
 c    = cc
 bbar = bbr2
 cbar = ccr2
 r    = rr2
 230 IF (tt1 == 1000000. .AND. tt2 == 1000000.) GO TO 275
 ib = b
 ic = c
 ibbar = bbar
 icbar = cbar
 ir    = r
 ix(1) = c
 ix(2) = r
 CALL page2 (4)
 WRITE  (nout,240) uim,b,bbar,c,cbar,r
 240 FORMAT (a29,' 3028',6X,3HB =,i5,5X,6HBBAR =,i5, /40X,3HC =,i5,5X,  &
     6HCBAR =,i5, /40X,3HR =,i5)
 CALL tfin (FLOAT(b),FLOAT(bbar),FLOAT(c),FLOAT(cbar),FLOAT(r),  &
     ientry,FLOAT(ncol),time)
 ix(1) = time
 CALL page2 (3)
 WRITE  (nout,250) uim,namin(1,ientry),namin(2,ientry),dbname,ncol, ix(1)
 250 FORMAT (a29,' 3027, UNSYMMETRIC ',2A4,' DECOMPOSITION OF DATA ',  &
     'BLOCK ',2A4,6H (n = ,i5,1H), /5X,'TIME ESTIMATE = ',i8, 8H seconds)
 CALL tmtogo (ixy)
 IF (ixy < ix(1)) CALL mesage (-50,ix(1),NAME)
 RETURN
 
!     TRY TO FIND POSSIBLE SOLUTION WITHIN FEASIBLE RANGE BY VARYING  BB
 
 275 i  = i + 2
 IF (i > ilast) GO TO 300
 bb = ix(i+1)
 cc = ix(i) + 1
 IF (bb > bmax) GO TO 300
 GO TO 159
 280 no = -1
 GO TO 310
 290 no = -2
 GO TO 310
 300 no = -8
 ifile = icrq
 310 CALL mesage (no,ifile,NAME)
 RETURN
 
!     NULL COLUMN DISCOVERED
 
 320 WRITE  (nout,325) ufm,i,namin(1,ientry),namin(2,ientry)
 325 FORMAT (a23,' 3097, COLUMN',i7,' IS SINGULAR.  UNSYMMETRIC ',2A4,  &
     'DECOMP ABORTED.')
 RETURN 1
 
END SUBROUTINE genvec
