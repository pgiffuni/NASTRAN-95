SUBROUTINE snpdf (sl,cl,tl,sgs,cgs,sgr,cgr,x0,y0,z0,ee,dij,beta,  &
        cv)
     
!     SNPDF CALCULATES THE STEADY PART OF THE INFLUENCE COEFFICIENT
!     MATRIX ELEMENTS
 
 test1 = 0.9999
 test2 = 0.0001*ee
 
!     ***  TEST1 AND TEST2  SERVE AS A MEASURE OF  'NEARNESS'  WITH
!     RESPECT TO THE BOUND-  AND TRAILING VORTICES RESPECTIVELY - SEE
!     TESTS BELOW
!     NOTE THAT THE MACH NUMBER EFFECT IS ACCOUNTED FOR BY STRETCHING
!     THE  X-COORDINATES AND THE  SWEEP ANGLE OF THE BOUND VORTEX LINE
 
 tlb   = tl/beta
 sqtlb = SQRT(1.0+tlb**2)
 slb   = tlb/sqtlb
 clb   = 1.0/sqtlb
 cave  = cv
 clsgs = clb*sgs
 clcgs = clb*cgs
 ex    = ee*tlb
 ey    = ee*cgs
 ez    = ee*sgs
 x0b   = x0/beta
 rix   = x0b+ ex
 riy   = y0 + ey
 riz   = z0 + ez
 rimag = SQRT(rix**2 + riy**2 + riz**2)
 rox   = x0b- ex
 roy   = y0 - ey
 roz   = z0 - ez
 romag = SQRT(rox**2 + roy**2 + roz**2)
 cab   = (rix*slb+ riy*clcgs + riz*clsgs)/rimag
 cbb   = (rox*slb+ roy*clcgs + roz*clsgs)/romag
 cbi   =-rix/rimag
 cao   = rox/romag
 ricab = rimag*cab
 dbx   = rix - ricab*slb
 dby   = riy - ricab*clcgs
 dbz   = riz - ricab*clsgs
 db2   = dbx**2 + dby**2 + dbz**2
 di2   = riy**2 + riz**2
 do2   = roy**2 + roz**2
 acab  = ABS(cab)
 acbb  = ABS(cbb)
 
!     ***  THE FOLLOWING IS A TEST TO SEE IF THE RECEIVING POINT LIES ON
!     OR NEAR THE BOUND VORTEX  --  IF SO, THE CONTRIBUTION OF THE BOUND
!     VORTEX IS SET TO ZERO
 
 IF (acab > test1) GO TO 30
 IF (acbb > test1) GO TO 30
 cacb = (cab-cbb)/db2
 GO TO  60
 30 IF (cab*cbb < 0.0) THEN
   GO TO    40
 ELSE
   GO TO    50
 END IF
 40 cacb = 0.
 GO TO 60
 50 cacb = 0.5*ABS((1./rimag**2)-(1./romag**2))
 60 CONTINUE
 vby = cacb*(dbx*clsgs - dbz*slb)
 vbz = cacb*(dby*slb - dbx*clcgs)
 
!     ***  TEST TO SEE IF THE RECEIVING POINT LIES ON OR NEAR THE
!     INBOARD TRAILING VORTEX  --  IF SO, THE CONTRIBUTION OF THE
!     INBOARD TRAILING VORTEX IS SET TO ZERO
 
 IF (di2 > test2) GO TO 62
 viy  = 0.0
 viz  = 0.0
 GO TO  64
 62 CONTINUE
 onecbi = (1.0-cbi)/di2
 viy =  onecbi*riz
 viz = -onecbi*riy
 64 CONTINUE
 
!     ***  TEST TO SEE IF THE RECEIVING POINT LIES ON OR NEAR THE
!     OUTBOARD TRAILING VORTEX  --  IF SO, THE CONTRIBUTION OF THE
!     OUTBOARD TRAILING VORTEX IS SET TO ZERO
 
 IF (do2 > test2) GO TO 66
 voy  = 0.0
 voz  = 0.0
 GO TO  68
 66 CONTINUE
 caoone = (1.0+cao)/do2
 voy = -caoone*roz
 voz =  caoone*roy
 68 CONTINUE
 vy  = vby + viy + voy
 vz  = vbz + viz + voz
 ww  = vy*sgr - vz*cgr
 dij = ww*cave/25.132741
 RETURN
END SUBROUTINE snpdf
