SUBROUTINE ddrmms (buf,eltype,buf4,buf6)
     
 
 REAL, INTENT(OUT)                        :: buf(16)
 INTEGER, INTENT(IN OUT)                  :: eltype
 INTEGER, INTENT(IN)                      :: buf4
 INTEGER, INTENT(IN)                      :: buf6
 EXTERNAL andf
 INTEGER :: andf   , dit   ,elm(4),elt   ,bufa(100)    ,  &
     est    ,ielid ,ieltmp,int1  ,z     ,FILE  , matflg,  &
     matid  ,mpt   ,mtd(4),n     ,nelt  ,nwords,n1mat ,n2mat , tmp(4),wrd(4)
 REAL :: eltemp,stress,sinth ,coth  ,e     ,g     ,nu    ,  &
     rho    ,alpha ,t0    ,gsube ,sigt  ,sigc  ,sigs  ,fint1 , temp   ,cprim
 COMMON /matin /  matid ,matflg,eltemp,stress,sinth ,coth
 COMMON /matout/  e,g   ,nu,rho,alpha ,t0    ,gsube ,sigt   ,sigc , sigs
 COMMON /zzzzzz/  z(1)
 COMMON /system/  isys(61)
 EQUIVALENCE     (int1,fint1) ,(ieltmp,eltemp)
 DATA    int1/1/ ,est/109/ ,mpt/110/ ,dit/111/ ,elm/1,3,10,34  / ,  &
     mtd   / 4,4,4,16/ ,tmp/ 17,16,17,42 / ,wrd/17,16,17,42/
 
 
 DO  i = 1,4
   IF (elm(i) == eltype) GO TO 215
 END DO
 GO TO 230
 215 nelt = i
 CALL OPEN (*240,est,z(buf4),0)
 220 CALL fwdrec (*250,est)
 CALL fread (est,elt,1,0)
 IF (elt /= eltype) GO TO 220
 nwords = wrd(nelt)
 CALL fread (est,bufa,nwords,0)
 CALL CLOSE (est,1)
 n1mat  = buf4 - buf6
 CALL premat (z(buf6),z(buf6),z(buf4),n1mat,n2mat,mpt,dit)
 matflg = 1
 itemp  = tmp(nelt)
 imatid = mtd(nelt)
 ieltmp = bufa(itemp)
 matid  = bufa(imatid)
 ielid  = bufa(1)
 CALL mat (ielid)
 230 CONTINUE
 IF (eltype > 0.0) THEN
   GO TO    10
 ELSE
   GO TO   200
 END IF
 10 IF (eltype > 100) GO TO 200
!              ROD       BEAM      TUBE      SHEAR     TWIST
 GO TO (  21       ,1        ,21       ,1        ,1  &
!              TRIA1     TRBSC     TRPLT     TRMEM     CONROD  &
 ,20       ,20       ,20       ,40       ,21  &
!              ELAS1     ELAS2     ELAS3     ELAS4     QDPLT  &
 ,1        ,1        ,1        ,1        ,20  &
!              QDMEM     TRIA2     QUAD2     QUAD1     DAMP1  &
 ,40       ,20       ,20       ,20       ,1  &
!              DAMP2     DAMP3     DAMP4     VISC      MASS1  &
 ,1        ,1        ,1        ,1        ,1  &
!              MASS2     MASS3     MASS4     CONM1     CONM2  &
 ,1        ,1        ,1        ,1        ,1  &
!              PLOTEL    REACT     QUAD3     BAR       CONE  &
 ,1        ,1        ,1        ,22       ,1  &
!              TRIARG    TRAPRG    TORDRG    TETRA     WEDGE  &
 ,1        ,1        ,1        ,150      ,150  &
!              HEXA1     HEXA2     FLUID2    FLUID3    FLUID4  &
 ,150      ,150      ,1        ,1        ,1  &
!              FLMASS    AXIF2     AXIF3     AXIF4     SLOT3  &
 ,1        ,1        ,1        ,1        ,1  &
!              SLOT4     HBDY      DUM1      DUM2      DUM3  &
 ,1        ,1        ,1        ,1        ,1  &
!              DUM4      DUM5      DUM6      DUM7      DUM8  &
 ,1        ,1        ,1        ,1        ,1  &
!              DUM9      QDMEM1    QDMEM2    QUAD4     IHEX1  &
 ,1        ,40       ,40       ,20       ,1  &
!              IHEX2     IHEX3     QUADTS    TRIATS    TRIAAX  &
 ,1        ,1        ,1        ,1        ,1  &
!              TRAPAX    AERO1     TRIM6     TRPLT1    TRSHL  &
 ,1        ,1        ,1        ,1        ,1  &
!              FHEX1     FHEX2     FTETRA    FWEDGE    IS2D8  &
 ,1        ,1        ,1        ,1        ,40  &
!              ELBOW     FTUBE     TRIA3     -----     -----  &
 ,22       ,1        ,20       ,1        ,1  &
!              -----     -----     -----     -----     -----  &
 ,1        ,1        ,1        ,1        ,1  &
!              -----     -----     -----     -----     -----  &
 ,1        ,1        ,1        ,1        ,1  &
!              -----     -----     -----     -----     -----  &
 ,1        ,1        ,1        ,1        ,1    ), eltype
 
!     ROD  CONROD  TUBE
 
 21 buf(3) = fint1
 buf(5) = fint1
 
!     M. S. IN TENSION OR COMPRESSION
 
 IF (buf(2) >= 0.0) GO TO 300
 IF (sigc   == 0.0) GO TO 301
 buf(3) = (-ABS(sigc)/buf(2))-1.0
 GO TO 301
 300 IF (sigt <= 0.0 .OR. buf(2) == 0.0) GO TO 301
 buf(3) = sigt/buf(2)-1.0
 
!     M. S. IN TORSION
 
 301 IF (buf(4) == 0.0 .OR. sigs <= 0.0) GO TO 200
 buf(3) = sigs/ABS(buf(4))-1.0
 GO TO 200
 
!     BAR  ELBOW
 
 22 buf( 7) = buf(6) + AMAX1(buf(2),buf(3),buf(4),buf(5))
 buf( 8) = buf(6) + AMIN1(buf(2),buf(3),buf(4),buf(5))
 buf( 9) = fint1
 buf(14) = buf(6) + AMAX1(buf(10),buf(11),buf(12),buf(13))
 buf(15) = buf(6) + AMIN1(buf(10),buf(11),buf(12),buf(13))
 buf(16) = fint1
 
!     M. S. IN TENSION
 
 IF (sigt <= 0.0) GO TO 302
 temp = buf(7)
 IF (buf(7) < buf(14)) temp = buf(14)
 IF (temp   <= 0.0) GO TO 302
 buf(9) = sigt/temp-1.0
 
!     M. S. IN COMPRESSION
 
 302 IF (sigc == 0.0) GO TO 200
 temp = buf(8)
 IF (buf(8) > buf(15)) temp = buf(15)
 IF (temp   >= 0.0) GO TO 200
 cprim   =-ABS(sigc)
 buf(16) = cprim/temp - 1.0
 GO TO 200
 
!     TRIA1  TRIA2  TRIA3  QUAD1  QUAD2  QUAD4  TRBSC  TRPLT  QDPLT
 
 20 i = 2
 ASSIGN 30 TO iretrn
 GO TO 100
 30 i = 10
 ASSIGN 200 TO iretrn
 GO TO 100
 
!     TRMEM  QDMEM  QDMEM1  QDMEM2  IS2D8
 
 40 i = 1
 ASSIGN 200 TO iretrn
 GO TO 100
 
!     PRINCIPAL STRESS EQUATIONS FOR 2-DIMENSIONAL ELEMENTS
 
 100 temp     = buf(i+1) - buf(i+2)
 buf(i+7) = SQRT((temp/2.0)**2 + buf(i+3)**2)
 delta    = (buf(i+1) + buf(i+2)) / 2.0
 buf(i+5) = delta + buf(i+7)
 buf(i+6) = delta - buf(i+7)
 
 IF (andf(isys(61),1) > 0.0) THEN
   GO TO   110
 ELSE
   GO TO   120
 END IF
 110 buf(i+7) = SQRT(buf(i+1)**2 + buf(i+2)**2 - buf(i+1)*buf(i+2)  &
     +  3.0*buf(i+3)**2)
 
 120 delta = 2.0*buf(i+3)
 IF (ABS(delta) < 1.0E-15 .AND. ABS(temp) < 1.0E-15) GO TO 121
 buf(i+4) = ATAN2(delta,temp)*28.6478898
 GO TO iretrn, (30,200)
 
 121 buf(i+4) = 0.0
 GO TO iretrn, (30,200)
 
!     TETRA  WEDGE  HEXA1  HEXA2
 
 150 buf(8) = SQRT(buf(2)*(buf(2)-buf(3)-buf(4))*2.0  &
     +  2.0*buf(3)*(buf(3)-buf(4)) + 2.0*buf(4)**2  &
     +  6.0*(buf(5)**2 + buf(6)**2 + buf(7)**2)) / 3.0
 GO TO 200
 
 1 CONTINUE
 200 RETURN
 
!     ERROR PROCESSING FOR DDRMMS
 
 240 n = -1
 FILE = est
 GO TO 260
 250 n = -2
 FILE = est
 260 CALL mesage (n,FILE,nam)
 RETURN
END SUBROUTINE ddrmms
