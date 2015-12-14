SUBROUTINE gust2(fol,wj,acpt,x0,v,cstm,qhjl)
     
!     GUST2 MAKE  WJ(W) MATRIX FOR GUST
 
 
 INTEGER, INTENT(IN)                      :: fol
 INTEGER, INTENT(IN)                      :: wj
 INTEGER, INTENT(IN)                      :: acpt
 REAL, INTENT(IN OUT)                     :: x0
 REAL, INTENT(IN OUT)                     :: v
 INTEGER, INTENT(IN OUT)                  :: cstm
 INTEGER, INTENT(IN)                      :: qhjl
 INTEGER :: buf1,FILE
 INTEGER :: sysbuf,iz(1),trl(7),acdr(13),nam(2)
 
 COMMON /condas/ pi,twopi
 COMMON /system/ sysbuf
 COMMON /zzzzzz/ z(1)
 COMMON /zblpkx/ a(4),irn
 
 EQUIVALENCE (z(1),iz(1))
 
 DATA  nam /4HGUST,1H2 /
 DATA nhnju,nhacj /4HNJU ,4HACJ /
 
 icore = korsz(iz) - sysbuf-2
 buf1 = icore+1
 
!     READ IN FREQUENCYS AND CONVERT TO OMEGA
 
 FILE = fol
 CALL OPEN(*999,fol,z(buf1),0)
 CALL fread(fol,z,-2,0)
 CALL READ(*998,*10,fol,z,icore,0,nfreq)
 GO TO 997
 10 DO  i=1,nfreq
   z(i) = z(i) * twopi
 END DO
 CALL CLOSE(fol,1)
 
!     SPACE FOR COLUMN OF W - 2 * J  LONG  1 J FOR A  1 J FOR COEF.
 
 FILE = qhjl
 trl(1) =  qhjl
 CALL rdtrl(trl)
 IF(trl(1) < 0) GO TO 999
 nj = trl(3)
 jap= nfreq
 jcp = jap + nj
 iacpt = jcp + nj + 1
 IF(iacpt > icore) GO TO 997
 DO  i=1,nj
   z(jap+i) = 0.0
 END DO
 
!     SET UP WJ
 
 trl(1) = wj
 trl(2) = 0
 trl(3) = nj
 trl(4) = 2
 trl(5) = 3
 trl(6) = 0
 trl(7) = 0
 
!     READ ACPT RECORDS BY METHOD AND FILL IN THE TWO COLUMNS
!     A =  COS G (CG) FOR DLB  1 FOR Z BODIES  0 FOR ALL ELSE
!     COEF =   XM  FOR PANELS AND BODIES
 
 CALL gopen(acpt,z(buf1),0)
 nju = 0
 FILE = acpt
 40 CALL READ(*100,*100,acpt,meth,1,0,nwr)
 SELECT CASE ( meth )
   CASE (    1)
     GO TO 50
   CASE (    2)
     GO TO 60
   CASE (    3)
     GO TO 90
   CASE (    4)
     GO TO 90
   CASE (    5)
     GO TO 90
 END SELECT
 
!     DOUBLET LATTICE WITHOUT BODIES
 
 50 CALL READ(*998,*995,acpt,acdr,4,0,nwr)
 np = acdr(1)
 nstrip = acdr(2)
 njg = acdr(3)
 nr = 2*np + 5*nstrip + 2*njg
 IF(iacpt+nr > icore) GO TO 997
 CALL READ(*998,*995,acpt,z(iacpt),nr,1,nwr)
 ixic  =  iacpt + 2*np + 5*nstrip - 1
 idelx =  ixic + njg
 icg   =  iacpt + 2*np + 4*nstrip
 k = 0
 ks= 0
 nbxr = iz(iacpt)
 DO  i = 1,njg
   z(jap+nju+i)=z(icg+ks)
   z(jcp+nju+i)=z(ixic+i) + .5* z(idelx+i)
   IF(i == njg) CYCLE
   IF(i == iz(iacpt+np+k)) k=k+1
   IF(i /= nbxr) CYCLE
   ks = ks+1
   nbxr = nbxr + iz(iacpt+k)
 END DO
 nju = nju+njg
 GO TO 40
 
!     DOUBLET LATTICE WITH BODIES
 
 60 CALL READ(*998,*995,acpt,acdr,13,0,nwr)
 njg = acdr(1)
 np  = acdr(3)
 nb = acdr(4)
 ntp = acdr(5)
 nto = acdr(10)
 ntzs= acdr(11)
 ntys = acdr(12)
 nstrip = acdr(13)
 ic = iacpt
 ib = ic + np
 ib1= ib + 2*np
 ibs= ib1+ 2*nb
 nr = 3*np + 3*nb
 CALL READ(*998,*995,acpt,z(iacpt),nr,0,nwr)
 nbei = 0
 nbes = 0
 DO  i=1,nb
   nbei= nbei+ iz(ib1+i-1)
   nbes= nbes+ iz(ibs+i-1)
 END DO
 icg = ib+ np
 ix  = icg + nstrip  -1
 ixs1= ix  + 4*ntp + 2*nbei + nbes
 ixs2= ixs1+ nbes
 nr = 11*nb + 4*nstrip
 CALL READ(*998,*995,acpt,z(icg),-nr,0,nwr)
 nr =  nstrip + 4*ntp + 2*nbei + 3* nbes
 IF(icg+nr > icore) GO TO 997
 CALL READ(*998,*995,acpt,z(icg),nr,1,nwr)
 IF(ntp == 0) GO TO 65
 k= 0
 ks=0
 nbxr = iz(ic)
 DO  i=1,ntp
   z(jap+nju+i)  =  z(icg+ks)
   z(jcp+nju+i)  =  z(ix+i)
   IF(i == ntp) CYCLE
   IF(i == iz(ib+k)) k=k+1
   IF(i /= nbxr) CYCLE
   ks = ks + 1
   nbxr = nbxr +  iz(ic+k)
 END DO
 65 nju = nju + nto
 IF(ntzs == 0) GO TO 80
 DO  i=1,ntzs
   z(jap+nju+i) = 1.0
   z(jcp+nju+i) =  .5 * (z(ixs1+i) + z(ixs2+i))
 END DO
 80 nju = nju + ntzs + ntys
 GO TO 40
 
!     MACH BOX  STRIP  PISTON  THEORIES
 
 90 CALL READ(*998,*995,acpt,njg,1,1,nwr)
 nju= nju + njg
 GO TO 40
 100 CALL CLOSE(acpt,1)
 CALL bug(nhnju ,100,nju,1)
 CALL bug(nhacj ,100,z(jap+1),2*nj)
 IF(nju /= nj) GO TO 996
 
!     BUILD WJ LOOP OVER ALL FREQUENCIES WITH AN INNER LOOP ON NJ
 
 CALL gopen(wj,z(buf1),1)
 DO  i=1,nfreq
   freq = z(i)
   CALL bldpk(3,3,wj,0,0)
   DO  j=1,nj
     am = z(jap+j)
     IF( am == 0.0 ) CYCLE
     irn = j
     temp   =   freq *((z(jcp+j)-x0)/v)
     a(1) = COS(temp)*am
     a(2) = -SIN(temp)*am
     CALL zblpki
   END DO
   CALL bldpkn(wj,0,trl)
 END DO
 CALL CLOSE(wj,1)
 CALL wrttrl(trl)
 CALL dmpfil(-wj,z,icore)
 1000 RETURN
 
!     ERROR MESSAGES
 
 995 CALL mesage(-3,FILE,nam)
 996 CALL mesage(-7,0,nam)
 997 CALL mesage(-8,0,nam)
 998 CALL mesage(-2,FILE,nam)
 999 CALL mesage(-1,FILE,nam)
 GO TO 1000
END SUBROUTINE gust2
