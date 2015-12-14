SUBROUTINE eqout1 (ia,len1,ns,len2,isil)
     
!     THIS ROUTINE GENERATES OUTPUT ENTRIES FOR CONNECTION TRACE
 
 
 INTEGER, INTENT(IN OUT)                  :: ia(1)
 INTEGER, INTENT(IN OUT)                  :: len1
 INTEGER, INTENT(IN)                      :: ns(1)
 INTEGER, INTENT(IN OUT)                  :: len2
 INTEGER, INTENT(OUT)                     :: isil
 EXTERNAL        lshift,rshift
 INTEGER :: ibits(2),n1(17),n2(14),rshift,outt
 COMMON /cmb002/ junk(8),outt
 COMMON /system/ junk1(8),nlpp,junk2(2),nline
 COMMON /cmb003/ icomb(7,5),conset,iauto,toler,npsub
 COMMON /machin/ mach,ihalf
 DATA    iblank/ 4H      /
 
!     SORT ON PSEUDOSTRUCTURE NUMBER
 
 ifirst = 1
 DO  k = 1,17
   n1(k) = iblank
 END DO
 DO  k = 1,14
   n2(k) = iblank
 END DO
 CALL sort (0,0,4,1,ia(1),len1)
 j = 1
 n1(1) = ia(j+2)
 icode = ia(j+3)
 CALL bitpat (icode,ibits)
 DO  i = 1,2
   n1(i+1) = ibits(i)
 END DO
 13 ips   = rshift(ia(j),ihalf)
 isub  = 2*(ips-1) + 4
 idbas = ia(j) - lshift(ips,ihalf)
 n1(isub  ) = ns(2*idbas-1)
 n1(isub+1) = ns(2*idbas  )
 ia(j) = -ia(j)
 CALL push (ia(j+1),n2(2*ips-1),1,8,1)
 jj = j
 12 IF (jj+4 > len1) GO TO 14
 IF (ia(jj+4) > 0) THEN
   GO TO    50
 ELSE
   GO TO    11
 END IF
 50  IF (rshift(IABS(ia(j)),ihalf) - rshift(ia(jj+4),ihalf) == 0.0) THEN
   GO TO    11
 ELSE
   GO TO    10
 END IF
 11 jj = jj + 4
 GO TO 12
 10 j = jj + 4
 GO TO 13
 
!     WRITE OUTPUT
 
 14 nline = nline + 3
 IF (nline <= nlpp) GO TO 20
 CALL page
 nline = nline + 3
 20 CONTINUE
 j = 3 + 2*npsub
 IF (ifirst == 1) WRITE(outt,1000) n1(1),isil,(n1(k),k=2,j)
 IF (ifirst == 0) WRITE(outt,1003) (n1(k),k=4,j)
 WRITE (outt,1001) (n2(k),k=1,14)
 ifirst = 0
 j = -3
 15 j = j + 4
 IF (j > len1) GO TO 17
 IF (ia(j) > 0) THEN
   GO TO    16
 ELSE
   GO TO    15
 END IF
 16 DO  k = 1,17
   n1(k) = iblank
 END DO
 DO  k = 1,14
   n2(k) = iblank
 END DO
 GO TO 13
 17 WRITE  (outt,1002)
 1000 FORMAT (8X,i6,6X,i6,8X,a4,a2,7(3X,2A4))
 1001 FORMAT (40X,7(3X,2A4) )
 1002 FORMAT (7X,4H  --,27(4H----),4H-    )
 1003 FORMAT (/40X,7(3X,2A4) )
 RETURN
END SUBROUTINE eqout1
