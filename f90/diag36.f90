SUBROUTINE diag36 (z,buf,gpl,sil,eqexin)
     
!     THIS ROUTINE PRINTS THE INTERNAL-EXTERNAL-SIL NOS. OF THE GRID
!     POINTS AND SCALAR POINTS, AS REQUESTED BY DIAG 36
 
 
 INTEGER, INTENT(OUT)                     :: z(2)
 INTEGER, INTENT(IN OUT)                  :: buf
 INTEGER, INTENT(IN)                      :: gpl
 INTEGER, INTENT(IN)                      :: sil
 INTEGER, INTENT(IN)                      :: eqexin
 INTEGER :: FILE,     nam(2)
 COMMON /system/ ibuf,     l,        dummy(6), nlpp
 COMMON /names / rd,       rdrew,    skip(2),  rew
 DATA            nam /     4HDIAG,   4H34      /
 
 FILE = gpl
 z(1) = gpl
 CALL rdtrl (z(1))
 n1 = z(2)
 n2 = n1 + n1
 n3 = n2 + n1 + 1
 IF (n1 <= 0) GO TO 150
 
 n = 1
 DO  i = 1,2
   CALL OPEN (*150,FILE,z(buf),rdrew)
   CALL fwdrec (*160,FILE)
   CALL READ (*150,*170,FILE,z(n),n1,1,j)
   CALL CLOSE (FILE,rew)
   FILE = sil
   n = n + n1
 END DO
 
!     HERE WE HAVE, IN INTERNAL NUMBER ORDER,
!        Z(   1 THRU N1) = EXTERNAL NOS.
!        Z(N1+1 THRU N2) = SIL NOS.
 
 nlpx = nlpp - 8
 n = nlpx*3
 DO  i = 1,n1,n
   CALL page1
   WRITE  (l,30)
   30   FORMAT (/46X,38HTABLE of internal-EXTERNAL-sil numbers,  &
       //10X,3(6X,30HINTERNAL  EXTERNAL      sil   ), /10X,3(6X,3(10H--------  )))
   im1 = i - 1
   DO  j = 1,nlpx
     j1 = im1 + j
     j2 = j1  + nlpx
     j3 = j2  + nlpx
     IF (j3 <= n1) WRITE (l,40)  &
         j1,z(j1),z(j1+n1), j2,z(j2),z(j2+n1), j3,z(j3),z(j3+n1)
     IF (j3 > n1 .AND. j2 <= n1) WRITE (l,40)  &
         j1,z(j1),z(j1+n1), j2,z(j2),z(j2+n1)
     IF (j2 > n1 .AND. j1 <= n1) WRITE (l,40) j1,z(j1),z(j1+n1)
     40   FORMAT (10X,3(4X,3I10,2X))
   END DO
 END DO
 
 CALL sswtch (20,j)
 IF (j == 0) RETURN
 
 FILE = eqexin
 CALL OPEN (*150,FILE,z(buf),rdrew)
 CALL fwdrec (*160,FILE)
 CALL READ (*150,*170,FILE,z( 1),n2,1,j)
 CALL READ (*150,*170,FILE,z(n3),n2,1,j)
 CALL CLOSE (FILE,rew)
 i = n3 - 1
 j = n2
 k = n3 + n2 - 1
 DO  n = 1,n1
   z(i  ) = z(k  )
   z(i-1) = z(j  )
   z(i-2) = z(j-1)
   i = i - 3
   j = j - 2
   k = k - 2
 END DO
 
!     HERE WE HAVE AN ARRAY OF EXTERNAL-INTERNAL-CODED SIL. PRINT IT OUT
 
 nlpx = nlpx*3
 n    = nlpx*3
 n3   = n3 - 1
 DO  i = 1,n3,n
   CALL page1
   WRITE  (l,80)
   80   FORMAT (/44X,44HTABLE of EXTERNAL-internal-coded sil numbers,  &
       //10X,3(6X,30HEXTERNAL  internal coded sil  ),  &
       /10X,3(5X,3(10H--------- ),1X))
   im1 = i - 1
   DO  j = 1,nlpx,3
     j1 = im1 + j
     j2 = j1  + nlpx
     j3 = j2  + nlpx
     IF (j3 <= n3) WRITE (l,40)  &
         z(j1),z(j1+1),z(j1+2), z(j2),z(j2+1),z(j2+2), z(j3),z(j3+1),z(j3+2)
     IF (j3 > n3 .AND. j2 <= n3) WRITE (l,40)  &
         z(j1),z(j1+1),z(j1+2), z(j2),z(j2+1),z(j2+2)
     IF (j2 > n3 .AND. j1 <= n3) WRITE (l,40) z(j1),z(j1+1),z(j1+2)
   END DO
 END DO
 
 WRITE  (l,110)
 110  FORMAT (//10X,33H*** job terminated by diag 20 ***)
 CALL pexit
 
 150  n = -1
 GO TO 180
 160  n = -2
 GO TO 180
 170  n = -7
 180  CALL mesage (n,FILE,nam)
 RETURN
END SUBROUTINE diag36
