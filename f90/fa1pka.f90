SUBROUTINE fa1pka(a,m1k,m1b,eiv,ncore,n)
     
!     FA1PKA BUILDS THE MATRIX FOR ALLMAT
 
 
 REAL, INTENT(OUT)                        :: a(1)
 REAL, INTENT(IN)                         :: m1k(1)
 REAL, INTENT(IN)                         :: m1b(1)
 REAL, INTENT(OUT)                        :: eiv(1)
 INTEGER, INTENT(IN OUT)                  :: ncore
 INTEGER, INTENT(IN)                      :: n
 INTEGER :: NAME(2)
 
 
 DATA NAME /4HFA1P,4HKA  /
 DATA nheigs,nheige /4HEIGS,4HEIGE/
 
 n2 = n*2
 iz = 0
 imk = n
 imi = n*n*2
 imb = imi + n
 k = 0
 DO  i = 1,n
   DO  j = 1,n
     k = k +1
     a(iz+j) = 0.0
     a(imk+j) = m1k(k)
     a(imb+j) = m1b(k)
     a(imi+j) = 0.0
     IF(i == j) a(imi+j) = 1.0
   END DO
   iz = iz + n2
   imk = imk + n2
   imi = imi + n2
   imb = imb + n2
 END DO
 
!     CALL HSBG AND ATEIG FOR EIVENVALUES
 
 n4=n2*2
 il = 1
 ih = il + n2
 im=ih+n4
 ii=im+n4
 IF(ii     > ncore) CALL mesage(-8,0,NAME)
 CALL sswtch(39,l39)
 IF(l39 /= 0) CALL conmsg(nheigs,1,0)
 CALL hsbg(n2,a,n2,a)
 CALL ateig(n2,a,eiv(ih),eiv(im),eiv(il),n2, a,eiv(ih),eiv(im))
 il = 0
 DO  i=1,n2
   eiv(i+il) = eiv(i+ih-1)
   eiv(i+il+1) = eiv(i + im -1)
   il = il +1
 END DO
 IF(l39 /= 0) CALL conmsg(nheige,1,0)
 RETURN
END SUBROUTINE fa1pka
