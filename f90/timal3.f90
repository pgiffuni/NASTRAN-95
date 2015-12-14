SUBROUTINE timal3
     
!     ASSEMBLY LANGUAGE ROUTINE FOR TIMTS3
 
!     NOTES FROM G.CHAN/UNISYS, 10/1989
!     THIS ROUTINE IS NOT USED IN NASTRAN. IF IT IS ACTIVATED, MAKE SURE
!     THE EQUIVALENCE OF A=B=C=D=AC=BC=CC=DC=AD=BD=CD=DD IS REMOVED
 
 REAL :: b(1),  c(1),  d(1)
 DOUBLE PRECISION :: ad(2), bd(1), cd(1), dd(1)
 COMPLEX :: ac(1), bc(1), cc(1), dc(1)
 
 COMMON /BLANK /  n,m
 COMMON /zzzzzz/  a(1)
 
 EQUIVALENCE     (a(1),ac(1),ad(1), b(1),bc(1),bd(1),  &
     c(1),cc(1),cd(1), d(1),dc(1),dd(1))
 
 
 ENTRY tmtrsp
!     ============
!     REAL SINGLE PRECISION - TIGHT LOOP
 
 DO  i = 1,n
   DO  j = 1,m
     d(j) = a(j)*b(j) + c(j)
   END DO
 END DO
 GO TO 999
 
 
 ENTRY tmmrsp
!     ============
!     REAL SINGLE PRECISION - MEDIUM LOOP
 
 DO  i = 1,n
   DO  j = 1,m
     d(j) = a(i)*b(j) + c(j)
   END DO
 END DO
 GO TO 999
 
 
 ENTRY tmlrsp
!     ============
!     REAL SINGLE PRECISION - LOOSE LOOP
 
 DO  i = 1,n
   DO  j = 1,m
     l = i+j-1
     d(j) = a(i)*b(l) + c(j)
   END DO
 END DO
 GO TO 999
 
 
 ENTRY tmtrdp
!     ============
!     REAL DOUBLE PRECISION - TIGHT LOOP
 
 DO  i = 1,n
   DO  j = 1,m
     dd(j) = ad(j)*bd(j) + cd(j)
   END DO
 END DO
 GO TO 999
 
 
 ENTRY tmmrdp
!     ============
!     REAL DOUBLE PRECISION - MEDIUM LOOP
 
 DO  i = 1,n
   DO  j = 1,m
     dd(j) = ad(i)*bd(j) + cd(j)
   END DO
 END DO
 GO TO 999
 
 
 ENTRY tmlrdp
!     ============
!     REAL DOUBLE PRECISION - LOOSE LOOP
 
 DO  i = 1,n
   DO  j = 1,m
     l = i + j - 1
     dd(j) = ad(i)*bd(l) + cd(j)
   END DO
 END DO
 GO TO 999
 
 
 ENTRY tmtcsp
!     ============
!     COMPLEX SINGLE PRECISION - TIGHT LOOP
 
 DO  i = 1,n
   DO  j = 1,m
     dc(j) = ac(j)*bc(j) + cc(j)
   END DO
 END DO
 GO TO 999
 
 
 ENTRY tmmcsp
!     ============
!     COMPLEX SINGLE PRECISION - MEDIUM LOOP
 
 DO  i = 1,n
   DO  j = 1,m
     dc(j) = ac(i)*bc(j) + cc(j)
   END DO
 END DO
 GO TO 999
 
 
 ENTRY tmlcsp
!     ============
!     COMPLEX SINGLE PRECISION - LOOSE LOOP
 
 DO  i = 1,n
   DO  j = 1,m
     l = i + j - 1
     dc(j) = ac(i)*bc(l) + cc(j)
   END DO
 END DO
 GO TO 999
 
 
 ENTRY tmtcdp
!     ============
!     COMPLEX DOUBLE PRECISION - TIGHT LOOP
 
 DO  i = 1,n
   DO  j = 1,m
     
!     D(J) AND D(J+1) CALCULATIONS WERE REVERSED
!     IN ORDER TO COUNTERACT THE ITERATIVE BUILD UP
     
     dd(j+1) = ad(j) * bd(j  ) - ad(j+1) * bd(j+1) + cd(j  )
     dd(j  ) = ad(j) * bd(j+1) + ad(j+1) * bd(j  ) + cd(j+1)
   END DO
 END DO
 GO TO 999
 
 
 ENTRY tmmcdp
!     ============
!     COMPLEX DOUBLE PRECISION - MEDIUM LOOP
 
 DO  i = 1,n
   DO  j = 1,m
     dd(j  ) = ad(i)*bd(j  ) - ad(i+1)*bd(j+1) + cd(j  )
     dd(j+1) = ad(i)*bd(j+1) + ad(i+1)*bd(j  ) + cd(j+1)
   END DO
 END DO
 GO TO 999
 
 
 ENTRY tmlcdp
!     ============
!     COMPLEX DOUBLE PRECISION - LOOSE LOOP
 
 DO  i = 1,n
   DO  j = 1,m
     l = i + j - 1
     dd(j  ) = ad(i)*bd(l  ) - ad(i+1)*bd(l+1) + cd(j  )
     dd(j+1) = ad(i)*bd(l+1) + ad(i+1)*bd(l  ) + cd(j+1)
   END DO
 END DO
 
 999 RETURN
END SUBROUTINE timal3
