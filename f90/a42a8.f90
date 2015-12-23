SUBROUTINE a42a8 (a,b,c)

    !     MERGES TWO A4 BCD WORDS (A AND B) TO ONE A8 BCD WORD (C)

    IMPLICIT NONE
    REAL, INTENT(IN OUT)     :: a,    b,    c(2)
    REAL                     :: dummy
    INTEGER                  :: ncpw
    CHARACTER (LEN=4)        :: ka,    kb
    CHARACTER (LEN=8)        :: kc,    d

    COMMON /system/ dummy(40), ncpw
 
    WRITE (d,10) a,b
    IF (ncpw < 8) READ (d,10) c(1),c(2)
    IF (ncpw >= 8) READ (d,20) c(1)
10  FORMAT (2A4)
20  FORMAT ( a8)

    RETURN
 
 
    ENTRY a42k8 (a,b,kc)
    !     ======================
 
    !     MERGES TWO A4 BCD WORDS (A AND B) TO ONE A8 CHARACTER WORD (KC)
 
    WRITE (kc,10) a,b
    RETURN
 
 
    ENTRY a42k4 (a,ka,notuse)
    !     ===========================
 
    !     CONVERTS ONE A4 BCD WORD (A) TO ONE A4 CHARACTER WORD (KA)
 
    WRITE  (ka,30) a
30  FORMAT (a4)
    RETURN
 
 
    ENTRY a82k8 (c,kc,notuse)
    !     ===========================
 
    !     CONVERTS ONE A8 BCD WORD (C) TO ONE A4 CHARACTER WORD (KC)
 
    IF (ncpw < 8) WRITE (kc,10) c(1),c(2)
    IF (ncpw >= 8) WRITE (kc,20) c(1)
    RETURN
 
 
    ENTRY k42k8 (ka,kb,kc)
    !     ========================
 
    !     MERGES TWO A4 CHARACTER WORDS (KA AND KB) TO ONE A8 CHARACTER
    !     WORD (KC)
 
    !     NOTE - SOME MACHINES, SUCH AS UNIVAC, HANDLE BCD WORD AND
    !            CHARACTER WORD QUIT DIFFERENTLY
 
    WRITE (kc,10) ka,kb
    RETURN
 
 
    ENTRY k42a8 (ka,kb,c)
    !     =======================
 
    !     MERGES TWO A4 CHARACTER WORDS (KA AND KB) TO ONE A8 BCD WORD (C)
 
    WRITE (d,10) ka,kb
    IF (ncpw < 8) READ (d,10) c(1),c(2)
    IF (ncpw >= 8) READ (d,20) c(1)
    RETURN
 
 
    ENTRY k42a4 (ka,a,notuse)
    !     ===========================
 
    !     CONVERTS ONE A4 CHARACTER WORD (KA) TO ONE A4 BCD WORD (A)
 
    READ (ka,30) a
    RETURN
 
 
    ENTRY k82a8 (kc,c,notuse)
    !     ===========================
 
    !     CONVERTS ONE A8 CHARACTER WORD (KC) TO ONE A8 BCD WORD (C)
 
    IF (ncpw < 8) READ (kc,10) c(1),c(2)
    IF (ncpw >= 8) READ (kc,20) c(1)
    RETURN
 
END SUBROUTINE a42a8
