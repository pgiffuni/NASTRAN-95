SUBROUTINE upcase (byte,n)
     
!     THIS ROUTINE CHANGES ALL LOWER CASE CHARACTERS INTO UPPER CASE.
!     IT ALSO CONVERTS BCD INPUT CODE TO EBCDIC FOR IBM MACHINE
 
 
 CHARACTER (LEN=1), INTENT(OUT)           :: byte(1)
 INTEGER, INTENT(IN)                      :: n
 LOGICAL :: flag
 INTEGER :: tab(20),   ffflag
 CHARACTER (LEN=1) :: bk1,      la,       lz,       il,  &
     ic,       ip,       lc(256)
 CHARACTER (LEN=56) :: kc(5)
 COMMON /machin/ machx
 COMMON /upcasx/ flag,     id,       ia,       iz
 COMMON /xechox/ ffflag
 EQUIVALENCE     (kc(1),lc(1))
 
!                     TAB = UPPER CASE 'A' TO LOWER CASE 'a' SPAN
 
 DATA            tab / +32, -64, +32, +3968, +32, +32, +32, +32 ,  &
     +32, +32, +32, +32,   +32, +32, +32, +32 , +32, +32, +32, +32     /
 DATA            bk1,      la,       lz,       il,       ic     /  &
     ' ',      'A',      'Z',      '(',      ','    /
 DATA            ip /      '%'       /
 
!     TAB IS DECIMAL VALUE BETWEEN UPPER CASE 'A' AND LOWER CASE 'a'
!     TAB IS POSITIVE IF LOWER CASE 'a' COMES AFTER UPPER CASE 'A' IN
!     MACHINE ASCII CHARACTER SET; OTHERWISE TAB IS NEGATIVE.
 
!     THE FOLLOWING KC TABLE MUST BE PUNCHED IN EBCDIC CODE (FOR IBM
!     ONLY)                          =======    ===========
 
 DATA            kc /  &
     '                                                        ',  &
     '                   .)(+ +          $*)  -/         ,(%  ',  &
     '           =''''=  ABCDEFGHI       JKLMNOPQR        STUVWX',  &
     'YZ                       ABCDEFGHI       JKLMNOPQR      ',  &
     '  STUVWXYZ      0123456789      WRITTEN BY G.CHAN/UNISYS'/
 
 IF (machx == 2) GO TO 30
 IF (flag) GO TO 10
 flag =.true.
 id = tab(machx)
 ia = ICHAR(la) + id
 iz = ICHAR(lz) + id
 
 10   DO  i = 1,n
   IF (byte(i) == bk1) CYCLE
   j = ICHAR(byte(i))
   IF (j < ia .OR. j > iz) CYCLE
   byte(i) = CHAR(j-id)
 END DO
 RETURN
 
!     IBM MACHINE ONLY, WHICH USES EBCDIC CODE
 
 30   DO  i = 1,n
   j = ICHAR(byte(i))
   byte(i) = lc(j+1)
 END DO
 
!     THE % SIGN MAY BE CHANGED TO ( IN BCD-EBCDIC CONVERSION,
!     CHANGE IT BACK TO %
 
 IF (ffflag /= 1234 .OR. n < 5) RETURN
 DO  i = 5,n
   IF (byte(i) == il .AND. byte(i+1) == il .AND. (byte(i-1) == ic  &
       .OR. byte(i-1) == bk1)) byte(i) = ip
 END DO
 RETURN
END SUBROUTINE upcase
