INTEGER FUNCTION khrfn5 (word1,i,word2,j)
     
    !     THIS FUNCTION IS SAME AS KHRFN1 EXECPT THAT THE BYTE (WORD1(1) ONL
    !     IS IN REVERSE ORDER.  (THIS FUNCTION IS MAINLY USED BY VAX)
 
 
    INTEGER, INTENT(IN OUT)                  :: word1(1)
    INTEGER, INTENT(IN)                      :: i
    INTEGER, INTENT(IN)                      :: word2(1)
    INTEGER, INTENT(IN OUT)                  :: j
 
    COMMON /system/ dummy(40), ncpw
 
    khrfn5=khrfn1(word1(1),ncpw-i+1,word2(1),j)

    RETURN
END FUNCTION khrfn5
