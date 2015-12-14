INTEGER FUNCTION isft (bf,sft,j)
     
 
 INTEGER, INTENT(IN OUT)                  :: bf
 INTEGER, INTENT(IN OUT)                  :: sft
 INTEGER, INTENT(IN OUT)                  :: j
 EXTERNAL lshift,rshift
 INTEGER :: rshift
 
 IF (j == 4) GO TO 10
 isft = rshift(bf,sft)
 RETURN
 10 isft = lshift(bf,sft)
 RETURN
END FUNCTION isft
