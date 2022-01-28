INTEGER FUNCTION iTabRow(i)
   INTEGER :: i
   !
   iTabRow=0
   If (i.gt. 0 .and. i.le. 2) Then
      iTabRow=1
   Else If (i.gt. 2 .and. i.le.10) Then
      iTabRow=2
   Else If (i.gt.10 .and. i.le.18) Then
      iTabRow=3
   Else If (i.gt.18 .and. i.le.36) Then
      iTabRow=3
   Else If (i.gt.36 .and. i.le.54) Then
      iTabRow=3
   Else If (i.gt.54 .and. i.le.86) Then
      iTabRow=3
   Else If (i.gt.86) Then
      iTabRow=3
   End If
   !
   Return
End function
