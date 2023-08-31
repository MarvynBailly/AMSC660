VERSION 1.0 CLASS
BEGIN
  MultiUse = -1  'True
END
Attribute VB_Name = "Sheet1"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = True
Sub SimpleForLoop()
    Dim i As Integer
    Dim runs As Integer
    
    runs = 100
    
    For j = 2 To 21
        Dim xDouble As Double
        xDouble = Rnd
        Cells(2, j).Value = xDouble
        For i = 2 To 101
            If xDouble < 1 / 2 Then
                xDouble = 2 * xDouble
            Else
                xDouble = 2 - 2 * xDouble
            End If
            
            Cells(i, j).Value = xDouble
        Next i
    Next j
    For j = 23 To 42
        Dim xSingle As Single
        xSingle = Rnd
        Cells(2, j).Value = xSingle
        For i = 2 To 101
            If xSingle < 1 / 2 Then
                xSingle = 2 * xSingle
            Else
                xSingle = 2 - 2 * xSingle
            End If
            
            Cells(i, j).Value = xSingle
        Next i
    Next j
    
End Sub
