' gratingAdjustment

Option Explicit

Sub Main()

	' Optimization variables
	Dim GratingGap As Double
	Dim WireDiameter As Double

	' Read-only variables
	Dim WaveguideSide As Double
	Dim WireCount As Double

	' Infinite loop indicator
	Dim safetyCounter As Integer

    WaveguideSide = GetParameter("WaveguideSide")
    GratingGap = GetParameter("GratingGap")
    WireDiameter = GetParameter("WireDiameter")
    WireCount = GetParameter("WireCount")

    ' Check if the grating is too narrow
    safetyCounter = 0
    Do While isGratingTooNarrow(WaveguideSide, GratingGap, WireDiameter, WireCount)
        WireCount = WireCount + 2

        safetyCounter = safetyCounter + 1
        If safetyCounter > 100 Then
            Err.Raise vbObjectError + 111, "Main", "Infinite loop detected!"
            Exit Sub
        End If
    Loop

    ' Check if the grating is too wide
    safetyCounter = 0
    Do While isIntersecting(WaveguideSide, GratingGap, WireDiameter, WireCount)
        WireCount = WireCount - 2
        If WireCount < 1 Then
            Err.Raise vbObjectError + 100, "Main", "No feasible grating geometry found."
            Exit Sub
        End If

        safetyCounter = safetyCounter + 1
        If safetyCounter > 100 Then
            Err.Raise vbObjectError + 111, "Main", "Infinite loop detected!"
            Exit Sub
        End If
    Loop

    ' Update the model parameters and rebuild structure
    StoreParameter "WireCount", WireCount
    Call RebuildOnParametricChange(True, False)

End Sub

Private Function GetParameter(ParameterName As String) As Double

	Dim i As Integer
	Dim ParameterCount As Double
	ParameterCount = GetNumberOfParameters()

	For i = 0 To ParameterCount - 1
		If GetParameterName(i) = ParameterName Then
			GetParameter = GetParameterSValue(i)
			Exit For
		End If
	Next i

	If i = ParameterCount Then
		Err.Raise vbObjectError + 200, "GetParameter", "Parameter '" & ParameterName & "' not found."
	End If

End Function

' Check if 
Private Function isGratingTooNarrow( _
    WaveguideSide As Double, _
    GratingGap As Double, _
    WireDiameter As Double, _
    WireCount As Integer _
) As Boolean

    Dim GapAroundGrating As Double
    GapAroundGrating = (WaveguideSide - WireCount * GratingGap - WireDiameter) / 2

    If GapAroundGrating > GratingGap Then
        isGratingTooNarrow = True
    Else
        isGratingTooNarrow = False
    End If

End Function

' Check if the grating is being pushed out into the waveguide walls
Private Function isIntersecting( _
    WaveguideSide As Double, _
    GratingGap As Double, _
    WireDiameter As Double, _
    WireCount As Integer _
) As Boolean

    Dim GratingEdge As Double
    GratingEdge = (WireCount - 1) / 2 * GratingGap + WireDiameter / 2

    If GratingEdge > WaveguideSide / 2 Then
        isIntersecting = True
    Else
        isIntersecting = False
    End If

End Function
