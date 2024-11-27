' gratingAdjustment

Option Explicit

Sub Main()

	' Optimization variables
	Dim gratingGap As Double
	Dim wireDiameter As Double

	' Read-only variables
	Dim waveguideSide As Double
	Dim wireCount As Double

	' Infinite loop indicator
	Dim safetyCounter As Integer

    waveguideSide = GetParameter("waveguideSide")
    gratingGap = GetParameter("gratingGap")
    wireDiameter = GetParameter("wireDiameter")
    wireCount = GetParameter("wireCount")

    ' Check if the grating is too narrow
    safetyCounter = 0
    Do While isGratingTooNarrow(waveguideSide, gratingGap, wireDiameter, wireCount)
        wireCount = wireCount + 2

        safetyCounter = safetyCounter + 1
        If safetyCounter > 100 Then
            Err.Raise vbObjectError + 111, "Main", "Infinite loop detected!"
            Exit Sub
        End If
    Loop

    ' Check if the grating is too wide
    safetyCounter = 0
    Do While isIntersecting(waveguideSide, gratingGap, wireDiameter, wireCount)
        wireCount = wireCount - 2
        If wireCount < 1 Then
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
    If wireCount <> GetParameter("wireCount") Then
        StoreParameter "wireCount", wireCount
        Call RebuildOnParametricChange(True, False)
    End If

End Sub

Private Function GetParameter(ParameterName As String) As Double

	Dim i As Integer
	Dim parameterCount As Double
	parameterCount = GetNumberOfParameters()

	For i = 0 To parameterCount - 1
		If GetParameterName(i) = ParameterName Then
			GetParameter = GetParameterSValue(i)
			Exit For
		End If
	Next i

	If i = parameterCount Then
		Err.Raise vbObjectError + 200, "GetParameter", "Parameter '" & ParameterName & "' not found."
	End If

End Function

' Check if
Private Function isGratingTooNarrow( _
    waveguideSide As Double, _
    gratingGap As Double, _
    wireDiameter As Double, _
    wireCount As Integer _
) As Boolean

    Dim gapAroundGrating As Double
    gapAroundGrating = (waveguideSide - wireCount * gratingGap - wireDiameter) / 2

    If gapAroundGrating > gratingGap Then
        isGratingTooNarrow = True
    Else
        isGratingTooNarrow = False
    End If

End Function

' Check if the grating is being pushed out into the waveguide walls
Private Function isIntersecting( _
    waveguideSide As Double, _
    gratingGap As Double, _
    wireDiameter As Double, _
    wireCount As Integer _
) As Boolean

    Dim gratingEdge As Double
    gratingEdge = (wireCount - 1) / 2 * gratingGap + wireDiameter / 2

    If gratingEdge > waveguideSide / 2 Then
        isIntersecting = True
    Else
        isIntersecting = False
    End If

End Function
