# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active source.
foamfoam = GetActiveSource()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1051, 547]

# show data in view
foamfoamDisplay = Show(foamfoam, renderView1)

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')

# trace defaults for the display properties.
foamfoamDisplay.Representation = 'Surface'
foamfoamDisplay.ColorArrayName = ['POINTS', 'p']
foamfoamDisplay.LookupTable = pLUT
foamfoamDisplay.OSPRayScaleArray = 'p'
foamfoamDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
foamfoamDisplay.SelectOrientationVectors = 'U'
foamfoamDisplay.ScaleFactor = 0.031059999205172065
foamfoamDisplay.SelectScaleArray = 'p'
foamfoamDisplay.GlyphType = 'Arrow'
foamfoamDisplay.GlyphTableIndexArray = 'p'
foamfoamDisplay.GaussianRadius = 0.001552999960258603
foamfoamDisplay.SetScaleArray = ['POINTS', 'p']
foamfoamDisplay.ScaleTransferFunction = 'PiecewiseFunction'
foamfoamDisplay.OpacityArray = ['POINTS', 'p']
foamfoamDisplay.OpacityTransferFunction = 'PiecewiseFunction'
foamfoamDisplay.DataAxesGrid = 'GridAxesRepresentation'
foamfoamDisplay.SelectionCellLabelFontFile = ''
foamfoamDisplay.SelectionPointLabelFontFile = ''
foamfoamDisplay.PolarAxes = 'PolarAxesRepresentation'
foamfoamDisplay.ScalarOpacityFunction = pPWF
foamfoamDisplay.ScalarOpacityUnitDistance = 0.013662170746235226

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
foamfoamDisplay.DataAxesGrid.XTitleFontFile = ''
foamfoamDisplay.DataAxesGrid.YTitleFontFile = ''
foamfoamDisplay.DataAxesGrid.ZTitleFontFile = ''
foamfoamDisplay.DataAxesGrid.XLabelFontFile = ''
foamfoamDisplay.DataAxesGrid.YLabelFontFile = ''
foamfoamDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
foamfoamDisplay.PolarAxes.PolarAxisTitleFontFile = ''
foamfoamDisplay.PolarAxes.PolarAxisLabelFontFile = ''
foamfoamDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
foamfoamDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
foamfoamDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Extract Edges'
extractEdges1 = ExtractEdges(Input=foamfoam)

# show data in view
extractEdges1Display = Show(extractEdges1, renderView1)

# trace defaults for the display properties.
extractEdges1Display.Representation = 'Surface'
extractEdges1Display.ColorArrayName = ['POINTS', 'p']
extractEdges1Display.LookupTable = pLUT
extractEdges1Display.OSPRayScaleArray = 'p'
extractEdges1Display.OSPRayScaleFunction = 'PiecewiseFunction'
extractEdges1Display.SelectOrientationVectors = 'U'
extractEdges1Display.ScaleFactor = 0.031059999205172065
extractEdges1Display.SelectScaleArray = 'p'
extractEdges1Display.GlyphType = 'Arrow'
extractEdges1Display.GlyphTableIndexArray = 'p'
extractEdges1Display.GaussianRadius = 0.001552999960258603
extractEdges1Display.SetScaleArray = ['POINTS', 'p']
extractEdges1Display.ScaleTransferFunction = 'PiecewiseFunction'
extractEdges1Display.OpacityArray = ['POINTS', 'p']
extractEdges1Display.OpacityTransferFunction = 'PiecewiseFunction'
extractEdges1Display.DataAxesGrid = 'GridAxesRepresentation'
extractEdges1Display.SelectionCellLabelFontFile = ''
extractEdges1Display.SelectionPointLabelFontFile = ''
extractEdges1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
extractEdges1Display.DataAxesGrid.XTitleFontFile = ''
extractEdges1Display.DataAxesGrid.YTitleFontFile = ''
extractEdges1Display.DataAxesGrid.ZTitleFontFile = ''
extractEdges1Display.DataAxesGrid.XLabelFontFile = ''
extractEdges1Display.DataAxesGrid.YLabelFontFile = ''
extractEdges1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
extractEdges1Display.PolarAxes.PolarAxisTitleFontFile = ''
extractEdges1Display.PolarAxes.PolarAxisLabelFontFile = ''
extractEdges1Display.PolarAxes.LastRadialAxisTextFontFile = ''
extractEdges1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# hide data in view
Hide(foamfoam, renderView1)

# show color bar/color legend
extractEdges1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Cell Size'
cellSize1 = CellSize(Input=extractEdges1)

# show data in view
cellSize1Display = Show(cellSize1, renderView1)

# trace defaults for the display properties.
cellSize1Display.Representation = 'Surface'
cellSize1Display.ColorArrayName = ['POINTS', 'p']
cellSize1Display.LookupTable = pLUT
cellSize1Display.OSPRayScaleArray = 'p'
cellSize1Display.OSPRayScaleFunction = 'PiecewiseFunction'
cellSize1Display.SelectOrientationVectors = 'U'
cellSize1Display.ScaleFactor = 0.031059999205172065
cellSize1Display.SelectScaleArray = 'p'
cellSize1Display.GlyphType = 'Arrow'
cellSize1Display.GlyphTableIndexArray = 'p'
cellSize1Display.GaussianRadius = 0.001552999960258603
cellSize1Display.SetScaleArray = ['POINTS', 'p']
cellSize1Display.ScaleTransferFunction = 'PiecewiseFunction'
cellSize1Display.OpacityArray = ['POINTS', 'p']
cellSize1Display.OpacityTransferFunction = 'PiecewiseFunction'
cellSize1Display.DataAxesGrid = 'GridAxesRepresentation'
cellSize1Display.SelectionCellLabelFontFile = ''
cellSize1Display.SelectionPointLabelFontFile = ''
cellSize1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
cellSize1Display.DataAxesGrid.XTitleFontFile = ''
cellSize1Display.DataAxesGrid.YTitleFontFile = ''
cellSize1Display.DataAxesGrid.ZTitleFontFile = ''
cellSize1Display.DataAxesGrid.XLabelFontFile = ''
cellSize1Display.DataAxesGrid.YLabelFontFile = ''
cellSize1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
cellSize1Display.PolarAxes.PolarAxisTitleFontFile = ''
cellSize1Display.PolarAxes.PolarAxisLabelFontFile = ''
cellSize1Display.PolarAxes.LastRadialAxisTextFontFile = ''
cellSize1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# hide data in view
Hide(extractEdges1, renderView1)

# show color bar/color legend
cellSize1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# save data
SaveData('/home/ckpark/OpenFOAM/ckpark-9/run/pitzDaily_kEpsilon/lengthData.csv', proxy=cellSize1, Precision=8,
    WriteTimeSteps=1,
    FieldAssociation='Cells')

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [0.13469999562948942, 0.0, 0.6080086693737206]
renderView1.CameraFocalPoint = [0.13469999562948942, 0.0, 0.0]
renderView1.CameraParallelScale = 0.1573642232213606

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).