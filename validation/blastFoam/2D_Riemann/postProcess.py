# trace generated using paraview version 5.11.2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'OpenFOAMReader'
a2D_Riemannfoam = OpenFOAMReader(registrationName='2D_Riemann.foam', FileName='./2D_Riemann.foam')
a2D_Riemannfoam.SkipZeroTime = 0
a2D_Riemannfoam.CaseType = 'Decomposed Case'
a2D_Riemannfoam.MeshRegions = ['internalMesh']
a2D_Riemannfoam.CellArrays = ['p', 'rho']



# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
a2D_RiemannfoamDisplay = Show(a2D_Riemannfoam, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
a2D_RiemannfoamDisplay.Representation = 'Surface'

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# update the view to ensure updated data information
renderView1.Update()
# reset view to fit data

# set scalar coloring
rhoLUT = GetColorTransferFunction('rho')
ColorBy(a2D_RiemannfoamDisplay, ('POINTS', 'rho'))
a2D_RiemannfoamDisplay.RescaleTransferFunctionToDataRange(True, False)
a2D_RiemannfoamDisplay.SetScalarBarVisibility(renderView1, True)

# get 2D transfer function for 'rho'
rhoTF2D = GetTransferFunction2D('rho')

# get color transfer function/color map for 'rho'
rhoLUT.TransferFunction2D = rhoTF2D
rhoLUT.RGBPoints = [0.5313000082969666, 0.239484, 0.00545035, 0.614821, 0.6048384043326974, 0.239484, 0.00545035, 0.614821, 0.6783768003684282, 0.220593, 0.0617459, 0.863547, 0.7519157824141681, 0.17509, 0.278988, 0.97794, 0.825454178449899, 0.143526, 0.576069, 0.998553, 0.8989925744856297, 0.166456, 0.871883, 0.96594, 0.9725309705213606, 0.376202, 0.993555, 0.981833, 1.0460693665570915, 0.681996, 0.991297, 0.999239, 1.1196080966185273, 0.954172, 0.952734, 0.94374, 1.1931467446385624, 0.999735, 0.99301, 0.662896, 1.266685140674293, 0.979399, 0.991466, 0.357973, 1.340223536710024, 0.968771, 0.854967, 0.162659, 1.413761932745755, 0.999245, 0.556697, 0.144323, 1.4873009147914946, 0.973959, 0.26223, 0.177946, 1.5608393108272254, 0.852358, 0.0526707, 0.222974, 1.6343777068629564, 0.593889, 0.00912724, 0.238855, 1.703320026397705, 0.593889, 0.00912724, 0.238855]
rhoLUT.ColorSpace = 'RGB'
rhoLUT.NanColor = [1.0, 0.0, 0.0]
rhoLUT.NumberOfTableValues = 21
rhoLUT.ScalarRangeInitialized = 1.0

# # get opacity transfer function/opacity map for 'rho'
rhoPWF = GetOpacityTransferFunction('rho')
rhoPWF.Points = [0.5313000082969666, 0.0, 0.5, 0.0, 1.703320026397705, 1.0, 0.5, 0.0]
rhoPWF.ScalarRangeInitialized = 1

# get color legend/bar for rhoLUT in view renderView1
rhoLUTColorBar = GetScalarBar(rhoLUT, renderView1)
rhoLUTColorBar.AutoOrient = 0
# rhoLUTColorBar.WindowLocation = 'Upper Right Corner'
rhoLUTColorBar.Orientation = 'Vertical'
rhoLUTColorBar.Title = 'Density [kg/m^3]'
rhoLUTColorBar.TitleFontFamily = 'Arial'
rhoLUTColorBar.LabelFontFamily = 'Arial'
rhoLUTColorBar.ScalarBarThickness = 25
rhoLUTColorBar.ScalarBarLength = 0.9
rhoLUTColorBar.Position = [0.8, 0.05]
rhoLUTColorBar.ScalarBarLength = 0.900000000000001




# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=a2D_Riemannfoam)
contour1.PointMergeMethod = 'Uniform Binning'
contour1.ContourBy = ['POINTS', 'rho']
contour1.Isosurfaces = [0.5313000082969666, 0.5871104853493827, 0.6429209624017987, 0.6987314394542149, 0.754541916506631, 0.8103523935590471, 0.8661628706114632, 0.9219733476638794, 0.9777838247162955, 1.0335943017687117, 1.0894047788211276, 1.145215255873544, 1.20102573292596, 1.256836209978376, 1.3126466870307922, 1.3684571640832082, 1.4242676411356245, 1.4800781181880405, 1.5358885952404566, 1.5916990722928728, 1.647509549345289, 1.703320026397705]


# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour1Display.Representation = 'Wireframe'
contour1Display.ColorArrayName = ['POINTS', 'None']
contour1Display.PointSize = 10.0
# contour1Display.SetScalarBarVisibility(renderView1, False)

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(600, 500)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

#change interaction mode for render view
renderView1.InteractionMode = '2D'

renderView1.ResetActiveCameraToNegativeZ()

# current camera placement for renderView1
# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.63, 0.5, 3.8460652149512318]
renderView1.CameraFocalPoint = [0.63, 0.5, 0.5]
renderView1.CameraParallelScale = 0.55

# Properties modified on rhoLUTColorBar


animationScene1.GoToLast()

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).

# save screenshot
SaveScreenshot('./2D_Riemann_blastFoam.png', renderView1,
ImageResolution=[1200, 1000])
