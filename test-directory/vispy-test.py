import numpy as np
from vispy import app, scene

# Create a canvas to plot on
canvas = scene.SceneCanvas(keys='interactive', bgcolor='white')
view = canvas.central_widget.add_view()

# Create some data to plot
x = np.linspace(-1, 1, 1000)
y = np.sin(x * 2 * np.pi)

# Create a Line visual
line = scene.Line(np.column_stack((x, y)), color='blue', width=2)

# Set the camera
view.camera = scene.PanZoomCamera(aspect=1)

# Show the plot
canvas.show()

# Run the Vispy app
app.run()

