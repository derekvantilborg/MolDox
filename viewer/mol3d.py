"""
Author: Derek van Tilborg
Data: dec 20 2021

Wrapper for py3Dmol and 3Dmol. Works with layers.

Example:
    model = Viewer()
    model.add_layer('1AZ8_clean.pdb')
    model.add_surface(color='white', opacity=0.75)
    model.cartoon()
    model.add_layer('1AZ8_lig.mol2')
    model.sticks(radius=0.2, colorscheme='orangeCarbon')
    model.add_outline()

    model.show()
"""

import py3Dmol


class Viewer:
    def __init__(self):
        self.view = py3Dmol.view()
        self.view.removeAllModels()
        self.layers = {}

    def add_outline(self, color='black', width=0.1):
        self.view.setViewStyle({'style': 'outline', color: 'black', 'width': width})

    def add_surface(self, opacity=1, color='spectrum', colorscheme=None):

        settings = {'opacity': opacity, 'color': color}
        if colorscheme is not None:
            settings['colorscheme'] = colorscheme

        self.view.addSurface(py3Dmol.VDW, settings)

    def add_layer(self, path, layername=None):

        layername = self.__layername(layername)

        # Load and add the layer
        self.view.addModel(open(path, 'r').read(), format=path.split('.')[-1])
        self.layers[layername] = self.view.getModel()

    def add_layer_from_object(self, obj, format='mol',  layername=None):

        layername = self.__layername(layername)

        # Load and add the layer
        self.view.addModel(obj, format=format)
        self.layers[layername] = self.view.getModel()

    def cartoon(self, layername=None, arrows=False, tubes=False, style='rectangle', color='white', opacity=1):
        # If no layername is provided, take the last added layer
        if layername is None:
            layername = list(self.layers.keys())[-1]

        self.layers[layername].setStyle({'cartoon': {'arrows': arrows,
                                                     'tubes': tubes,
                                                     'opacity': opacity,
                                                     'style': style,  # trace, oval, rectangle (def), parabola, edged
                                                     'color': color}})

    def sticks(self, layername=None, color='green', colorscheme='greenCarbon', radius=0.2):
        # If no layername is provided, take the last added layer
        if layername is None:
            layername = list(self.layers.keys())[-1]

        settings = {'color': color, 'radius': radius}
        if colorscheme is not None:
            settings['colorscheme'] = colorscheme

        self.layers[layername].setStyle({}, {'stick': settings})

    def lines(self, layername=None, color='green', colorscheme='greenCarbon'):
        # If no layername is provided, take the last added layer
        if layername is None:
            layername = list(self.layers.keys())[-1]

        settings = {'color': color}
        if colorscheme is not None:
            settings['colorscheme'] = colorscheme

        self.layers[layername].setStyle({}, {'line': settings})

    def spheres(self, layername=None, color='green', colorscheme='greenCarbon', scale=1):
        # If no layername is provided, take the last added layer
        if layername is None:
            layername = list(self.layers.keys())[-1]

        settings = {'color': color, 'scale': scale}
        if colorscheme is not None:
            settings['colorscheme'] = colorscheme

        self.layers[layername].setStyle({}, {'sphere': settings})

    def __layername(self, layername):
        if layername is None:
            layernames = list(self.layers.keys())
            layer_nr = [int(lay[-1]) for lay in layernames if lay.startswith('layer_') and lay[-1].isdigit()]
            if len(layer_nr) == 0:
                layername = 'layer_0'
            else:
                layername = f'layer_{max(layer_nr)+1}'
        return layername

    def show(self):
        self.view.zoomTo()
        self.view.show()
