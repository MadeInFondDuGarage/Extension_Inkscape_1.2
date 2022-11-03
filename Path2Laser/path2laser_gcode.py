#!/usr/bin/env python

################################################################################
#
# Engraving
#
# LT Notes to self: See wiki.inkscape.org/wiki/index.php/PythonEffectTutorial
# To create anything in the Inkscape document, look at the XML editor for
# details of how such an element looks in XML, then follow this model.
# layer number n appears in XML as <svg:g id="layern" inkscape:label="layername">
#
# to create it, use
# Mylayer = self.svg.add(Layer.new('layername'))
#
# group appears in XML as <svg:g id="gnnnnn"> where nnnnn is a number
#
# to create it, use
# Mygroup = parent.add(Group(gcodetools="My group label")
# where parent may be the layer or a parent group. To get the parent group, you can use
# parent = self.selected_paths[layer][0].getparent()
################################################################################
import cmath
import math
import os
import re
import sys

import numpy

import inkex
from inkex.bezier import bezierlength, bezierparameterize, beziertatlength
from inkex import Transform, PathElement, TextElement, Tspan, Group, Layer, CubicSuperPath


################################################################################
#
# Styles and additional parameters
#
################################################################################

TAU = math.pi * 2
STRAIGHT_TOLERANCE = 0.01 #0.0001 #tolérance pour passer en ligne
STRAIGHT_DISTANCE_TOLERANCE = 0.01 #0.0001 

EMC_TOLERANCE_EQUAL = 0.1 #0.1 #0.00001 #lier a la tolérance d'arc minimum a voir si il faudrait pas le remonter

options = {}


INTERSECTION_RECURSION_DEPTH = 10
INTERSECTION_TOLERANCE = 0.00001

MARKER_STYLE = {
        "biarc_style" : {
                'biarc0': { 'stroke': '#88f', 'fill': 'none',  'stroke-width':'1' },
                'biarc1':{ 'stroke': '#8f8', 'fill': 'none', 'stroke-width':'1' },
                'line':{ 'stroke': '#f88', 'fill': 'none','stroke-width':'1' },
                'area':{ 'stroke': '#777', 'fill': 'none', 'stroke-width':'0.1' },
            },
        "dxf_points":{ "stroke": "#ff0000", "fill": "#ff0000"},

    }

################################################################################
###
###        Class d'entrée du module
###
################################################################################

class Path2Laser(inkex.EffectExtension):

    def __init__(self):
        inkex.Effect.__init__(self)
        #Gravure
        self.arg_parser.add_argument("--grav-travel-speed",dest="grav_travel_speed", default="3000", help="Travel speed (mm/min)")
        self.arg_parser.add_argument("--grav-laser-speed", type=int,dest="grav_laser_speed", default="750",  help="Laser speed (mm/min)")
        self.arg_parser.add_argument("--grav-laser-power",  type=int,dest="grav_laser_power",  default="255",   help="S# is 256 or 10000 for full power")
        self.arg_parser.add_argument("--grav-passes",  type=int,dest="grav_passes",  default="1", help="Quantity of passes")
        self.arg_parser.add_argument("--grav-pass-depth", type=float ,dest="grav_pass_depth", default="0",  help="Depth of laser cut")
        self.arg_parser.add_argument("--grav-color",  dest="grav_color",  default="0x000000ff",   help="Color for engraving")

        #Decoupe
        self.arg_parser.add_argument("--cut-travel-speed",dest="cut_travel_speed", default="3000", help="Travel speed (mm/min)")
        self.arg_parser.add_argument("--cut-laser-speed", type=int,dest="cut_laser_speed", default="750",  help="Laser speed (mm/min)")
        self.arg_parser.add_argument("--cut-laser-power",  type=int,dest="cut_laser_power",  default="255",   help="S# is 256 or 10000 for full power")
        self.arg_parser.add_argument("--cut-passes",  type=int,dest="cut_passes",  default="3", help="Quantity of passes")
        self.arg_parser.add_argument("--cut-pass-depth", type=float ,dest="cut_pass_depth", default="0.5",  help="Depth of laser cut")
        self.arg_parser.add_argument("--cut-color",  dest="cut_color",  default="0xff0000ff",   help="Color for cuting")
        
        #Gcode
        self.arg_parser.add_argument("--Gcode-start", dest="Gcode_start",default="",help="start gcode command")
        self.arg_parser.add_argument("--Gcode-end", dest="Gcode_end",default="",help="end gcode command")
        self.arg_parser.add_argument("--laser-on-command", dest="laser_on_command",default="M4",help="Laser gcode command")
        self.arg_parser.add_argument("--laser-off-command", dest="laser_off_command", default="M5", help="Laser gcode end command")
        self.arg_parser.add_argument("--unit", dest="unit", default="G21 (All units in mm)",        help="Units either mm or inches")       
        self.arg_parser.add_argument("--power-delay",type=float , dest="power_delay", default=0,help="Laser power-on delay (ms)")

        #Fichier
        self.arg_parser.add_argument("-d","--directory",    dest="directory",   default="/home/david/",     help="Output directory")
        self.arg_parser.add_argument("-f","--filename",     dest="file",        default="output.gcode",     help="File name")
        self.arg_parser.add_argument("--add-numeric-suffix-to-filename",        type=inkex.Boolean,         dest="add_numeric_suffix_to_filename",      default=False,                          help="Add numeric suffix to file name")
        self.arg_parser.add_argument("--only-one-file",        type=inkex.Boolean,         dest="only_one_file",      default=True,                          help="Add numeric suffix to file name")

        self.arg_parser.add_argument("--tab", dest="tab",   help="The selected UI-tab when OK was pressed")
        self.arg_parser.add_argument("--path-to-gcode-depth-function", dest="path_to_gcode_depth_function", default="zd", help="Path to gcode depth function.")
        self.arg_parser.add_argument("--biarc-tolerance", type=float, dest="biarc_tolerance",  default="1", help="Tolerance used when calculating biarc interpolation.")
        self.arg_parser.add_argument("--biarc-max-split-depth", type=int, dest="biarc_max_split_depth", default="5", help="Defines maximum depth of splitting while approximating using biarcs.")
        self.arg_parser.add_argument("--min-arc-radius", type=float, dest="min_arc_radius", default=".5", help="All arc having radius less than minimum will be considered as straight line")

       
################################################################################
###
###        Effect
###
###        Main function of Gcodetools class
###
################################################################################
    def effect(self) :

        global options
        options = self.options
        options.self = self
        options.doc_root = self.document.getroot()

        # This automatically calls any `tab_{tab_name_in_inx}` which in this
        # extension is A LOT of different functions. So see all method prefixed
        # with tab_ to find out what's supported here.
        
        #defnition de l'outils pour la gravure laser, version par défaut dans un dict
        if self.options.power_delay < 0.1:
            Pause=""
        else:
            Pause="\nG4 P{}".format(round(self.options.power_delay,1))
        
        self.grav_tools = {
            "name": "Laser Engraver",
            "id": "Laser Engraver",
            "penetration feed": self.options.grav_laser_speed,
            "feed": self.options.grav_laser_speed,
            "gcode before path": self.options.laser_on_command + " S" + str(int(self.options.grav_laser_power)) + Pause+ "(A)", 
            "gcode after path": self.options.laser_off_command + " S0" + "(B)\n" + "G0 F" + self.options.grav_travel_speed + "(B)",
            "Pass_in":self.options.grav_passes,
            "Pass_depth":self.options.grav_pass_depth
                    }
        
        #defnition de l'outils pour la découpe laser, version par défaut dans un dict        
        self.cut_tools = {
            "name": "Laser cutter",
            "id": "Laser cutter",
            "penetration feed": self.options.cut_laser_speed,
            "feed": self.options.cut_laser_speed,
            "gcode before path":  self.options.laser_on_command + " S" + str(int(self.options.cut_laser_power)) + Pause +"(A)",  
            "gcode after path": self.options.laser_off_command + " S0" + "(B)\n" + "G0 F" + self.options.cut_travel_speed + "(B)",
            "Pass_in":self.options.cut_passes,
            "Pass_depth":self.options.cut_pass_depth                       

                    }

        self.Zone_name={'Definition':'Path2Laser'}
        
        #apres la vérification de la préparation, si tout est ok on lance le precessus de creation
        if self.Preparation():
            self.tools=self.grav_tools #tools par défaut
            self.Go_Gcode()
        

################################################################################
###
###        go gcode
###
###        Main function of Gcodetools class
###TODO
################################################################################
    def Preparation(self):
        #"Orientation points have not been defined! A default set of orientation points has been automatically added."),"warning")
        self.get_info()
        self.cut_rgb=self.options.cut_color# [self.options.cut_red,self.options.cut_green,self.options.cut_blue]
        self.grav_rgb=self.options.grav_color #[self.options.grav_red,self.options.grav_green,self.options.grav_blue]
    
        if self.orientation_points == {}:
            self.orientation(self.layers[min(0,len(self.layers)-1)] )
            self.error("Aucun point d'orientation present. Creation des éléments.","warning")
            return False
        return True
        
    def Go_Gcode(self):
        self.path_to_gcode()

################################################################################
#
# Path to Gcode
#
################################################################################
    def path_to_gcode(self):

        def point_to_point_d2(a, b):
            return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2

        def sort_lines(lines):
            if len(lines) == 0:
                return []
            lines = [[key] + lines[key] for key in range(len(lines))]
            keys = [0]
            end_point = lines[0][3:]
            del lines[0]
            while len(lines) > 0:
                dist = [[point_to_point_d2(end_point, lines[i][1:3]), i] for i in range(len(lines))]
                i = min(dist)[1]
                keys.append(lines[i][0])
                end_point = lines[i][3:]
                del lines[i]
            return keys
    
        def sort_curves(curves):
            lines = []
            for curve in curves:
                lines += [curve[0][0][0] + curve[-1][-1][0]]
            return sort_lines(lines)
         

        if self.selected_paths == {} :
            paths = self.paths
            # self.error("No paths are selected! Trying to work on all available paths.")
        else:
            paths = self.selected_paths
        self.check_dir()
        gcode = ""
        comment="()"
        
        attr = {self.Zone_name['Definition']: "Preview groups"}
        PreviewGroups = self.svg.add(Group(**attr))
        PreviewGroups.label="Preview groups"
        
        for layer in self.layers:
            #je ne traite que les calques/groupes visible qui sont récuperer dans la partie get_infos()
            comment_layer=""
 
            if layer in paths:# la selection de l'état de visibilité des calques et groups ce fait sur la gestion des layers
                # transform simple path to get all var about orientation
                self.transform_csp([[[[0, 0], [0, 0], [0, 0]], [[0, 0], [0, 0], [0, 0]]]], layer)

                #je fait un commentaire spécifique a bCNC pour chaque calque/groupe 
                curves = []
                if layer.label != None:
                    comment_layer= "\n(Block-name:" + layer.label +")\n"  #pour chaque calque qui font partie des ensemble de chemin
                else:
                    comment_layer= "\n(Block-name:" + layer.get("id") +")\n"
                    
                for path in paths[layer]:
                    if "d" not in path.keys():
                        self.error("Warning: One or more paths do not have 'd' parameter, try to Ungroup (Ctrl+Shift+G) and Object to Path (Ctrl+Shift+C)!")
                        continue
                    csp = path.path.to_superpath()
                    csp = self.apply_transforms(path, csp)
                    id_ = path.get("id")

                    stroke = path.style('stroke')
                    display=path.style('display')
                       
                    comment_path="\n("+path.get("id")+")\n" #pour chaque chemin
                                        
                    comment=comment_layer + comment_path
                    
                    #ici je definie un set d'outils en fonction du calque ou je me trouve.
                    Ok=False
                    if inkex.Color(stroke).to_rgb() ==inkex.Color(self.cut_rgb).to_rgb():  
                        self.tools=self.cut_tools
                        Ok=True
                    elif inkex.Color(stroke).to_rgb() == inkex.Color(self.grav_rgb).to_rgb():
                        self.tools=self.grav_tools
                        Ok=True

                    #gestion des chemin qui sont rendue invisible  ou avec la couleur verte
                    if Ok and display=='inline':
                        curves += [
                            [
                                [id_, self.tools["Pass_in"],self.tools["Pass_depth"], comment],
                                [self.parse_curve([subpath], layer) for subpath in csp] #ceci découpe les chemins fait en plusieur partie de courbe simple
                            ]
                        ]
                        comment_layer= ""  #je reinitialise le commentaire du claque/group une fois au moins un chemin a été pris en compte
   
                #maintenant que les chemin sont répertorier je les dessine
                for curve in curves:
                    for subcurve in curve[1]:
                        self.draw_curve(subcurve, layer, PreviewGroups,name=curve[0][0])

                #passage en ordonnancement direct
                keys = sort_curves([curve[1] for curve in curves])
                for key in keys:
                    #attention ici on est sur un ensemble avec une profondeur de passe
                    for step in range(0, int(curves[key][0][1])):
                        z =   -(curves[key][0][2]* (step))
                        if curves[key][0][3] != "()" and step ==0:
                            gcode += curves[key][0][3]  # add comment

                        for curve in curves[key][1]:
                            gcode += self.generate_gcode(curve, layer, z)
 
            
        self.export_gcode(gcode)
################################################################################
#
# Generate Gcode
# Generates Gcode on given curve.
#
# Curve definition [start point, type = {'arc','line','move','end'}, arc center, arc angle, end point, [zstart, zend]]
#TODO
################################################################################
    def generate_gcode(self, curve, layer, depth):
        tool = self.tools
        g = ""
        
        def c(c):
            c = [c[i] if i < len(c) else None for i in range(6)]
            if c[5] == 0:
                c[5] = None
            s = [" X", " Y", " Z", " I", " J", " K"]
            s1 = ["", "", "", "", "", ""]
            m = [1, 1, 1, 1, 1, 1]
            a = [0, 0, 0.0, 0, 0, 0]
            r = ''
            for i in range(6):
                if c[i] is not None:
                    r += s[i] + ("{:f}".format(c[i] * m[i] + a[i])) + s1[i]
            return r

        def calculate_angle(a, current_a):
            return min(
                    [abs(a - current_a % TAU + TAU), a + current_a - current_a % TAU + TAU],
                    [abs(a - current_a % TAU - TAU), a + current_a - current_a % TAU - TAU],
                    [abs(a - current_a % TAU), a + current_a - current_a % TAU])[1]

        if len(curve) == 0:
            return ""

        try:
            self.last_used_tool is None
        except:
            self.last_used_tool = None

        lg = 'G00'
        #TODO a supprimer
        f = " F{}".format(tool['feed'])
        for i in range(1, len(curve)):
            #    Creating Gcode for curve between s=curve[i-1] and si=curve[i] start at s[0] end at s[4]=si[0]
            s = curve[i - 1]
            si = curve[i]
            feed = f if lg not in ['G01', 'G02', 'G03'] else ''
            if s[1] == 'move':
                g += "G00" + c(si[0]+[depth]) + "\n" + tool['gcode before path'] + "\n"
                lg = 'G00'
            elif s[1] == 'end':
                g += tool['gcode after path'] + "\n"
                lg = 'G00'
            elif s[1] == 'line':
                # if lg == "G00":
                #     g += "G00" + c([None, None, s[5][0] + depth]) + feed + "(2)\n"
                g += "G01" + c(si[0] ) + feed + "(3)\n"
                lg = 'G01'
            elif s[1] == 'arc':
                r = [(s[2][0] - s[0][0]), (s[2][1] - s[0][1])]
                # if lg == "G00":
                    # g += "G00" + c([None, None, s[5][0] + depth]) + feed + "(4)\n"
                if (r[0] ** 2 + r[1] ** 2) > self.options.min_arc_radius ** 2:
                    r1 = (P(s[0]) - P(s[2]))
                    r2 = (P(si[0]) - P(s[2]))
                    if abs(r1.mag() - r2.mag()) < 0.001:
                        g += ("G02" if s[3] < 0 else "G03") + c(si[0] + [s[5][1] + depth, (s[2][0] - s[0][0]), (s[2][1] - s[0][1])]) + feed + "(5)\n"
                    else:
                        r = (r1.mag() + r2.mag()) / 2
                        g += ("G02" if s[3] < 0 else "G03") + c(si[0] + [s[5][1] + depth]) + " R{:f}".format(r) + feed + "(6)\n"
                    lg = 'G02'
                else:
                    g += "G00" + c(si[0] + [s[5][1] + depth]) + feed + "(7)\n"
                    lg = 'G00'
        if si[1] == 'end':
            g += tool['gcode after path'] + "(8)\n"
        return g

################################################################################
###        Errors handling function, notes are just printed into Logfile, 
###        warnings are printed into log file and warning message is displayed but
###        extension continues working, errors causes log and execution is halted
###        Notes, warnings adn errors could be assigned to space or comma or dot 
###        sepparated strings (case is ignoreg).
################################################################################
    def error(self, s, type_= "Warning"):
        errors = """
                        Error     
                        wrong_orientation_points    
                        area_tools_diameter_error
                        no_tool_error
                        active_layer_already_has_tool
                        active_layer_already_has_orientation_points
                    """
        s = str(s)
        #version simplifier sans le log, si c'est pas une erreur on met juste un message d'info, sinon j'interompt le process
        if type_.lower() in re.split("[\s\n,\.]+", errors.lower()) :
            inkex.errormsg(s+"\n")        
            sys.exit()
        else :
            inkex.errormsg(s)        

################################################################################
###         
###        Zone Get_transform reverse transform et apply_transform
################################################################################
    def get_transforms(self, g):
        root = self.document.getroot()
        trans = []
        while g != root:
            if 'transform' in g.keys():
                t = g.get('transform')
                t = Transform(t).matrix
                trans = (Transform(t) @ Transform(trans)).matrix if trans != [] else t

            g = g.getparent()
        return trans

    def reverse_transform(self, transform):
        trans = numpy.array(transform + ([0, 0, 1],))
        if numpy.linalg.det(trans) != 0:
            trans = numpy.linalg.inv(trans).tolist()[:2]
            return trans
        else:
            return transform

    def apply_transforms(self, g, csp, reverse=False):
        trans = self.get_transforms(g)
        if trans:
            if not reverse:
                for comp in csp:
                    for ctl in comp:
                        for pt in ctl:
                            pt[0], pt[1] = Transform(trans).apply_to_point(pt)

            else:
                for comp in csp:
                    for ctl in comp:
                        for pt in ctl:
                            pt[0], pt[1] = Transform(self.reverse_transform(trans)).apply_to_point(pt)
        return csp

    def transform_scalar(self, x, layer, reverse=False):
        return self.transform([x, 0], layer, reverse)[0] - self.transform([0, 0], layer, reverse)[0]

    def transform(self, source_point, layer, reverse=False):
        if layer not in self.transform_matrix:
            for i in range(self.layers.index(layer), -1, -1):
                if self.layers[i] in self.orientation_points:
                    break
            if self.layers[i] not in self.orientation_points:
                self.error("Orientation points have not been found!", "error")
            elif self.layers[i] in self.transform_matrix:
                self.transform_matrix[layer] = self.transform_matrix[self.layers[i]]
                self.Zcoordinates[layer] = self.Zcoordinates[self.layers[i]]
            else:
                orientation_layer = self.layers[i]
                if len(self.orientation_points[orientation_layer]) > 1:
                    self.error(f"There are more than one orientation point groups in '{orientation_layer.label}' layer")
                points = self.orientation_points[orientation_layer][0]
                if len(points) == 2:
                    points += [[[(points[1][0][1] - points[0][0][1]) + points[0][0][0], -(points[1][0][0] - points[0][0][0]) + points[0][0][1]], [-(points[1][1][1] - points[0][1][1]) + points[0][1][0], points[1][1][0] - points[0][1][0] + points[0][1][1]]]]
                if len(points) == 3:

                    # for point in points:
                    #    Zcoordinates definition taken from Orientatnion point 1 and 2
                    self.Zcoordinates[layer] = [max(points[0][1][2], points[1][1][2]), min(points[0][1][2], points[1][1][2])]
                    matrix = numpy.array([
                        [points[0][0][0], points[0][0][1], 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, points[0][0][0], points[0][0][1], 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, points[0][0][0], points[0][0][1], 1],
                        [points[1][0][0], points[1][0][1], 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, points[1][0][0], points[1][0][1], 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, points[1][0][0], points[1][0][1], 1],
                        [points[2][0][0], points[2][0][1], 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, points[2][0][0], points[2][0][1], 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, points[2][0][0], points[2][0][1], 1]
                    ])

                    if numpy.linalg.det(matrix) != 0:
                        m = numpy.linalg.solve(matrix,
                                               numpy.array(
                                                       [[points[0][1][0]], [points[0][1][1]], [1], [points[1][1][0]], [points[1][1][1]], [1], [points[2][1][0]], [points[2][1][1]], [1]]
                                               )
                                               ).tolist()
                        self.transform_matrix[layer] = [[m[j * 3 + i][0] for i in range(3)] for j in range(3)]

                    else:
                        self.error("Orientation points are wrong! (if there are two orientation points they should not be the same. If there are three orientation points they should not be in a straight line.)", "error")
                else:
                    self.error("Orientation points are wrong! (if there are two orientation points they should not be the same. If there are three orientation points they should not be in a straight line.)", "error")

            self.transform_matrix_reverse[layer] = numpy.linalg.inv(self.transform_matrix[layer]).tolist()


            # Zautoscale is obsolete
            self.Zauto_scale[layer] = 1

        x = source_point[0]
        y = source_point[1]
        if not reverse:
            t = self.transform_matrix[layer]
        else:
            t = self.transform_matrix_reverse[layer]
        return [t[0][0] * x + t[0][1] * y + t[0][2], t[1][0] * x + t[1][1] * y + t[1][2]]

    def transform_csp(self, csp_, layer, reverse=False):
        csp = [[[csp_[i][j][0][:], csp_[i][j][1][:], csp_[i][j][2][:]] for j in range(len(csp_[i]))] for i in range(len(csp_))]
        for i in range(len(csp)):
            for j in range(len(csp[i])):
                for k in range(len(csp[i][j])):
                    csp[i][j][k] = self.transform(csp[i][j][k], layer, reverse)
        return csp

        
###############################################################################
###
###        Get Gcodetools info from the svg
###
################################################################################
    def get_info(self):
        """Get Gcodetools info from the svg"""
        self.selected_paths = {}
        self.paths = {}
        self.tools = {}
        self.orientation_points = {}
        self.graffiti_reference_points = {}
        self.layers = [self.document.getroot()]
        self.Zcoordinates = {}
        self.transform_matrix = {}
        self.transform_matrix_reverse = {}
        self.Zauto_scale = {}
        self.in_out_reference_points = []
        self.DocName=""
              
        def recursive_search(g, layer, selected=False,level=""):

            items = g.getchildren()
            items.reverse()
            for i in items:
                if selected:
                    self.svg.selected[i.get("id")] = i

                if i.get(self.Zone_name['Definition']) == "Orientation group":
                    points = self.get_orientation_points(i)

                    if points is not None:
                        self.orientation_points[layer] = self.orientation_points[layer] + [points[:]] if layer in self.orientation_points else [points[:]]
                    else:
                        self.error(f"Warning! Found bad orientation points in '{layer.label}' layer. Resulting Gcode could be corrupt!")
                elif isinstance(i, inkex.Group) or isinstance(i, Layer):
                    invisible=i.style('display')
                    if invisible=='inline':
                        self.layers += [i]
                        recursive_search(i, i, (i.get("id") in self.svg.selected))#la partie selected permet de concerver le récursif sur le group
                        
                elif isinstance(i, inkex.PathElement):
                    if self.Zone_name['Definition'] not in i.keys():
                        self.paths[layer] = self.paths[layer] + [i] if layer in self.paths else [i]

                        if i.get("id") in self.svg.selected.ids:
                            self.selected_paths[layer] = self.selected_paths[layer] + [i] if layer in self.selected_paths else [i]
        
                elif i.get("id") in self.svg.selected:
                    # xgettext:no-pango-format
                    if i.label != None:
                        objet=i.label
                    else:
                        objet=i.get("id")
                        
                    self.error("This extension works with Paths and Dynamic Offsets and groups of them only! "
                               "All other objects will be ignored!\n"
                               "Solution 1: press Path->Object to path or Shift+Ctrl+C.\n"
                               "Solution 2: Path->Dynamic offset or Ctrl+J.\n"
                               "Solution 3: export all contours to PostScript level 2 (File->Save As->.ps) and File->Import this file.\n{}".format(objet))
    
        recursive_search(self.document.getroot(), self.document.getroot())
    
        if len(self.layers) == 1:
            self.error("Document has no layers! Add at least one layer using layers panel (Ctrl+Shift+L)", "error")
        root = self.document.getroot()
    
        #récupération du nom
        self.DocName=root.get("sodipodi:docname")
        
        if root in self.selected_paths or root in self.paths:
            self.error("Warning! There are some paths in the root of the document, but not in any layer! Using bottom-most layer for them.")
    
        if root in self.selected_paths:
            if self.layers[-1] in self.selected_paths:
                self.selected_paths[self.layers[-1]] += self.selected_paths[root][:]
            else:
                self.selected_paths[self.layers[-1]] = self.selected_paths[root][:]
            del self.selected_paths[root]
    
        if root in self.paths:
            if self.layers[-1] in self.paths:
                self.paths[self.layers[-1]] += self.paths[root][:]
            else:
                self.paths[self.layers[-1]] = self.paths[root][:]
            del self.paths[root]
    
    def get_orientation_points(self, g):
        items = g.getchildren()
        items.reverse()
        p2 = []
        p = None
        for i in items:
            if isinstance(i, inkex.Group):
                p2 += [i]
        p = p2
        if p is None:
            return None
        points = []
        for i in p:
            point = [[], []]
            for node in i:
                if node.get(self.Zone_name['Definition']) == "Orientation point arrow":
                    csp = node.path.transform(node.composed_transform()).to_superpath()
                    point[0] = csp[0][0][1]
                if node.get(self.Zone_name['Definition']) == "Orientation point text":
                    r = re.match(r'(?i)\s*\(\s*(-?\s*\d*(?:,|\.)*\d*)\s*;\s*(-?\s*\d*(?:,|\.)*\d*)\s*;\s*(-?\s*\d*(?:,|\.)*\d*)\s*\)\s*', node.get_text())
                    point[1] = [float(r.group(1)), float(r.group(2)), float(r.group(3))]
            if point[0] != [] and point[1] != []:
                points += [point]
        if len(points) == len(p2) == 2 :
            return points
        else:
            return None


################################################################################
###
###        Orientation
###
################################################################################
    def orientation(self, layer=None):
        if layer is None:
            layer = self.svg.get_current_layer() if self.svg.get_current_layer() is not None else self.document.getroot()

        transform = self.get_transforms(layer)
        if transform:
            transform = self.reverse_transform(transform)
            transform = str(Transform(transform))

        if layer in self.orientation_points:
            None
            #"Active layer already has orientation points! Remove them or select another layer!", "error")
        else:
            attr = {self.Zone_name['Definition']: "Orientation group"}
            
            if transform:
                attr["transform"] = transform
    
            orientation_group = layer.add(Group(**attr))
            orientation_group.label="Orientation group"
            doc_height = self.svg.unittouu(self.document.getroot().get('height'))
            if self.document.getroot().get('height') == "100%":
                doc_height = 1052.3622047
            if self.options.unit == "G21 (All units in mm)":
                points = [[0., 0., 0.], [100., 0., 0.], [0., 100., 0.]]
            elif self.options.unit == "G20 (All units in inches)":
                points = [[0., 0., 0.], [5., 0., 0.], [0., 5., 0.]]
            points = points[:2]
            for i in points:
                name = {self.Zone_name['Definition']:"Orientation point (2 points)"}
                grp = orientation_group.add(Group(**name))
                grp.label="({}; {}; {})".format(i[0], i[1], i[2])
                elem = grp.add(PathElement(style="stroke:none;fill:#000000;"))
                elem.set(self.Zone_name['Definition'], "Orientation point arrow")
                elem.path = 'm {},{} 2.9375,-6.343750000001 0.8125,1.90625 6.843748640396,'\
                    '-6.84374864039 0,0 0.6875,0.6875 -6.84375,6.84375 1.90625,0.812500000'\
                    '001 z'.format(i[0], -i[1] + doc_height)
                elem.label="Orientation points arrow"
                draw_text("({}; {}; {})".format(i[0], i[1], i[2]), (i[0] + 10), (-i[1] - 10 + doc_height), group=grp, gcodetools_tag="Orientation point text",label_text="Orientation point text")


################################################################################
###
###       Export_gcode
###
###        Main function of Gcodetools class
###TODO faire les export en 2  ou 1 fichier
################################################################################
    def check_dir(self):
        #TODO voir si je prépare les différent type de fichier ici
        if os.path.isdir(self.options.directory):
            if os.path.isfile(os.path.join(self.options.directory, 'header')):
                with open(os.path.join(self.options.directory, 'header')) as f:
                    self.header = f.read()
            else:
                self.header ="\n(GCODE realiser avec l'extension inskape Path2Laser_GCODE de madeinfonddugarage)\n{}\n".format(self.options.Gcode_start).replace("\\n","\n")
            if os.path.isfile(os.path.join(self.options.directory, 'footer')):
                with open(os.path.join(self.options.directory, 'footer')) as f:
                    self.footer = f.read()
            else:
                self.footer = self.options.Gcode_end.replace("\\n","\n")
            self.header += self.options.unit + "\n"
        else:
            self.error("Directory does not exist! Please specify existing directory at Preferences tab!", "error")
            return False
    
        if self.options.add_numeric_suffix_to_filename:
            dir_list = os.listdir(self.options.directory)
            if "." in self.options.file:
                r = re.match(r"^(.*)(\..*)$", self.options.file)
                ext = r.group(2)
                name = r.group(1)
            else:
                ext = ""
                name = self.options.file
            max_n = 0
            for s in dir_list:
                r = re.match(r"^{}_0*(\d+){}$".format(re.escape(name), re.escape(ext)), s)
                if r:
                    max_n = max(max_n, int(r.group(1)))
            filename = name + "_" + ("0" * (4 - len(str(max_n + 1))) + str(max_n + 1)) + ext
            self.options.file = filename
    
        try:
            with open(os.path.join(self.options.directory, self.options.file), "w") as f:
                pass
        except:
            self.error("Can not write to specified file!\n{}".format(os.path.join(self.options.directory, self.options.file)), "error")
            return False
        return True              

    def parse_curve(self, p, layer, w=None, f=None):
        c = []
        if len(p) == 0:
            return []
        p = self.transform_csp(p, layer)

        # Sort to reduce Rapid distance
        k = list(range(1, len(p)))
        keys = [0]
        while len(k) > 0:
            end = p[keys[-1]][-1][1]
            dist = None
            for i in range(len(k)):
                start = p[k[i]][0][1]
                dist = max((-((end[0] - start[0]) ** 2 + (end[1] - start[1]) ** 2), i), dist)
            keys += [k[dist[1]]]
            del k[dist[1]]
        for k in keys:
            subpath = p[k]
            c += [[[subpath[0][1][0], subpath[0][1][1]], 'move', 0, 0]]
            for i in range(1, len(subpath)):
                sp1 = [[subpath[i - 1][j][0], subpath[i - 1][j][1]] for j in range(3)]
                sp2 = [[subpath[i][j][0], subpath[i][j][1]] for j in range(3)]
                c += biarc(sp1, sp2, 0, 0) if w is None else biarc(sp1, sp2, -f(w[k][i - 1]), -f(w[k][i]))
            c += [[[subpath[-1][1][0], subpath[-1][1][1]], 'end', 0, 0]]
        return c
       
    def export_gcode(self,gcode):
        #export du gcode un fois finie
        f = open(self.options.directory+self.options.file, "w")
        f.write( "(Block-name:Start)\n"+"\n(Gcode creer a partir du fichier " + self.DocName + ")" +self.header +self.options.laser_off_command + " S0" + "\n" +  "G0 F" + self.options.grav_travel_speed + "\n" + gcode + "\n(Block-name:End)\n" + self.footer) 
        f.close()


 

################################################################################
#
# draw csp
#
################################################################################
    def draw_csp(self, csp, layer=None, group=None, fill='none', stroke='#178ade', width=0.354, style=None):
        if layer is not None:
            csp = self.transform_csp(csp, layer, reverse=True)
        if group is None and layer is None:
            group = self.document.getroot()
        elif group is None and layer is not None:
            group = layer
        csp = self.apply_transforms(group, csp, reverse=True)
        if style is not None:
            return draw_csp(csp, group=group, style=style)
        else:
            return draw_csp(csp, group=group, fill=fill, stroke=stroke, width=width)

    def draw_curve(self, curve, layer, groupbase=None, style=MARKER_STYLE["biarc_style"],name=""):
        
        attr={self.Zone_name['Definition']:'Preview group'}
        group= groupbase.add(Group(**attr))
        group.label='Preview de {}'.format(name)
            
        s = ''
        arcn = 0

        transform = self.get_transforms(group)
        if transform:
            transform = self.reverse_transform(transform)
            transform = str(Transform(transform))

        a = [0., 0.]
        b = [1., 0.]
        c = [0., 1.]
        k = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])
        a = self.transform(a, layer, True)
        b = self.transform(b, layer, True)
        c = self.transform(c, layer, True)
        if ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) * k > 0:
            reverse_angle = 1
        else:
            reverse_angle = -1

        for sk in curve:
            si = sk[:] #sk contient les movement et les position relative
            si[0] = self.transform(si[0], layer, True)
            si[2] = self.transform(si[2], layer, True) if type(si[2]) == type([]) and len(si[2]) == 2 else si[2]

            if s != '':
                if s[1] == 'line':
                    attr={self.Zone_name['Definition']:'Preview line'}
                    elem = group.add(PathElement(**attr))
                    elem.label='Preview line {}'.format(name)
                    elem.transform = transform
                    elem.style = style['line']
                    elem.path = 'M {},{} L {},{}'.format(s[0][0], s[0][1], si[0][0], si[0][1])
                elif s[1] == 'arc':
                    arcn += 1
                    sp = s[0]
                    c = s[2]
                    s[3] = s[3] * reverse_angle

                    a = ((P(si[0]) - P(c)).angle() - (P(s[0]) - P(c)).angle()) % TAU  # s[3]
                    if s[3] * a < 0:
                        if a > 0:
                            a = a - TAU
                        else:
                            a = TAU + a
                    r = math.sqrt((sp[0] - c[0]) ** 2 + (sp[1] - c[1]) ** 2)
                    a_st = (math.atan2(sp[0] - c[0], - (sp[1] - c[1])) - math.pi / 2) % (math.pi * 2)
                    if a > 0:
                        a_end = a_st + a
                        st = style['biarc{}'.format(arcn % 2)]
                    else:
                        a_end = a_st * 1
                        a_st = a_st + a
                        st = style['biarc{}'.format(arcn % 2)]
                    attr={self.Zone_name['Definition']:'Preview arc'}
                    elem = group.add(PathElement.arc(c, r, start=a_st, end=a_end,
                                                     open=True, **attr))
                    elem.label='Preview arc {}'.format(name)
                    elem.transform = transform
                    elem.style = st

            s = si

################################################################################
#
# Biarc function
#
# Calculates biarc approximation of cubic super path segment
#  splits segment if needed or approximates it with straight line
#
################################################################################
def biarc(sp1, sp2, z1, z2, depth=0):
    def biarc_split(sp1, sp2, z1, z2, depth):
        if depth < options.biarc_max_split_depth:
            sp1, sp2, sp3 = csp_split(sp1, sp2)
            l1 = cspseglength(sp1, sp2)
            l2 = cspseglength(sp2, sp3)
            if l1 + l2 == 0:
                zm = z1
            else:
                zm = z1 + (z2 - z1) * l1 / (l1 + l2)
            return biarc(sp1, sp2, z1, zm, depth + 1) + biarc(sp2, sp3, zm, z2, depth + 1)
        else:
            return [[sp1[1], 'line', 0, 0, sp2[1], [z1, z2]]]

    P0 = P(sp1[1])
    P4 = P(sp2[1])
    TS = (P(sp1[2]) - P0)
    TE = -(P(sp2[0]) - P4)
    v = P0 - P4
    tsa = TS.angle()
    tea = TE.angle()
    if TE.mag() < STRAIGHT_DISTANCE_TOLERANCE and TS.mag() < STRAIGHT_DISTANCE_TOLERANCE:
        # Both tangents are zero - line straight
        return [[sp1[1], 'line', 0, 0, sp2[1], [z1, z2]]]
    if TE.mag() < STRAIGHT_DISTANCE_TOLERANCE:
        TE = -(TS + v).unit()
        r = TS.mag() / v.mag() * 2
    elif TS.mag() < STRAIGHT_DISTANCE_TOLERANCE:
        TS = -(TE + v).unit()
        r = 1 / (TE.mag() / v.mag() * 2)
    else:
        r = TS.mag() / TE.mag()
    TS = TS.unit()
    TE = TE.unit()
    tang_are_parallel = ((tsa - tea) % math.pi < STRAIGHT_TOLERANCE or math.pi - (tsa - tea) % math.pi < STRAIGHT_TOLERANCE)
    if (tang_are_parallel and
            ((v.mag() < STRAIGHT_DISTANCE_TOLERANCE or TE.mag() < STRAIGHT_DISTANCE_TOLERANCE or TS.mag() < STRAIGHT_DISTANCE_TOLERANCE) or
             1 - abs(TS * v / (TS.mag() * v.mag())) < STRAIGHT_TOLERANCE)):
        # Both tangents are parallel and start and end are the same - line straight
        # or one of tangents still smaller then tolerance

        # Both tangents and v are parallel - line straight
        return [[sp1[1], 'line', 0, 0, sp2[1], [z1, z2]]]

    c = v * v
    b = 2 * v * (r * TS + TE)
    a = 2 * r * (TS * TE - 1)
    if v.mag() == 0:
        return biarc_split(sp1, sp2, z1, z2, depth)
    asmall = abs(a) < 10 ** -10
    bsmall = abs(b) < 10 ** -10
    csmall = abs(c) < 10 ** -10
    if asmall and b != 0:
        beta = -c / b
    elif csmall and a != 0:
        beta = -b / a
    elif not asmall:
        discr = b * b - 4 * a * c
        if discr < 0:
            raise ValueError(a, b, c, discr)
        disq = discr ** .5
        beta1 = (-b - disq) / 2 / a
        beta2 = (-b + disq) / 2 / a
        if beta1 * beta2 > 0:
            raise ValueError(a, b, c, disq, beta1, beta2)
        beta = max(beta1, beta2)
    elif asmall and bsmall:
        return biarc_split(sp1, sp2, z1, z2, depth)
    alpha = beta * r
    ab = alpha + beta
    P1 = P0 + alpha * TS
    P3 = P4 - beta * TE
    P2 = (beta / ab) * P1 + (alpha / ab) * P3

    def calculate_arc_params(P0, P1, P2):
        D = (P0 + P2) / 2
        if (D - P1).mag() == 0:
            return None, None
        R = D - ((D - P0).mag() ** 2 / (D - P1).mag()) * (P1 - D).unit()
        p0a = (P0 - R).angle() % (2 * math.pi)
        p1a = (P1 - R).angle() % (2 * math.pi)
        p2a = (P2 - R).angle() % (2 * math.pi)
        alpha = (p2a - p0a) % (2 * math.pi)
        if (p0a < p2a and (p1a < p0a or p2a < p1a)) or (p2a < p1a < p0a):
            alpha = -2 * math.pi + alpha
        if abs(R.x) > 1000000 or abs(R.y) > 1000000 or (R - P0).mag() < options.min_arc_radius ** 2:
            return None, None
        else:
            return R, alpha

    R1, a1 = calculate_arc_params(P0, P1, P2)
    R2, a2 = calculate_arc_params(P2, P3, P4)
    if R1 is None or R2 is None or (R1 - P0).mag() < STRAIGHT_TOLERANCE or (R2 - P2).mag() < STRAIGHT_TOLERANCE:
        return [[sp1[1], 'line', 0, 0, sp2[1], [z1, z2]]]

    d = csp_to_arc_distance(sp1, sp2, [P0, P2, R1, a1], [P2, P4, R2, a2])
    if d > options.biarc_tolerance and depth < options.biarc_max_split_depth:
        return biarc_split(sp1, sp2, z1, z2, depth)
    else:
        if R2.mag() * a2 == 0:
            zm = z2
        else:
            zm = z1 + (z2 - z1) * (abs(R1.mag() * a1)) / (abs(R2.mag() * a2) + abs(R1.mag() * a1))

        l = (P0 - P2).l2()
        if l < EMC_TOLERANCE_EQUAL ** 2 or l < EMC_TOLERANCE_EQUAL ** 2 * R1.l2() / 100:
            # arc should be straight otherwise it could be treated as full circle
            arc1 = [sp1[1], 'line', 0, 0, [P2.x, P2.y], [z1, zm]]
        else:
            arc1 = [sp1[1], 'arc', [R1.x, R1.y], a1, [P2.x, P2.y], [z1, zm]]

        l = (P4 - P2).l2()
        if l < EMC_TOLERANCE_EQUAL ** 2 or l < EMC_TOLERANCE_EQUAL ** 2 * R2.l2() / 100:
            # arc should be straight otherwise it could be treated as full circle
            arc2 = [[P2.x, P2.y], 'line', 0, 0, [P4.x, P4.y], [zm, z2]]
        else:
            arc2 = [[P2.x, P2.y], 'arc', [R2.x, R2.y], a2, [P4.x, P4.y], [zm, z2]]

        return [arc1, arc2]
    

################################################################################
# Cubic Super Path additional functions
################################################################################


def csp_from_polyline(line):
    return [[[point[:] for _ in range(3)] for point in subline] for subline in line]


def csp_remove_zero_segments(csp, tolerance=1e-7):
    res = []
    for subpath in csp:
        if len(subpath) > 0:
            res.append([subpath[0]])
            for sp1, sp2 in zip(subpath, subpath[1:]):
                if point_to_point_d2(sp1[1], sp2[1]) <= tolerance and point_to_point_d2(sp1[2], sp2[1]) <= tolerance and point_to_point_d2(sp1[1], sp2[0]) <= tolerance:
                    res[-1][-1][2] = sp2[2]
                else:
                    res[-1].append(sp2)
    return res


def point_inside_csp(p, csp, on_the_path=True):
    # we'll do the raytracing and see how many intersections are there on the ray's way.
    # if number of intersections is even then point is outside.
    # ray will be x=p.x and y=>p.y
    # you can assign any value to on_the_path, by default if point is on the path
    # function will return thai it's inside the path.
    x, y = p
    ray_intersections_count = 0
    for subpath in csp:

        for i in range(1, len(subpath)):
            sp1 = subpath[i - 1]
            sp2 = subpath[i]
            ax, ay, bx, by, cx, cy, dx, dy = csp_parameterize(sp1, sp2)
            if ax == 0 and bx == 0 and cx == 0 and dx == x:
                # we've got a special case here
                b = csp_true_bounds([[sp1, sp2]])
                if b[1][1] <= y <= b[3][1]:
                    # points is on the path
                    return on_the_path
                else:
                    # we can skip this segment because it won't influence the answer.
                    pass
            else:
                for t in csp_line_intersection([x, y], [x, y + 5], sp1, sp2):
                    if t == 0 or t == 1:
                        # we've got another special case here
                        x1, y1 = csp_at_t(sp1, sp2, t)
                        if y1 == y:
                            # the point is on the path
                            return on_the_path
                        # if t == 0 we should have considered this case previously.
                        if t == 1:
                            # we have to check the next segment if it is on the same side of the ray
                            st_d = csp_normalized_slope(sp1, sp2, 1)[0]
                            if st_d == 0:
                                st_d = csp_normalized_slope(sp1, sp2, 0.99)[0]

                            for j in range(1, len(subpath) + 1):
                                if (i + j) % len(subpath) == 0:
                                    continue  # skip the closing segment
                                sp11 = subpath[(i - 1 + j) % len(subpath)]
                                sp22 = subpath[(i + j) % len(subpath)]
                                ax1, ay1, bx1, by1, cx1, cy1, dx1, dy1 = csp_parameterize(sp1, sp2)
                                if ax1 == 0 and bx1 == 0 and cx1 == 0 and dx1 == x:
                                    continue  # this segment parallel to the ray, so skip it
                                en_d = csp_normalized_slope(sp11, sp22, 0)[0]
                                if en_d == 0:
                                    en_d = csp_normalized_slope(sp11, sp22, 0.01)[0]
                                if st_d * en_d <= 0:
                                    ray_intersections_count += 1
                                    break
                    else:
                        x1, y1 = csp_at_t(sp1, sp2, t)
                        if y1 == y:
                            # the point is on the path
                            return on_the_path
                        else:
                            if y1 > y and 3 * ax * t ** 2 + 2 * bx * t + cx != 0:  # if it's 0 the path only touches the ray
                                ray_intersections_count += 1
    return ray_intersections_count % 2 == 1


def csp_close_all_subpaths(csp, tolerance=0.000001):
    for i in range(len(csp)):
        if point_to_point_d2(csp[i][0][1], csp[i][-1][1]) > tolerance ** 2:
            csp[i][-1][2] = csp[i][-1][1][:]
            csp[i] += [[csp[i][0][1][:] for _ in range(3)]]
        else:
            if csp[i][0][1] != csp[i][-1][1]:
                csp[i][-1][1] = csp[i][0][1][:]
    return csp


def csp_simple_bound(csp):
    minx = None
    miny = None
    maxx = None
    maxy = None

    for subpath in csp:
        for sp in subpath:
            for p in sp:
                minx = min(minx, p[0]) if minx is not None else p[0]
                miny = min(miny, p[1]) if miny is not None else p[1]
                maxx = max(maxx, p[0]) if maxx is not None else p[0]
                maxy = max(maxy, p[1]) if maxy is not None else p[1]
    return minx, miny, maxx, maxy


def csp_segment_to_bez(sp1, sp2):
    return sp1[1:] + sp2[:2]


def csp_to_point_distance(csp, p, dist_bounds=(0, 1e100)):
    min_dist = [1e100, 0, 0, 0]
    for j in range(len(csp)):
        for i in range(1, len(csp[j])):
            d = csp_seg_to_point_distance(csp[j][i - 1], csp[j][i], p, sample_points=5)
            if d[0] < dist_bounds[0]:
                return [d[0], j, i, d[1]]
            else:
                if d[0] < min_dist[0]:
                    min_dist = [d[0], j, i, d[1]]
    return min_dist

def csp_seg_to_point_distance(sp1, sp2, p, sample_points=5):
    ax, ay, bx, by, cx, cy, dx, dy = csp_parameterize(sp1, sp2)
    f=0
    dx = dx - p[0]
    dy = dy - p[1]
    if sample_points < 2:
        sample_points = 2
    d = min([(p[0] - sp1[1][0]) ** 2 + (p[1] - sp1[1][1]) ** 2, 0.], [(p[0] - sp2[1][0]) ** 2 + (p[1] - sp2[1][1]) ** 2, 1.])
    for k in range(sample_points):
        t = float(k) / (sample_points - 1)
        i = 0
        while i == 0 or abs(f) > 0.000001 and i < 20:
            t2 = t ** 2
            t3 = t ** 3
            f = (ax * t3 + bx * t2 + cx * t + dx) * (3 * ax * t2 + 2 * bx * t + cx) + (ay * t3 + by * t2 + cy * t + dy) * (3 * ay * t2 + 2 * by * t + cy)
            df = (6 * ax * t + 2 * bx) * (ax * t3 + bx * t2 + cx * t + dx) + (3 * ax * t2 + 2 * bx * t + cx) ** 2 + (6 * ay * t + 2 * by) * (ay * t3 + by * t2 + cy * t + dy) + (3 * ay * t2 + 2 * by * t + cy) ** 2
            if df != 0:
                t = t - f / df
            else:
                break
            i += 1
        if 0 <= t <= 1:
            p1 = csp_at_t(sp1, sp2, t)
            d1 = (p1[0] - p[0]) ** 2 + (p1[1] - p[1]) ** 2
            if d1 < d[0]:
                d = [d1, t]
    return d


def csp_seg_to_csp_seg_distance(sp1, sp2, sp3, sp4, dist_bounds=(0, 1e100), sample_points=5, tolerance=.01):
    # check the ending points first
    Flast=0
    dist = csp_seg_to_point_distance(sp1, sp2, sp3[1], sample_points)
    dist += [0.]
    if dist[0] <= dist_bounds[0]:
        return dist
    d = csp_seg_to_point_distance(sp1, sp2, sp4[1], sample_points)
    if d[0] < dist[0]:
        dist = d + [1.]
        if dist[0] <= dist_bounds[0]:
            return dist
    d = csp_seg_to_point_distance(sp3, sp4, sp1[1], sample_points)
    if d[0] < dist[0]:
        dist = [d[0], 0., d[1]]
        if dist[0] <= dist_bounds[0]:
            return dist
    d = csp_seg_to_point_distance(sp3, sp4, sp2[1], sample_points)
    if d[0] < dist[0]:
        dist = [d[0], 1., d[1]]
        if dist[0] <= dist_bounds[0]:
            return dist
    sample_points -= 2
    if sample_points < 1:
        sample_points = 1
    ax1, ay1, bx1, by1, cx1, cy1, dx1, dy1 = csp_parameterize(sp1, sp2)
    ax2, ay2, bx2, by2, cx2, cy2, dx2, dy2 = csp_parameterize(sp3, sp4)
    #    try to find closes points using Newtons method
    for k in range(sample_points):
        for j in range(sample_points):
            t1 = float(k + 1) / (sample_points + 1)
            t2 = float(j) / (sample_points + 1)

            t12 = t1 * t1
            t13 = t1 * t1 * t1
            t22 = t2 * t2
            t23 = t2 * t2 * t2
            i = 0

            F1 = [0, 0]
            F2 = [[0, 0], [0, 0]]
            F = 1e100
            x = ax1 * t13 + bx1 * t12 + cx1 * t1 + dx1 - (ax2 * t23 + bx2 * t22 + cx2 * t2 + dx2)
            y = ay1 * t13 + by1 * t12 + cy1 * t1 + dy1 - (ay2 * t23 + by2 * t22 + cy2 * t2 + dy2)
            while i < 2 or abs(F - Flast) > tolerance and i < 30:
                f1x = 3 * ax1 * t12 + 2 * bx1 * t1 + cx1
                f1y = 3 * ay1 * t12 + 2 * by1 * t1 + cy1
                f2x = 3 * ax2 * t22 + 2 * bx2 * t2 + cx2
                f2y = 3 * ay2 * t22 + 2 * by2 * t2 + cy2
                F1[0] = 2 * f1x * x + 2 * f1y * y
                F1[1] = -2 * f2x * x - 2 * f2y * y
                F2[0][0] = 2 * (6 * ax1 * t1 + 2 * bx1) * x + 2 * f1x * f1x + 2 * (6 * ay1 * t1 + 2 * by1) * y + 2 * f1y * f1y
                F2[0][1] = -2 * f1x * f2x - 2 * f1y * f2y
                F2[1][0] = -2 * f2x * f1x - 2 * f2y * f1y
                F2[1][1] = -2 * (6 * ax2 * t2 + 2 * bx2) * x + 2 * f2x * f2x - 2 * (6 * ay2 * t2 + 2 * by2) * y + 2 * f2y * f2y
                F2 = inv_2x2(F2)
                if F2 is not None:
                    t1 -= (F2[0][0] * F1[0] + F2[0][1] * F1[1])
                    t2 -= (F2[1][0] * F1[0] + F2[1][1] * F1[1])
                    t12 = t1 * t1
                    t13 = t1 * t1 * t1
                    t22 = t2 * t2
                    t23 = t2 * t2 * t2
                    x = ax1 * t13 + bx1 * t12 + cx1 * t1 + dx1 - (ax2 * t23 + bx2 * t22 + cx2 * t2 + dx2)
                    y = ay1 * t13 + by1 * t12 + cy1 * t1 + dy1 - (ay2 * t23 + by2 * t22 + cy2 * t2 + dy2)
                    Flast = F
                    F = x * x + y * y
                else:
                    break
                i += 1
            if F < dist[0] and 0 <= t1 <= 1 and 0 <= t2 <= 1:
                dist = [F, t1, t2]
                if dist[0] <= dist_bounds[0]:
                    return dist
    return dist


def csp_to_csp_distance(csp1, csp2, dist_bounds=(0, 1e100), tolerance=.01):
    dist = [1e100, 0, 0, 0, 0, 0, 0]
    for i1 in range(len(csp1)):
        for j1 in range(1, len(csp1[i1])):
            for i2 in range(len(csp2)):
                for j2 in range(1, len(csp2[i2])):
                    d = csp_seg_bound_to_csp_seg_bound_max_min_distance(csp1[i1][j1 - 1], csp1[i1][j1], csp2[i2][j2 - 1], csp2[i2][j2])
                    if d[0] >= dist_bounds[1]:
                        continue
                    if d[1] < dist_bounds[0]:
                        return [d[1], i1, j1, 1, i2, j2, 1]
                    d = csp_seg_to_csp_seg_distance(csp1[i1][j1 - 1], csp1[i1][j1], csp2[i2][j2 - 1], csp2[i2][j2], dist_bounds, tolerance=tolerance)
                    if d[0] < dist[0]:
                        dist = [d[0], i1, j1, d[1], i2, j2, d[2]]
                    if dist[0] <= dist_bounds[0]:
                        return dist
            if dist[0] >= dist_bounds[1]:
                return dist
    return dist


def csp_split(sp1, sp2, t=.5):
    [x1, y1] = sp1[1]
    [x2, y2] = sp1[2]
    [x3, y3] = sp2[0]
    [x4, y4] = sp2[1]
    x12 = x1 + (x2 - x1) * t
    y12 = y1 + (y2 - y1) * t
    x23 = x2 + (x3 - x2) * t
    y23 = y2 + (y3 - y2) * t
    x34 = x3 + (x4 - x3) * t
    y34 = y3 + (y4 - y3) * t
    x1223 = x12 + (x23 - x12) * t
    y1223 = y12 + (y23 - y12) * t
    x2334 = x23 + (x34 - x23) * t
    y2334 = y23 + (y34 - y23) * t
    x = x1223 + (x2334 - x1223) * t
    y = y1223 + (y2334 - y1223) * t
    return [sp1[0], sp1[1], [x12, y12]], [[x1223, y1223], [x, y], [x2334, y2334]], [[x34, y34], sp2[1], sp2[2]]


def csp_true_bounds(csp):
    # Finds minx,miny,maxx,maxy of the csp and return their (x,y,i,j,t)
    minx = [float("inf"), 0, 0, 0]
    maxx = [float("-inf"), 0, 0, 0]
    miny = [float("inf"), 0, 0, 0]
    maxy = [float("-inf"), 0, 0, 0]
    for i in range(len(csp)):
        for j in range(1, len(csp[i])):
            ax, ay, bx, by, cx, cy, x0, y0 = bezierparameterize((csp[i][j - 1][1], csp[i][j - 1][2], csp[i][j][0], csp[i][j][1]))
            roots = cubic_solver(0, 3 * ax, 2 * bx, cx) + [0, 1]
            for root in roots:
                if type(root) is complex and abs(root.imag) < 1e-10:
                    root = root.real
                if type(root) is not complex and 0 <= root <= 1:
                    y = ay * (root ** 3) + by * (root ** 2) + cy * root + y0
                    x = ax * (root ** 3) + bx * (root ** 2) + cx * root + x0
                    maxx = max([x, y, i, j, root], maxx)
                    minx = min([x, y, i, j, root], minx)

            roots = cubic_solver(0, 3 * ay, 2 * by, cy) + [0, 1]
            for root in roots:
                if type(root) is complex and root.imag == 0:
                    root = root.real
                if type(root) is not complex and 0 <= root <= 1:
                    y = ay * (root ** 3) + by * (root ** 2) + cy * root + y0
                    x = ax * (root ** 3) + bx * (root ** 2) + cx * root + x0
                    maxy = max([y, x, i, j, root], maxy)
                    miny = min([y, x, i, j, root], miny)
    maxy[0], maxy[1] = maxy[1], maxy[0]
    miny[0], miny[1] = miny[1], miny[0]

    return minx, miny, maxx, maxy


############################################################################
# csp_segments_intersection(sp1,sp2,sp3,sp4)
#
# Returns array containing all intersections between two segments of cubic
# super path. Results are [ta,tb], or [ta0, ta1, tb0, tb1, "Overlap"]
# where ta, tb are values of t for the intersection point.
############################################################################
def csp_segments_intersection(sp1, sp2, sp3, sp4):
    a = csp_segment_to_bez(sp1, sp2)
    b = csp_segment_to_bez(sp3, sp4)

    def polish_intersection(a, b, ta, tb, tolerance=INTERSECTION_TOLERANCE):
        ax, ay, bx, by, cx, cy, dx, dy = bezierparameterize(a)
        ax1, ay1, bx1, by1, cx1, cy1, dx1, dy1 = bezierparameterize(b)
        i = 0
        F = [.0, .0]
        F1 = [[.0, .0], [.0, .0]]
        while i == 0 or (abs(F[0]) ** 2 + abs(F[1]) ** 2 > tolerance and i < 10):
            ta3 = ta ** 3
            ta2 = ta ** 2
            tb3 = tb ** 3
            tb2 = tb ** 2
            F[0] = ax * ta3 + bx * ta2 + cx * ta + dx - ax1 * tb3 - bx1 * tb2 - cx1 * tb - dx1
            F[1] = ay * ta3 + by * ta2 + cy * ta + dy - ay1 * tb3 - by1 * tb2 - cy1 * tb - dy1
            F1[0][0] = 3 * ax * ta2 + 2 * bx * ta + cx
            F1[0][1] = -3 * ax1 * tb2 - 2 * bx1 * tb - cx1
            F1[1][0] = 3 * ay * ta2 + 2 * by * ta + cy
            F1[1][1] = -3 * ay1 * tb2 - 2 * by1 * tb - cy1
            det = F1[0][0] * F1[1][1] - F1[0][1] * F1[1][0]
            if det != 0:
                F1 = [[F1[1][1] / det, -F1[0][1] / det], [-F1[1][0] / det, F1[0][0] / det]]
                ta = ta - (F1[0][0] * F[0] + F1[0][1] * F[1])
                tb = tb - (F1[1][0] * F[0] + F1[1][1] * F[1])
            else:
                break
            i += 1

        return ta, tb

    def recursion(a, b, ta0, ta1, tb0, tb1, depth_a, depth_b):
        global bezier_intersection_recursive_result
        if a == b:
            bezier_intersection_recursive_result += [[ta0, tb0, ta1, tb1, "Overlap"]]
            return
        tam = (ta0 + ta1) / 2
        tbm = (tb0 + tb1) / 2
        if depth_a > 0 and depth_b > 0:
            a1, a2 = bez_split(a, 0.5)
            b1, b2 = bez_split(b, 0.5)
            if bez_bounds_intersect(a1, b1):
                recursion(a1, b1, ta0, tam, tb0, tbm, depth_a - 1, depth_b - 1)
            if bez_bounds_intersect(a2, b1):
                recursion(a2, b1, tam, ta1, tb0, tbm, depth_a - 1, depth_b - 1)
            if bez_bounds_intersect(a1, b2):
                recursion(a1, b2, ta0, tam, tbm, tb1, depth_a - 1, depth_b - 1)
            if bez_bounds_intersect(a2, b2):
                recursion(a2, b2, tam, ta1, tbm, tb1, depth_a - 1, depth_b - 1)
        elif depth_a > 0:
            a1, a2 = bez_split(a, 0.5)
            if bez_bounds_intersect(a1, b):
                recursion(a1, b, ta0, tam, tb0, tb1, depth_a - 1, depth_b)
            if bez_bounds_intersect(a2, b):
                recursion(a2, b, tam, ta1, tb0, tb1, depth_a - 1, depth_b)
        elif depth_b > 0:
            b1, b2 = bez_split(b, 0.5)
            if bez_bounds_intersect(a, b1):
                recursion(a, b1, ta0, ta1, tb0, tbm, depth_a, depth_b - 1)
            if bez_bounds_intersect(a, b2):
                recursion(a, b2, ta0, ta1, tbm, tb1, depth_a, depth_b - 1)
        else:  # Both segments have been subdivided enough. Let's get some intersections :).
            intersection, t1, t2 = straight_segments_intersection([a[0]] + [a[3]], [b[0]] + [b[3]])
            if intersection:
                if intersection == "Overlap":
                    t1 = (max(0, min(1, t1[0])) + max(0, min(1, t1[1]))) / 2
                    t2 = (max(0, min(1, t2[0])) + max(0, min(1, t2[1]))) / 2
                bezier_intersection_recursive_result += [[ta0 + t1 * (ta1 - ta0), tb0 + t2 * (tb1 - tb0)]]

    global bezier_intersection_recursive_result
    bezier_intersection_recursive_result = []
    recursion(a, b, 0., 1., 0., 1., INTERSECTION_RECURSION_DEPTH, INTERSECTION_RECURSION_DEPTH)
    intersections = bezier_intersection_recursive_result
    for i in range(len(intersections)):
        if len(intersections[i]) < 5 or intersections[i][4] != "Overlap":
            intersections[i] = polish_intersection(a, b, intersections[i][0], intersections[i][1])
    return intersections


def csp_segments_true_intersection(sp1, sp2, sp3, sp4):
    intersections = csp_segments_intersection(sp1, sp2, sp3, sp4)
    res = []
    for intersection in intersections:
        if (
                (len(intersection) == 5 and intersection[4] == "Overlap" and (0 <= intersection[0] <= 1 or 0 <= intersection[1] <= 1) and (0 <= intersection[2] <= 1 or 0 <= intersection[3] <= 1))
                or (0 <= intersection[0] <= 1 and 0 <= intersection[1] <= 1)
        ):
            res += [intersection]
    return res


def csp_get_t_at_curvature(sp1, sp2, c, sample_points=16):
    # returns a list containing [t1,t2,t3,...,tn],  0<=ti<=1...
    if sample_points < 2:
        sample_points = 2
    tolerance = .0000000001
    res = []
    ax, ay, bx, by, cx, cy, dx, dy = csp_parameterize(sp1, sp2)
    for k in range(sample_points):
        t = float(k) / (sample_points - 1)
        i = 0
        F = 1e100
        while i < 2 or abs(F) > tolerance and i < 17:
            try:  # some numerical calculation could exceed the limits
                t2 = t * t
                # slopes...
                f1x = 3 * ax * t2 + 2 * bx * t + cx
                f1y = 3 * ay * t2 + 2 * by * t + cy
                f2x = 6 * ax * t + 2 * bx
                f2y = 6 * ay * t + 2 * by
                f3x = 6 * ax
                f3y = 6 * ay
                d = (f1x ** 2 + f1y ** 2) ** 1.5
                F1 = (
                        ((f1x * f3y - f3x * f1y) * d - (f1x * f2y - f2x * f1y) * 3. * (f2x * f1x + f2y * f1y) * ((f1x ** 2 + f1y ** 2) ** .5)) /
                        ((f1x ** 2 + f1y ** 2) ** 3)
                )
                F = (f1x * f2y - f1y * f2x) / d - c
                t -= F / F1
            except:
                break
            i += 1
        if 0 <= t <= 1 and F <= tolerance:
            if len(res) == 0:
                res.append(t)
            for i in res:
                if abs(t - i) <= 0.001:
                    break
            if not abs(t - i) <= 0.001:
                res.append(t)
    return res

def csp_max_curvature(sp1, sp2):
    ax, ay, bx, by, cx, cy, dx, dy = csp_parameterize(sp1, sp2)
    Flast=0
    tolerance = .0001
    F = 0.
    i = 0
    while i < 2 or F - Flast < tolerance and i < 10:
        t = .5
        f1x = 3 * ax * t ** 2 + 2 * bx * t + cx
        f1y = 3 * ay * t ** 2 + 2 * by * t + cy
        f2x = 6 * ax * t + 2 * bx
        f2y = 6 * ay * t + 2 * by
        f3x = 6 * ax
        f3y = 6 * ay
        d = pow(f1x ** 2 + f1y ** 2, 1.5)
        if d != 0:
            Flast = F
            F = (f1x * f2y - f1y * f2x) / d
            F1 = (
                    (d * (f1x * f3y - f3x * f1y) - (f1x * f2y - f2x * f1y) * 3. * (f2x * f1x + f2y * f1y) * pow(f1x ** 2 + f1y ** 2, .5)) /
                    (f1x ** 2 + f1y ** 2) ** 3
            )
            i += 1
            if F1 != 0:
                t -= F / F1
            else:
                break
        else:
            break
    return t


def csp_curvature_at_t(sp1, sp2, t, depth=3):
    ax, ay, bx, by, cx, cy, dx, dy = bezierparameterize(csp_segment_to_bez(sp1, sp2))

    # curvature = (x'y''-y'x'') / (x'^2+y'^2)^1.5

    f1x = 3 * ax * t ** 2 + 2 * bx * t + cx
    f1y = 3 * ay * t ** 2 + 2 * by * t + cy
    f2x = 6 * ax * t + 2 * bx
    f2y = 6 * ay * t + 2 * by
    d = (f1x ** 2 + f1y ** 2) ** 1.5
    if d != 0:
        return (f1x * f2y - f1y * f2x) / d
    else:
        t1 = f1x * f2y - f1y * f2x
        if t1 > 0:
            return 1e100
        if t1 < 0:
            return -1e100
        # Use the Lapitals rule to solve 0/0 problem for 2 times...
        t1 = 2 * (bx * ay - ax * by) * t + (ay * cx - ax * cy)
        if t1 > 0:
            return 1e100
        if t1 < 0:
            return -1e100
        t1 = bx * ay - ax * by
        if t1 > 0:
            return 1e100
        if t1 < 0:
            return -1e100
        if depth > 0:
            # little hack ;^) hope it won't influence anything...
            return csp_curvature_at_t(sp1, sp2, t * 1.004, depth - 1)
        return 1e100


def csp_subpath_ccw(subpath):
    # Remove all zero length segments
    s = 0
    if (P(subpath[-1][1]) - P(subpath[0][1])).l2() > 1e-10:
        subpath[-1][2] = subpath[-1][1]
        subpath[0][0] = subpath[0][1]
        subpath += [[subpath[0][1], subpath[0][1], subpath[0][1]]]
    pl = subpath[-1][2]
    for sp1 in subpath:
        for p in sp1:
            s += (p[0] - pl[0]) * (p[1] + pl[1])
            pl = p
    return s < 0


def csp_at_t(sp1, sp2, t):
    ax = sp1[1][0]
    bx = sp1[2][0]
    cx = sp2[0][0]
    dx = sp2[1][0]

    ay = sp1[1][1]
    by = sp1[2][1]
    cy = sp2[0][1]
    dy = sp2[1][1]

    x1 = ax + (bx - ax) * t
    y1 = ay + (by - ay) * t

    x2 = bx + (cx - bx) * t
    y2 = by + (cy - by) * t

    x3 = cx + (dx - cx) * t
    y3 = cy + (dy - cy) * t

    x4 = x1 + (x2 - x1) * t
    y4 = y1 + (y2 - y1) * t

    x5 = x2 + (x3 - x2) * t
    y5 = y2 + (y3 - y2) * t

    x = x4 + (x5 - x4) * t
    y = y4 + (y5 - y4) * t

    return [x, y]


def csp_at_length(sp1, sp2, l=0.5, tolerance=0.01):
    bez = (sp1[1][:], sp1[2][:], sp2[0][:], sp2[1][:])
    t = beziertatlength(bez, l, tolerance)
    return csp_at_t(sp1, sp2, t)


def cspseglength(sp1, sp2, tolerance=0.01):
    bez = (sp1[1][:], sp1[2][:], sp2[0][:], sp2[1][:])
    return bezierlength(bez, tolerance)


def csp_line_intersection(l1, l2, sp1, sp2):
    dd = l1[0]
    cc = l2[0] - l1[0]
    bb = l1[1]
    aa = l2[1] - l1[1]
    if aa == cc == 0:
        return []
    if aa:
        coef1 = cc / aa
        coef2 = 1
    else:
        coef1 = 1
        coef2 = aa / cc
    bez = (sp1[1][:], sp1[2][:], sp2[0][:], sp2[1][:])
    ax, ay, bx, by, cx, cy, x0, y0 = bezierparameterize(bez)
    a = coef1 * ay - coef2 * ax
    b = coef1 * by - coef2 * bx
    c = coef1 * cy - coef2 * cx
    d = coef1 * (y0 - bb) - coef2 * (x0 - dd)
    roots = cubic_solver(a, b, c, d)
    retval = []
    for i in roots:
        if type(i) is complex and abs(i.imag) < 1e-7:
            i = i.real
        if type(i) is not complex and -1e-10 <= i <= 1. + 1e-10:
            retval.append(i)
    return retval


def csp_split_by_two_points(sp1, sp2, t1, t2):
    if t1 > t2:
        t1, t2 = t2, t1
    if t1 == t2:
        sp1, sp2, sp3 = csp_split(sp1, sp2, t1)
        return [sp1, sp2, sp2, sp3]
    elif t1 <= 1e-10 and t2 >= 1. - 1e-10:
        return [sp1, sp1, sp2, sp2]
    elif t1 <= 1e-10:
        sp1, sp2, sp3 = csp_split(sp1, sp2, t2)
        return [sp1, sp1, sp2, sp3]
    elif t2 >= 1. - 1e-10:
        sp1, sp2, sp3 = csp_split(sp1, sp2, t1)
        return [sp1, sp2, sp3, sp3]
    else:
        sp1, sp2, sp3 = csp_split(sp1, sp2, t1)
        sp2, sp3, sp4 = csp_split(sp2, sp3, (t2 - t1) / (1 - t1))
        return [sp1, sp2, sp3, sp4]


def csp_seg_split(sp1, sp2, points):
    # points is float=t or list [t1, t2, ..., tn]
    if type(points) is float:
        points = [points]
    points.sort()
    res = [sp1, sp2]
    last_t = 0
    for t in points:
        if 1e-10 < t < 1. - 1e-10:
            sp3, sp4, sp5 = csp_split(res[-2], res[-1], (t - last_t) / (1 - last_t))
            last_t = t
            res[-2:] = [sp3, sp4, sp5]
    return res


def csp_subpath_split_by_points(subpath, points):
    # points are [[i,t]...] where i-segment's number
    points.sort()
    points = [[1, 0.]] + points + [[len(subpath) - 1, 1.]]
    parts = []
    for int1, int2 in zip(points, points[1:]):
        if int1 == int2:
            continue
        if int1[1] == 1.:
            int1[0] += 1
            int1[1] = 0.
        if int1 == int2:
            continue
        if int2[1] == 0.:
            int2[0] -= 1
            int2[1] = 1.
        if int1[0] == 0 and int2[0] == len(subpath) - 1:  # and small(int1[1]) and small(int2[1]-1) :
            continue
        if int1[0] == int2[0]:  # same segment
            sp = csp_split_by_two_points(subpath[int1[0] - 1], subpath[int1[0]], int1[1], int2[1])
            if sp[1] != sp[2]:
                parts += [[sp[1], sp[2]]]
        else:
            sp5, sp1, sp2 = csp_split(subpath[int1[0] - 1], subpath[int1[0]], int1[1])
            sp3, sp4, sp5 = csp_split(subpath[int2[0] - 1], subpath[int2[0]], int2[1])
            if int1[0] == int2[0] - 1:
                parts += [[sp1, [sp2[0], sp2[1], sp3[2]], sp4]]
            else:
                parts += [[sp1, sp2] + subpath[int1[0] + 1:int2[0] - 1] + [sp3, sp4]]
    return parts


def arc_from_s_r_n_l(s, r, n, l):
    if abs(n[0] ** 2 + n[1] ** 2 - 1) > 1e-10:
        n = normalize(n)
    return arc_from_c_s_l([s[0] + n[0] * r, s[1] + n[1] * r], s, l)


def arc_from_c_s_l(c, s, l):
    r = point_to_point_d(c, s)
    if r == 0:
        return []
    alpha = l / r
    cos_ = math.cos(alpha)
    sin_ = math.sin(alpha)
    e = [c[0] + (s[0] - c[0]) * cos_ - (s[1] - c[1]) * sin_, c[1] + (s[0] - c[0]) * sin_ + (s[1] - c[1]) * cos_]
    n = [c[0] - s[0], c[1] - s[1]]
    slope = rotate_cw(n) if l > 0 else rotate_ccw(n)
    return csp_from_arc(s, e, c, r, slope)


def csp_from_arc(start, end, center, r, slope_st):
    # Creates csp that approximise specified arc
    r = abs(r)
    alpha = (atan2(end[0] - center[0], end[1] - center[1]) - atan2(start[0] - center[0], start[1] - center[1])) % TAU

    sectors = int(abs(alpha) * 2 / math.pi) + 1
    alpha_start = atan2(start[0] - center[0], start[1] - center[1])
    cos_ = math.cos(alpha_start)
    sin_ = math.sin(alpha_start)
    k = (4. * math.tan(alpha / sectors / 4.) / 3.)
    if dot(slope_st, [- sin_ * k * r, cos_ * k * r]) < 0:
        if alpha > 0:
            alpha -= TAU
        else:
            alpha += TAU
    if abs(alpha * r) < 0.001:
        return []

    sectors = int(abs(alpha) * 2 / math.pi) + 1
    k = (4. * math.tan(alpha / sectors / 4.) / 3.)
    result = []
    for i in range(sectors + 1):
        cos_ = math.cos(alpha_start + alpha * i / sectors)
        sin_ = math.sin(alpha_start + alpha * i / sectors)
        sp = [[], [center[0] + cos_ * r, center[1] + sin_ * r], []]
        sp[0] = [sp[1][0] + sin_ * k * r, sp[1][1] - cos_ * k * r]
        sp[2] = [sp[1][0] - sin_ * k * r, sp[1][1] + cos_ * k * r]
        result += [sp]
    result[0][0] = result[0][1][:]
    result[-1][2] = result[-1][1]

    return result


def point_to_arc_distance(p, arc):
    # Distance calculattion from point to arc
    P0, P2, c, a = arc
    p = P(p)
    r = (P0 - c).mag()
    if r > 0:
        i = c + (p - c).unit() * r
        alpha = ((i - c).angle() - (P0 - c).angle())
        if a * alpha < 0:
            if alpha > 0:
                alpha = alpha - TAU
            else:
                alpha = TAU + alpha
        if between(alpha, 0, a) or min(abs(alpha), abs(alpha - a)) < STRAIGHT_TOLERANCE:
            return (p - i).mag(), [i.x, i.y]
        else:
            d1 = (p - P0).mag()
            d2 = (p - P2).mag()
            if d1 < d2:
                return d1, [P0.x, P0.y]
            else:
                return d2, [P2.x, P2.y]


def csp_to_arc_distance(sp1, sp2, arc1, arc2, tolerance=0.01):  # arc = [start,end,center,alpha]
    n = 10
    i = 0
    d = (0, [0, 0])
    d1 = (0, [0, 0])
    dl = 0
    while i < 1 or (abs(d1[0] - dl[0]) > tolerance and i < 4):
        i += 1
        dl = d1 * 1
        for j in range(n + 1):
            t = float(j) / n
            p = csp_at_t(sp1, sp2, t)
            d = min(point_to_arc_distance(p, arc1), point_to_arc_distance(p, arc2))
            d1 = max(d1, d)
        n = n * 2
    return d1[0]


def csp_point_inside_bound(sp1, sp2, p):
    bez = [sp1[1], sp1[2], sp2[0], sp2[1]]
    x, y = p
    c = 0
    # CLT added test of x in range
    xmin = 1e100
    xmax = -1e100
    for i in range(4):
        [x0, y0] = bez[i - 1]
        [x1, y1] = bez[i]
        xmin = min(xmin, x0)
        xmax = max(xmax, x0)
        if x0 - x1 != 0 and (y - y0) * (x1 - x0) >= (x - x0) * (y1 - y0) and x > min(x0, x1) and x <= max(x0, x1):
            c += 1
    return xmin <= x <= xmax and c % 2 == 0

def line_line_intersect(p1, p2, p3, p4):  # Return only true intersection.
    if (p1[0] == p2[0] and p1[1] == p2[1]) or (p3[0] == p4[0] and p3[1] == p4[1]):
        return False
    x = (p2[0] - p1[0]) * (p4[1] - p3[1]) - (p2[1] - p1[1]) * (p4[0] - p3[0])
    if x == 0:  # Lines are parallel
        if (p3[0] - p1[0]) * (p2[1] - p1[1]) == (p3[1] - p1[1]) * (p2[0] - p1[0]):
            if p3[0] != p4[0]:
                t11 = (p1[0] - p3[0]) / (p4[0] - p3[0])
                t12 = (p2[0] - p3[0]) / (p4[0] - p3[0])
                t21 = (p3[0] - p1[0]) / (p2[0] - p1[0])
                t22 = (p4[0] - p1[0]) / (p2[0] - p1[0])
            else:
                t11 = (p1[1] - p3[1]) / (p4[1] - p3[1])
                t12 = (p2[1] - p3[1]) / (p4[1] - p3[1])
                t21 = (p3[1] - p1[1]) / (p2[1] - p1[1])
                t22 = (p4[1] - p1[1]) / (p2[1] - p1[1])
            return "Overlap" if (0 <= t11 <= 1 or 0 <= t12 <= 1) and (0 <= t21 <= 1 or 0 <= t22 <= 1) else False
        else:
            return False
    else:
        return (
                0 <= ((p4[0] - p3[0]) * (p1[1] - p3[1]) - (p4[1] - p3[1]) * (p1[0] - p3[0])) / x <= 1 and
                0 <= ((p2[0] - p1[0]) * (p1[1] - p3[1]) - (p2[1] - p1[1]) * (p1[0] - p3[0])) / x <= 1)


def line_line_intersection_points(p1, p2, p3, p4):  # Return only points [ (x,y) ]
    if (p1[0] == p2[0] and p1[1] == p2[1]) or (p3[0] == p4[0] and p3[1] == p4[1]):
        return []
    x = (p2[0] - p1[0]) * (p4[1] - p3[1]) - (p2[1] - p1[1]) * (p4[0] - p3[0])
    if x == 0:  # Lines are parallel
        if (p3[0] - p1[0]) * (p2[1] - p1[1]) == (p3[1] - p1[1]) * (p2[0] - p1[0]):
            if p3[0] != p4[0]:
                t11 = (p1[0] - p3[0]) / (p4[0] - p3[0])
                t12 = (p2[0] - p3[0]) / (p4[0] - p3[0])
                t21 = (p3[0] - p1[0]) / (p2[0] - p1[0])
                t22 = (p4[0] - p1[0]) / (p2[0] - p1[0])
            else:
                t11 = (p1[1] - p3[1]) / (p4[1] - p3[1])
                t12 = (p2[1] - p3[1]) / (p4[1] - p3[1])
                t21 = (p3[1] - p1[1]) / (p2[1] - p1[1])
                t22 = (p4[1] - p1[1]) / (p2[1] - p1[1])
            res = []
            if (0 <= t11 <= 1 or 0 <= t12 <= 1) and (0 <= t21 <= 1 or 0 <= t22 <= 1):
                if 0 <= t11 <= 1:
                    res += [p1]
                if 0 <= t12 <= 1:
                    res += [p2]
                if 0 <= t21 <= 1:
                    res += [p3]
                if 0 <= t22 <= 1:
                    res += [p4]
            return res
        else:
            return []
    else:
        t1 = ((p4[0] - p3[0]) * (p1[1] - p3[1]) - (p4[1] - p3[1]) * (p1[0] - p3[0])) / x
        t2 = ((p2[0] - p1[0]) * (p1[1] - p3[1]) - (p2[1] - p1[1]) * (p1[0] - p3[0])) / x
        if 0 <= t1 <= 1 and 0 <= t2 <= 1:
            return [[p1[0] * (1 - t1) + p2[0] * t1, p1[1] * (1 - t1) + p2[1] * t1]]
        else:
            return []


def point_to_point_d2(a, b):
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2


def point_to_point_d(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)


def point_to_line_segment_distance_2(p1, p2, p3):
    # p1 - point, p2,p3 - line segment
    # draw_pointer(p1)
    w0 = [p1[0] - p2[0], p1[1] - p2[1]]
    v = [p3[0] - p2[0], p3[1] - p2[1]]
    c1 = w0[0] * v[0] + w0[1] * v[1]
    if c1 <= 0:
        return w0[0] * w0[0] + w0[1] * w0[1]
    c2 = v[0] * v[0] + v[1] * v[1]
    if c2 <= c1:
        return (p1[0] - p3[0]) ** 2 + (p1[1] - p3[1]) ** 2
    return (p1[0] - p2[0] - v[0] * c1 / c2) ** 2 + (p1[1] - p2[1] - v[1] * c1 / c2)


def line_to_line_distance_2(p1, p2, p3, p4):
    if line_line_intersect(p1, p2, p3, p4):
        return 0
    return min(
            point_to_line_segment_distance_2(p1, p3, p4),
            point_to_line_segment_distance_2(p2, p3, p4),
            point_to_line_segment_distance_2(p3, p1, p2),
            point_to_line_segment_distance_2(p4, p1, p2))


def csp_seg_bound_to_csp_seg_bound_max_min_distance(sp1, sp2, sp3, sp4):
    bez1 = csp_segment_to_bez(sp1, sp2)
    bez2 = csp_segment_to_bez(sp3, sp4)
    min_dist = 1e100
    max_dist = 0.
    for i in range(4):
        if csp_point_inside_bound(sp1, sp2, bez2[i]) or csp_point_inside_bound(sp3, sp4, bez1[i]):
            min_dist = 0.
            break
    for i in range(4):
        for j in range(4):
            d = line_to_line_distance_2(bez1[i - 1], bez1[i], bez2[j - 1], bez2[j])
            if d < min_dist:
                min_dist = d
            d = (bez2[j][0] - bez1[i][0]) ** 2 + (bez2[j][1] - bez1[i][1]) ** 2
            if max_dist < d:
                max_dist = d
    return min_dist, max_dist


def csp_reverse(csp):
    for i in range(len(csp)):
        n = []
        for j in csp[i]:
            n = [[j[2][:], j[1][:], j[0][:]]] + n
        csp[i] = n[:]
    return csp


def csp_normalized_slope(sp1, sp2, t):
    ax, ay, bx, by, cx, cy, dx, dy = bezierparameterize((sp1[1][:], sp1[2][:], sp2[0][:], sp2[1][:]))
    if sp1[1] == sp2[1] == sp1[2] == sp2[0]:
        return [1., 0.]
    f1x = 3 * ax * t * t + 2 * bx * t + cx
    f1y = 3 * ay * t * t + 2 * by * t + cy
    if abs(f1x * f1x + f1y * f1y) > 1e-9:  # LT changed this from 1e-20, which caused problems
        l = math.sqrt(f1x * f1x + f1y * f1y)
        return [f1x / l, f1y / l]

    if t == 0:
        f1x = sp2[0][0] - sp1[1][0]
        f1y = sp2[0][1] - sp1[1][1]
        if abs(f1x * f1x + f1y * f1y) > 1e-9:  # LT changed this from 1e-20, which caused problems
            l = math.sqrt(f1x * f1x + f1y * f1y)
            return [f1x / l, f1y / l]
        else:
            f1x = sp2[1][0] - sp1[1][0]
            f1y = sp2[1][1] - sp1[1][1]
            if f1x * f1x + f1y * f1y != 0:
                l = math.sqrt(f1x * f1x + f1y * f1y)
                return [f1x / l, f1y / l]
    elif t == 1:
        f1x = sp2[1][0] - sp1[2][0]
        f1y = sp2[1][1] - sp1[2][1]
        if abs(f1x * f1x + f1y * f1y) > 1e-9:
            l = math.sqrt(f1x * f1x + f1y * f1y)
            return [f1x / l, f1y / l]
        else:
            f1x = sp2[1][0] - sp1[1][0]
            f1y = sp2[1][1] - sp1[1][1]
            if f1x * f1x + f1y * f1y != 0:
                l = math.sqrt(f1x * f1x + f1y * f1y)
                return [f1x / l, f1y / l]
    else:
        return [1., 0.]


def csp_normalized_normal(sp1, sp2, t):
    nx, ny = csp_normalized_slope(sp1, sp2, t)
    return [-ny, nx]


def csp_parameterize(sp1, sp2):
    return bezierparameterize(csp_segment_to_bez(sp1, sp2))


def csp_concat_subpaths(*s):
    def concat(s1, s2):
        if not s1:
            return s2
        if not s2:
            return s1
        if (s1[-1][1][0] - s2[0][1][0]) ** 2 + (s1[-1][1][1] - s2[0][1][1]) ** 2 > 0.00001:
            return s1[:-1] + [[s1[-1][0], s1[-1][1], s1[-1][1]], [s2[0][1], s2[0][1], s2[0][2]]] + s2[1:]
        else:
            return s1[:-1] + [[s1[-1][0], s2[0][1], s2[0][2]]] + s2[1:]

    if len(s) == 0:
        return []
    if len(s) == 1:
        return s[0]
    result = s[0]
    for s1 in s[1:]:
        result = concat(result, s1)
    return result


def csp_subpaths_end_to_start_distance2(s1, s2):
    return (s1[-1][1][0] - s2[0][1][0]) ** 2 + (s1[-1][1][1] - s2[0][1][1]) ** 2


def csp_clip_by_line(csp, l1, l2):
    result = []
    for i in range(len(csp)):
        s = csp[i]
        intersections = []
        for j in range(1, len(s)):
            intersections += [[j, int_] for int_ in csp_line_intersection(l1, l2, s[j - 1], s[j])]
        splitted_s = csp_subpath_split_by_points(s, intersections)
        for s in splitted_s[:]:
            clip = False
            for p in csp_true_bounds([s]):
                if (l1[1] - l2[1]) * p[0] + (l2[0] - l1[0]) * p[1] + (l1[0] * l2[1] - l2[0] * l1[1]) < -0.01:
                    clip = True
                    break
            if clip:
                splitted_s.remove(s)
        result += splitted_s
    return result


def csp_subpath_line_to(subpath, points, prepend=False):
    # Appends subpath with line or polyline.
    if len(points) > 0:
        if not prepend:
            if len(subpath) > 0:
                subpath[-1][2] = subpath[-1][1][:]
            if type(points[0]) == type([1, 1]):
                for p in points:
                    subpath += [[p[:], p[:], p[:]]]
            else:
                subpath += [[points, points, points]]
        else:
            if len(subpath) > 0:
                subpath[0][0] = subpath[0][1][:]
            if type(points[0]) == type([1, 1]):
                for p in points:
                    subpath = [[p[:], p[:], p[:]]] + subpath
            else:
                subpath = [[points, points, points]] + subpath
    return subpath


def csp_join_subpaths(csp):
    result = csp[:]
    done_smf = True
    joined_result = []
    while done_smf:
        done_smf = False
        while len(result) > 0:
            s1 = result[-1][:]
            del (result[-1])
            j = 0
            joined_smf = False
            while j < len(joined_result):
                if csp_subpaths_end_to_start_distance2(joined_result[j], s1) < 0.000001:
                    joined_result[j] = csp_concat_subpaths(joined_result[j], s1)
                    done_smf = True
                    joined_smf = True
                    break
                if csp_subpaths_end_to_start_distance2(s1, joined_result[j]) < 0.000001:
                    joined_result[j] = csp_concat_subpaths(s1, joined_result[j])
                    done_smf = True
                    joined_smf = True
                    break
                j += 1
            if not joined_smf:
                joined_result += [s1[:]]
        if done_smf:
            result = joined_result[:]
            joined_result = []
    return joined_result


def triangle_cross(a, b, c):
    return (a[0] - b[0]) * (c[1] - b[1]) - (c[0] - b[0]) * (a[1] - b[1])


def csp_segment_convex_hull(sp1, sp2):
    a = sp1[1][:]
    b = sp1[2][:]
    c = sp2[0][:]
    d = sp2[1][:]

    abc = triangle_cross(a, b, c)
    abd = triangle_cross(a, b, d)
    bcd = triangle_cross(b, c, d)
    cad = triangle_cross(c, a, d)
    if abc == 0 and abd == 0:
        return [min(a, b, c, d), max(a, b, c, d)]
    if abc == 0:
        return [d, min(a, b, c), max(a, b, c)]
    if abd == 0:
        return [c, min(a, b, d), max(a, b, d)]
    if bcd == 0:
        return [a, min(b, c, d), max(b, c, d)]
    if cad == 0:
        return [b, min(c, a, d), max(c, a, d)]

    m1 = abc * abd > 0
    m2 = abc * bcd > 0
    m3 = abc * cad > 0

    if m1 and m2 and m3:
        return [a, b, c]
    if m1 and m2 and not m3:
        return [a, b, c, d]
    if m1 and not m2 and m3:
        return [a, b, d, c]
    if not m1 and m2 and m3:
        return [a, d, b, c]
    if m1 and not (m2 and m3):
        return [a, b, d]
    if not (m1 and m2) and m3:
        return [c, a, d]
    if not (m1 and m3) and m2:
        return [b, c, d]

    raise ValueError("csp_segment_convex_hull happened which is something that shouldn't happen!")


################################################################################
# Bezier additional functions
################################################################################

def bez_bounds_intersect(bez1, bez2):
    return bounds_intersect(bez_bound(bez2), bez_bound(bez1))


def bez_bound(bez):
    return [
        min(bez[0][0], bez[1][0], bez[2][0], bez[3][0]),
        min(bez[0][1], bez[1][1], bez[2][1], bez[3][1]),
        max(bez[0][0], bez[1][0], bez[2][0], bez[3][0]),
        max(bez[0][1], bez[1][1], bez[2][1], bez[3][1]),
    ]


def bounds_intersect(a, b):
    return not ((a[0] > b[2]) or (b[0] > a[2]) or (a[1] > b[3]) or (b[1] > a[3]))


def tpoint(xy1, xy2, t):
    (x1, y1) = xy1
    (x2, y2) = xy2
    return [x1 + t * (x2 - x1), y1 + t * (y2 - y1)]


def bez_split(a, t=0.5):
    a1 = tpoint(a[0], a[1], t)
    at = tpoint(a[1], a[2], t)
    b2 = tpoint(a[2], a[3], t)
    a2 = tpoint(a1, at, t)
    b1 = tpoint(b2, at, t)
    a3 = tpoint(a2, b1, t)
    return [a[0], a1, a2, a3], [a3, b1, b2, a[3]]


################################################################################
# Some vector functions
################################################################################

def normalize(xy):
    (x, y) = xy
    l = math.sqrt(x ** 2 + y ** 2)
    if l == 0:
        return [0., 0.]
    else:
        return [x / l, y / l]


def cross(a, b):
    return a[1] * b[0] - a[0] * b[1]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def rotate_ccw(d):
    return [-d[1], d[0]]


def rotate_cw(d):
    return [d[1], -d[0]]


def vectors_ccw(a, b):
    return a[0] * b[1] - b[0] * a[1] < 0


################################################################################
# Common functions
################################################################################

def inv_2x2(a):  # invert matrix 2x2
    det = a[0][0] * a[1][1] - a[1][0] * a[0][1]
    if det == 0:
        return None
    return [
        [a[1][1] / det, -a[0][1] / det],
        [-a[1][0] / det, a[0][0] / det]
    ]


def small(a):
    global small_tolerance
    return abs(a) < small_tolerance


def atan2(*arg):
    if len(arg) == 1 and (type(arg[0]) == type([0., 0.]) or type(arg[0]) == type((0., 0.))):
        return (math.pi / 2 - math.atan2(arg[0][0], arg[0][1])) % TAU
    elif len(arg) == 2:
        return (math.pi / 2 - math.atan2(arg[0], arg[1])) % TAU
    else:
        raise ValueError("Bad argumets for atan! ({})".format(*arg))


def draw_text(text, x, y, group=None, style=None, font_size=10, gcodetools_tag=None,label_text=None):
    if style is None:
        style = "font-family:DejaVu Sans;font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;font-family:DejaVu Sans;fill:#000000;fill-opacity:1;stroke:none;"
    style += "font-size:{:f}px;".format(font_size)
    attributes = {'x': str(x), 'y': str(y), 'style': style}
    if gcodetools_tag is not None:
        attributes[options.self.Zone_name['Definition']] = str(gcodetools_tag)

    if group is None:
        group = options.doc_root

    text_elem = group.add(TextElement(**attributes))
    if label_text != None:
        text_elem.label=label_text
    text_elem.set("xml:space", "preserve")
    text = str(text).split("\n")
    for string in text:
        span = text_elem.add(Tspan(x=str(x), y=str(y)))
        span.set('sodipodi:role', 'line')
        y += font_size
        span.text = str(string)



def draw_csp(csp, stroke="#f00", fill="none", comment="", width=0.354, group=None, style=None):
    if group is None:
        group = options.doc_root
    node = group.add(PathElement())

    node.style = style if style is not None else \
        {'fill': fill, 'fill-opacity': 1, 'stroke': stroke, 'stroke-width': width}

    node.path = CubicSuperPath(csp)

    if comment != '':
        node.set('comment', comment)

    return node


def draw_pointer(x, color="#f00", figure="cross", group=None, comment="", fill=None, width=.1, size=10., text=None, font_size=None, pointer_type=None, attrib=None):
    size = size / 2
    if attrib is None:
        attrib = {}
    if pointer_type is None:
        pointer_type = "Pointer"
    attrib[options.self.Zone_name['Definition']] = pointer_type
    if group is None:
        group = options.self.svg.get_current_layer()
    if text is not None:
        if font_size is None:
            font_size = 7
        group = group.add(Group(gcodetools=pointer_type + " group"))
        draw_text(text, x[0] + size * 2.2, x[1] - size, group=group, font_size=font_size)
    if figure == "line":
        s = ""
        for i in range(1, len(x) / 2):
            s += " {}, {} ".format(x[i * 2], x[i * 2 + 1])
        attrib.update({"d": "M {},{} L {}".format(x[0], x[1], s), "style": "fill:none;stroke:{};stroke-width:{:f};".format(color, width), "comment": str(comment)})
    elif figure == "arrow":
        if fill is None:
            fill = "#12b3ff"
        fill_opacity = "0.8"
        d = "m {},{} ".format(x[0], x[1]) + re.sub("([0-9\\-.e]+)", (lambda match: str(float(match.group(1)) * size * 2.)), "0.88464,-0.40404 c -0.0987,-0.0162 -0.186549,-0.0589 -0.26147,-0.1173 l 0.357342,-0.35625 c 0.04631,-0.039 0.0031,-0.13174 -0.05665,-0.12164 -0.0029,-1.4e-4 -0.0058,-1.4e-4 -0.0087,0 l -2.2e-5,2e-5 c -0.01189,0.004 -0.02257,0.0119 -0.0305,0.0217 l -0.357342,0.35625 c -0.05818,-0.0743 -0.102813,-0.16338 -0.117662,-0.26067 l -0.409636,0.88193 z")
        attrib.update({"d": d, "style": "fill:{};stroke:none;fill-opacity:{};".format(fill, fill_opacity), "comment": str(comment)})
    else:
        attrib.update({"d": "m {},{} l {:f},{:f} {:f},{:f} {:f},{:f} {:f},{:f} , {:f},{:f}".format(x[0], x[1], size, size, -2 * size, -2 * size, size, size, size, -size, -2 * size, 2 * size), "style": "fill:none;stroke:{};stroke-width:{:f};".format(color, width), "comment": str(comment)})
    group.add(PathElement(**attrib))


def straight_segments_intersection(a, b, true_intersection=True):  # (True intersection means check ta and tb are in [0,1])
    ax = a[0][0]
    bx = a[1][0]
    cx = b[0][0]
    dx = b[1][0]
    ay = a[0][1]
    by = a[1][1]
    cy = b[0][1]
    dy = b[1][1]
    if (ax == bx and ay == by) or (cx == dx and cy == dy):
        return False, 0, 0
    if (bx - ax) * (dy - cy) - (by - ay) * (dx - cx) == 0:  # Lines are parallel
        ta = (ax - cx) / (dx - cx) if cx != dx else (ay - cy) / (dy - cy)
        tb = (bx - cx) / (dx - cx) if cx != dx else (by - cy) / (dy - cy)
        tc = (cx - ax) / (bx - ax) if ax != bx else (cy - ay) / (by - ay)
        td = (dx - ax) / (bx - ax) if ax != bx else (dy - ay) / (by - ay)
        return ("Overlap" if 0 <= ta <= 1 or 0 <= tb <= 1 or 0 <= tc <= 1 or 0 <= td <= 1 or not true_intersection else False), (ta, tb), (tc, td)
    else:
        ta = ((ay - cy) * (dx - cx) - (ax - cx) * (dy - cy)) / ((bx - ax) * (dy - cy) - (by - ay) * (dx - cx))
        tb = (ax - cx + ta * (bx - ax)) / (dx - cx) if dx != cx else (ay - cy + ta * (by - ay)) / (dy - cy)
        return (0 <= ta <= 1 and 0 <= tb <= 1 or not true_intersection), ta, tb


def between(c, x, y):
    return x - STRAIGHT_TOLERANCE <= c <= y + STRAIGHT_TOLERANCE or y - STRAIGHT_TOLERANCE <= c <= x + STRAIGHT_TOLERANCE


def cubic_solver_real(a, b, c, d):
    # returns only real roots of a cubic equation.
    roots = cubic_solver(a, b, c, d)
    res = []
    for root in roots:
        if type(root) is complex:
            if -1e-10 < root.imag < 1e-10:
                res.append(root.real)
        else:
            res.append(root)
    return res


def cubic_solver(a, b, c, d):
    if a != 0:
        #    Monics formula see http://en.wikipedia.org/wiki/Cubic_function#Monic_formula_of_roots
        a, b, c = (b / a, c / a, d / a)
        m = 2 * a ** 3 - 9 * a * b + 27 * c
        k = a ** 2 - 3 * b
        n = m ** 2 - 4 * k ** 3
        w1 = -.5 + .5 * cmath.sqrt(3) * 1j
        w2 = -.5 - .5 * cmath.sqrt(3) * 1j
        if n >= 0:
            t = m + math.sqrt(n)
            m1 = pow(t / 2, 1. / 3) if t >= 0 else -pow(-t / 2, 1. / 3)
            t = m - math.sqrt(n)
            n1 = pow(t / 2, 1. / 3) if t >= 0 else -pow(-t / 2, 1. / 3)
        else:
            m1 = pow(complex((m + cmath.sqrt(n)) / 2), 1. / 3)
            n1 = pow(complex((m - cmath.sqrt(n)) / 2), 1. / 3)
        x1 = -1. / 3 * (a + m1 + n1)
        x2 = -1. / 3 * (a + w1 * m1 + w2 * n1)
        x3 = -1. / 3 * (a + w2 * m1 + w1 * n1)
        return [x1, x2, x3]
    elif b != 0:
        det = c ** 2 - 4 * b * d
        if det > 0:
            return [(-c + math.sqrt(det)) / (2 * b), (-c - math.sqrt(det)) / (2 * b)]
        elif d == 0:
            return [-c / (b * b)]
        else:
            return [(-c + cmath.sqrt(det)) / (2 * b), (-c - cmath.sqrt(det)) / (2 * b)]
    elif c != 0:
        return [-d / c]
    else:
        return []


################################################################################
# Point (x,y) operations
################################################################################
class P(object):
    def __init__(self, x, y=None):
        if not y is None:
            self.x = float(x)
            self.y = float(y)
        else:
            self.x = float(x[0])
            self.y = float(x[1])

    def __add__(self, other):
        return P(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return P(self.x - other.x, self.y - other.y)

    def __neg__(self):
        return P(-self.x, -self.y)

    def __mul__(self, other):
        if isinstance(other, P):
            return self.x * other.x + self.y * other.y
        return P(self.x * other, self.y * other)

    __rmul__ = __mul__

    def __div__(self, other):
        return P(self.x / other, self.y / other)

    def __truediv__(self, other):
        return self.__div__(other)

    def mag(self):
        return math.hypot(self.x, self.y)

    def unit(self):
        h_mag = self.mag()
        if h_mag:
            return self / h_mag
        return P(0, 0)

    def dot(self, other):
        return self.x * other.x + self.y * other.y

    def rot(self, theta):
        c = math.cos(theta)
        s = math.sin(theta)
        return P(self.x * c - self.y * s, self.x * s + self.y * c)

    def angle(self):
        return math.atan2(self.y, self.x)

    def __repr__(self):
        return '{:f},{:f}'.format(self.x, self.y)

    def pr(self):
        return "{:.2f},{:.2f}".format(self.x, self.y)

    def to_list(self):
        return [self.x, self.y]

    def ccw(self):
        return P(-self.y, self.x)

    def l2(self):
        return self.x * self.x + self.y * self.y


#point de départ    
if __name__ == '__main__':
    Path2Laser().run()
    
#a traiter    
#     actif Decoupe
# actif Dessin
# /home/David/.var/app/org.inkscape.Inkscape/config/inkscape/extensions/path2laser_gcode.py:441: DeprecationWarning: inkex.deprecated.main.transform_mul -> Use @ operator instead
#   trans = (Transform(t) * Transform(trans)).matrix if trans != [] else t
# actif Gravure
