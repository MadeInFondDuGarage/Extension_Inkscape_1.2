#!/usr/bin/env python
'''

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
'''
'''
ce programme a été réaliser dans le but de creer des gabarits pour la realisation de patron
fait par D.Vantieghem aka madeinfonddugarage sur youtube
'''

import inkex,  math
from inkex import Circle,Rectangle, bezier
import re
import cmath
from lxml import etree


class P:
    def __init__(self, x, y=None):
        if not y==None:
            self.x, self.y = float(x), float(y)
        else:
            self.x, self.y = float(x[0]), float(x[1])
    def __sub__(self, other): return P(self.x - other.x, self.y - other.y)
    def l2(self): return self.x*self.x + self.y*self.y

def csp_subpath_ccw(subpath):
    # Remove all zerro length segments
    s = 0
    #subpath = subpath[:]
    if (P(subpath[-1][1])-P(subpath[0][1])).l2() > 1e-10 :
        subpath[-1][2] = subpath[-1][1]
        subpath[0][0] = subpath[0][1]
        subpath += [ [subpath[0][1],subpath[0][1],subpath[0][1]]  ]
    pl = subpath[-1][2]
    for sp1 in subpath:
        for p in sp1 :
            s += (p[0]-pl[0])*(p[1]+pl[1])
            pl = p
    return s<0

def cubic_solver(a,b,c,d):
    if a!=0:
        #   Monics formula see http://en.wikipedia.org/wiki/Cubic_function#Monic_formula_of_roots
        a,b,c = (b/a, c/a, d/a)
        m = 2*a**3 - 9*a*b + 27*c
        k = a**2 - 3*b
        n = m**2 - 4*k**3
        w1 = -.5 + .5*cmath.sqrt(3)*1j
        w2 = -.5 - .5*cmath.sqrt(3)*1j
        if n>=0 :
            t = m+math.sqrt(n)
            m1 = pow(t/2,1./3) if t>=0 else -pow(-t/2,1./3)
            t = m-math.sqrt(n)
            n1 = pow(t/2,1./3) if t>=0 else -pow(-t/2,1./3)
        else :
            m1 = pow(complex((m+cmath.sqrt(n))/2),1./3)
            n1 = pow(complex((m-cmath.sqrt(n))/2),1./3)
        x1 = -1./3 * (a + m1 + n1)
        x2 = -1./3 * (a + w1*m1 + w2*n1)
        x3 = -1./3 * (a + w2*m1 + w1*n1)
        return [x1,x2,x3]
    elif b!=0:
        det = c**2-4*b*d
        if det>0 :
            return [(-c+math.sqrt(det))/(2*b),(-c-math.sqrt(det))/(2*b)]
        elif d == 0 :
            return [-c/(b*b)]
        else :
            return [(-c+cmath.sqrt(det))/(2*b),(-c-cmath.sqrt(det))/(2*b)]
    elif c!=0 :
        return [-d/c]
    else : return []

def csp_true_bounds(csp) :
    # Finds minx,miny,maxx,maxy of the csp and return their (x,y,i,j,t)
    minx = [float("inf"), 0, 0, 0]
    maxx = [float("-inf"), 0, 0, 0]
    miny = [float("inf"), 0, 0, 0]
    maxy = [float("-inf"), 0, 0, 0]
    for i in range(len(csp)):
        for j in range(1,len(csp[i])):
            ax,ay,bx,by,cx,cy,x0,y0 = bezier.bezierparameterize((csp[i][j-1][1],csp[i][j-1][2],csp[i][j][0],csp[i][j][1]))

            roots = cubic_solver(0, 3*ax, 2*bx, cx)  + [0,1]
            for root in roots :
                if type(root) is complex and abs(root.imag)<1e-10:
                    root = root.real
                if type(root) is not complex and 0<=root<=1:
                    y = ay*(root**3)+by*(root**2)+cy*root+y0
                    x = ax*(root**3)+bx*(root**2)+cx*root+x0
                    maxx = max([x,y,i,j,root],maxx)
                    minx = min([x,y,i,j,root],minx)
            roots = cubic_solver(0, 3*ay, 2*by, cy)  + [0,1]
            for root in roots :
                if type(root) is complex and root.imag==0:
                    root = root.real
                if type(root) is not complex and 0<=root<=1:
                    y = ay*(root**3)+by*(root**2)+cy*root+y0
                    x = ax*(root**3)+bx*(root**2)+cx*root+x0
                    maxy = max([y,x,i,j,root],maxy)
                    miny = min([y,x,i,j,root],miny)
    maxy[0],maxy[1] = maxy[1],maxy[0]
    miny[0],miny[1] = miny[1],miny[0]

    return minx,miny,maxx,maxy

def csp_split(sp1,sp2,t=.5) :
    [x1,y1],[x2,y2],[x3,y3],[x4,y4] = sp1[1], sp1[2], sp2[0], sp2[1]
    x12 = x1+(x2-x1)*t
    y12 = y1+(y2-y1)*t
    x23 = x2+(x3-x2)*t
    y23 = y2+(y3-y2)*t
    x34 = x3+(x4-x3)*t
    y34 = y3+(y4-y3)*t
    x1223 = x12+(x23-x12)*t
    y1223 = y12+(y23-y12)*t
    x2334 = x23+(x34-x23)*t
    y2334 = y23+(y34-y23)*t
    x = x1223+(x2334-x1223)*t
    y = y1223+(y2334-y1223)*t
    return [sp1[0],sp1[1],[x12,y12]], [[x1223,y1223],[x,y],[x2334,y2334]], [[x34,y34],sp2[1],sp2[2]]

class Gabarits(inkex.Effect):
    def __init__(self):
        inkex.Effect.__init__(self)
        self.arg_parser.add_argument("-d", "--dotsize",dest="dotsize", default="10mm",help="Size of the dots placed at path nodes")
        self.arg_parser.add_argument("--stylegab", dest="stylegab", default="Gabcarrer",help="Type de Gabarit")
        self.arg_parser.add_argument("--tab", dest="tab",help="The selected UI-tab when OK was pressed")

    def addDot(self, MonChemin, carrer_cercle):
        #je cre le groupe qui acceuille les point (cercle ou carrer)
        NouveauGroup = MonChemin.getparent().add(inkex.Group())
        #je le nomme
        NouveauGroup.label = 'gabarit: ' + self.options.dotsize +" " + carrer_cercle 
        # je prédéfinie la taille a partir des infos du formulaire
        Dot_Size=self.svg.unittouu(self.options.dotsize)
        try:
            t = MonChemin.get('transform') #le node est la courbe selectionner, je récupère le 'transform' si il existe
            MonChemin.set('transform', t) #et je l'affecte au nouveau groupe
        except:
            pass
        
        #je definie le style
        style = inkex.Style({ 'stroke': '#00ff00', 'fill': 'none' ,'stroke-width': str(self.svg.unittouu('1px'))})
        #je récupère le path et le transforme pour pouvoir le modifier 
        path_trans_applied = MonChemin.path.transform(MonChemin.composed_transform())
        
        #je cree un liste vide
        MonCheminList=list() #init d'une list
        
        #je lance la creation des points
        for step, (x, y) in enumerate(path_trans_applied.end_points):
            MonCheminList.append((x,y)) #jenregistre les points 
            if (MonCheminList[0] != MonCheminList[step]) or step ==0: 
            #je verifie que le point courant n'est pas identique au premier point. sinon je ne fait le dessin
                if carrer_cercle=="Carrer":
                    rectangle = NouveauGroup.add(Rectangle(x=str(x-Dot_Size),y=str(y-Dot_Size),width=str(Dot_Size*2),height=str(Dot_Size*2)))
                    rectangle.style=style
                else:
                    circle = NouveauGroup.add(Circle(cx=str(x), cy=str(y),r=str(Dot_Size)))
                    circle.style = style
                    
    def effect(self):
        #je verifie que j'ai au moins 1 path de sélectionner
        if len(self.options.ids)<=0:
            raise inkex.AbortExtension("This extension requires at least one selected path.")
            
        #je récupère le chemin a travers un filtre pour etre sur que c'est un path
        MaSelection = self.svg.selection.filter(inkex.PathElement)
        
        #si il est vide j'envoie une erreur
        if not MaSelection:
            raise inkex.AbortExtension("This extension need a path, not groups.")
                    
        #si c'est bon j'envoie les transformation
        if self.options.stylegab == "Carrer" or self.options.stylegab =="Cercle":
            for MonChemin in MaSelection:    
                self.addDot(MonChemin, self.options.stylegab)
                    
        else:
            for MonChemin in MaSelection:
                #je defini si j'ai choisi l'intérieur ou l'extérieur
                if self.options.stylegab == "Interieur":
                    in_out=-1
                else:
                    in_out=1
                
                #je récupère la structure du chemin
                d = MonChemin.get('d')
                
                csp = inkex.paths.CubicSuperPath(d)
                
                #le code qui suis est issu d'un script pour réaliser un décallage, puis pour faire une courbe en mode offset dynamique
                min_x,min_y,min_i,min_j,min_t = csp_true_bounds(csp)[1]

                if min_y!=float("-inf"):
                    subp = csp[min_i]
                    del csp[min_i]
                    j = min_j
                    if min_t in [0,1]:
                        if min_t == 0: j=j-1
                        subp[-1][2], subp[0][0] = subp[-1][1], subp[0][1]
                        subp = [ [subp[j][1], subp[j][1], subp[j][2]] ] + subp[j+1:] + subp[:j] + [ [subp[j][0], subp[j][1], subp[j][1]] ]
                    else:
                        sp1,sp2,sp3 = csp_split(subp[j-1],subp[j],min_t)
                        subp[-1][2], subp[0][0] = subp[-1][1], subp[0][1]
                        subp = [ [ sp2[1], sp2[1],sp2[2] ] ] + [sp3] + subp[j+1:] + subp[:j-1] + [sp1] + [[ sp2[0], sp2[1],sp2[1] ]]
                    csp = [subp] + csp
                    if csp_subpath_ccw(csp[0]) :
                        for i in range(len(csp)):
                            n = []
                            for j in csp[i]:
                                n = [  [j[2][:],j[1][:],j[0][:]]  ] + n
                            csp[i] = n[:]

                    d = str(inkex.paths.CubicSuperPath(csp).to_path())

                    d = re.sub(r'(?i)(m[^mz]+)',r'\1 Z ',d)
                    d = re.sub(r'(?i)\s*z\s*z\s*',r' Z ',d)
                    d = re.sub(r'(?i)\s*([A-Za-z])\s*',r' \1 ',d)
                
                Rayon = self.svg.unittouu(self.options.dotsize)*in_out
                
                           
                NouvelleCourbe =etree.SubElement(MonChemin.getparent(), inkex.addNS('path','svg'),
                        {   inkex.addNS('label','inkscape'): 'gabarit: ' + self.options.dotsize + " " + self.options.stylegab,
                            inkex.addNS('type','sodipodi'): 'inkscape:offset',
                            inkex.addNS('radius','inkscape'):   str(Rayon),
                            inkex.addNS('original','inkscape'): d,
                        })
                NouvelleCourbe.style=inkex.Style({ 'stroke': '#00ff00', 'fill': 'none' ,'stroke-width': str(self.svg.unittouu('1px'))})


if __name__ == '__main__':

    Gabarits().run()





