#!/usr/bin/env/python
'''
Copyright (C) 2021 Vantieghem David, 
Code partiellement repris de l'excellent path along path
Mod pour creation de point couture par Vantieghem David 2018-2019
'''

#library
import copy
import math
import random
from inkex import TextElement,bezier
from inkex.bezier import csplength
import inkex
import pathmodifier
from inkex.paths import CubicSuperPath

inkex.localization.localize()

def offset(pathcomp,dx,dy):
    for ctl in pathcomp:
        for pt in ctl:
            pt[0]+=dx
            pt[1]+=dy

def linearize(p,tolerance=0.001):
    '''
    This function recieves a component of a 'cubicsuperpath' and returns two things:
    The path subdivided in many straight segments, and an array containing the length of each segment.

    We could work with bezier path as well, but bezier arc lengths are (re)computed for each point
    in the deformed object. For complex paths, this might take a while.
    '''
    zero=0.000001
    i=0
    d=0
    lengths=[]
    while i<len(p)-1:
        box  = bezier.pointdistance(p[i  ][1],p[i  ][2])
        box += bezier.pointdistance(p[i  ][2],p[i+1][0])
        box += bezier.pointdistance(p[i+1][0],p[i+1][1])
        chord = bezier.pointdistance(p[i][1], p[i+1][1])
        if (box - chord) > tolerance:
            b1, b2 = bezier.beziersplitatt([p[i][1],p[i][2],p[i+1][0],p[i+1][1]], 0.5)
            p[i  ][2][0],p[i  ][2][1]=b1[1]
            p[i+1][0][0],p[i+1][0][1]=b2[2]
            p.insert(i+1,[[b1[2][0],b1[2][1]],[b1[3][0],b1[3][1]],[b2[1][0],b2[1][1]]])
        else:
            d=(box+chord)/2
            lengths.append(d)
            i+=1
    new=[p[i][1] for i in range(0,len(p)-1) if lengths[i]>zero]
    new.append(p[-1][1])
    lengths=[l for l in lengths if l>zero]
    return(new,lengths)

class Pointsellier(pathmodifier.Diffeo):
    def __init__(self):
        pathmodifier.Diffeo.__init__(self)
        
        self.arg_parser.add_argument("--title")

        self.arg_parser.add_argument("--diamlong",dest="diamlong", default="1.0mm")

        self.arg_parser.add_argument("--typePoint",dest="typePoint", default="LigneH")

        self.arg_parser.add_argument("--textInfos",type=inkex.Boolean, dest="textInfos", default=False)

        self.arg_parser.add_argument("-t", "--toffset",dest="toffset", default="0.1mm")

        self.arg_parser.add_argument("--autoSpace",type=inkex.Boolean, dest="autoSpace", default=False)

        self.arg_parser.add_argument("-p", "--space", dest="space", default="3.0mm")

        self.arg_parser.add_argument("--autoOffset",  type=inkex.Boolean,dest="autoOffset", default=False)

        self.arg_parser.add_argument("-r","--nrepeat", type=int, dest="nrepeat", default=1,help="nombre d'objets")

        self.arg_parser.add_argument("--autoRepeat", type=inkex.Boolean, dest="autoRepeat", default=False)

        self.arg_parser.add_argument("--autoMask", type=inkex.Boolean,dest="autoMask", default=False)

        self.arg_parser.add_argument("--autoMark", type=inkex.Boolean,dest="autoMark", default=False)

        self.arg_parser.add_argument("--typeMark",dest="typeMark", default="markX")

        self.arg_parser.add_argument( "--nrepeat2", type=int,dest="nrepeat2", default=1,help="nombre d'objets")

        self.arg_parser.add_argument("--tab",dest="tab",help="The selected UI-tab when OK was pressed")


    def lengthtotime(self,l):
        '''
        Recieves an arc length l, and returns the index of the segment in self.skelcomp
        containing the coresponding point, to gether with the position of the point on this segment.

        If the deformer is closed, do computations modulo the toal length.
        '''
        if self.skelcompIsClosed:
            l=l % sum(self.lengths)
        if l<=0:
            return 0,l/self.lengths[0]
        i=0
        while (i<len(self.lengths)) and (self.lengths[i]<=l):
            l-=self.lengths[i]
            i+=1
        t=l/self.lengths[min(i,len(self.lengths)-1)]
        return i, t

    def applyDiffeo(self,bpt,vects=()):
        '''
        The kernel of this stuff:
        bpt is a base point and for v in vectors, v'=v-p is a tangent vector at bpt.
        '''
        s=bpt[0]-self.skelcomp[0][0]
        i,t=self.lengthtotime(s)
        if i==len(self.skelcomp)-1:#je regarde si je suis au debut du skelete car sinon j'ai pas de vecteur
            x,y=bezier.tpoint(self.skelcomp[i-1],self.skelcomp[i],1+t)
            dx=(self.skelcomp[i][0]-self.skelcomp[i-1][0])/self.lengths[-1]
            dy=(self.skelcomp[i][1]-self.skelcomp[i-1][1])/self.lengths[-1]
        else:
            x,y=bezier.tpoint(self.skelcomp[i],self.skelcomp[i+1],t)
            dx=(self.skelcomp[i+1][0]-self.skelcomp[i][0])/self.lengths[i]
            dy=(self.skelcomp[i+1][1]-self.skelcomp[i][1])/self.lengths[i]

        vx=0
        vy=bpt[1]-self.skelcomp[0][1]
        bpt[0]=x+vx*dx-vy*dy
        bpt[1]=y+vx*dy+vy*dx

        for v in vects:
            vx=v[0]-self.skelcomp[0][0]-s
            vy=v[1]-self.skelcomp[0][1]
            v[0]=x+vx*dx-vy*dy
            v[1]=y+vx*dy+vy*dx

    def CreatePatern(self,node,idPoint,labelPoint,diametre,typepoint):
        #creation des paterns par défauts

        Dimensions=self.svg.unittouu(diametre)
        Patern = inkex.PathElement()
        
        if typepoint=="Cercle":
            d = 'M dim,0 A dim,dim 0 0 1 0,dim dim,dim 0 0 1 -dim,0 dim,dim 0 0 1 0,-dim dim,dim 0 0 1 dim,0 Z'.replace('dim',str(Dimensions/2))
        elif typepoint=="LigneV":
            d= 'M 0,0 V dim'.replace('dim',str(Dimensions))
        elif typepoint=="LigneH":
            d= 'M 0,0 H dim'.replace('dim',str(Dimensions))
        elif typepoint=="LigneG45":
            d= 'M 0,dim dim,0'.replace('dim',str(Dimensions))
        else:
            d= 'M dim,dim 0,0'.replace('dim',str(Dimensions))
        
        Patern.path=d  
        Patern.set('id',idPoint)
        Style= { 'stroke': '#000000', 'fill': 'none','stroke-opacity':'1', 'stroke-width': str(self.svg.unittouu('1px')) }
        Patern.set('style', Style)
        Patern.label=labelPoint
        
        return Patern 
    
    def addText(self,node,x,y,text):
        MonText=node.getparent().add(inkex.TextElement())
        MonText.set('style', "font-style:normal;font-weight:normal;font-size:3px;line-height:100%;font-family:sans-serif;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1")#simplestyle.formatStyle(s))
        MonText.set('x', str(x))
        MonText.set('y', str(y))
        MonText.text = str(text)
        MonText.label="Texte - " + text

    def CreateMark(self,node,idPoint,labelPoint,diametre, Couleur):
        #creation des paterns par défauts
        Dimensions=self.svg.unittouu(diametre)
        Mark = inkex.PathElement()
        d = 'M 0,0 V dim'.replace('dim',str(Dimensions))
        Mark.path=d  
        Mark.set('id',idPoint)
        Style= { 'stroke': Couleur, 'fill': 'none','stroke-opacity':'1', 'stroke-width': str(self.svg.unittouu('1px')) }
        Mark.set('style', Style)
        Mark.label=labelPoint
        
        return Mark

    def effect(self):

        if len(self.options.ids)<1 and len(self.options.ids)>1:
            raise inkex.AbortExtension("This extension requires only one selected paths.")

        if self.options.autoSpace and self.options.autoRepeat:
            raise inkex.AbortExtension("Vous ne pouvez pas avoir l autoSpace et l autoRepeat en auto en meme temps.")
        
        #liste des chemins, preparation
        MaSelection = self.svg.selection.filter(inkex.PathElement)
        
        if not MaSelection:
            raise inkex.AbortExtension("This extension need a path, not groups.")
  
        for MonChemin in MaSelection:

            #je recupère l'id du chemin pour crer les autres
            idUnique = MonChemin.get('id')
            idpoint=idUnique+'-'+ str(random.randint(1, 99)) #id du paterns creer a partir du chemin selectionner
            idpointMark=idUnique+'-'+ str(random.randint(1, 99))
        
            MonStyle = MonChemin.style #je recupere l'ancien style pour avoir un sructure complete
            #je modifie la valeur je modifie la couleur du chemin selectionner
            MonStyle['stroke']='#00ff00' 
            if self.options.autoMask==True:
                MonStyle['display']='none'
            MonChemin.style= MonStyle  #j'applique la modif

            #calcul de la longeur du chemin
            csp = MonChemin.path.transform(MonChemin.composed_transform()).to_superpath()
            Longeur_ligne, Longeur = csplength(csp)

            #génération du patern de base
            MonPattern=self.CreatePatern(MonChemin,idpoint,"Temp",self.options.diamlong,self.options.typePoint) #creation du patern sélectioner de base de base
            MonChemin.getparent().add(MonPattern)
            MonPatern_box=MonPattern.bounding_box()
            
            #je calcul les différente informations qui me sont nécéssaire
            
            if not self.options.autoOffset: #gestion du decallage automatique
                tOffset=self.svg.unittouu(self.options.toffset)
                Longeur=Longeur-tOffset # si je fais un décallaged'offset j'enlève d'autant la taille restante
                
            Taille= self.svg.unittouu(self.options.diamlong)
            Distance=self.svg.unittouu(self.options.space) #pour prendre en compte la taille du point ajout de ceci +bbox.x.size #calcul de la distance en prenant en compte la taille en x du patern
            MaxCopies=max(1,int(round((Longeur+Distance)/Distance)))
            NbCopies= self.options.nrepeat #nombre de copie desirer a integrer dans les choix a modifier pour ne pas depasser les valeurs maxi

            if NbCopies > MaxCopies:
                NbCopies=MaxCopies #on limitte le nombre de copie au maxi possible sur le chemin

            if self.options.autoRepeat and not self.options.autoSpace: #auto pour le nombre de marque avec un espace constant
                NbCopies=MaxCopies
            
            if self.options.autoSpace and not self.options.autoRepeat: #auto espacement avec un nombre de marque constantes (pour les cercles)
                NbCopies= self.options.nrepeat
                Distance=Longeur/NbCopies

            if Distance < 0.01:
                raise inkex.AbortExtension("The total length of the pattern is too small :\nPlease choose a larger object or set 'Space between copies' > 0")
                
            if self.options.autoOffset: #gestion du decallage automatique
                tOffset=((Longeur-(NbCopies-1)*Distance)/2)-Taille/2
            else:
                tOffset=self.svg.unittouu(self.options.toffset)
                
            #gestion du paterns
            labelpoint='Point: '+ idpoint+ ' Nbr:' + str(NbCopies)+' longueur:'+str(round(self.svg.uutounit(Longeur,'mm'),2))+'mm'            
            MonPattern.label=labelpoint
            #je récupère le chemin en mode inkscape            
            d = MonPattern.get('d')
            #et je le transforme en mode point
            p0 = CubicSuperPath(d)
            
            #creation de la liste pour le positionnement
            newp=[]

            #je récupère le chemin en mode points
            self.curSekeleton=CubicSuperPath(MonChemin.get('d'))
            
            # et je le parcours
            for comp in self.curSekeleton:
                p=copy.deepcopy(p0)
                #je crer un vecteur a partir des point
                self.skelcomp,self.lengths=linearize(comp)
                #test simple pour savoir si il est clot ou pas
                self.skelcompIsClosed = (self.skelcomp[0]==self.skelcomp[-1])

                #j'applique des offset de décallage
                xoffset = self.skelcomp[0][0] - MonPatern_box.x.minimum + tOffset
                yoffset = self.skelcomp[0][1] - MonPatern_box.y.center
                
                #si j'ai l'options de texte
                if self.options.textInfos:
                    self.addText(MonChemin,xoffset,yoffset,labelpoint)

                #une nouvelle liste pour chaque points
                new=[]
                for sub in p: #creation du nombre de patern
                    for i in range(0,NbCopies,1):
                        new.append(copy.deepcopy(sub)) #realise une copie de sub pour chaque nouveau element du patern                                        
                        offset(sub,Distance,0)
                        
                #je refait un passe en mappant sur le chemin
                p=new
                for sub in p:
                    offset(sub,xoffset,yoffset)
                for sub in p: #une fois tous creer, on les mets en place
                    for ctlpt in sub:#pose le patern sur le chemin
                        self.applyDiffeo(ctlpt[1],(ctlpt[0],ctlpt[2]))
                #je met le patern dans la liste     
                newp+=p
            #je retransforme le mode vectoriel/point en mode inkscape
            MonPattern.set('d', CubicSuperPath(newp).to_path())

        # je fait les calcul pour le marquage
        if self.options.autoMark:
            if self.options.typeMark=="markFraction":
                Fraction= self.options.nrepeat2 #en mode fraction 1= au debut et a la fin, 2= un demi, 3= 1/3 etc 
                Distance=Longeur/Fraction #distance inter point
                NbrMark=max(1,int(round((Longeur+Distance)/Distance))) 
                infos= " Marquage 1/"+ str(Fraction)
                couleur= '#ff0000'
            else:
                Repeat= self.options.nrepeat2  #en mode fraction 1= au debut et a la fin, 2= 1 sur 2, 3= 1 sur 3 etc 
                NbrMark=max(1,int(round((NbCopies/Repeat)))) 
                Distance=Distance*Repeat #distance inter point
                infos=" Marquage tous les " + str(Repeat) + " points"
                couleur= '#ffaa00'

        #Puis je l'applique avec la même méthode que au dessus
            labelMark="Mark: "+idpoint + infos
            MaMark=self.CreateMark(MonChemin,idpoint,labelMark,self.options.diamlong,couleur) #creation du patern sélectioner de base de base
            MonChemin.getparent().add(MaMark)
            MaMark_box=MaMark.bounding_box()
            d = MaMark.get('d')
            p0 = CubicSuperPath(d)
                
            newp=[]

            self.curSekeleton=CubicSuperPath(MonChemin.get('d'))
            for comp in self.curSekeleton:
                p=copy.deepcopy(p0)
                self.skelcomp,self.lengths=linearize(comp)
                #test simple pour savoir si il est clot ou pas
                self.skelcompIsClosed = (self.skelcomp[0]==self.skelcomp[-1])

                xoffset = self.skelcomp[0][0] - MaMark_box.x.minimum + tOffset + MonPatern_box.x.size/2 #la formule permet de prendre en compte la taille en x du point :(Longeur-(NbrMark-1)*Distance)/2 -Taille/2
                yoffset = self.skelcomp[0][1] - MaMark_box.y.center
                
                if self.options.textInfos:
                    self.addText(MonChemin,xoffset,yoffset,labelpoint)

                new=[]
                for sub in p: #creation du nombre de patern
                    for i in range(0,NbrMark,1):
                        new.append(copy.deepcopy(sub)) #realise une copie de sub pour chaque nouveau element du patern                                        
                        offset(sub,Distance,0)
                p=new
                for sub in p:
                    offset(sub,xoffset,yoffset)
                for sub in p: #une fois tous creer, on les mets en place
                    for ctlpt in sub:#pose le patern sur le chemin
                        self.applyDiffeo(ctlpt[1],(ctlpt[0],ctlpt[2]))
                        
                newp+=p
            MaMark.set('d', CubicSuperPath(newp).to_path())
                
                
if __name__ == '__main__':
    
    Pointsellier().run()

