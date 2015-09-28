#!/usr/bin/env python
'''
read data files for quantum yield measurements
20150921 

'''
import sys
import os
import math
import ROOT 
from ROOT import TH1D, TFile, gROOT, TCanvas, TLegend, TGraph, TDatime, TMultiGraph, gStyle, TGraphErrors, TLine
from array import array
import numpy
from scipy.interpolate import interp1d
from scipy.optimize    import curve_fit

class reader():
    def __init__(self):
        self.inputData = {}
        
        self.figDir = 'Figures/'
        self.histFile = 'Output/banana.root'

        self.filelist = []
        self.Graphs = []
        self.Graphdict = {}
        self.Hists  = []
        self.MultiGraphs = {}

    # parametrizations describing baseline
        self.baselinePar = {}
        

    # use in color()
        self.goodColors = [x for x in range(1,10)]
        self.goodColors.extend( [11, 12, 18] )
        self.goodColors.extend( [x for x in range(28,50)] )
        self.goodMarkers = [x for x in range(20,31) ]

        return
    def read2Col(self,filename):
        '''
        read a two column file, unless it has already been read. In that case, get from dict
        line0 = skip
        line1 = number of measurements
        line2 = description
        line3 = column headers
        line4:line4+N = data
        Return Nmeas, Descrip, ColumnHeaders, Data

        Requires some monkey business to avoid running into problems with control characters
        '''
        if filename not in self.inputData: 
            f = open(filename,'r')
            Nmeas = None
            Descrip = None
            ColH = None
            Contents = []
            for line in f:
                if 'Trace' not in line:
                    lm1 = line[:-1]
                    if Nmeas is None:
                        Nmeas = int(lm1)
                    elif Descrip is None:
                        Descrip = ''
                        for a in line.split():
                            Descrip += a + ' ' 
                    elif ColH is None:
                        ColH = line.split()
                    else:
                        words = line.split()
                        a = []
                        for v in words:
                            a.append(float(v))
                        Contents.append(a)
            self.inputData[filename] = Nmeas,Descrip,ColH,Contents
        return self.inputData[filename]
    def readOrdinates(self,filename):
        N,D,H,C = self.read2Col(filename)
        y = []
        for pair in C: y.append(pair[1])
        return y
    def readPairs(self,filename):
        N,D,H,C = self.read2Col(filename)
        x,y = [],[]
        for pair in C:
            x.append(pair[0])
            y.append(pair[1])
        return x,y
    def compareFiles(self,file1,file2):
        '''
        compare contents of two files
        '''
        N1,D1,H1,C1 = self.read2Col(file1)
        N2,D2,H2,C2 = self.read2Col(file2)

        print 'Compare file1',file1,'and file2',file2
        print 'Nmeas',self.printComp(N1,N2)
        print 'Descrip',self.printComp(D1,D2)
        print 'ColHdrs',self.printComp(H1,H2)
        if N1!=N2:
            print 'Different number of measurements. Don`t compare contents'
        else:
            for A,B in zip(C1,C2):
                if A[0]==B[0]:
                    r = None
                    if B[1]!=0.: r = A[1]/B[1]
                    print H1[0],A[0],H1[1],A[1],B[1],'ratio',r
                else:
                    print H1[0],A[0],H2[0],B[0]
        return
    def testMode(self,filename):
        '''
        test estimation of mode of data
        '''
        y = self.readOrdinates(filename)
        print filename
        f = -1.
        for nIter in range(4):
            mode,N,binWidth = self.getMode(y,f=f,nIter=nIter)
            print 'f',f,'nIter',nIter,'mode',mode,'N',N,'binWidth',binWidth,'len(list)',len(y)
        if 0:
            print 'list',y
            print 'list-mode',
            for x in y: print x-mode,
            print ''
        return
    def getModeF(self,filename,f=-1.,nIter=1):
        '''
        getMode using filename instead of list as input
        '''
        y = self.readOrdinates(filename)
        return self.getMode(y,f=f,nIter=nIter)
    def getMode(self,v,f=-1.,nIter=1):
        '''
        Given input list, method to determine bin width and # of iterations
        Return Mode, number of values used to compute Mode and bin width
        '''
        best = list(v)
        binWidth = -1. 
        for iter in range(nIter):
            best, binWidth = self.getModeList(best,f=f)
        mode = sum(best)/float(len(best))
        return mode,len(best),binWidth
    def getModeList(self,v,f=-1.):
        '''
        try to determine most probable value in input list using binning
        from spacing between values.
        Copy input list
        Sort new list
        Get minimum value in list
        if f>=0:   Find minimum difference between new list values and multiply by max(f,1.) to get bin width
        if f<0:    Use average difference between new list values as bin width
        Bin the data
        Get frequency of values in each bin
        Find highest frequency bin
        Check if next bin has same frequency
        Mode = mean value of entries in the highest frequency bin(s)
        Return values to be used to compute Mode and bin width
        '''
        y = list(v)
        y.sort()
        ymi = float(min(y))
        G = numpy.diff( numpy.unique(y) )
        if f<0:
            binWidth = sum(G)/float(len(G))
        else:
            binWidth = max(f,1.)*float(min(numpy.diff( numpy.unique(y) )))
        #print 'ymi',ymi,'f',f,'binWidth',binWidth,'d/f',d/f
        
        z = []
        for x in y: z.append(int((x-ymi)/binWidth))
        
        freq = {}
        for q in numpy.unique(z):
            freq[q] = z.count(q)
            #print 'q',q,'z.count(q)',z.count(q)
        km = max(freq,key=freq.get)
        
        L = 1
        if km+1 in freq:
            if freq[km+1]==freq[km]: L += 1
        #print 'getMode: f',f,'km',km,'freq[km]',freq[km]

                
        best = []
        ylo = ymi + binWidth*float(km)
        yhi = ymi + binWidth*float(km+L)
        for q in y:
            if ylo<=q and q<=yhi: best.append(q)
        return best, binWidth
    def printComp(self,A,B):
        '''
        string for printing A,B and note if identical or not
        '''
        if A==B:
            s = str(A) + ' identical'
        else:
            s = 'first: ' + str(A) + ' second: ' + str(B)
        return s
    def getNameFromFilename(self,filename):
        '''
        filename includes path = topdir/subdir/file.ext
        name = subdir_file with blanks replaced by underscores
        
        '''
        S = os.path.split(filename)
        subdir = S[0].split('/')[1]
        name = S[1]
        name = name.replace('.txt','')
        name = subdir + '_' + name
        name = name.replace(' ','_')
        return name
    def getTitleFromFilename(self,filename):
        title = filename.replace(' ','_')
        return title
    def getSubDirFromFilename(self,filename):
        S = os.path.split(filename)
        return S[0].split('/')[1]
    def makeGraph(self,filename,Yoffset=0.,suffix=''):
        '''
        make TGraph from input file, subtracting Yoffset from ordinate values and appending suffix to name and title
        name is taken from filename without file extension, with blanks replaced by underscores, prepended with a character to make root happy
        title is taken from description in file and filename with path. blanks replaced by underscores
        '''
        N,D,H,C = self.read2Col(filename)
        u,v = [],[]
        for pair in C:
            u.append(pair[0])
            v.append(pair[1] - Yoffset)
        name = self.getNameFromFilename(filename)
        if name[0] in '0123456789': name = 'g'+name # prepend character because root doesn't like variables that start with integer
        name += suffix    
        title = self.getTitleFromFilename(filename)
        title = D + '_' + title + '_' + suffix
        g = self.makeTGraph(u,v,title,name)
        g.SetMarkerStyle(20)
        return g
            
    def drawGraph(self,g):
        '''
        output graph to file
        '''
        title = g.GetTitle()
        name  = g.GetName()
        pdf   = self.figDir + name + '.pdf'
        
        xsize,ysize = 1100,850 # landscape style
        noPopUp = True
        if noPopUp : gROOT.ProcessLine("gROOT->SetBatch()")
        canvas = TCanvas(pdf,title,xsize,ysize)

        g.Draw()
        
        canvas.Draw()
        canvas.SetGrid(1)
        canvas.SetTicks(1)
        canvas.cd()
        canvas.Modified()
        canvas.Update()
        canvas.Print(pdf,'pdf')
        return

    def makeTGraph(self,u,v,title,name,ex=None,ey=None):
        if ex is None:
            g = TGraph(len(u),array('d',u), array('d',v))
        else:
            dy = ey
            if ey is None: dy = [0. for x in range(len(ex))]
            g = TGraphErrors(len(u),array('d',u),array('d',v),array('d',ex),array('d',dy))
        g.SetTitle(title)
        g.SetName(name)
        return g
    def makeHist(self,filename,prefix='Hist',nx=40,xmi=1.,xma=-1.):
        '''
        make TH1D from input file
        '''
        N,D,H,C = self.read2Col(filename)
        v = []
        for pair in C: v.append(pair[1])
        name = self.getNameFromFilename(filename)
        title= D + self.getTitleFromFilename(filename)
        if len(prefix)>0:
            name = prefix[0] + name
            title = prefix + title
        else:
            name = 'H'+name
            title = 'Hist ' + title
        h = self.makeTH1D(v,title,name,nx=nx,xmi=xmi,xma=xma)
        return h
    def makeTH1D(self,v,title,name,nx=100,xmi=1,xma=-1):
        if xmi>xma:
            xmi = min(v)
            xma = max(v)
        dx = (xma-xmi)/float(nx)
        xmi -= dx/2.
        xma += dx/2.
        h = TH1D(name,title,nx,xmi,xma)
        for y in v: h.Fill(y)
        return h
    def gangGraphs(self,keywords=[],vetos=[]):
        '''
        create multigraph containing all graphs with titles that match keywords and do not match vetos
        '''
        debug = False
        if len(keywords)==0: return # this bores me
        key = ''
        for k in keywords:
            key += k
            if keywords.index(k)<len(keywords)-1 : key += '__'
        n = 0
        if debug: print 'gangGraphs: key',key,'graphs:'
        for g in self.Graphs:
            title = g.GetTitle()
            OK = True
            for k in keywords: OK = OK and (k in title)
            for v in  vetos: OK = OK and not (v in title)
            if OK:
                n += 1
                if key not in self.MultiGraphs:
                    self.MultiGraphs[key] = TMultiGraph()
                    self.MultiGraphs[key].SetTitle(key)
                    self.MultiGraphs[key].SetName(key)
                self.color(g,n,n)
                self.MultiGraphs[key].Add(g)
                if debug: print 'add graph',g.GetName()
                
        return
    def color(self,obj,n,M):
        '''
        set line color and marker type for obj based on indices n and M
        if M=n then use M to set marker type, otherwise determine marker type from n
        '''
        debug = False
        LC = len(self.goodColors)
        LM = len(self.goodMarkers)
        c = n%LC
        obj.SetLineColor( self.goodColors[c] )
        if debug: print 'color: obj',obj,'n',n,'obj.IsA().GetName()',obj.IsA().GetName()
        if obj.IsA().GetName()=='TGraph':
            if M==n:
                m = M%LM
            else:
                m = int(float(n)/float(LC))%LM
            obj.SetMarkerStyle( self.goodMarkers[m] )
            if debug: print 'color: m',m,'self.goodMarkers[m]',self.goodMarkers[m]
        return
    def getListOfAllDataFiles(self,topdir='data',fext=None):
        '''
        get list of all data files with file extension matching fext
        assumes directory structure is topdir/subdir/file 
        '''
        filenames = []
        for dirname in os.listdir(topdir):
            dn = topdir + '/' + dirname
            if os.path.isdir(dn):
                for x in os.listdir(dn):
                    if (fext is None) or '.'+fext in x:
                        fn = dn + '/' + x
                        if os.path.isfile(fn): filenames.append(fn)
        return filenames
    def getBaseline(self,filename,BLguess=0.,BLwidth=1.e20,n=0,Axis='ordinate'):
        '''
        Return parameters to linear function that gives the baseline for the data in filename.
        Also return values used in baseline estimate
        if Axis=='ordinate', then
           Provide initial guess at baseline and range of ordinate values to include when estimating baseline.
        else
           Use the first and last n values of abscissa 
        '''
        x,y = self.readPairs(filename)

        if Axis=='ordinate':
            aX,aY = self.getConsistentPairs(x,y,BLguess,BLwidth,Axis=Axis)
            P = self.fitLinear(aX,aY,p0=BLguess,p1=0.)
        else:
            width=abs(x[n]-x[0])
            aX,aY = self.getConsistentPairs(x,y,x[0],width,Axis=Axis)
            bX,bY = self.getConsistentPairs(x,y,x[-1],width,Axis=Axis)
            aX.extend(bX)
            aY.extend(bY)
            P = self.fitLinear(aX,aY)
        
        return P, aX,aY
    def fitLinear(self,aX,aY, p0=None,p1=None):
        '''
        Return parameters to linear function y = p0 + p1*x that fits the input values.
        p0,p1 are initial guesses
        '''
        if p0 is None: p0 = sum(aX)/float(len(aX))
        if p1 is None: p1 = 0.
        P0 = numpy.array( [p0,p1] )
        X = numpy.array( aX )
        Y = numpy.array( aY )
        P,cov = curve_fit(self.linFun, X,Y,P0)
        return P
    def getConsistentPairs(self,x,y,z0,dz,Axis='ordinate'):
        '''
        return abscissa and ordinate values with Axis values consistent with z0+-dz
        If Axis=='ordinate', then use ordinate values, otherwise use abscissa
        '''
        aX,aY = [],[]
        for u,v in zip(x,y):
            z = u
            if Axis=='ordinate': z = v
            if (abs((z-z0)/dz)<=1.) :
                aX.append(u)
                aY.append(v)
        return aX,aY
    def graphBaseline(self,filename,Par,prefix='Baseline'):
        '''
        create and return graph of baseline defined by parameters Par of polynomial function
        '''
        #print 'graphBaseline: filename',filename,'Par',Par
        x,y = self.readPairs(filename)
        name = prefix[0] + self.getNameFromFilename(filename)
        title = prefix + ' ' + self.getTitleFromFilename(filename)
        Ax,Ay = array('d',[]),array('d',[])
        for u in x:
            v = self.linFun(u,Par[0],Par[1])
            Ax.append(u)
            Ay.append(v)
        g = self.makeTGraph(Ax,Ay,title,name)
        return g
    def linFun(self,X,a,b):
        return a+b*X
    def graphValues(self,filename,Ax,Ay,prefix='Showvalues'):
        '''
        graph values used to estimate baseline
        '''
        name = prefix[0] +  self.getNameFromFilename(filename)
        title = prefix + ' ' + self.getTitleFromFilename(filename)
        
        g = self.makeTGraph(Ax,Ay,title,name)
        return g
    def findGraph(self,name):
        '''
        return graph given name
        '''
        for g in self.Graphs:
            if g.GetName()==name: return g
        sys.exit('findGraph: ERROR Failed to find ' + name)
        return None
    def getGraphContent(self,g):
        '''
        return lists of contents of graph
        from Brett: a,b are memory locations, so must copy them with float 
        '''
        n = g.GetN()
        x,y = [],[]
        a,b = ROOT.Double(0.),ROOT.Double(0.)
        for i in range(n):
            g.GetPoint(i,a,b)
            #print 'g.GetName(),i,a,b',g.GetName(),i,a,b
            x.append(float(a))
            y.append(float(b))
        #print 'g.GetName(),x,y',x,y
        return x,y
    def compareGraphs(self,name):
        '''
        Create pdf file that compares graphs created from same raw data
        Plot raw data(points+line) and either points or lines for the following.
        X = [Line] baseline estimate using initial and final abscissa values
        A = [Points] initial and final abscissa values used for estimate
        P = [Line] baseline esimate using mode of ordinate values
        O = [Points] values used for mode-based estimate
        Two panels are created.
        Top panel shows full ordinate range.
        Lower panel concentrates on baseline.
        '''
        prefixes = {'X':'L', 'A':'P', 'P':'L', 'O':'P'}
        graw = self.findGraph(name)
        title = 'Comparison' + graw.GetTitle()
        purename = name.replace('raw','')
        pdf = self.figDir+'compare_'+purename+'.pdf'

        xsize,ysize = 1100,850
        noPopUp = True
        if noPopUp: gROOT.ProcessLine("gROOT->SetBatch()")
        canvas = TCanvas(pdf,title,xsize,ysize)
        canvas.Divide(1,2)
        x,y = self.getGraphContent(graw)
        xmi,xma = min(x),max(x)
        #print 'compareGraphs:graw.GetName(),x[0],x[1],x[-1],len(x),y[0],y[1],y[-1],len(y)',graw.GetName(),x[0],x[1],x[-1],len(x),y[0],y[1],y[-1],len(y)
        abit = (xma-xmi)/40.
        xmi -= abit
        xma += abit
        ymi,yma = 1.e20,-1.e20
        options = {1:"ALP",2:"LP"}
        for i in [1,2]:

            canvas.cd(i)
            if ymi<yma:
                abit = (yma-ymi)
                yma += abit
                ymi -= abit/10.
                #print 'compareGraphs: graw',graw.GetName(),'xmi,xma,ymi,yma',xmi,xma,ymi,yma
                canvas.DrawFrame(xmi,ymi,xma,yma)
            graw.Draw(options[i])
            for key in prefixes:
                n = key+purename
                g = self.findGraph(n)
                g.Draw(prefixes[key])
                x,y = self.getGraphContent(g)
                ylo = min(y)
                yhi = max(y)
                if ylo<ymi: ymi=ylo
                if yma<yhi: yma=yhi

            canvas.Draw()
            canvas.SetGrid(1)
            canvas.SetTicks(1)
            canvas.cd()
            canvas.Modified()
            canvas.Update()
        canvas.Print(pdf,'pdf')

        return
    def storeGraph(self,g):
        '''
        append graph to global list and dict and bleat if requested
        '''
        debug = True
        if debug:         print 'Created graph #bins,name,title',g.GetN(),g.GetName(),g.GetTitle()
        self.Graphs.append(g)
        if g.GetName() in self.Graphdict: # error?
            sys.exit('storeGraph: ERROR graph name '+g.GetName()+' already exists')
        else:
            self.Graphdict[g.GetName()] = g
        return
    def calcCorrFun(self,subdir):
        '''
        calculate correction function using different baseline subtractions using data from subdir
        Get graph of raw data
        Get list of parameter values and keywords
        Do baseline subtraction.
        Plot result.
        Integrate
        Store resulting integrals and method as function of wavelength
        '''
        Results = {}

        for filename in self.filelist:
            sd = self.getSubDirFromFilename(filename)
            if subdir==sd or subdir==sd.replace(' ','_'):
                name = self.getNameFromFilename(filename) 
                rawname = name + 'raw'
                graw = self.Graphdict[rawname]
                xraw,yraw = self.getGraphContent(graw)
                listOfPar = self.baselinePar[filename] 
                for pair in listOfPar:
                    Par,words = pair
                    y = self.createGraphMinusFunction(xraw,yraw,Par)
                    gname = 'blSub_'+words+name
                    title = 'Baseline-subtracted ' + words + ' ' + self.getTitleFromFilename(filename)
                    gnew = self.makeTGraph(xraw,y,title,gname)
                    self.storeGraph(gnew)
                    Integral = self.integrate(xraw,y)
                    print name,words,Integral
                    wl = float(name.split('_')[1])
                    if words not in Results:
                        Results[words] = [ [wl], [Integral] ]
                    else:
                        Results[words][0].append( wl )
                        Results[words][1].append( Integral )

        for words in Results:
            x,y = Results[words]
            Ry = []
            minRy = 1.e20
            for q in y:
                if q>0:
                    minRy = min(minRy,1./q)
                    Ry.append(1./q)
                else:
                    Ry.append(-1.)

            title = subdir + ' Baseline-subtracted Integral vs wavelength ' + words
            name = 'blSub_' + subdir + '_' + words
            ###print 'method,name,title,x,y',method,name,title,x,y
            gnew = self.makeTGraph(x,y,title,name)
            self.storeGraph(gnew)
            title = subdir + ' Reciprocal of baseline-sub integral vs wavelength ' + words
            name = 'blSubR_' + subdir + '_' + words
            gnew = self.makeTGraph(x,Ry,title,name)
            self.storeGraph(gnew)
            y = []
            for q in Ry:
                if q<0.:
                    y.append(-1.)
                else:
                    y.append(q/minRy)
            title = subdir + 'Normed, recip. of baseline-sub integral vs wavelength ' + words
            name = 'blSubRN_' + subdir + '_' + words
            gnew = self.makeTGraph(x,y,title,name)
            self.storeGraph(gnew)
                        
        return
    def graphLaurenCorr(self,subdir):
        '''
        make a graph of Lauren's correction function
        '''
        laurenDir = '../lcapelluto_WbLS/correction/'
        fn = laurenDir + subdir + '.txt'
        if os.path.isfile(fn):
            f = open(fn,'r')
            print 'graphLaurenCorr: Opened',fn
            x,y = [],[]
            for line in f:
#                print line
                a,b = line.split(' ')
                x.append(float(a))
                y.append(float(b))
            title = fn
            name  = 'lcapelluto_'+subdir
            gnew = self.makeTGraph(x,y,title,name)
            self.storeGraph(gnew)
        else:
            print 'graphLaurenCorr: ERROR non-existent file',fn
        return
    def createGraphMinusFunction(self,xIn,yIn,Par):
        '''
        yOut = yIn - f(xIn,Par)
        Get parameters
        Make copy of input so that output is same type
        loop over x values
        compute function value
        delete old entry
        insert new entry
        '''
        a,b = Par
        yOut = yIn[:]
        for i,x in enumerate(xIn):
            fx = self.linFun(x,a,b)
            y = yIn[i]
            del yOut[i]
            yOut.insert(i,y-fx)
        return yOut
    def integrate(self,xIn,yIn):
        '''
        compute integral of yIn, take bin widths from separation between xIn values
        Assumes Xin is monotonically increasing
        '''
        sY = 0.
        sX = 0.
        Lm1 = len(yIn)-1
        for i,y in enumerate(yIn):
            dxlo = 0.
            if i>1: dxlo = xIn[i]-xIn[i-1]
            dxhi = 0.
            if i<Lm1: dxhi = xIn[i+1]-xIn[i]
            sY += y*(dxhi+dxlo)
            sX += dxhi+dxlo
        if sX>0.: sY = sY/sX
        return sY
    def compareGraphContents(self,gnames,descrip):
        '''
        compare contents of graphs named in input list
        create multigraph
        create graphs of ratios
        '''
        glist = []
        gdict = {}
        goodnames = []
        finterp = {}
        nmax = 0
        xBig = None
        if descrip not in self.MultiGraphs:
            self.MultiGraphs[descrip] = TMultiGraph()
            self.MultiGraphs[descrip].SetTitle(descrip)
            self.MultiGraphs[descrip].SetName(descrip)
        for name in gnames:
            if name in self.Graphdict:
                goodnames.append(name)
                g = self.Graphdict[name]
                glist.append(g)
                N = goodnames.index(name)
                self.color(g,N,N)
                self.MultiGraphs[descrip].Add(g)
                x,y = self.getGraphContent(g)
                gdict[name] = [x,y]
                if len(x)>nmax:
                    nmax = len(x)
                    xBig = x
                finterp[name] = interp1d(x,y,bounds_error=False)
            else:
                print 'compareGraphContents: ERROR',name,'is not the name of a graph'

        print ''
        for name in goodnames:
            print 'Column',goodnames.index(name),'Data',name
        line = ''
        nl = 10
        sl = '{0:^'+str(nl)+'}'
        line += sl.format('x')
        for name in goodnames:
            line += sl.format('y'+str(goodnames.index(name)))
        for name in goodnames:
            if name!=goodnames[0]:
                line += sl.format('y'+str(goodnames.index(name))+'/y0')
        print line

        sf = '{0:>'+str(nl)+'.2f}'
        si = '{1:>'+str(nl-1)+'.2f}{0:1.1}'
        ratio = {}
        for ax in xBig:
            line = ''
            line += sf.format(ax)
            y0 = None
            r  = []
            for name in goodnames:
                x,y = gdict[name]
                if ax in x:
                    Y = y[x.index(ax)]
                    line += sf.format(Y)
                else:
                    Y = float(finterp[name](ax)) # otherwise output type is numpy.ndarray
                    line += si.format('I',Y)
                if y0 is None:
                    y0 = Y
                else:
                    if name not in ratio: ratio[name] = []
                    R = Y/y0
                    r.append(R)
                    ratio[name].append(R)
            for Y in r:
                line += sf.format(Y)
            print line
        print ''

        if len(ratio)>0:
            name0 = goodnames[0]
            mgname = descrip + '_Ratio'
            self.MultiGraphs[mgname] = TMultiGraph()
            self.MultiGraphs[mgname].SetName(mgname)
            self.MultiGraphs[mgname].SetTitle(mgname)
            for name in ratio:
                title = 'ratio of correction ' + name + ' to ' + name0
                gname = 'Ratio_'+name
                g = self.makeTGraph(xBig,ratio[name],title,gname)
                self.storeGraph(g)
                N = goodnames.index(name)
                self.color(g,N,N)
                self.MultiGraphs[mgname].Add(g)
        return
    def Main(self):
        '''
        evocatively titled.
        Get list of all input data files.
        Make a TGraph for each file
        Write TGraphs to output file
        '''
        self.filelist = self.getListOfAllDataFiles(fext='txt')
        print 'Created file list with',len(self.filelist),'files'
        
        for filename in self.filelist:
            #self.testMode(filename)

            # plot raw data
            g = self.makeGraph(filename,suffix='raw')
            rawname = g.GetName()
            self.storeGraph(g)

            # plot data with mode subtracted
            mode,Nval,binWidth = self.getModeF(filename,nIter=1)
            #print 'mode',mode,'Nval',Nval,'binWidth',binWidth
            g = self.makeGraph(filename,Yoffset=mode,suffix='subMode')
            self.storeGraph(g)

            # using values consistent with mode, estimate baseline with linear function and plot it
            Par, aX,aY = self.getBaseline(filename,BLguess=mode,BLwidth=binWidth)
            if filename not in self.baselinePar: self.baselinePar[filename] = [ [Par, 'modeEst'] ]
            g = self.graphBaseline(filename,Par,'Plot for mode')
            g.SetLineColor(ROOT.kBlue)
            self.storeGraph(g)

            # show raw data values used to estimate baseline
            g = self.graphValues(filename,aX,aY,prefix='Ordinate')
            g.SetLineColor(ROOT.kBlue)
            g.SetMarkerStyle(2) # +
            self.storeGraph(g)


            # using n initial and final values, estimate baseline with linear function and plot it
            Par, aX,aY = self.getBaseline(filename,Axis='Abscissa',n=5)
            self.baselinePar[filename].append( [Par, 'inifinEst'] )
            g = self.graphBaseline(filename,Par,'X-axis selection')
            g.SetLineColor(ROOT.kRed)
            self.storeGraph(g)

            # show raw data values used to estimate baseline
            g = self.graphValues(filename,aX,aY,prefix='Abscissa')
            g.SetLineColor(ROOT.kRed)
            g.SetMarkerStyle(5) # X
            self.storeGraph(g)

            # comparison of graphs and methods to estimate baseline
            self.compareGraphs(rawname)
            
            
            # histograms
            h = self.makeHist(filename)
            print 'Created hist',h.GetName(),h.GetTitle()
            self.Hists.append(h)

            h = self.makeHist(filename,prefix='ModeSub',xmi=mode-binWidth,xma=mode+binWidth)
            print 'Created hist',h.GetName(),h.GetTitle()
            self.Hists.append(h)
            
        print 'Created',len(self.Graphs),'graphs and',len(self.Hists),'hists'

        # plot multiple graphs
        self.gangGraphs(keywords=['water','raw'],vetos=['WbLS','Baseline'])
        self.gangGraphs(['water','subMode']     ,vetos=['WbLS','Baseline'])

        # calculation of correction function
        self.calcCorrFun('water')
        self.calcCorrFun('emcorr\ w\ cyhx')

        # Lauren's correction functions
        for name in ['water','normalizer','alignedDefault','alignedCyhx',]:
            self.graphLaurenCorr(name)

        # compare correction functions
        self.compareGraphContents(['lcapelluto_water','blSubRN_water_inifinEst','blSubRN_water_modeEst'],'waterBased_correction')
        
        outfile = TFile(self.histFile,'RECREATE')

        for g in self.Graphs:   outfile.WriteTObject(g)
        print 'Wrote',len(self.Graphs),'graphs to',self.histFile

        for h in self.Hists:
            outfile.WriteTObject(h)
        print 'Wrote',len(self.Hists),'hists to',self.histFile
        
        print 'MultiGraphs: '
        for k in self.MultiGraphs:
            g = self.MultiGraphs[k]
            print g.GetName()
            outfile.WriteTObject(g)
        print ''
        print 'Wrote',len(self.MultiGraphs),'multigraphs to',self.histFile



        outfile.Close()
        return
    
if __name__ == '__main__' :
    r = reader()
    r.Main()
    sys.exit()
    print fl
    if 1: sys.exit()
    fn1 = 'data/empty-sphere/17_370nm_excorr.txt'
    fn2 = 'data/empty-sphere/17_370nm_uncorr.txt'
    #r.compareFiles(fn1,fn2)
    for fn in [fn1,fn2]:
        if fn is not None:
            g = r.makeGraph(fn)
            r.drawGraph(g)
            if 0:
                N,D,H,C = r.read2Col(fn)
                print 'Nmeas',N
                print 'Descrip',D
                print 'ColH',H
                print 'Contents',C
                print 'len(Contents)',len(C)
    
