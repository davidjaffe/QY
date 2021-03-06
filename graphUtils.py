#!/usr/bin/env python
'''
basic modules to make, print tgraphs
20150921 

'''
import time
import datetime
import sys
import os
import math
import ROOT 
from ROOT import TH1D, TFile, gROOT, TCanvas, TLegend, TGraph, TDatime, TMultiGraph, gStyle, TGraphErrors, TLine
from array import array
import re # regular expression

class graphUtils():
    def __init__(self):
    # use in color()
        self.goodColors = [1,2,3,4, 6,7,8,9] # no yellow(5) or white (0,10)
        self.goodColors.extend( [11, 12, 18] )
        self.goodColors.extend( [x for x in range(28,50)] )
        self.goodMarkers = [x for x in range(20,31) ]
        
        return
    
    def t2dt(self,t):
        '''
        convert struct_time to datetime object
        '''
        return datetime.datetime.fromtimestamp(time.mktime(t))
    def dt2t(self,dt):
        '''
        convert datetime object to struct_time (time object)
        '''
        fmt = "%Y %M %d %H %m %S"
        return time.strptime(dt.strftime(fmt),fmt)
    def addSeconds(self,t,seconds=0):
        '''
        add seconds to struct_time t by converting to datetime object,
        using timedelta and converting back
        '''
        dt = self.t2dt(t)
        dt += datetime.timedelta(seconds=seconds)
        return self.dt2t(dt)
    def convertTime(self,day,fmt,text):
        c = text.count(":")
        if c==1: return time.strptime(day+text,fmt+"%H:%M")
        if c==2: return time.strptime(day+text,fmt+"%H:%M:%S")
        sys.exit("graphUtils.convertTime ERROR Unknown input " + str(text))
        return
    def fixTimeDisplay(self,g,showDate=False,maybeShowDate=True):
        '''
        set time axis to display nicely
        '''
        if g:
            g.GetXaxis().SetTimeDisplay(1)
            g.GetXaxis().SetTimeFormat("%H:%M")
            if showDate:
                g.GetXaxis().SetTimeFormat("#splitline{%H:%M}{%y/%m/%d}")
            else:
                if maybeShowDate:
                    x1 = g.GetXaxis().GetXmin()
                    x2 = g.GetXaxis().GetXmax()
                    if x2-x1>24.*60.*60.:
                        g.GetXaxis().SetTimeFormat("#splitline{%H:%M}{%y/%m/%d}")
                        #print 'graphUtils.fixTimeDisplay: >1 day, so use splitline in SetTimeFormat'
            g.GetXaxis().SetNdivisions(-409)
            g.GetXaxis().SetLabelSize(0.025) #0.5*lx)
            g.GetXaxis().SetTimeOffset(0,"local") # what does this do?
#            g.GetXaxis().SetTimeOffset(0,"gmt") # using gmt option gives times that are only off by 1 hour on tgraph
        else:
            print 'graphUtils.fixTimeDisplay: WARNING Null pointer passed to fixTimeDisplay?????'
        return

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
    def makeTH1Dwtd(self,x,y,title,Name='',NX=None,XMI=None,XMA=None):
        '''
        fill 1d hist with weights y
        given equal size, monotonically increasing bin centers x
        '''
        name = Name
        if Name=='': name = title.replace(' ','_').replace('.','_')
        nx = len(x)
        dx = x[1]-x[0]
        xmi = min(x)-dx/2.
        xma = max(x)+dx/2.
        if NX is not None: nx = NX
        if XMI is not None:xmi =XMI
        if XMA is not None:xma =XMA
        h = TH1D(name,title,nx,xmi,xma)
        for a,b in zip(x,y): h.Fill(a,b)
        ymi,yma = min(y),max(y)
        dy = (yma-ymi)/20.
        ymi,yma = ymi-dy/2.,yma+dy/2.
        h.SetMaximum(yma)
        h.SetMinimum(ymi)
        return h
    def printHistStats(self,h):
        '''
        print some stats for input hist
        '''
        N,mean,stddev,underflow,overflow = self.getHistStats(h)
        print h.GetTitle(),'mean',mean,'stddev',stddev,'Nentries',N,'uflow',underflow,'oflow',overflow
        return
    def getHistStats(self,h):
        '''
        return histogram stats
        '''
        axis = 1 # 1d hist only
        mean = h.GetMean(axis)
        stddev = h.GetStdDev(axis)
        N = h.GetEntries()
        underflow = h.GetBinContent(0)
        if axis==1: nbins = h.GetNbinsX()
        if axis==2: nbins = h.GetNbinsY()
        if axis==3: nbins = h.GetNbinsZ()
        overflow = h.GetBinContent(nbins+1)
        return N,mean,stddev,underflow,overflow
    def drawGraph(self,g,figDir="",SetLogx=False,SetLogy=False,option='APL'):
        '''
        output graph to file
        '''
        title = g.GetTitle()
        name  = g.GetName()
        if SetLogx: name += '_logx'
        if SetLogy: name += '_logy'

        pdf   = os.path.join( figDir,  name + '.pdf')
    
        xsize,ysize = 1100,850 # landscape style
        noPopUp = True
        if noPopUp : gROOT.ProcessLine("gROOT->SetBatch()")
        canvas = TCanvas(pdf,title,xsize,ysize)

        
        g.Draw(option)

        if SetLogy: canvas.SetLogy(1)
        if SetLogx: canvas.SetLogx(1)
    
        canvas.Draw()
        canvas.SetGrid(1)
        canvas.SetTicks(1)
        canvas.cd()
        canvas.Modified()
        canvas.Update()
        canvas.Print(pdf,'pdf')
        return
    def drawFit(self,h,figdir='',SetLogy=False,SetLogx=False,extraName=None):
        '''
        draw histogram with fit parameters
        '''
        name = h.GetName()
        if extraName is not None: name += '_' + extraName
        title = h.GetTitle()
        if SetLogx: name += '_logx'
        if SetLogy: name += '_logy'
        pdf = os.path.join( figdir ,  name + '.pdf' )
        ps  = os.path.join( figdir , name + '.ps' )
        xsize,ysize = 1100,850 # landscape style
        noPopUp = True
        if noPopUp : gROOT.ProcessLine("gROOT->SetBatch()")
        canvas = TCanvas(pdf,title,xsize,ysize)
        
        gStyle.SetOptFit(1111)
        h.Draw()
        if SetLogy: canvas.SetLogy(1)
        if SetLogx: canvas.SetLogx(1)
        
        canvas.Draw()
        canvas.SetGrid(1)
        canvas.SetTicks(1)
        canvas.cd()
        canvas.Modified()
        canvas.Update()
        
        canvas.Print(ps,'Landscape')
        os.system('ps2pdf ' + ps + ' ' + pdf)
        if os.path.exists(pdf): os.system('rm ' + ps)

        return
    def finishDraw(self,canvas,ps,pdf,setGrid=True,setTicks=True,ctitle=None):
        '''
        standard nonsense to finish drawing
        ctitle can be considered 'global' title
        '''
        canvas.Draw()
        canvas.SetGrid(setGrid)
        canvas.SetTicks(setTicks)
        canvas.cd()
        canvas.Modified()
        canvas.Update()
        ct = None
        if ctitle is not None:
            ct = ROOT.TText(0.5,0.975,ctitle)
            ct.SetTextAlign(20) # horizontally centered
            s = ct.GetTextSize()
            ct.SetTextSize(s/2.) 
            ct.Draw()
        
        canvas.Print(ps,'Landscape')
        os.system('ps2pdf ' + ps + ' ' + pdf)
        if os.path.exists(pdf): os.system('rm ' + ps)
        return
    def drawMultiHists(self,histlist,fname='',figdir='',statOpt=1111111,setLogy=False,setLogx=False,dopt='',abscissaIsTime=False):
        '''
        draw multiple histograms on single pdf output file
        '''
        nHist = len(histlist)
        if nHist<=0:
            print 'graphUtils.drawMultiHists: ERROR zero length histogram list'
            return
        if nHist==1:
            nX = nY = 1
        else:
            nX = 2
            nY = int(float(nHist)/float(nX) + 0.5)
            if nHist<4: nX,nY = 1,nHist

        #print 'nHist,nX,nY=',nHist,nX,nY
        # create output directory if it does not exist
        if len(figdir)>0:
            if os.path.isdir(figdir):
               pass
            else:
                try:
                    os.mkdir(figdir)
                except IOError,e:
                    print 'graphUtils.drawMultiHists',e
                else:
                    print 'graphUtils.drawMultiHists created',figdir
        # set output file name and canvas title

        ctitle = None
        if fname!='':
            pdf = os.path.join( figdir, fname )
            ctitle = fname
        else:
            words = ''
            for h in histlist:
                name = h.GetName()
                words += name
                if h!=histlist[-1]: words += '_'
            pdf = os.path.join( figdir, words )

        if setLogx: pdf += '_logx'
        if setLogy: pdf += '_logy'
        ps = pdf + '.ps'
        pdf += '.pdf'
        
        # open canvas, draw on it
        title = ''
        xsize,ysize = 1100,850 # landscape style
        noPopUp = True
        if noPopUp : gROOT.ProcessLine("gROOT->SetBatch()")
        canvas = TCanvas(pdf,title,xsize,ysize)

        gStyle.SetOptStat(statOpt)
        spaceBtwnPads = 0.01 / 1000.
        canvas.Divide(nX,nY,spaceBtwnPads,spaceBtwnPads)
        for i,h in enumerate(histlist):
            canvas.cd(i+1).SetLogy(setLogy)
            canvas.cd(i+1).SetLogx(setLogx)
            if abscissaIsTime : self.fixTimeDisplay(h)

            h.Draw(dopt)
            self.biggerLabels(h)
            if abscissaIsTime : self.fixTimeDisplay(h)
            #print i+1,h.GetName()

        self.finishDraw(canvas,ps,pdf,ctitle=ctitle)
        return
    def biggerLabels(self,h):
        '''
        increase axis label size

        '''
        factor = 2.0 # empirically determined
        sx = h.GetXaxis().GetLabelSize()
        h.GetXaxis().SetLabelSize(factor*sx)
        sy = h.GetYaxis().GetLabelSize()
        h.GetYaxis().SetLabelSize(factor*sy)

        return
        
    def drawMultiGraph(self,TMG,figdir='',SetLogy=False, SetLogx=False, abscissaIsTime = True, drawLines=True, xAxisLabel=None,yAxisLabel=None):
        '''
        draw TMultiGraph with legend and output as pdf
        Default is that abscissa is calendar time.
        
        '''
        debugMG = False
        if not TMG.GetListOfGraphs(): return  # empty
        title = TMG.GetTitle()
        name  = TMG.GetName()
        if SetLogx: name += '_logx'
        if SetLogy: name += '_logy'
        if debugMG: print 'graphUtils.drawMultiGraph',title,name,'TMG.GetListOfGraphs()',TMG.GetListOfGraphs(),'TMG.GetListOfGraphs().GetSize()',TMG.GetListOfGraphs().GetSize()
        nGraphs = TMG.GetListOfGraphs().GetSize()


        pdf = os.path.join( figdir , name  + '.pdf' )
        ps  = os.path.join( figdir , name  + '.ps')
        xsize,ysize = 1100,850 # landscape style
        noPopUp = True
        if noPopUp : gROOT.ProcessLine("gROOT->SetBatch()")
        canvas = TCanvas(pdf,title,xsize,ysize)
        canvas.SetLogy(SetLogy)
        canvas.SetLogx(SetLogx)

        # move title to left in order to put legend above plot
        gStyle.SetTitleX(0.3)
        x1 = 0.5
        x2 = x1 + .5
        y1 = 0.9
        y2 = y1 + .1
        lg = TLegend(x1,y1,x2,y2)
        NGraph = 0
        for g in TMG.GetListOfGraphs():
            NGraph += 1
            t = g.GetTitle()
            lg.AddEntry(g, t, "LP" )
            if abscissaIsTime : self.fixTimeDisplay(g)
        if NGraph>6: lg.SetNColumns(2)

        dOption = "AP"
        if drawLines: dOption += "L"

        # complicated monkey business because of idiotic way that logY is set
        if SetLogy:
            ymi,yma = 1.e20,1.e-20
            for g in TMG.GetListOfGraphs():
                x,y = self.getPoints(g)
                ymi = min(ymi,min(y))
                yma = max(yma,max(y))
            if ymi<=0: ymi = 0.1
            ymi = ymi/2.
            yma = 2.*yma
            TMG.SetMinimum(ymi)
            TMG.SetMaximum(yma)
            for g in TMG.GetListOfGraphs():
                g.SetMinimum(ymi)
                g.SetMaximum(yma)
                if "A" in dOption:
                    g.SetTitle( TMG.GetTitle() )
                    if xAxisLabel is not None: g.GetXaxis().SetTitle(xAxisLabel)
                    if yAxisLabel is not None: g.GetYaxis().SetTitle(yAxisLabel)
                g.Draw(dOption)
                dOption = dOption.replace("A","")
        else:
            TMG.Draw(dOption)
            if xAxisLabel is not None: TMG.GetXaxis().SetTitle(xAxisLabel)
            if yAxisLabel is not None: TMG.GetYaxis().SetTitle(yAxisLabel)
            if abscissaIsTime : self.fixTimeDisplay(TMG)


        self.labelTMultiGraph(TMG,debug=debugMG)
        lg.Draw()
        canvas.Draw()
        canvas.SetGrid(1)
        canvas.SetTicks(1)
        canvas.cd()
        canvas.Modified()
        canvas.Update()
        if 0:
            canvas.Print(pdf,'pdf')
        else:
            canvas.Print(ps,'Landscape')
            os.system('ps2pdf ' + ps + ' ' + pdf)
            if os.path.exists(pdf): os.system('rm ' + ps)

        if debugMG: print 'graphUtils.drawMultiGraph',title,'complete'
        return
    def makeTMultiGraph(self,name,tit=None):
        title = tit
        if tit is None:title = name.replace('_',' ')
        tmg = TMultiGraph()
        tmg.SetName(name)
        tmg.SetTitle(title)
        return tmg
    def labelTMultiGraph(self,tmg,debug=False):
        name = tmg.GetName()
        if 'vs' in name:
            s = name.split('_')
            xt = s[2]
            yt = s[0]
            xt = xt.replace('by','/')
            xt = xt.replace('BY','/')
            yt = yt.replace('by','/')
            yt = yt.replace('BY','/')
            if debug:
                print 'graphUtils.labelTMultiGraph: xt,yt',xt,yt,'tmg',tmg
                print 'tmg.GetXaxis()',tmg.GetXaxis(),'tmg.GetYaxis()',tmg.GetYaxis()
            if tmg.GetXaxis(): tmg.GetXaxis().SetTitle(xt)
            if tmg.GetYaxis(): tmg.GetYaxis().SetTitle(yt)
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
    def color(self,obj,n,M,setMarkerColor=False,setMarkerType=True):
        '''
        set line color and marker type for obj based on indices n and M
        if M=n then use M to set marker type, otherwise determine marker type from n
        unless setMarkerType is False
        '''
        debug = False
        LC = len(self.goodColors)
        LM = len(self.goodMarkers)
        c = n%LC
        obj.SetLineColor( self.goodColors[c] )
        if debug: print 'color: obj',obj,'n',n,'obj.IsA().GetName()',obj.IsA().GetName()
        if setMarkerType:
            oName = obj.IsA().GetName()
            if oName=='TGraph' or oName=='TGraphErrors':
                if M==n:
                    m = M%LM
                else:
                    m = int(float(n)/float(LC))%LM
                obj.SetMarkerStyle( self.goodMarkers[m] )
                if setMarkerColor: obj.SetMarkerColor( self.goodColors[c] )
                if debug: print 'color:',obj.GetName(),'m',m,'self.goodMarkers[m]',self.goodMarkers[m]
        return
    def getPoints(self,g,getErrors=False):
        '''
        return abscissa,ordinate values of input graph g
        also return errors if getErrors is True
        '''
        x,y = [],[]
        if getErrors: dx,dy = [],[]
        for i in range(g.GetN()):
            a,b = ROOT.Double(0),ROOT.Double(0)
            OK = g.GetPoint(i,a,b)
            if OK!=-1:
                x.append(a)
                y.append(b)
                if getErrors:
                    dx.append(g.GetErrorX(i))
                    dy.append(g.GetErrorY(i))
        if getErrors: return x,y,dx,dy
        return x,y
    def getTDatime(self,dt,fmt='%Y/%m/%d %H:%M:%S'):
        '''
        convert date/time text to TDatime object
        '''
        datetimeObj = self.getdatetime(dt,fmt=fmt)
        return TDatime( datetimeObj.strftime('%Y-%m-%d %H:%M:%S') ).Convert()
    def getdatetime(self,dt,fmt='%Y/%m/%d %H:%M:%S'):
        ''' convert timestamp dt to text '''
        return datetime.datetime.strptime(dt,fmt)
    def reportHist(self,h):
        '''
        write out some properties of hist h
        '''
        name = h.GetName()
        title = h.GetTitle()
        xa = h.GetXaxis()
        nx = xa.GetNbins()
        xmi= xa.GetXmin()
        xma= xa.GetXmax()
        xex= xa.CanExtend()
        nd = h.GetDimension()
        words = 'graphUtils.reportHist',name,title,'nx,xmi,xma',nx,xmi,xma
        if xex: words += 'can extend x-axis.'
        if nd>1:
            ya = h.GetYaxis()
            ny = ya.GetNbins()
            ymi= ya.GetXmin()
            yma= ya.GetXmax()
            yex= ya.CanExtend()
            words += 'ny,ymi,yma=',ny,ymi,yma
            if yex: words += 'can extend y-axis.'
        print words
        return
