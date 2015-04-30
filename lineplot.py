import matplotlib as mpl
import matplotlib.pyplot as plt
import mydictionaries as Mdict

colors = ['k','r','b','g','m','c','y']
linestyles = ['-','--','-.',':']

def linePlot(data,linedict=None,scale=[1,1],xlabel=None,ylabel=None,axis=None,legendlabels=None,legendlocation=2,fontsize1=21,fontsize2=18,ticksize=[8,2],figuresize=[8,6],font='Arial',fignum=1):
    plt.figure(fignum,figsize=figuresize)
    plt.clf()
    preProcessPlot(font,figuresize)
    if linedict is None:
        linedict = lineDictStandard()
    for i, datacurr in enumerate(data):
        linedictcurr = Mdict.writeNewDict(linedict,i)
        if len(datacurr) == 2: # x and y specified individually
            line, = plt.plot(scale[0]*datacurr[0],scale[1]*datacurr[1],**linedictcurr)
        else: # in single array
            line, = plt.plot(scale[0]*datacurr[:,0],scale[1]*datacurr[:,1],**linedictcurr)
        if legendlabels is not None:
            line.set_label(legendlabels[i])
    postProcessPlot(xlabel,ylabel,axis,legendlabels,legendlocation,fontsize1,fontsize2,ticksize)

def preProcessPlot(font,figuresize):
    mpl.rcParams['font.sans-serif'] = font
    mpl.rcParams['pdf.fonttype'] = 42
    
def postProcessPlot(xlabel,ylabel,axis,legendlabels,legendlocation,fontsize1,fontsize2,ticksize):
    if xlabel is not None:
        plt.xlabel(xlabel,fontsize=fontsize1)
    if ylabel is not None:
        plt.ylabel(ylabel,fontsize=fontsize1)
    if axis is not None:
        plt.axis(axis)
    plt.legend(loc=legendlocation,fontsize=fontsize2)
    plt.tick_params(axis='both',labelsize=fontsize2,width=ticksize[1],length=ticksize[0])
    plt.tight_layout()
    plt.show()

def genLabels(numvec,label):
    return [label + ' = ' + str(num) for num in numvec]
    
def alternatingList(stringlist,n,nuniques):
    return stringlist[0:nuniques]*(n//nuniques)
    
def alternatingList2(stringlist,n,nuniques):
    listfinal = []
    for i in range(nuniques):
        listfinal = listfinal + stringlist[i:i+1]*(n//nuniques)
    return listfinal

def lineDictStandard(color=colors*10,linestyle=linestyles*10,marker='',linewidth=2.5,markersize=10): # every line is a different color and linetype, no markers
    keys = ['color','linestyle','marker','linewidth','markersize']
    values = [color,linestyle,marker,linewidth,markersize]
    return dict(zip(keys,values))
    
def lineDictAlt1(n,ncolors): # cycle through colors and linestyles, type 1, no markers
    # e.g. color = ['k','k','r','r','b','b'], linestyle = ['-','--','-','--','-','--']
    linedict = lineDictStandard()
    linedict['color'] = alternatingList2(linedict['color'],n,ncolors)
    linedict['linestyle'] = alternatingList(linedict['linestyle'],n,n//ncolors)
    return linedict
    
def lineDictAlt2(n,ncolors): # cycle through colors and linestyles, type 2, no markers
    # e.g. color = ['k','r','b','k','r','b'], linestyle = ['-','-','-','--','--','--']
    linedict = lineDictStandard()
    linedict['color'] = alternatingList(linedict['color'],n,ncolors)
    linedict['linestyle'] = alternatingList2(linedict['linestyle'],n,n//ncolors)
    return linedict
