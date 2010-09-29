""" This module contains functions I use regularly"""
__author__='John Kitchin <jkitchin@udel.edu>'

import re, os, sys, getopt, os, operator
from Scientific.Functions.Interpolation import InterpolatingFunction
import Numeric
from Simulations.Dacapo import *
from Scientific.IO.NetCDF import NetCDFFile
from Simulations.Dacapo.AtomProjectedDOS import AtomProjectedDOS    # import AtomProjectedDOS class
from Simulations.Dacapo.ListOfEigenStates import ListOfEigenStates  # import ListOfEigenstates class
from GracePlot.GracePlot import *
from struct import unpack

def strip(*lists):
    """
    Used to strip None values out of a group of lists for analysis or plotting
    """
    ind=[]
    for i in range(len(lists[0])):
        for list in lists:
            if list[i]==None:
                ind.append(i)
                break
    count=0
    for i in ind:
        for list in lists:
            del list[i-count]
        count=count+1
        
def GetCalculationTime(filename):
    """
    Used to get the CPU time used in a calculation
    """
    file,ext=os.path.splitext(filename)
    file=file+'.txt'
    try:
        f=open(file,'r')
    except:
        print "that text file (%s) doesn't exist!" % file
        

    regexp=re.compile('User')
    number=re.compile('\d+\.\d+')
    for line in f.readlines():
        if regexp.search(line):
            time=number.search(line).group()
                
    return string.atof(time)

def PlotConvergence(filename):
    """
    this function only serves to plot the energy vs. time for a calculation
    to visualize convergence.
    """
    file,ext=os.path.splitext(filename)
    file=file+'.txt'
    try:
        f=open(file,'r')
    except:
        print "that text file (%s) doesn't exist!" % file
        sys.exit()

    regexp=re.compile('DFT:\s+\d+')
    times,energies=[],[]
    for line in f.readlines():
        if regexp.search(line):
            fields=line.split()
            times.append(string.atof(fields[1]))
            energies.append(string.atof(fields[3]))

    data1=Data(x=times,y=energies)
    p=GracePlot()
    p.plot(data1)


def GetVersion(filename):
    """
    Used to get the version of dacapo for a calculation.
    Actually, this assumes the ncfile and txt file are identical except for extension,
    which is how I do it.

    It would be really great if the version was stored somehow in the netcdf file.
    """
    file,ext=os.path.splitext(filename)
    file=file+'.txt'
    try:
        f=open(file,'r')
    except:
        print "that text file (%s) doesn't exist!" % file
        sys.exit()

    v27="dacapo-\d-\d"
    p26=re.compile('2\.6\.0')
    p27=re.compile(v27)
    version = "Unknown" # until it is redefined later
    for s in f.readlines():    
        if p26.search(s):
            version=p26.search(s).group()
        if p27.search(s):
            version=p27.search(s).group()
        
    return version

def IsVersion(filename,version):
    """
    tests if the filename was created with version, where version is one of the following
    IsVersion(ncfile,'2.6.0')
    IsVersion(ncfile,'dacapo-2.7')
    or it is another unique string that would be found in the corresponding text file.
    """
    return GetVersion(filename)==version

def Integrate(xdata,ydata):
    x_array=Numeric.array(xdata)
    y_array=Numeric.array(ydata)

    # Create the interpolating function
    f = InterpolatingFunction((x_array,), y_array)

    f_integral=f.definiteIntegral()

    return f_integral

def Integrate_Trapezoid(xdata,ydata):
    sum=0
    for i in range(len(xdata)-1):
        area_i=0.5*(xdata[i+1]-xdata[i])*(ydata[i+1]+ydata[i])
        sum=sum+area_i
    return sum

def Integrate_Simpson(xdata,ydata):
   n=len(xdata)
   h=(xdata[-1]-xdata[0])/(n-1)

   sum = ydata[0]

   for i in range(len(xdata)-1):
       if (i+1)%2 == 0:
           sum = sum+ydata[i+1]*2.0
       else:
           sum = sum+ydata[i+1]*4.0

   sum = sum + ydata[-1]
        
   integral = sum * h / 3.0

   return integral        


def Moments(xdata,ydata):
    """
    calculates the first,second, third and 4th moments of a distribution
    """
    xdata=Numeric.array(xdata)
    ydata=Numeric.array(ydata)

    xydata=xdata*ydata #first moment=Center of Gravity
    xxydata=xdata*xdata*ydata # Second moment = mean square width
    xxxydata=xdata*xdata*xdata*ydata # third moment = skewness
    xxxxydata=xdata*xdata*xdata*xdata*ydata # 4th moment = modalness

    normalization=Integrate(xdata,ydata)
    mom1=Integrate(xdata,xydata)/normalization
    mom2=Integrate(xdata,xxydata)/normalization
    mom3=Integrate(xdata,xxxydata)/normalization
    mom4=Integrate(xdata,xxxxydata)/normalization

    return mom1, mom2, mom3, mom4

def Moment(xdata,ydata,nth=1):
    """
    Calculate the nth moment  of the distribution xdata,ydata
 
    John Kitchin <jkitchin@udel.edu> 10-14-02
    """
    xdata=Numeric.array(xdata)
    ydata=Numeric.array(ydata)
    
    moment_data=(xdata**(nth))*ydata
    normalization=Integrate(xdata,ydata)
    
    return Integrate(xdata,moment_data)/normalization

from math import *
def Smooth(xdata,ydata,broadening=1.,xstep=0.1):
    xdata=N.array(xdata)
    ydata=N.array(ydata)
    
    norm_factor=sqrt(broadening/3.141592653589793)

    xmin=min(xdata)
    xmax=max(xdata)
    
    num_x_steps=int(ceil(fabs(xmax-xmin)/xstep)+1)

    new_xdata,new_ydata=[],[]
    for i in range(num_x_steps):
        new_xdata.append(xmin+float(i)*xstep)

        y_here=0.
        for k in range(len(ydata)):
            xdiff=new_xdata[i]-xdata[k]
            y_here=y_here+ydata[k]*exp(-broadening*xdiff**2.0)

        y_here=y_here/norm_factor
        new_ydata.append(y_here)

    return new_xdata,new_ydata


#Do Linear Least Squares Fit, returning the intercept, slope and R^2
from Scientific.Functions.LeastSquares import *
from Scientific.Functions.FirstDerivatives import DerivVar
from Scientific.Statistics import mean

def LinearLeastSquaresFit(xdata,ydata):
    """
    Special case of polynomialLeastSquaresFit,
    which returns the R^2 value

    intercept,slope,R_squared = LinearLeastSquaresFit(xdata,ydata)
    
    """
    xdata=Numeric.array(xdata)
    ydata=Numeric.array(ydata)
    # form that LeastSquaresFit wants the data in
    data=map(lambda x,y:(x,y),xdata,ydata)
    # I am not sure if the particular guess for parameters will result in a error in some cases
    lsq=polynomialLeastSquaresFit(parameters=(0,0), data=data)
    
    intercept=lsq[0][0]
    slope=lsq[0][1]

    avg_y=mean(ydata)
    
    predicted_data = intercept + slope*xdata
    numerator=Numeric.sum((predicted_data-ydata)**2)
    denominator=Numeric.sum((ydata-avg_y)**2)

    R_squared = 1 - numerator/denominator

    return intercept,slope, R_squared

##################################################################
### BANDS
###################################################################

def bands(ListOfAtoms,extra_bands=10):
    """
    This function is intended to calculate the number of bands you want in a dacapo calculation.
    It calculates the number of valence electrons in the ListOfAtoms, then adds the number of extra
    bands (default=10) to half of that number (making sure it is an integer)
    """
    v_electrons=ListOfAtoms.GetNumberOfValenceElectrons()
    return int(v_electrons/2+extra_bands)


def GetDistance(atom1,atom2):
    """
    Calculates the geometric distance between atom1 and atom2
    
    GetDistance(sim.loa[0],sim.loa[1]) #calculates distance between atom1 and atom2 in the ListOfAtoms
    """
    pos1=Numeric.array(atom1.GetCartesianPosition())
    pos2=Numeric.array(atom2.GetCartesianPosition())
    return Numeric.sqrt(reduce(operator.add,(pos1-pos2)*(pos1-pos2)))


def CheckRasmol(sim,repetitions=[1,1,1]):
    sim.WriteAsNetCDFFile('test.nc')        
    os.system('plotnetcdf -x %i -y %i -z %i test.nc' %(repetitions[0],repetitions[1],repetitions[2]))
    os.system('rm test.nc')

from Analysis.DensityOfStates import DensityOfStates,_LocalDensityOfStates_

def GetDOS(ListOfEigenStates,smoothfactor):
    Efermi = ListOfEigenStates.GetEFermi() # get the fermi level

    dos=DensityOfStates(energyinterval=[-15.+Efermi,Efermi])
    dos.Npoints=500  #number of desired points
    dos.SetListOfEigenStates(ListOfEigenStates)
    dos.SetSmoothFactor(smoothfactor)
    dos.SetSpin(None)
    dos.ReadFromListOfEigenStates()

    dos_energygrid  = dos.GetEnergies() # energy (x-axis) for DOS
    dos_energygrid  = dos_energygrid - Efermi       # shift energy zero point
    dos_gridvalues  = dos.GetArray()       # raw DOS (the y-axis)

    return dos_energygrid,dos_gridvalues


####################################################################################
class Calculation:
    """
    This class is intended to be a container for data extracted from dacapo
    calculations. It will continue to grow as new types of data are available,
    or i need to get different types of data from the ncfile.
    """
    def __init__(self,filename):
        """
        the class is initialized with a netcdf filename.
        There may be bugs if the netcdf filename and associated ascii filename differ by more
        than the extension.
        """
        self.filename=filename
        
        sim = Simulation()
        sim.loa = ListOfAtoms()
        sim.loe  = ListOfEigenStates()
        sim.ados = AtomProjectedDOS()
        sim.pdos = NetCDF.Entry("PrintAtomProjectedDOS") # used for retrieving short cutoff later
        sim.TotalEnergy=NetCDF.Entry("EvaluateTotalEnergy")
        
        sim.UpdateFromNetCDFFile(filename)

        self.sim = sim
        
        #self.energy=self.GetEnergy()
        self.version=GetVersion(filename)
        

    def PrintLOA(self):
        """
        Prints the LOA with less information than repr(LOA)
        """
        
        print " # %12s   %7s   %7s   %7s %12s" % ('species','x','y','z','label')
        print "---------------------------------------------------------------"
        
        for atomnumber in range(len(self.sim.loa)):
            type=self.sim.loa[atomnumber].GetType()
            pos=self.sim.loa[atomnumber].GetCartesianPosition()
            
            # I have modified my local CamposASE modules so I can put labels on the atoms. This should
            # allow it to work for unmodified modules and should be transparent for most users.
            try:
                print "%2d" % (atomnumber+1),"%12s" % type,"  %7f   %7f   %7f %12s" % (pos[0],pos[1],pos[2],self.sim.loa[atomnumber].GetLabel())
            except:
                print "%2d" % (atomnumber+1),"%12s" % type,"  %7f   %7f   %7f" % (pos[0],pos[1],pos[2])

    def GetNumberOfAtoms(self):
        """
        returns number of atoms in the unit cell
        """
        return len(self.sim.loa)

    def GetNumberOfValenceElectrons(self):
        """
        returns the number of Valence Electrons in the ListOfAtoms
        """
        return ValenceElectrons(self.sim.loa)
    
    def GetEnergy(self,*functional):
        """
        has been implemented in loa now
        """
        
        return self.loa.GetTotalPotentialEnergy(functional)

    def GetEFermi(self):
        """
        Returns the Fermi Level
        """
        return self.sim.loe.GetEFermi()

    def GetCutOffRadius(self):
        """
        Gets short cutoff radius used in AtomProjectedDOS. It returns the value it is set to
        if it was set, or 1.0 otherwise
        """
        try:
            return self.sim.pdos.CutoffRadius[0]
        except:
            return 1.0 #Default is 1 angstrom


##	def GetVariance(self):
##		"""
##		This function returns the variance of the ProjectedDOS
##		which is a measure of how spread out the band is about the
##		center. The standard deviation would be the square root of the variance.
		
##		John Kitchin <jkitchin@udel.edu>
##		"""
##		e1,d1=self._GetRelativeData()	       		
##		mean=self.GetMoments()[0]
		
##		sum=Numeric.sum((e1-mean)*(e1-mean))
		
##		variance=sum/(len(e1)-1)

##		return variance
