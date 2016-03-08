#!/usr/local/python2.1/bin/python
# this is an executable python script
import sys,os,math

import Tkinter
from Tkinter import *
import SimpleDialog,tkMessageBox

class QuestionTemplate2:

    def __init__(self,master,pile_number,filename,directory,max_height):
        
        Label(master, font="LucidaBright 18", text="Information for Pile Number "+str(pile_number+1)).grid(row=0,column=0,sticky=W,columnspan=2)
        Label(master, font="LucidaBright 14", text="Thickness of Initial Volume, h(x,y):").grid(row=1,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="P*(1-((x-xc)/xr)^2 - ((y-yc)/yr)^2)").grid(row=1,column=1,sticky=W)
        Label(master, font="LucidaBright 14", text="Maximum Initial Thickness, P (m):").grid(row=2,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="Center of Initial Volume, xc, yc (UTM E, UTM N):").grid(row=3,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="X and Y Extent, xr, yr (m, m):").grid(row=4,column=0,sticky=W)        
        self.pileheight = Entry(master,font="LucidaTypewriter 14")
        self.xpilecenter = Entry(master,font="LucidaTypewriter 14")
	self.ypilecenter = Entry(master,font="LucidaTypewriter 14")
        self.xradius = Entry(master,font="LucidaTypewriter 14")
        self.yradius = Entry(master,font="LucidaTypewriter 14")

        self.pileheight.grid(row=2,column=1)
        self.xpilecenter.grid(row=3,column=1)
	self.ypilecenter.grid(row=3,column=2)
        self.xradius.grid(row=4,column=1)
        self.yradius.grid(row=4,column=2)

        Button(master, font="LucidaBright 14", text="Done", command=self.done).grid(row=6,column=0)
        Button(master, font="LucidaBright 14", text="Quit", command=master.quit).grid(row=6,column=1)
        Button(master, font="LucidaBright 10", text="Calculate\n Volume", command=self.showVolume).grid(row=6,column=2)
        #passed in variables
        self.master = master
        self.top = master.winfo_toplevel()
        self.filename = filename
        self.max_height = max_height
        self.write_flag = 0 # used to make sure that the information is written out only once per call
        self.directory = directory
        self.input_flag = 0 #used to make sure that the values are entered
        
    def showVolume(self):
        if self.pileheight.get() == '':
            self.ph = 0.
        else:
            self.ph=float(self.pileheight.get())

        if self.xradius.get() == '':
            self.xr = 0.
        else:
            self.xr=float(self.xradius.get())
            
        if self.yradius.get() == '':
            self.yr = 0
        else:
            self.yr=float(self.yradius.get())
            
        self.volume=math.pi*self.ph*self.xr*self.yr/2.
        tkMessageBox.showinfo(self,message=self.volume)

    def done(self):
	# check values
	if self.pileheight.get() == '':
	    pileheight = 0.0
	else:
	    pileheight = float(self.pileheight.get())
	    if pileheight <= float(0.):
		pileheight = 0.

	if self.xpilecenter.get() == '':
	    xpilecenter = 1.0
	else:
	    xpilecenter = float(self.xpilecenter.get())

	if self.ypilecenter.get() == '':
	    ypilecenter = 1.0
	else:
	    ypilecenter = float(self.ypilecenter.get())

	if self.xradius.get() == '':
	    xradius = 1.0
	else:
	    xradius = float(self.xradius.get())
	    if xradius <= 0:
		xradius = 1

	if self.yradius.get() == '':
	    yradius = 1.0
	else:
	    yradius = float(self.yradius.get())
	    if yradius <= 0:
		yradius = 1

        if self.write_flag == 0:
            self.write_flag = 1
            file = open(self.filename, "a+", 0)
            file.write( str(pileheight) + '\n' + str(xpilecenter) + '\n' + str(ypilecenter) + '\n' + str(xradius) + '\n' +str(yradius) + '\n')
            file.close
            if pileheight > self.max_height:
                self.max_height = pileheight

            self.input_flag = 1
            
        
class QuestionTemplate:

    def __init__(self,master):

        Label(master, font="LucidaBright 14", text="GIS Information Main Directory:").grid(row=0,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="GIS Sub-Directory:").grid(row=1,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="GIS Map Set:").grid(row=2,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="GIS Map:").grid(row=3,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="Simulation Directory Location:").grid(row=4,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="Number of Processors:").grid(row=5,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="Computational Mesh Points in Y-Direction:").grid(row=6,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="Internal Friction Angle (deg):").grid(row=7,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="Bed Friction Angle (deg):").grid(row=8,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="Number of Piles:").grid(row=9,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="Scale Simulation? ").grid(row=10,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="If Scaled, Length Scale (m):").grid(row=11,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="Maximum Number of Time Steps:").grid(row=12,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="Maximum Time (s):").grid(row=13,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="Time Steps per Results Output:").grid(row=14,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="Adapt the Grid?").grid(row=15,column=0,sticky=W)	
        Label(master, font="LucidaBright 14", text="Visualization Output:").grid(row=16,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="First/Second Order Method:").grid(row=17,column=0,sticky=W)
        Label(master, font="LucidaBright 14", text="Email Address:").grid(row=18,column=0,sticky=W)

        self.topomain = Entry(master,font="LucidaTypewriter 14")
        self.toposub = Entry(master,font="LucidaTypewriter 14")
        self.topomapset = Entry(master,font="LucidaTypewriter 14")
        self.topomap = Entry(master,font="LucidaTypewriter 14")
        self.directory = Entry(master,font="LucidaTypewriter 14")
        self.numprocs = Entry(master,font="LucidaTypewriter 14")
        self.numygrid = Entry(master,font="LucidaTypewriter 14")
        self.intfrict = Entry(master,font="LucidaTypewriter 14")
        self.bedfrict = Entry(master,font="LucidaTypewriter 14")
        self.numpiles = Entry(master,font="LucidaTypewriter 14")
        self.scale = IntVar()
        self.scale_fn = Checkbutton(master,variable=self.scale,text="Yes",font="LucidaTypewriter 14")
        self.lengthscale = Entry(master,font="LucidaTypewriter 14")
        self.steps = Entry(master,font="LucidaTypewriter 14")
        self.maxtime = Entry(master,font="LucidaTypewriter 14")
        self.numoutput = Entry(master,font="LucidaTypewriter 14")
        self.adapt = IntVar()
        self.adapt_fn = Checkbutton(master,variable=self.adapt,text="Yes",font="LucidaTypewriter 14")
        self.order = IntVar()
        self.order_fn = Checkbutton(master,variable=self.order,text="Second",font="LucidaTypewriter 14")
        self.vizoutput = Menubutton(master, text="Choose Formats",relief=RAISED)
        self.emailaddress = Entry(master,font="LucidaTypewriter 14") 
	
        self.topomain.grid(row=0,column=1)
        self.toposub.grid(row=1,column=1)
        self.topomapset.grid(row=2,column=1)
        self.topomap.grid(row=3,column=1)
        self.directory.grid(row=4,column=1)
        self.numprocs.grid(row=5,column=1)
        self.numygrid.grid(row=6,column=1)
        self.intfrict.grid(row=7,column=1)
        self.bedfrict.grid(row=8,column=1)
        self.numpiles.grid(row=9,column=1)
        self.scale_fn.grid(row=10,column=1)
        self.lengthscale.grid(row=11,column=1)
        self.steps.grid(row=12,column=1)
        self.maxtime.grid(row=13,column=1)
        self.numoutput.grid(row=14,column=1)
        self.adapt_fn.grid(row=15,column=1)
        self.vizoutput.grid(row=16,column=1)
        self.order_fn.grid(row=17,column=1)
        #self.vizoutput.grid()
        self.emailaddress.grid(row=18,column=1)

        Button(master, font="LucidaBright 14", text="Run", command=self.run).grid(row=19,column=1)
        Button(master, font="LucidaBright 14", text="Quit", command=master.quit).grid(row=19,column=2)
        Button(master, font="LucidaBright 14", bitmap="question", command=self.getHelp).grid(row=19,column=4)
        Button(master, font="LucidaBright 14", text="Clean", command=self.clean).grid(row=19,column=0)

        self.vizoutput.menu = Menu(self.vizoutput, tearoff=0)
        self.vizoutput["menu"] = self.vizoutput.menu
        self.tecplotVar = IntVar()
        self.mshplotVar = IntVar()
        self.hdfVar = IntVar()
        self.padyVar = IntVar()
        
        
        self.vizoutput.menu.add_checkbutton(label="tecplotxxxx.plt",variable=self.tecplotVar)
        self.vizoutput.menu.add_checkbutton(label="mshplotxxxx.plt",variable=self.mshplotVar)
        self.vizoutput.menu.add_checkbutton(label="GMFG Viz",variable=self.padyVar)
        self.vizoutput.menu.add_checkbutton(label="HDF",variable=self.hdfVar)
        

    def getHelp(self):
        a=open('../README').read()
        tk = Tkinter.Tk()
        frame = Tkinter.Frame(tk, relief=FLAT, borderwidth=0)
        frame.pack(fill=BOTH,expand=1)
        label = Tkinter.Label(frame, text="Instructions")
        label.pack(fill=X, expand=1)
        text=Tkinter.Text(frame,relief=RIDGE,width=75)
        text.insert(END,a)
        text.pack(side=LEFT)
        button=Tkinter.Button(frame,text="OK",command=tk.destroy)
        button.pack(side=RIGHT,anchor=S)
        scroll = Tkinter.Scrollbar(frame,command=text.yview,width=20)
        scroll.pack(side=RIGHT,ipady=80,anchor=W)
        text["yscrollcommand"] = scroll.set
        tk.mainloop()

    def clean(self):
        os.system('gmake clean') 
        os.system('cd PRE; gmake -f Makefile_C clean; gmake clean')
        

    def run(self):
	os.system('echo "Please email acbauer@eng.buffalo.edu with any problems"')
        #get system information so it is known which system the script is running on
        machine = os.uname()
        #print 'machine is ' + machine[0] + ' ' + machine[1] + ' '+machine[2]+ ' ' + machine[3] + ' '+machine[4]
        print 'Trying to submit a job on ' + machine[1]
	# if there is no topo file, quit
	if self.topomain.get() == '' or self.toposub.get() == '' or self.topomapset.get() == '' or self.topomap.get() == '':
            print 'Missing GIS information.  No job will be run.'
            return
	#elif os.access(self.topofile.get(), os.F_OK) == 0:
	#    print 'cannot find '+self.topofile.get()
	#    return
		
	str_numprocs = self.numprocs.get()
	if str_numprocs == '':
	    numprocs = 1
	else:
	    numprocs = int(str_numprocs)
	    
	if numprocs <= 0:
	    print 'numprocs must be greater than 0, it is ' + str(numprocs)
	    numprocs = 1

        if numprocs != 1 and numprocs != 2 and numprocs != 4 and numprocs != 8 and numprocs != 16 and numprocs != 32 and numprocs != 64:
            print 'wrong amount of processors!'
            return
            
        if numprocs > 8 and machine[1] == 'crosby':
            print 'you need to have access to the mp queues for that many processors on crosby, email acbauer@eng.buffalo.edu for help'
            return

        if numprocs > 32 and machine[1] == 'stills':
            print 'must use 32 or less processors on stills'
            return
        
	#create proper directory
	if os.access(self.directory.get(), os.F_OK) == 1:
	    print self.directory.get() + ' already exists, will not overwrite it'
	    return

	dir1 = self.directory.get()
        if dir1 == '':
            print 'Must specify a simulation directory.'
            return
        
	dirlen = len(dir1)
	if dir1[dirlen-1] != '/':
	    directory = dir1 + '/'
	else:
	    directory = dir1 

	os.system('mkdir ' + directory)

        #check values
        if self.numpiles.get() == '':
            numpiles = 0
	    print 'Number of piles not set.  Setting it to ' + str(numpiles)
        else:
            numpiles = int(self.numpiles.get())
            if numpiles < 0:
                print 'Number of piles cannot be a negative number.  Setting it to 0'
                numpiles = 0
                
	if self.intfrict.get() == '':
	    intfrict = 30.
	    print 'Internal friction angle not set.  Setting it to ' + str(intfrict)
	else:
	    intfrict = float(self.intfrict.get())
	    if intfrict <= 0:
		intfrict = 30.
                    
	if self.bedfrict.get() == '':
	    bedfrict = 15.
	    print 'Bed friction angle not set.  Setting it to ' + str(bedfrict)
	else:
	    bedfrict = float(self.bedfrict.get())
	    if bedfrict <= 0:
		bedfrict = 15

        # get all of the pile information
        f_p = directory+'simulation.data'
        f_p2=open(f_p, "w", 0)
        f_p2.write(str(numpiles) + '\n')
        f_p2.close
        counter = 0
        max_height = 0.0
        while counter < numpiles:
            root2=Tk()
            app=QuestionTemplate2(root2,counter,f_p,directory,max_height)
            root2.mainloop()
            #print 'counter is ' + str(counter) + ' and numpiles is ' + str(numpiles)
            counter = counter + app.input_flag
            max_height = app.max_height
            app.top.withdraw()
            del app
            
        print 'max height is ' + str(max_height)
        f=open("PRE/frict.data", "w", 0)
        f.write(str(intfrict) + '\n')
        f.write(str(bedfrict))
        f.close

	numygrid = 20
	if self.numygrid.get() != '':
	    numygrid = int(self.numygrid.get())
	    if numygrid <= 0:
		numygrid = 20

	if self.steps.get() == '':
	    steps = 100
	else:
	    steps = self.steps.get()
	    if steps < 1:
		steps = 100

	if self.numoutput.get() == '':
	    numoutput = 100
	else:
	    numoutput = self.numoutput.get()
	    if numoutput < 1:
		numoutput = 100

	if self.maxtime.get() == '':
	    maxtime = 1.5
	else:
	    maxtime = float(self.maxtime.get())
	    if maxtime <= 0:
		maxtime = 1.5
	
        #scaling stuff
        scale = self.scale.get()
        lengthscale = 1.
        if scale == 1:
	    if self.lengthscale.get() == '':
		lengthscale = 1.
	    else:
		lengthscale = float(self.lengthscale.get())
		if lengthscale <= 0.:
		    lengthscale = 1.
            if max_height <= 0.:
                max_height = 1.
                
	    output1 = str(lengthscale) + "\n" + str(max_height) + "\n9.80" 
	    
	else:
	    output1 = "1 \n1 \n1"
	    
	f2 = directory+ 'scale.data'
	f=open(f2, "w", 0)
	f.write(output1)
	f.close
	os.system('chmod a+wrx '+f2)

	#run the makefile (clean out all c and c++ object files first)
        if os.access('hpfem',os.X_OK) == 0:
            os.system('gmake cleanC')
            os.system('gmake')
            os.system('chmod a+wrx hpfem')

        #pile geometry stuff, friction coefficients, etc.
	output2 = str(steps) + '\n' + str(maxtime) + '\n' + str(numoutput) + '\n' + str(self.adapt.get())
        f_p2=open(f_p, "a+", 0)        
	f_p2.write(output2)
        #put in the stuff for viz the idea is to use prime numbers and the remainder function to determine which formats to output in
        viz_num = 1
        #print 'tecplot var is ' + str(self.tecplotVar.get())
        if self.tecplotVar.get() == 1:
            viz_num = viz_num * 2

        #print 'meshplot var is ' + str(self.mshplotVar.get())
        if self.mshplotVar.get() == 1:
            viz_num = viz_num * 3

        #print 'gmfg viz ' + str(self.padyVar.get())
        if self.padyVar.get() == 1:
            viz_num = viz_num * 5

        #print 'hdf var is ' + str(self.hdfVar.get())
        if self.hdfVar.get() == 1:
            viz_num = viz_num * 7

        f_p2.write('\n' + str(viz_num) + '\n' + str(self.order.get()+1))
        #GIS stuff
        f_p2.write('\n' + self.topomain.get() + '\n' + self.toposub.get() + '\n' + self.topomapset.get() + '\n' + self.topomap.get())
	f_p2.close
	os.system('chmod a+wrx '+f_p)

        # make a copy of the input file for reference
        if self.scale.get() == 1:
            str_scale = 'yes'
        else:
            str_scale = 'no'

        if self.adapt.get() == 1:
            str_adapt = 'yes'
        else:
            str_adapt = 'no'

        run_file = directory+'python_input.data'
        run_info = """NOTE:  This file outputs the data used in the runs, NOT the data input into the Python script
GIS Information Main Directory:  """ + self.topomain.get() + """
GIS Sub-Directory:  """ + self.toposub.get() + """
GIS Map Set:  """ + self.topomapset.get() + """
GIS Map:  """ + self.topomap.get() + """
Simulation Directory Location:  """ + self.directory.get() + """
Number of Processors:  """ + str(numprocs) + """
Computational Mesh Points in Y-Direction:  """ + str(numygrid) + """
Internal Friction Angle (deg):  """ + str(intfrict) + """
Bed Friction Angle (deg):  """+ str(bedfrict) + """
Number of Piles:  """ + str(numpiles) + """
Scale Simulation? """ + str_scale + """
If Scaled, Length Scale (m):  """ +str(lengthscale) + """
Maximum Number of Time Steps:  """ + str(steps) + """
Maximum Time (s):  """ + str(maxtime) +"""
Time Steps per Results Output:  """ + str(numoutput) + """
Adapt the Grid?  """ + str_adapt + """
"""
	f=open(run_file, "w", 0)
        f.write(run_info)
	f.close

	# run preproc.x to create the fem grid, if it is not already there
	if os.access('PRE/preproc.x',os.X_OK)==0:
            os.system('cd PRE;gmake')
            os.system('chmod a+wrx PRE/preproc.x')
	    
        #os.system('cd PRE;./preproc.x '+str(numygrid)+' /volcano/elchichon1/acbauer/grass.data/grass5 Colima ColimaR ColimaR')
        os.system('cd PRE;./preproc.x '+str(numygrid)+' '+self.topomain.get() +' '+self.toposub.get()+' '+self.topomapset.get()+' '+self.topomap.get())
	os.system('chmod a+wrx PRE/funky.dat')
	if os.access('PRE/preprocess',os.X_OK) != 1:
	    os.system('cd PRE;gmake -f Makefile_C')
	    os.system('chmod a+wrx PRE/preprocess')

	os.system('cd PRE;./preprocess '+str(numprocs))
	#os.system('cp -f ' + self.topofile.get() +' '+ directory +'/topo.data')
	os.system('mv PRE/funky*inp '+directory)
	os.system('cp hpfem ' + directory)
        #	os.system('mpirun -np ' + str(numprocs) + directory+'hpfem')
############################################################
############################################################
###  machine specific stuff here...
############################################################
############################################################
        if machine[1] == 'crosby' or machine[1] == 'nash.ccr.buffalo.edu' or machine[1] == 'joplin.ccr.buffalo.edu':
            iam = self.emailaddress.get()
            print 'sending email to : '+iam
        elif machine[1] == 'stills':
            print 'no email notification on stills'
        elif machine[0] == 'Linux':
            print 'no email notification on Linux machines'
            
        if machine[1] == 'crosby':            
            batchscript = """#!/bin/csh -f 
#PBS -m e
#PBS -l ncpus="""+ str(numprocs)+ """
#PBS -l walltime=50:00:00
#PBS -l mem=128mb
#PBS -M """+iam +"""
#PBS -o volcano.out.crosby
#PBS -j oe
#PBS -N volc."""+str(numprocs)+ """.crosby
#PBS
limit coredumpsize 0
#
# stage input files and executable to scratch - $PBS_O_WORKDIR is
#  the directory from which the job was submitted.
#
cp $PBS_O_WORKDIR/*.inp $PBSTMPDIR
cp $PBS_O_WORKDIR/hpfem $PBSTMPDIR
cp $PBS_O_WORKDIR/*.data $PBSTMPDIR
# goto scratch dir and run job
cd $PBSTMPDIR
date
mpirun -cpr -np """+str(numprocs)+"""  ./hpfem
date
#
# stage output files back
## cp *.out $PBS_O_WORKDIR/
gzip *out *gz *h5
# remove scratch directory
#rm -r $PBSTMPDIR
"""
            f6 = directory+'pbs_script'
            f=open(f6, "w", 0)
            f.write(batchscript)
            f.close
            os.system('cd '+directory+';qsub pbs_script')


#################################
# nash.ccr.buffalo.edu
#################################
        if machine[1] == 'nash.ccr.buffalo.edu':
            if numprocs ==1:
                b2 = "#PBS -l nodes=1:ppn=1"
            else:
                b2 = "#PBS -l nodes="+str(numprocs/2)+":ppn=2"
                
            batchscript = """#!/bin/csh
#PBS -m e
""" + b2 + """
#PBS -l walltime=30:00:00
#PBS -M """+iam +"""
#PBS -o volcano.out.nash
#PBS -j oe
#PBS -N volcano."""+str(numprocs)+ """.nash
#PBS
limit coredumpsize 0
#
cd $PBS_O_WORKDIR
set NP = `cat $PBS_NODEFILE | wc -l`
set NN = `cat $PBS_NODEFILE | wc -l`
echo "NN = "$NN
set plist = `cat $PBS_NODEFILE | sed "s/\.ccr\.buffalo\.edu//"`
set uplist = `cat $PBS_NODEFILE | uniq | sed "s/\.ccr\.buffalo\.edu//"`
echo "unique nodes list = "$uplist
echo "full process list = "$plist

@ NN = -1
foreach p ($plist)
 @ NN ++
  set NUM = `echo $NN | awk '{printf "%.4d", $1}'`
  rsh $p "cp $PBS_O_WORKDIR/funky$NUM.inp $PBSTMPDIR/."
end

foreach p ($uplist)
  rsh $p "cp $PBS_O_WORKDIR/hpfem $PBSTMPDIR/"
  rsh $p "cp $PBS_O_WORKDIR/*.data $PBSTMPDIR/"
  echo "Contents of "$PBSTMPDIR" on node "$p
  rsh $p "ls -l $PBSTMPDIR/"
end

source /util/tag-gm-gnu.csh
cd $PBSTMPDIR
#date
#/util/mpich-gm/gnu/ch_gm/bin/mpirun.ch_gm -v -machinefile $PBS_NODEFILE  -np """+str(numprocs)+"""  ./hpfem
#date
#
foreach p ($uplist)
  rsh $p "cd $PBSTMPDIR; gzip *out *plt *h5"
  rsh $p "cp $PBSTMPDIR/*gz $PBS_O_WORKDIR/"
end
"""
            f6 = directory+'pbs_script'
            f=open(f6, "w", 0)
            f.write(batchscript)
            f.close
            os.system('cd '+directory+';qsub pbs_script')
#################################
# joplin.ccr.buffalo.edu -- added by abani 02/27
################################
        elif machine[1] == 'joplin.ccr.buffalo.edu':
            if numprocs ==1:
                b2 = "#PBS -l nodes=1:ppn=1"
            else:
                b2 = "#PBS -l nodes="+str(numprocs/2)+":ppn=2"
                
            batchscript = """#!/bin/csh -f 
""" + b2 + """
#PBS -l walltime=00:59:00
#PBS -M """ + iam + """
#PBS -m e
#PBS -j oe
#PBS -o volcano."""+str(numprocs)+""".joplin
#
# Set executable name & mpiexec path
#  EXE = Executable name
#  INITDIR = dir holding input & executable
#  SAVEDIR = dir to place output files
#
source /util/tag-gm-gnu.csh
set EXE = hpfem
set INITDIR = $PBS_O_WORKDIR
set SAVEDIR = $PBS_O_WORKDIR
#
# Get list of processors from PBS

set plist = `cat $PBS_NODEFILE | sed "s/\.ccr\.buffalo\.edu//"`
set uplist = `cat $PBS_NODEFILE | uniq | sed "s/\.ccr\.buffalo\.edu//"`
echo "unique nodes list = "$uplist
echo "full process list = "$plist
#
# Stage out indexed input files (and common ones too) to each
#  processor
# (also build an mpiexec config file with the ranks to correspond
#  to the staged inputs, or so we hope - I still don't see that
#  as guaranteed!)
#
cd $PBS_O_WORKDIR
set conf = $PBSTMPDIR/conf.$$
if (-e $conf) rm -f $conf
@ NN = -1
foreach p ($plist)
 @ NN ++
  set NUM = `echo $NN | awk '{printf "%.4d", $1}'`
  rsh $p "cp $INITDIR/funky$NUM.inp $PBSTMPDIR"
  echo "$p : $EXE" >> $conf
end
#
#  files that only need to appear once on a node, no matter
#  how many procs
#
foreach p ($uplist)
  rsh $p "cp $INITDIR/$EXE $PBSTMPDIR/"
  rsh $p "cp $INITDIR/*.data $PBSTMPDIR/"
  echo "Contents of "$PBSTMPDIR" on node "$p
  rsh $p "ls -l $PBSTMPDIR/"
end
echo "Contents of mpiexec config file:"
cat $conf
#
# Run code from $PBSTMPDIR
#
cd $PBSTMPDIR
/usr/bin/time mpiexec -config $conf
#
# Stage back output files
#
foreach p ($uplist)
  rsh $p "gzip $PBSTMPDIR/*out; gzip $PBSTMPDIR/*plt; gzip $PBSTMPDIR/*h5; cp $PBSTMPDIR/*gz $SAVEDIR/"
end
rm -f $conf
"""
            f6 = directory+'pbs_script'
            f=open(f6, "w", 0)
            f.write(batchscript)
            f.close
            os.system('cd '+directory+';qsub pbs_script')

#################################
# stills.ccr.buffalo.edu - not working as of 5/7/03 but left in for future use/installation of proper libraries
#################################
        elif machine[1] == 'stills':
            #figure out the distribution of processors
            if numprocs <=4:
                nodes = 1
                procs_per_node = numprocs
                requirements = '(Memory>256)'
            else:
                nodes = numprocs/4
                procs_per_node = 4
                requirements = '(Pool==2)'

            batchscript = """###########################################################################
## A sample Loadl file. 
## Only run jobs in $LOADL_SCRATCH or /gpfs.
###########################################################################
## classes are:  Short, Medium, Long, V.Long
#!/bin/csh -xf 
# @ class =  V.Long
# @ environment = COPY_ALL
# @ error       = volcano"""+str(numprocs)+""".err..$(jobid)
# @ output      = volcano"""+str(numprocs)+""".out.$(jobid)
# @ network.MPI = css0,not_shared,US
# @ job_type = parallel
# @ requirements = """+requirements+"""
# @ tasks_per_node="""+str(procs_per_node)+"""
# @ node ="""+str(nodes)+"""
# @ notification = never
# @ queue
#
##########################################################################
##
##  Change to suit your needs

set EXECDIR=`pwd`
set EXECNAME=hpfem
set WORKDIR=`pwd`

if (-e $EXECDIR/EXECNAME) then
   echo ' Executible not found'
   echo ' Directory checked was $EXECDIR'
   echo ' EXE filename checked was $EXECNAME'
   exit1
else
   echo ' Found the EXE file $EXECDIR/$EXECNAME'
endif
set SCRATCH = $LOADL_SCRATCHDIR

############################################################################
##
## No need to change this. 
## Automatically sets up a proper scratch directory on all nodes. 
##

set SCRATCH = $LOADL_SCRATCHDIR
cd $SCRATCH

set LL_hosts = ($LOADL_PROCESSOR_LIST)
echo 'Number of processes is' $#LL_hosts
echo 'The scheduled processors are on nodes ' $LL_hosts
set noclobber


#############################################################################
## Change to suit
## Add to the copy list as necessary - never remove $EXECNAME copy
##      
mcp $EXECDIR/$EXECNAME $SCRATCH/.
mcp $EXECDIR/scale.data $SCRATCH/.
mcp $EXECDIR/topo.data $SCRATCH/.
mcp $EXECDIR/simulation.data $SCRATCH/.
"""
            f6 = directory+'loadleveler_script'
            f=open(f6, "w", 0)
            f.write(batchscript)
            for i in range(numprocs):
                if i < 10:
                    f.write("mcp $EXECDIR/funky0"+str(i)+".inp $SCRATCH/.\n")
                else:
                    f.write("mcp $EXECDIR/funky"+str(i)+".inp $SCRATCH/.\n")

            f.write("""############################################################################
## Change to suit
## Invoke POE
time poe << EOF
    $SCRATCH/$EXECNAME 
    quit
EOF

############################################################################
## Change to suit
## Copy output files back to $WORKDIR
## need to do some tricky stuff here because the output files are
## in scratch space on all of the nodes and only 1 node runs this script

foreach node ($LL_hosts)
	rcp "$node":"$SCRATCH/viz_output*.out" $WORKDIR/.
end
gzip viz*out

# ----------------------------------------------------------------
# Last second information before giving up the nodes/disks
#
date
""")
            f.close
            os.system('cd '+directory+';llsubmit loadleveler_script')
            
        elif machine[0] == 'Linux':  
            f8 = directory+'machines'
            f=open(f8, "w",0)
            f.write(machine[1])
            f.close
            i = 0
            while i == 0:
                i =  os.access(f_p, os.F_OK)
                                        
#            os.system('cd '+directory+';/usr/local/mpich-1.2.4/ch_p4/bin/mpirun -machinefile machines -np '+str(numprocs)+' hpfem')

        
root=Tk()
app=QuestionTemplate(root)
root.mainloop()
            




