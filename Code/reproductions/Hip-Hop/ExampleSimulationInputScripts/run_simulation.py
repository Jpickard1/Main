
import random
import sys, getopt
import numpy
import types
from lammps import lammps

def idx( id ):
    return 3*id-3

def idy( id ):
    return 3*id-2

def idz( id ):
    return 3*id-1

def sep(l,r,coords,box) :
    dx = coords[idx(l)]-coords[idx(r)]
    if ( abs(dx)>0.5*box[0] ) : 
        dx=box[0]-abs(dx)
    dy = coords[idy(l)]-coords[idy(r)]
    if ( abs(dy)>0.5*box[1] ) : 
        dy=box[1]-abs(dy)
    dz = coords[idz(l)]-coords[idz(r)]
    if ( abs(dz)>0.5*box[2] ) : 
        dz=box[2]-abs(dz)
    if ( numpy.sqrt( dx*dx + dy*dy + dz*dz )>50.0 ) :
        print "********************",dx,dy,dz
    return numpy.sqrt( dx*dx + dy*dy + dz*dz )

quiet = False                          # flag to switch off screen output
runnumber = 1                          # an id number for input and output files

try:
    opts, args = getopt.getopt(sys.argv[1:],"qn:",)
except getopt.GetoptError:
    print 'test.py [-q] [-n <runnumber>]'
    print '      -q runs quietly without output to screen'
    print '      -n <runnumber> specifies a run number, default 1.'
    sys.exit(2)
for opt, arg in opts:
    if (opt=="-q"):
        quiet = True
    elif (opt=="-n"):
        runnumber = int(arg)


# Brownian time for a DNA bead is 2 simulation time units

# parameter values
dt = 0.01                                             # time step in simulation time units
septhresh=6.0                                         # max distance for bond creation
runtime = 50e4                                        # total tun time in simulation time units
Nbeads = 5000                                         # number of DNA beads
extruders_per_bead = 30.0/5000.0                      # to set total number of extruders
extrude_rate = 1.0/500.0                              # rate at which extruders move in inverse simulation time units
update_interval = 10.0                                # how often to update extruders in simulation time units
onrate = 1.0/1000.0                                   # rate at which extruders attach in inverse simulation time units
offrate = 1.0/40000.0                                 # rate at which extruders detach in inverse simulation time units
mapfile = "extruders_locations_%i.dat"%runnumber      # output file for extruders map
switch_rate = 1.0/8000.0                              # in inverse LJ times


# set some variables
runtime_steps = int(runtime/dt)                            # runtime in timesteps
Nextruders = int(Nbeads*extruders_per_bead)                # number of extruders in system
update_interval_steps = int(update_interval/dt)            # how often to do an extruder update in timesteps
Tintervals = int(runtime_steps/update_interval_steps)-1    # number of times extruders will be updated
extrude_rate_steps = extrude_rate*dt                       # rate at which extruders move in inverse timesteps
onrate_steps = onrate*dt                                   # rate at which extruders attach in inverse timesteps
offrate_steps = offrate*dt                                 # rate at which extruders detach in inverse timesteps
switch_steps = switch_rate*dt                              # rate at which prote

step_prob = update_interval_steps*extrude_rate_steps    # probability an extruder advances
add_prob = update_interval_steps*onrate_steps           # probability an extruder attaches
remove_prob = update_interval_steps*offrate_steps       # probability an extruder detaches
switch_prob = update_interval_steps*switch_steps        # probability a protein 

# set up arrays for extruders
left = []
right = []

# random numbers for extruders
seedfile = open("seed.%i"%runnumber,"r")  # load a seed from an external file, different for each run
seedline = seedfile.read()
seedfile.close()
seedline = seedline.split(" ")
seed = int(seedline[3])
random.seed(seed)

# set up arrays for beads, and where CTCFs are
occupied=[0 for i in range(Nbeads+1)]

# load CTCF information from external file
sitesCTCF=["0" for i in range(Nbeads+1)]
scoresCTCF=[0.0 for i in range(Nbeads+1)]
execfile('setCTCFscores.py')
# set up CTCFs
CTCF=["0" for i in range(Nbeads+1)]
for i in range(0,len(scoresCTCF)):
    if (random.random() < scoresCTCF[i]):
        CTCF[i]=sitesCTCF[i]
# if a peak is BOTH, then select F,R or B with 1:1:2 ratio 
for i in range(Nbeads+1):
    if CTCF[i]=="B" and random.random() < 0.5:
        if random.random() < 0.5:
            CTCF[i]="F"
        else:
            CTCF[i]="R"
ctcfoutput = open("CTFCmap_%i.dat"%runnumber,"w")
for i in range(1,len(CTCF)):
    ctcfoutput.write("%i %s\n"%(i,CTCF[i]))
ctcfoutput.close()

# set up lammps and load initial input file
lmp = lammps(cmdargs=["-screen","screen.%i.log"%runnumber,"-log","log.%i.log"%runnumber])
# pass 'runnumber' to lammps for output names
lmp.command("variable runnumber equal %i"%runnumber)
lmp.file("in.initial")
lmp.command("timestep %f"%dt)
boxsize = [ lmp.extract_global("boxxhi",1)-lmp.extract_global("boxxlo",1) , lmp.extract_global("boxyhi",1)-lmp.extract_global("boxylo",1), lmp.extract_global("boxzhi",1)-lmp.extract_global("boxzlo",1) ]

# some messages
if (quiet == False) :
    ostreams = [sys.stdout, open("extruders_%i.log"%runnumber,"w")]
else :
    ostreams = [open("extruders_%i.log"%runnumber,"w")]

for messages in ostreams:
    messages.write("#******************************************************************************************************************************************\n")
    messages.write("#  Running simulation %i with %i polymer beads and %i extruders for %i time steps\n"%(runnumber,Nbeads,Nextruders, runtime_steps))
    messages.write("#\n")
    messages.write("#  extruders are updated every %i timesteps\n"%update_interval_steps)
    messages.write("#  extruders advance every %i timesteps\n"%(1.0/extrude_rate_steps))
    messages.write("#  extruders detach after %f timesteps\n"%(1.0/offrate_steps))
    messages.write("#  extruders detach after travelling %f beads\n"%(extrude_rate/offrate))
    messages.write("#  free extruders attach after %f timesteps\n"%(1.0/onrate_steps))
    messages.write("#\n")
    messages.write("#  step_prob = %f\n"%step_prob)
    messages.write("#  add_prob = %f\n"%add_prob)
    messages.write("#  remove_prob = %f\n"%remove_prob)
    messages.write("#  switch_prob = %f\n"%switch_prob)
    messages.write("#  random seed = %i\n"%seed)
    messages.write("#******************************************************************************************************************************************\n\n")

messages.write("#\n# step, number of extruders\n")

# set up output file for extruders
allextrudefile = open(mapfile,"w")
allextrudefile.write("# timestep, left side, right side")


# set up switching rate
lmp.command("variable k equal %f "%switch_prob)


# do an initial run without extruders
lmp.command("run %i"%update_interval_steps)

messages.write("%i %i\n"%(0,len(left)))

# run subsequent steps, updating according to step_prob
for step in range(0,Tintervals):
    if (quiet == False):
        print "running ",step*update_interval_steps,"   ( ", 100.0*(step)/float(Tintervals)," % done )     \r",
        sys.stdout.flush()

    coords = lmp.gather_atoms("x",1,3)

    # output extruders positions
    for i in range(0,len(left)): 
        allextrudefile.write("%i %i %i\n" %(step*update_interval_steps,left[i],right[i]) )

    # remove extruders with rate remove_prob
    for j in range(0,len(left)):
        if ( random.random() < remove_prob ):
            lmp.command( "group newbond id %i %i" % (left[j],right[j]) )   # delete the old bonds
            lmp.command("delete_bonds newbond bond 3 remove special")
            lmp.command("group newbond delete")
            occupied[left[j]]=0
            occupied[right[j]]=0
            left[j]=-1
            right[j]=-1
    left[:] = (value for value in left if value != -1)
    right[:] = (value for value in right if value != -1)


    # update extruders
    for j in range(0,len(left)): # left and right get updated here, but the length and order do not change
        lastleft = list(left)
        lastright = list(right)
        if ( left[j]-1>0 and occupied[left[j]-1]==0 and CTCF[left[j]-1]!="F" and CTCF[left[j]-1]!="B" and random.random() < step_prob ):
            occupied[left[j]] = 0
            left[j] -= 1
            occupied[left[j]] = 1
        if ( right[j]+1<=Nbeads and occupied[right[j]+1]==0 and CTCF[right[j]+1]!="R" and CTCF[right[j]+1]!="B" and random.random() < step_prob):
            occupied[right[j]] = 0
            right[j] += 1
            occupied[right[j]] = 1
        if ( lastleft[j]!=left[j] or lastright[j]!=right[j] ):
            if ( sep(left[j],right[j],coords,boxsize) > septhresh ):
                print "ERROR - bond cannot be made. Separation is ",sep(left[j],right[j],coords,boxsize),boxsize[0],boxsize[1],boxsize[2]
                print left[j],right[j],sep(left[j],right[j],coords,boxsize),"update"
                sys.exit()
            lmp.command("group newbond id %i %i" % (lastleft[j],lastright[j]) ) # delete the old bonds
            lmp.command("delete_bonds newbond bond 3 remove special") # need special so that bond lists are recalculated properly
            lmp.command("group newbond delete")
            lmp.command("group newbond id %i %i" % (left[j],right[j]) ) # add the new bonds
            lmp.command("create_bonds many newbond newbond 3 0.0 %f"%septhresh)
            lmp.command("group newbond delete")

    # add extruders at rate add_prob
    for i in range(0,Nextruders-len(left)):
        if ( random.random() < add_prob ):
            where = random.randrange(1,Nbeads-2-1)
            if ( occupied[ where ]==0 and occupied[ where+1 ]==0 and occupied[ where+2 ]==0 and occupied[ where+3 ]==0 and CTCF[ where ]=="0" and CTCF[ where+1 ]=="0" and CTCF[ where+2 ]=="0" and CTCF[ where+3 ]=="0" ):
                left.append(where)
                right.append(where+3)
                occupied[ left[-1] ]=1
                occupied[ right[-1] ]=1
                lmp.command("group newbond id %i %i" % (left[-1],right[-1]) ) # add the new bonds
                lmp.command("create_bonds many newbond newbond 3 0.0 %f"%septhresh)
                lmp.command("group newbond delete")

    if ( len(left)>Nextruders ) :
        print "ERROR - too many extruders"
        sys.exit()

    # do switching of proteins with lammps commands
    lmp.command("group atac_off type 6")
    lmp.command("group atac_on type 5")
    lmp.command("set group atac_off type/fraction 5 ${k} %i"%random.randrange(1000,99999))
    lmp.command("set group atac_on type/fraction 6 ${k} %i"%random.randrange(1000,99999))
    lmp.command("group atac_off delete")
    lmp.command("group atac_on delete")

    # do the run
    lmp.command("run %i"%update_interval_steps)
    messages.write("%i %i\n"%((step+1)*update_interval_steps,len(left)))

 
if ( quiet == False ):
    print ""
    print "Sucsess!"

allextrudefile.close()

lmp.command("write_restart end.%i.restart"%runnumber)
lmp.close()  




