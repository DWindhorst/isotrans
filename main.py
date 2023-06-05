# This is a sample Python script.
import os
import cmf

os.system('cls')
p = cmf.project()
# Create W1 in project p
W1 = p.NewStorage(name="W1")
# Create W2 in project p without any volume as an initial state
W2 = p.NewStorage(name="W2")
# Create a linear storage equation from W1 to W2 with a residence time tr of one day
q = cmf.LinearStorageConnection(source=W1,target=W2,residencetime=1.0)
# Set the initial state of w1 to 1m³ of water.
W1.volume = 1.0


# Create an integrator (solver) for the ODE represented by project p,
# with an error tolerance of 1e-9
solver = cmf.RKFIntegrator(p, 1e-9)
# Import Python's datetime module
import datetime
# Set the intitial time of the solver
solver.t = datetime.datetime(2012,1,1)

# Iterate the solver hourly through the time range and return for each time step the volume in W1 and W2
result = [[W1.volume,W2.volume] for t in solver.run(datetime.datetime(2012,1,1),datetime.datetime(2012,1,7),datetime.timedelta(hours=1))]
import pylab as plt
plt.plot(result)
plt.xlabel('hours')
plt.ylabel('Volume in $m^3$')
plt.legend(('W1','W2'))
plt.show()
