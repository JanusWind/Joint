from numpy import mean, std

anglefile = open("theta_phi_reduced_2.txt", "r")
data = anglefile.read()
data = data.split(" ")
data[-1]=data[-1].split("\n")[0]

angles = [float(datum) for datum in data]

fc_theta = []
fc_phi = []
pl_theta = []
pl_phi = []

for i in range(len(angles)/4) :

	fc_theta += [angles[0]]
	fc_phi   += [angles[1]]
	pl_theta += [angles[2]]
	pl_phi   += [angles[3]]

	angles = angles[4::]

for i in range(len(fc_theta)-1):

	print i, fc_theta[i+1]-fc_theta[i]

dtheta = [round(fc_theta[i]-pl_theta[i], 2) for i in range(len(fc_theta))]
dphi   = [round(fc_phi[i]-pl_phi[i], 2) for i in range(len(fc_phi))]

print "dtheta = {}".format(dtheta)

print "mean = {}, std = {} ".format(mean(dtheta), std(dtheta))

print "dphi = {}".format(dphi)

print "mean = {}, std = {}".format(mean(dphi), std(dphi))

