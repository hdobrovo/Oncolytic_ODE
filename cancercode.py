from optparse import OptionParser
from pylab import *
import numpy as np

parser = OptionParser()
parser.add_option("-f", "--betavslambdadata.dat", metavar="FILE", default="combined", help="read data file FILE [combined]")
(options, args) = parser.parse_args()
del parser, args

#[rl,rb,max(vl),(dur[-1]-dur[0])/100,ucl[-1]/1e6,uncl[-1]/1e6]
data = np.loadtxt("bvsde2.dat")
print(shape(data))
rd = data[:,1] # Find y data
rk = data[:,0] # Find x data
dst = data[:,5] # Find z data

dose = compress( rk == rk[0], rd ) # Find unique x values
time = compress( rd == rd[0], rk ) # Find unique y values
z = reshape( dst, (len(dose), len(time)) ) # Reshape z into matrix

fig = plt.figure()
ax = fig.add_axes([0.15, 0.15, 0.7, 0.7])

#plt.plot(log10(time), log10(176*10**6*4.26*10**(-4)*49/(0.13*(1+dose))-1),'k',linewidth=3)
#print(dose)
#print(0.046*3.2*10**(-5)*4*10**8/5.2-dose-5.2)
#plt.imshow( np.transpose(abs(z)/10000), aspect='auto',vmin=min(abs(dst)/10000), vmax=max(abs(dst)/10000) ,origin='lower', cmap=cm.jet, extent=(-3,3,-3,3), axes=ax)
plt.imshow( z, aspect='auto',vmin=0, vmax=950 ,origin='lower', cmap=cm.jet, extent=(-1,5,-1,5), axes=ax)

yticks( arange(-1,6,1), ('$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'), fontsize=20 )
xticks( arange(-1,6,1), ('$10^{-6}$','$10^{-5}$','$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$'), fontsize=20)

xlabel(r'Ratio of infection rates, $r_\beta$', fontsize=20)
ylabel(r'Ratio of death rates, $r_\delta$', fontsize=20)

#text(-1.75, 10.25, r'Peak viral titer (upper)', fontsize=24)

ax = fig.add_axes([0.89, 0.15, 0.03, 0.7])
cb =colorbar(cax=ax ,orientation='vertical')

#labels = ['$ $', '50', '$ $', '100', '$ $', '150', '$ $',  '200' ]
#cb.ax.set_yticks([0, 2, 4, 6, 8])
#ax.set_yticklabels(labels)
text(-1.3, 200, r'Non-cancer cells', fontsize=20, rotation=90)
yticks(fontsize=20)
savefig("figures/bvsde2_noncancer.pdf", bbox_inches='tight', pad_inches=0.2)
show()
