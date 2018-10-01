from pylab import *
from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

N = [50,100,150,200,250,300,350,400,450]
J = [82*1e-3,1104*1e-3,5355.34*1e-3,16703*1e-3, 42353*1e-3,85811*1e-3,161355*1e-3,274910*1e-3, 485301*1e-3]
A = [0.82*1e-3,4.1*1e-3,9.522*1e-3,16.253*1e-3,24.231*1e-3,33.61*1e-3,45.888*1e-3,61.244*1e-3, 85.813*1e-3]
sims = [4112,16401,36570,64509,100077,143659,194687, 252e3, 319572]

f = 19

plot(N,J)
xlabel(r'$N$', fontsize=f)
ylabel('Time [s]', fontsize=f)
legend(['Jacobi Method'], fontsize=f)
plt.xticks(size=f)
plt.yticks(size=f)
tight_layout()
plt.savefig("J_time.png")
show()

plot(N,A)
xlabel(r'$N$', fontsize=f)
ylabel(r'Time [s]', fontsize=f)
legend([r'Armadillo'], fontsize=f)
plt.xticks(size=f)
plt.yticks(size=f)
tight_layout()
plt.savefig("A_time.png")

show()


plot(N,sims)
xlabel(r'$N$', fontsize=f)
ylabel(r'Number', fontsize=f)
legend(['Number of similarity transformations'], fontsize=f)
plt.xticks(size=f)
plt.yticks(size=f)
tight_layout()
plt.savefig("Sims.png")
show()
