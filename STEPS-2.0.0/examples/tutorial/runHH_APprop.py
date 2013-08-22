# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# Hodgkin-Huxley Action Potential propagation model, plotting functions
# Author Iain Hepburn

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

from pylab import *

from HH_APprop import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

tpnt = arange(0.0, N_timepoints*DT_sim, DT_sim)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

def plotVz(tidx):
    if (tidx >= tpnt.size): 
        print 'Time index out of range'
        return
    plot(results[1]*1e6, results[0][tidx], \
         label=str(1e3*tidx*DT_sim)+'ms', linewidth=3)
    legend(numpoints=1)
    xlim(0, 1000)
    ylim(-80,40)
    xlabel('Z-axis (um)')
    ylabel('Membrane potential (mV)')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

def plotIz(tidx, plotstyles = ['-', '--']):
    if (tidx >= tpnt.size): 
        print 'Time index out of range'
        return
    plot(results[4]*1e6, results[2][tidx], plotstyles[0],\
         label = 'Na: '+str(1e3*tidx*DT_sim)+'ms', linewidth=3)
    plot(results[4]*1e6, results[3][tidx], plotstyles[1],\
         label = 'K: '+str(1e3*tidx*DT_sim)+'ms', linewidth=3)
    legend(loc='best')
    xlim(0, 1000)
    ylim(-10, 15)
    xlabel('Z-axis (um)')
    ylabel('Current  (pA/um^2)')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# END
