#### IMPORT MODULES and Function definition
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

pi=np.pi
def shift(d,lambda_0):
    return np.exp(-1j*4*pi*d/lambda_0)

if __name__ == "__main__":

    ### Set the parameters
    lambda_0 = 1  # does not matter as we set w and dx based on lambda_0
    R = 2 # does not matter (just for visualization)
    w  = lambda_0/4
    dx = lambda_0/5
    sigma_1 = 1/2
    sigma_2 = 1/2
    theta_1 = -pi/3 # -60 deg
    theta_2 = +pi/6 # 30 deg

    ##Compute Ideal Signals
    sin_th1=np.sin(theta_1)
    sin_th2=np.sin(theta_2)

    rC1 = sigma_1+  sigma_2
    rH1 = sigma_1*shift(d=w*sin_th1,lambda_0=lambda_0)      +  sigma_2*shift(d=w*sin_th2,lambda_0=lambda_0)
    rC2 = sigma_1*shift(d=dx*sin_th1    ,lambda_0=lambda_0) +  sigma_2*shift(d=dx*sin_th2    ,lambda_0=lambda_0)
    rH2 = sigma_1*shift(d=(w+dx)*sin_th1,lambda_0=lambda_0) +  sigma_2*shift(d=(w+dx)*sin_th2,lambda_0=lambda_0) # not used for this analysis

    ## Estimation
    ux = (np.angle(rC1*np.conj(rH1)))*(lambda_0/(4*np.pi))*(1/w)
    uy = np.sqrt(1-ux**2)

    dr = (np.angle(rC1*np.conj(rC2)))*(lambda_0/(4*np.pi))
    dx_est=ux*dr/(ux**2)

    ## Visualize the resuts
    plt.figure(figsize=(10,5))
    plt.scatter([R*np.sin(theta_1),R*np.sin(theta_2)],[R*np.cos(theta_1),R*np.cos(theta_2)],color="red",marker="o")
    plt.text(x=R*np.sin(theta_1),y=R*np.cos(theta_1),s="target 1\n$(\sigma_1={:.2f})$".format(sigma_1),ha="left",va="bottom",color="red")
    plt.text(x=R*np.sin(theta_2),y=R*np.cos(theta_2),s="target 2\n$(\sigma_2={:.2f})$".format(sigma_2),ha="right",va="bottom",color="red")

    plt.quiver(0, 0, np.sin(theta_1), np.cos(theta_1), angles='xy', scale_units='xy', scale=1, color='red')

    plt.quiver(0, 0,  np.sin(theta_2), np.cos(theta_2), angles='xy', scale_units='xy', scale=1, color='red')

    plt.plot([R*np.sin(theta) for theta in np.linspace(-pi/2,pi/2,100)],[R*np.cos(theta) for theta in np.linspace(-pi/2,pi/2,100)], color="k",linestyle="--",alpha=0.5)
    plt.plot([r*np.sin(theta_1) for r in np.linspace(0,R,100)],[r*np.cos(theta_1) for r in np.linspace(0,R,100)], color="k",linestyle="--",alpha=0.5)
    plt.plot([r*np.sin(theta_2) for r in np.linspace(0,R,100)],[r*np.cos(theta_2) for r in np.linspace(0,R,100)], color="k",linestyle="--",alpha=0.5)
    plt.plot([r*np.sin(0) for r in np.linspace(0,R,100)],[r*np.cos(0) for r in np.linspace(0,R,100)], color="k",linestyle="--",alpha=0.5)
    plt.plot([.6*R*np.sin(theta) for theta in np.linspace(min(0,theta_1),max(0,theta_1),100)],
            [.6*R*np.cos(theta) for theta in np.linspace(min(0,theta_1),max(0,theta_1),100)], color="k",linestyle="-",alpha=0.5)
    plt.text(x=.5*R*np.sin(theta_1/2),y=.5*R*np.cos(theta_1/2),s="{:.2f}$^\circ$".format(theta_1*180/pi),ha="center",va="bottom",color="black",alpha=.5)

    plt.plot([.55*R*np.sin(theta) for theta in np.linspace(min(0,theta_2),max(0,theta_2),100)],
            [.55*R*np.cos(theta) for theta in np.linspace(min(0,theta_2),max(0,theta_2),100)], color="k",linestyle="-",alpha=0.5)
    plt.text(x=.45*R*np.sin(theta_2/2),y=.45*R*np.cos(theta_2/2),s="{:.2f}$^\circ$".format(theta_2*180/pi),ha="center",va="bottom",color="black",alpha=.5)

    arrow = FancyArrowPatch((0, 0),
                        (dx, 0),
                        arrowstyle='->',   # try 'simple', 'fancy', '<|-|>', '-[' ...
                        color='green',
                        mutation_scale=15, # overall size scaling
                        linewidth=1.0)
    plt.text(x=dx,y=0,s="$d_x$="+"{:.2f}".format(dx),ha="left",va="bottom",color="green")

    plt.gca().add_patch(arrow)
    plt.gca().set_aspect(1)
    plt.quiver(0, 0,  ux, uy, angles='xy', scale_units='xy', scale=1, color='purple')
    plt.text(x=ux,y=uy,s="$\hat{\mathbf{u}}$="+"[{:.2f},{:.2f}]$^T$".format(ux,uy),ha="left",va="bottom",color="purple")

    plt.scatter([R*ux,],[R*uy,],marker="o",color="purple")
    plt.plot([r*ux for r in np.linspace(0,R,100)],[r*uy for r in np.linspace(0,R,100)],linestyle="--",alpha=.5,color="purple")
    plt.text(x=R*ux,y=R*uy,s="detected\ntarget".format(ux,uy),ha="left",va="bottom",color="purple")

    plt.quiver(0, 0,  ux, uy, angles='xy', scale_units='xy', scale=1, color='purple')
    plt.text(x=ux,y=uy,s="$\hat{\mathbf{u}}$="+"[{:.2f},{:.2f}]$^T$".format(ux,uy),ha="left",va="bottom",color="purple")


    arrow = FancyArrowPatch((R*ux, R*uy),
                            (R*ux - dr*ux, R*uy - dr*uy),
                            arrowstyle='->',   # try 'simple', 'fancy', '<|-|>', '-[' ...
                            color='purple',
                            mutation_scale=15, # overall size scaling
                            linewidth=1.0)

    plt.gca().add_patch(arrow)
    plt.text(x=R*ux - dr*ux,y=R*uy - dr*uy,s="$d_r$="+"{:.2f}".format(dr),ha="left",va="top",color="purple")

    arrow = FancyArrowPatch((0, -.2),
                            (dx_est, -.2),
                            arrowstyle='->',   # try 'simple', 'fancy', '<|-|>', '-[' ...
                            linestyle="--",
                            color='green',
                            mutation_scale=15, # overall size scaling
                            linewidth=1.0)

    plt.gca().add_patch(arrow)
    plt.text(x=dx_est,y=-.2,s="$\hat{d}_x$="+"{:.2f}".format(dx_est),ha="left",va="center",color="green")
    plt.show()