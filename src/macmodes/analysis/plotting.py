import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm_y

def non_uniform_first_deriv(f, grid_params):
    dr, r_1, r_2, dr1a, dr1b, dr1c, dr2a, dr2b, dr2c = grid_params
    fp = np.empty_like(f)
    for ii in range(1, len(f)-1):
        fp[ii] = dr1a[ii]*f[ii-1] + dr1b[ii]*f[ii] + dr1c[ii]*f[ii+1]
    fp[0] = dr1a[0]*f[0] + dr1b[0]*f[1] + dr1c[0]*f[2]
    fp[-1] = dr1a[-1]*f[-1] + dr1b[-1]*f[-2] + dr1c[-1]*f[-2]
    return fp

# Scipy spherical harmonics
def my_sph_harm(l, theta, m=0, phi=0):
    # Returns spherical harmonic of degree l, order m and first derivative
    return sph_harm_y(l, m, theta, phi, diff_n=1)

# Go from potential scalars to velocity field
def mode_to_spatial_MAC(mode_ind, r, grid_params, n_rad, n_theta, ellmax, evecs):
    
    ellmax_h = ellmax//2 # number of harmonics
    mode = evecs[:, mode_ind-1]

    theta = np.linspace(0, np.pi, n_theta) # colatitudinal grid
    dtheta = np.pi/(n_theta-1) # grid spacing
    
    vr = np.zeros((n_rad, n_theta), dtype=complex) # radial velocity 
    vtheta = np.zeros((n_rad, n_theta), dtype=complex) # meridional velocity 
    vphi = np.zeros((n_rad, n_theta), dtype=complex) # zonal velocity
    bphi = np.zeros((n_rad, n_theta), dtype=complex) # zonal magnetic field 
    
    NW = ellmax_h*n_rad # index where W part of eigenvector begins
    NZ = 3*ellmax_h*n_rad # index where Z part of eigenvector begins
    Nt = 4*ellmax_h*n_rad # index where tau part of eigenvector begins
    
    for nl in range(ellmax_h):
        ello = 2*nl+1
        ellpo = ello*(ello+1)
        elle = 2*nl+2
        ellpe = elle*(elle+1)
        
        W = mode[NW+nl*n_rad:NW+(nl+1)*n_rad]
        dW = non_uniform_first_deriv(W, grid_params)
        Z = mode[NZ+nl*n_rad:NZ+(nl+1)*n_rad]
        tau = mode[Nt+nl*n_rad:Nt+(nl+1)*n_rad]
        
        Ye, dYe_ = my_sph_harm(elle, theta) # even degree spherical harmonic and derivative
        Yo, dYo_ = my_sph_harm(ello, theta) # odd degree spherical harmonic and derivative
        dYe = dYe_[:,0]
        dYo = dYo_[:,0]

        for n in range(n_theta): 
            for nr in range(n_rad):
                vr[nr, n] += ellpe/r[nr]**2 * W[nr] * Ye[n]
                vtheta[nr, n] += 1/r[nr]**2 * dW[nr] * dYe[n]
                vphi[nr, n] += -1/r[nr] * Z[nr] * dYo[n]
                bphi[nr, n] += -1/r[nr] * tau[nr] * dYo[n]

    return theta, vr, vtheta, vphi, bphi

def contour_plot_mac(mode_ind, r, grid_params, n_rad, n_theta, ellmax, evecs):
    theta, vr, vtheta, vphi, bphi = mode_to_spatial_MAC(mode_ind, r, grid_params, n_rad, n_theta, ellmax, evecs)
    deg = theta*180/np.pi
    levels=50

    ellmax_h = ellmax//2 # number of harmonics
    mode = evecs[:, mode_ind-1]
    kd = 2*ellmax_h*n_rad
    kmid = kd + n_rad//2
    phase_d = np.atan2(mode[kmid].imag, mode[kmid].real)
    amp_d = np.abs(mode[kmid])

    fig1, axs1 = plt.subplots(3,1,sharex=True, figsize=(6, 6))
    axs1[0].set_xticks([0, 30, 60, 90, 120, 150, 180])
    c1a = axs1[0].contourf(deg, r, np.real(vr), levels=levels, cmap='jet')
    axs1[0].set_ylabel(r'$Re(v_r)$')
    cbar1a = fig1.colorbar(c1a)
    c1b = axs1[1].contourf(deg, r, np.real(vtheta), levels=levels, cmap='jet')
    axs1[1].set_ylabel(r'$Re(v_{\theta})$')
    cbar1b = fig1.colorbar(c1b)
    c1c = axs1[2].contourf(deg, r, np.real(vphi), levels=levels, cmap='jet')
    axs1[2].set_ylabel(r'$Re(v_{\phi})$')
    cbar1c = fig1.colorbar(c1c)
    fig1.supxlabel('Colatitude (degree)')
    fig1.supylabel('Radius ($r/r_f$)')

    fig2, axs2 = plt.subplots(3,1,sharex=True, figsize=(6, 6))
    axs2[0].set_xticks([0, 30, 60, 90, 120, 150, 180])
    c2a = axs2[0].contourf(deg, r, vr.imag, levels=levels, cmap='jet')
    axs2[0].set_ylabel(r'$Im(v_{r})$')
    cbar2a = fig2.colorbar(c2a)
    c2b = axs2[1].contourf(deg, r, vtheta.imag, levels=levels, cmap='jet')
    axs2[1].set_ylabel(r'$Im(v_{\theta})$')
    cbar2b = fig2.colorbar(c2b)
    c2c = axs2[2].contourf(deg, r, vphi.imag, levels=levels, cmap='jet')
    axs2[2].set_ylabel(r'$Im(v_{\phi})$')
    cbar2c = fig2.colorbar(c2c)
    fig2.supxlabel('Colatitude (degree)')
    fig2.supylabel('Radius ($r/r_f$)')

    fig3, axs3 = plt.subplots(2,1,sharex=True, figsize=(6, 6))
    c3a = axs3[0].contourf(deg, r, bphi.real, levels=levels, cmap='jet')
    axs3[0].set_ylabel(r'$Re(b_{\phi})$')
    cbar3a = fig3.colorbar(c3a)
    c3b = axs3[1].contourf(deg, r, bphi.imag, levels=levels, cmap='jet')
    axs3[1].set_ylabel(r'$Im(b_{\phi})$')
    cbar3b = fig3.colorbar(c3b)
    fig3.supxlabel('Colatitude (degree)')
    fig3.supylabel('Radius ($r/r_f$)')

    fig4, [ax1, ax2, ax3] = plt.subplots(3, 1, figsize=(6,6))
    c1 = ax1.contourf(deg, r, np.real(vr*np.exp(-1j*phase_d)), levels=levels, cmap='jet')
    ax1.set_xticks([0, 30, 60, 90, 120, 150, 180])
    ax1.set_ylabel(r'$v_r$')
    cbar1 = fig4.colorbar(c1)
    c2 = ax2.contourf(deg, r, np.real(vphi*np.exp(-1j*phase_d)), levels=levels, cmap='jet')
    ax2.set_xticks([0, 30, 60, 90, 120, 150, 180])
    ax2.set_ylabel(r'$v_{\phi}$')
    cbar2 = fig4.colorbar(c2)
    c3 = ax3.contourf(deg, r, np.real(bphi*np.exp(-1j*phase_d)), levels=levels, cmap='jet')
    ax3.set_xticks([0, 30, 60, 90, 120, 150, 180])
    ax3.set_ylabel(r'$b_{\phi}$')
    cbar3 = fig4.colorbar(c3)
    fig4.supxlabel('Colatitude (degree)')
    fig4.supylabel('Radius ($r/r_f$)')

    return fig1, fig2, fig3, fig4

# In directory mode_folder, saves images of velocty fields and different points in phase
def vr_phases(mode_ind, mode_folder, evecs, n_rad, r, n_theta, ellmax, grid_params):

    theta, vr, vtheta, vphi, bphi = mode_to_spatial_MAC(mode_ind, evecs, n_rad, r, n_theta, ellmax, grid_params)
    deg = theta*180/np.pi
    levels=30

    num_frames = 10
    phases = np.linspace(0, 2*np.pi, num_frames)
    filename = os.path.join(mode_folder, f"mode_{mode_ind}_flow")
    os.makedirs(filename, exist_ok=True)

    for ii in range(len(phases)):
        fig1, axs1 = plt.subplots(3,1,sharex=True, figsize=(6, 6))
        axs1[0].set_xticks([0,   30, 60, 90, 120, 150, 180])
        c1a = axs1[0].contourf(deg, r, np.real(vr*np.exp(-1j*phases[ii])), levels=levels, cmap='jet')
        axs1[0].set_ylabel(r'$v_r$')
        cbar1a = fig1.colorbar(c1a)
        c1b = axs1[1].contourf(deg, r, np.real(vtheta*np.exp(-1j*phases[ii])), levels=levels, cmap='jet')
        axs1[1].set_ylabel(r'$v_{\theta}$')
        cbar1b = fig1.colorbar(c1b)
        c1c = axs1[2].contourf(deg, r, np.real(vphi*np.exp(-1j*phases[ii])), levels=levels, cmap='jet')
        axs1[2].set_ylabel(r'$v_{\phi}$')
        cbar1c = fig1.colorbar(c1c)
        fig1.supxlabel('Colatitude (degree)')
        fig1.supylabel('Radius ($r/r_f$)')
        fig1.suptitle(rf'$\omega t = {phases[ii]:.3f}$')
        path = os.path.join(filename, f"vr_{ii+1:05d}")
        fig1.savefig(path)
        plt.close()

