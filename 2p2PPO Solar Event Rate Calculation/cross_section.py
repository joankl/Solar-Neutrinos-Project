'''
Python Module build to compute the cross-section values
of the elastic scattering of neutrino with electron. It
contain the cross-section calculation functions with and
without QED corrections.
'''

import numpy as np
from scipy.special import spence
import scipy.integrate as integrate


def cross_sec_ES_nu_e(E_nu, flavor = 'nue'):
    
    """
    Compute the cross-section of neutrino-electron ES without radiative 
    corrections by analytically integrating dσ/dT between the recoil electron
    energies Tmin and Tmax for a given neutrino energy.

    Parameters:
    -E_nu: scalar or array of neutrino energies (MeV)
    -flavor: neutrino flavor to compute the interactions cross-section.
                Default is electron neutrino flavor.

    returns: sigma(E) in cm2
    """

    E_nu = np.asarray(E_nu, dtype=float)

    # -------- constants --------

    pi = np.pi
    G_F = 1.16639 * 10 **(-11)               # Fermi Constant in MeV^(-2)
    h_bar = (4.135668 * 10**(-21))/(2*pi)   # Plank bar constant in MeV.s
    m_e = 0.5109989461                                          # electron mass energy in MeV/c^2 
    c = 3 * 10 ** (10)                                          # speed of light cm/s
    sigma_0 = (2 * (G_F **2) * (m_e**2)/pi) * (h_bar * c)**2    # cm^2
    g_R = 0.23116                                               # sin^2(Weinberg angle)
    
    if flavor == 'nue':
        g_L = g_R + 1/2
    
    elif flavor == 'numu' or nu_flavor == 'nutau':
        g_L = g_R - 1/2 
        
    else:
        print('Unknown neutrino flavor')
    
    
    # Recoil Electron Energy values:
    Tmax = 2.0 * E_nu**2 / (m_e + 2.0 * E_nu) # Maximum Kinetic energy of recoil electron
    #Tmin = integrate.quad(Tmin_func, 0, 2*np.pi, args=(E_nu))[0] # Minimum kinetic energy of recoil electron (Integration over theta)
    Tmin = 0
    # -------- Cross-section Computation --------

    # coefficients from algebra of the cross-section
    C0 = (g_L**2 + g_R**2)
    A1 = (-2.0 * g_R**2 / E_nu) - (g_L * g_R * m_e / (E_nu**2 + 1e-300))  # avoid div by 0
    A2 = g_R**2 / (E_nu**2 + 1e-300)

    # compute differences of integration
    dT = Tmax - Tmin
    dT2 = (Tmax**2 - Tmin**2)
    dT3 = (Tmax**3 - Tmin**3)

    integral = C0 * dT + 0.5 * A1 * dT2 + (1.0/3.0) * A2 * dT3
    
    # Final result
    sigma = (sigma_0 / m_e) * integral
    
    return sigma


def dSigma_dT_corr(T, E_nu, flavor = 'nue'):
    
    '''
    Function that defines the differential cross section equation
    as a function of the recoil electron energy with QED radiative
    corrections.
    
    Parameters:
    -E_nu: scalar of neutrino energies (MeV)  
    -flavor: neutrino flavor to compute the interactions cross-section.
                Default is electron neutrino flavor.
    
    return: The value of the differential equation at each point
    
    NOTE: The Spence function function in scipy is different than 
          the convention in arXiv:astro-ph/9502003, so L(x) must be 
          written as L(1-x)
    '''
    
    E_nu = np.asarray(E_nu, dtype=float)
    
    # ---------- Constants ----------
    pi = np.pi
    G_F = 1.16639 * 10 **(-11)               # Fermi Constant in MeV^(-2)
    h_bar = (4.135668 * 10**(-21))/(2*pi)    # Plank bar constant in MeV.s
    c = 3 * 10 ** (10)                       # speed of light cm/s
    m_e = 0.5109989461                       # electron mass energy in MeV/c^2 
    alpha = 1/137                            # Fine structure constant
    sin_theta_w_sqr = 0.23116                # sin^2(Weinberg angle)
    sigma_0 = (2 * (G_F **2) * (m_e**2)/pi) * (h_bar * c)**2
    
    # ---------- Recoil Electron Energy values ----------
    
    #Tmax = 2.0 * E_nu**2 / (m_e + 2.0 * E_nu)                    # Maximum Kinetic energy of recoil electron
    #Tmin = integrate.quad(Tmin_func, 0, 2*np.pi, args=(E_nu))[0] # Minimum kinetic energy of recoil electron (Integration over theta)
    #T = np.linspace(Tmin, Tmax, 200)
    
    # -------- g_L(T) and g_R(T) functions for nu_e --------
    
    rho_NC = 1.0126
    x = np.sqrt(1 + 2*m_e/T)
    I = (1/6) * (1/3 + (3 - x**2) * (0.5 * x * np.log((x + 1)/(x - 1)) - 1))

    if flavor == 'nue':
        k = 0.9791 + 0.0097*I 
        g_L = rho_NC * (0.5 - k * sin_theta_w_sqr) - 1
        g_R = -rho_NC * k * sin_theta_w_sqr
    
    elif flavor == 'numu' or nu_flavor == 'nutau':
        k = 0.9970 - 0.00037*I 
        g_L = rho_NC * (0.5 - k * sin_theta_w_sqr)
        g_R = -rho_NC * k * sin_theta_w_sqr

    else:
        print('Unknown neutrino flavor')
    
    
    # -------- QED functions f+, f-, f+/- --------
    
    z = T / (E_nu + 1e-300)
    E_e = T + m_e                 #Electron energy = kinetic energy + rest mass energy
    p_e = np.sqrt(E_e**2 - m_e**2)  # Electron 3-momentum modulo
    beta = p_e / E_e
    
    epsilon = 1e-100  #Logaritmic regularizer
    
    # f-
    term1_minus = ((E_e/p_e)*np.log((E_e + p_e)/m_e) - 1)*(2*np.log(np.maximum(1 - z - m_e/(E_e + p_e), epsilon)) - np.log(np.maximum(1 - z, epsilon)) - 0.5*np.log(np.maximum(z, epsilon)) - (5/12))
    term2_minus = 0.5*(spence(1 - z) - spence(1 - beta)) - 0.5 * (np.log(np.maximum(1 - z, epsilon)))**(2) - (11/12 + z/2) * np.log(np.maximum(1 - z, epsilon))
    term3_minus = z * (np.log(np.maximum(z, epsilon)) + 0.5*np.log(2*E_nu/m_e)) - (31/18 + (1/12)*np.log(np.maximum(z, epsilon)))*beta - (11/12)*z + (z**2)/24
    
    f_minus = term1_minus + term2_minus + term3_minus
    
    # (1-z)^2 * f+
    term1_plus = ((E_e/p_e)*np.log((E_e + p_e)/m_e) - 1)*(((1 - z)**2) * (2*np.log(np.maximum(1 - z - m_e/(E_e + p_e), epsilon)) - np.log(np.maximum(1 - z, epsilon)) - 0.5*np.log(np.maximum(z, epsilon)) - 2/3) - 0.5*((z**2)*np.log(np.maximum(z, epsilon)) + 1 - z))
    term2_plus = -(0.5*(1 - z)**2)*((np.log(np.maximum(1 - z, epsilon)))**2 + beta * (spence(z) - np.log(np.maximum(z, epsilon))*np.log(np.maximum(1 - z, epsilon))))
    term3_plus = np.log(np.maximum(1 - z, epsilon))*(0.5*(z**2) * np.log(np.maximum(z, epsilon)) + ((1 - z)/3) * (2*z - 0.5))
    term4_plus = -0.5*(z**2)*spence(z) - (z*(1 - 2*z)/3) * np.log(np.maximum(z, epsilon)) - z*(1-z)/6
    term5_plus = -(beta/12) * (np.log(np.maximum(z, epsilon)) + (1 - z)*((115 - 109*z)/6))
    
    f_plus_prod = term1_plus + term2_plus + term3_plus + term4_plus + term5_plus
    f_plus = f_plus_prod / ((1-z)**2)
    
    # f+/-
    f_pm = ((E_e/p_e)*np.log((E_e + p_e) / m_e) - 1) * 2*np.log(np.maximum(1 - z - m_e/(E_e + p_e), epsilon))
    
    # ------ Differential Cross Section Calculation in cm^2 / MeV ------

    term1 = (g_L**2) * (1 + (alpha/pi) * f_minus)
    term2 = ((g_R * (1 - z))**2) * (1 + (alpha/pi) * (f_plus))
    term3 = -g_R * g_L * ((m_e * T)/(E_nu**2)) * (1 + (alpha/pi) * (f_pm))
    
    result = (sigma_0 / m_e) * (term1 + term2 + term3)  
    
    return result


# ========== Numerical integration of dSigma_dT_corr ==========

def integrated_sigma(E_nu, flavor = 'nue'):

    '''
    Function that takes the dSigma_dT_corr function and integrate with respect to the
    electron kinetic energy. The integration is done numerically with the quad method 
    of scipy. The function computes the maximum kinetic energy of the recoil electron 
    for a given neutrino energy and takes the values as the integral limits.

    Parameters:
    - E_nu: Energy of the incident neutrino
    - flavor: neutrino flavor to compute the interactions cross-section.
                Default is electron neutrino flavor.

    return: Values of the cross-section for a neutrino energy
    
    '''

    m_e = 0.5109989461                      # electron mass energy in MeV/c^2 
    pi = np.pi
    
    Tmax = 2.0 * E_nu**2 / (m_e + 2.0 * E_nu)
    Tmin = 0
    #Tmin, _ = integrate.quad(Tmin_func, 0, 2*pi, args=(E_nu))
    
    #print(f"--- Integrating for E_nu = {E_nu:.2f} MeV ---")
    #print(f"Limits of T_e: [{Tmin:.4f}, {Tmax:.4f}] MeV")
    
    sigma, abs_error = integrate.quad(
        dSigma_dT_corr, 
        Tmin, 
        Tmax, 
        args=(E_nu, flavor)
    )
    
    #print(f"result of integration: {sigma:.6e} cm^2")
    
    return sigma