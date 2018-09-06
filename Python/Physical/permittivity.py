#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Collection of subfunctions related to the computation of permittivities for 
a stack of ice and snow layers. 
Based on the MEMLS model (copyright (c) 1997 by the Institute of Applied 
Physics, University of Bern, Switzerland), revised and extended to sea ice 
applications.  

@fork from: mwi, 2018

"""
import numpy as np
import brine 

def permittivity(rho,T,W,sal,sitype,freq):
    """   
    Compute the complex permittivity profile for snow and ice pack, with
    profiles of changing density, salinity, temperature, wetness with depth, 
    for a given frequency.
    
    eps = permittivity(freq,rho,T,W,sal,sitype)
        epsi:   relative permittivity (real part of dielectric constant)
        epsii:  dielectric loss factor (imaginary part of dielectric constant)
        rho:    density [g/cm^3]
        T:      temperature [K]
        W:      wetness of snow
        sal:    salinity profile
        sitype: type of snow/ice
        freq:   frequency [GHz]
        
    Uses: drysnow, wetsnow, salineicewithairbubbles
    """
    rho = np.array(rho)
    T = np.array(T)
    W = np.array(W)
    sal = np.array(sal)
    
    # Initialize:
    eps = np.complex(len(rho))

    # For snow on the sea ice: 
    # Assuming the same permittivity for either of the two snow types
    mask = (sitype == 1)|(sitype == 2)
    
    # Dry, clean!, snow permittivity:    
    eps[mask]= drysnow(rho[mask],T[mask],freq,'ImprovedBorn')    
    # Adjust epsi and epsii if wet snow:
    #epsi[mask], epsii[mask] = wetsnow(epsi[mask],epsii[mask],W[mask],rho[mask],freq,'Debye-like') # Probably better
    eps[mask] = wetsnow(eps[mask],W[mask],rho[mask],freq,'original')

    # For the ice: Permittivity depends on whether first-year (FY, type 3) 
    # or multi-year (MY, type 4) sea ice.
    # First-year sea ice: Sea ice with brine inclusions
    mask = (sitype == 3)
    # Setting rho = [] to indicate that there are no air inclusions.
    eps[mask] = salineicewithairbubbles(sal[mask],T[mask],[],freq)
        
    # Multi-year sea ice: Sea ice with brine and air inclusions
    mask = (sitype==4)
    eps[mask] = salineicewithairbubbles(sal[mask],T[mask],rho[mask],freq)
    return eps

#%% Dry snow dielectric constant:
def drysnow(rho,T,freq,method='Tinga-Voss-Blossey'):
    """
    Dielectric constant for clean, dry snow at various frequencies. 
    
    Real part only depends on density (and very weakly on temperature), while 
    imaginary part also depends on temperature and frequency.
    
    eps = drysnow(rho,T,freq,method)
        epsi:   relative permittivity (real part of dielectric constant)
        epsii:  dielectric loss factor (imaginary part of dielectric constant)
        rho:    density [g/cm^3]
        T:      snow temperature [K]
        freq:   frequency [GHz]
        method: Available models: 'Tinga-Voss-Blossey' (recommended), 
                'Tinga-Voss-Blossey-red', 'empirical', 'original' (NOT recommended). 
    
    Uses: epureice
    """
    # Ice volume factor:
    rho_pureice = 0.9167
    vi = rho/rho_pureice

    # Dielectric constant of pure ice: 
    #??? Assuming clean snow, which is not realistic!! (for epsii_ice) - salinity will change this, but not dust
    epsi_ice, epsii_ice = epureice(T,freq)
   
    if method == 'Tinga-Voss-Blossey':
        # Using the Tinga-Voss-Blossey model for spherical inclusions: 
        # Full dielectric constant for pure ice (real and imaginary parts)
        eps_ice = epsi_ice-1j*epsii_ice
    
        # Tinga-Voss_Blossey model (with dielectric constant of air: eps_air = 1-1j*0)
        eps = 1 + 3*vi*(eps_ice-1)/(2+eps_ice-vi*(eps_ice-1))
        epsi = eps.real # Relative permittivity
        epsii = -eps.imag # Loss factor
        
    if method == 'Tinga-Voss-Blossey-red':
        # Using the Tinga-Voss-Blossey mixture model for spherical inclusions, 
        # when removing small terms: 
        epsi = (1+0.84*vi)/(1-0.42*vi)
        #epsii = 9*vi*epsii_ice/(2+vi+epsi_ice*(1-vi))**2 #, or:
        epsii = 0.34*vi*epsii_ice/(1-0.42*vi)**2

    elif method == 'ImprovedBorn':
        # Using improved born approximation with spherical inclusions
        # (air is host material, with ice inclusions): 
        eps_air = 1
        eps = sphericalinclusions(eps_air,eps_ice,vi,'ImprovedBorn')
        
    elif method == 'ImprovedBorn_shells':
        # Using improved born approximation with shell inclusions
        # (air is host material, with ice inclusions): 
        eps_air = 1
        eps = sphericalinclusions(eps_air,eps_ice,vi,'ImprovedBorn_shells')
        
    elif method == 'empirical':
        # Empirical relationship, from Matzler (2006)  (Ulaby & Long, eq. 4.55, p.141)
        # Relative permittivity:
        epsi = np.zeros(len(rho))
        mask = (rho<=0.45) #??? Check grænse-værdi. Andet steds: 0.4, men 0.45 er mere kontinuert
        epsi[mask] =  1 + 1.4667*vi[mask] + 1.435*vi[mask]**3
        epsi[~mask] = (1 + 0.4759*vi[~mask])**3
    
        # Dielectric loss factor (epsii): 
        # After Tiuri (1984):
        epsii = epsii_ice*(0.52*rho + 0.62*(rho**2))     
        
    elif method == 'original':
        # This approach below is NOT recommended as it gives rise to 
        # discontinuities in epsii. 
        
        # Relative permittivity: 
        epsi = np.zeros(len(rho))
        mask = rho<=0.4
        epsi[mask] = 1 + 1.5995*rho[mask] + 1.861*rho[mask]**3
        epsi[~mask] = ((1-vi[~mask])*0.99913 + vi[~mask]*1.4759)**3
    
        # Loss factor:
        # After Polder and van Santen 1946 (Effective-Medium Approx):
        A = np.ones((len(rho)))*0.3 # Value for f>=0.55
        mask = (vi<0.55) & (vi>0.333)
        A[mask] = 0.476 - 0.64*vi[mask]
        mask = (vi<=0.333)
        A[mask] = 0.1 + 0.5*vi[mask]
    
        ei = 3.185  
        A3 = 1.-2.*A
        ea = (epsi*(1-A))+A
        ea3 = epsi*(1-A3)+A3
        K1 = (ea/(ea+A*(ei-1)))**2
        K3 = (ea3/(ea3+A3*(ei-1)))**2
        Ksq = (2*K1+K3)/3.
        # Imaginary part of dielectric permittivity of dry snow:
        epsii = np.sqrt(epsi)*epsii_ice*Ksq*vi
    
    eps=epsi-1j*epsii
    return eps

#%% Wet snow dielectric constant:
def wetsnow(Wetness,rho,freq,method='Debye-like'):
    """
    The dielectric constant for wet snow (W>0). 
    
    The Debye-like model is based on empirical constants appropriate for the 
    intervals: 3<f<37 GHz, 0.09<rho<0.38 g/cm³, 0.01<W<0.12
        
    The original model is a physical mixing model (Weise 97) after Matzler 1987 (corrected).
    Water temperature is assumed constant at 273.15 K. 

        eps  =  wetsnow(epsi,epsii,W,freq,method)
        eps:    dielectric constant
        W:      liquid water content of snow or "snow wetness"; the volume 
                fraction of liquid water in the snow mixture (values: 0-1)
        freq:   frequency [GHz]
        method: Available models: 'Debye-like' (recommended), 'original'
    """
    W=np.clip(Wetness,0.01,0.12)
    if method == 'Debye-like':
        # Using amodified Debye-like model with empirical constants:
        # Constants are found from fit to data in ranges: 
        # 3<f<37 GHz, 0.09<rho<0.38 g/cm³, 0.01<W<0.12
        # Hallikainen (1986), from Ulaby and Long, E4.60,4.61, p. 143
    
        A1 = 0.78+0.03*freq-0.58e-3*freq**2
        A2 = 0.97-0.39e-2*freq+0.39e-3*freq**2
        B1 = 0.31-0.05*freq+0.87e-3*freq**2
        A = A1*(1.0+1.83*rho+0.02*(W*100)**1.015)+B1
        B = 0.073*A1
        C = 0.073*A2
        x = 1.31
        f0 = 9.07 # Effective relaxation frequency of wet snow
        
        epsi = A + B*(W*100)**x/(1.0+(freq/f0)**2)
        epsii = C*(freq/f0)*(W*100.0)**x/(1.0+(freq/f0)**2)

    elif method == 'original':
        Aa = 0.005
        Ab = 0.4975
        Ac = 0.4975
        euw = 4.9
        esw = 88.045 
        frw = 0.11109 # inverse relaxation frequency of water
  
        esa = (esw - epsi)/(3*(1+Aa*(esw/epsi-1)))
        esb = (esw - epsi)/(3*(1+Ab*(esw/epsi-1)))
        esc = (esw - epsi)/(3*(1+Ac*(esw/epsi-1)))
        eua = (euw - epsi)/(3*(1+Aa*(euw/epsi-1)))
        eub = (euw - epsi)/(3*(1+Ab*(euw/epsi-1)))
        euc = (euw - epsi)/(3*(1+Ac*(euw/epsi-1)))
  
        fa = 1. + Aa * (esw-euw)/(epsi+Aa*(euw-epsi))
        fb = 1. + Ab * (esw-euw)/(epsi+Ab*(euw-epsi))
        fc = 1. + Ac * (esw-euw)/(epsi+Ac*(euw-epsi))
  
        eea = esa - eua
        eeb = esb - eub
        eec = esc - euc
  
        fwa = frw/fa
        fwb = frw/fb
        fwc = frw/fc
        
        # Adjustment to relative permittivity (real part):
        depsia = eua + eea / (1+(fwa*freq)**2)
        depsib = eub + eeb / (1+(fwb*freq)**2)
        depsic = euc + eec / (1+(fwc*freq)**2)
        depsi = W * (depsia+depsib+depsic)

        # Adjustment to loss factor (imaginary part):
        depsiia = fwa*freq*eea / (1+(fwa*freq)**2)
        depsiib = fwb*freq*eeb / (1+(fwb*freq)**2)
        depsiic = fwc*freq*eec / (1+(fwc*freq)**2)
        depsii = W * (depsiia+depsiib+depsiic)

        # Adjusted values:
        epsi = epsi + depsi
        epsii = epsii + depsii
        
    eps=epsi-1j* epsii
    return eps

def epureice(T,freq,method='Hufford_corrected'):
    """
    Dielectric permittivity of pure water ice.
    The relative permittivity is only very weakly temperature dependent. 
    The dielectric loss factor is temperature and frequency dependent. Several 
    options are available for calculating the dielectric loss factor. 
    
    Note: There is some confusion where the various models are coming from.
    
    eps_ice = epureice(T,freq,method)
        epsice:  real part of dielectric permittivity
        epsiice: imaginary part of dielectric permittivity
        T:       temperature [K]
        freq:    frequency [GHz]
        method:  'Hufford_corrected' (default), 'Hufford1991', 
                 'Hufford_Mishima' (not recommended), 'Hufford_Matzler', 'Matzler1998' (not recommended)
        
        Unsure about: 
            A correction proposed by Pulliainen (1997) (??)
            Slightly adjusted parameters (from original files, quoting Hufford, Mitzima, Matzler)??

    """
    
    # Temperature in oC:
    Tc = T-273.15

    # The relative permittivity of pure water ice (real part):
    # Derived by Matzler and Wegmuller (1987); Valid for frequency range: 
    # 10MHz to 300 GHz, and temperatures: -40<T<0oC
    epsi_ice=3.1884+9.1e-4*Tc 
    
    # Dielectric loss of pure ice:
    if method == 'Hufford_corrected':
        # From Hufford (1991), and explained in Matzler (2006), p. 456-460. 
        # This is the form used in Ulaby and Long, 2014. 
        theta = (300.0/T)-1
        B1 = 0.0207
        b = 335.25 # or 335.0?
        B2 = 1.16e-11
        alpha = (0.00504 + 0.0062*theta)*np.exp(-22.1*theta)        
        beta = B1/T * np.exp(b/T)/((np.exp(b/T)-1)**2) + B2*freq**2 + np.exp(-9.963+0.0372*Tc)

    elif method == 'Hufford1991':
        # The Hufford 1991 model (according to the Matzlers 2004 note on MEMLS (eq 3-4, p. 25))
        theta = (300.0/T)-1
        alpha = (0.00504 + 0.0062*theta)*np.exp(-22.1*theta)
        beta = (0.502-0.131*theta)/(1+theta)*10**-4 + 0.542e-6*((1+theta)/(theta+0.0073))**2 
    
    elif method == 'Hufford_Mishima':
        # The Hufford 1991 model, updated by Mishima (1983):  
        # However, Mishima may include a third term, depending on freq**3
        theta = (300.0/T)-1
        B1 = 0.0207
        b = 335
        B2 = 1.16e-11
        alpha = (0.00504 + 0.0062*theta)*np.exp(-22.1*theta)
        beta = B1/T * np.exp(b/T) / ((np.exp(b/T)-1)**2) + B2*freq**2 
        
    elif method == 'Hufford_Matzler':
        # The Hufford 1991 model, updated by Mishima (1983) and with additional 
        # correction by Matzler (MEMLS notes, note 5) for temperatures T>200K:
        theta = (300.0/T)-1
        B1 = 0.0207
        b = 335
        B2 = 1.16e-11
        alpha = (0.00504 + 0.0062*theta)*np.exp(-22.1*theta)
        beta = B1/T * np.exp(b/T) / ((np.exp(b/T)-1)**2) + B2*freq**2 + np.exp(-10.02+0.0364*(T-273))
        
    elif method == 'Matzler1998':
        # From: Mätzler, 1998, Microwave properties of ice and snow, in:
        # B. Smitht et al. (eds.) solar system ices, 241-257, Kluwer.
        # Gives rather different values! Not recommended.    
        theta=300.0/T
        alpha=(0.00504+0.0062*theta)*np.exp(-22.1*theta)
        beta=(0.502-0.131*theta/(1+theta))*1e-4 +(0.542e-6*((1+theta)/(theta+0.0073))**2)   

    epsii_ice = alpha/freq + beta*freq
    eps_ice=epsi_ice-1j*epsii_ice
    return eps_ice

def ebrine(T,freq,method='StogrynDesargant'):
    """
    Permittivity of brine.
    Equations from Ulaby et al, 1986.
    Coefficients are derived for a fit between -2 - -43.2oC.
    Alternative equations from Stogryn and Desargant, 1985, Eq.1, are also 
    provided; coefficients are here derived for a fit down to -25°C. 
    The main difference between the two lies in the calculation of brine 
    conductivity. 
    For small salinity values (and temperatures close to 0oC), the saline-water 
    Double-Debye model (D3M) may be used for the permittivity calculations. 
    
    eps_brine = ebrine(T, freq, method)
        epsbrine:     dielectric constant
        freq:         frequency [GHz]
        T:            temperature [K]
        method:       'Stogryn','StogrynDesargant' (default)
                
    Uses: brine.conductivity, brine.salinity, brine.normality
    """
    
    Tc = T-273.15 # Temperature in Celcius
    freqHz = freq*1e9 # Frequency in Hz
    e0 = 8.85419e-12 # Vaccuum electric permittivity

    # Salinity of brine:
    # Version from Ulaby et al, 1986 (E63), originally by Assur. 
    # Valid range: -43.2oC--1.8oC
    sal = brine.salinity(T,'Assur')    # unit: psu = o/oo
    
    if method == 'Stogryn':
        # Ionic conductivity of brine:
        condbrine, N = brine.conductivity(T,sal,'Stogryn1971')
    
        # Relaxation time:
        # From Ulaby et al. 1986 vol. III E67, E17
        rel=(1.1109e-10-3.824e-12*Tc+6.938e-14*Tc**2-5.096e-16*Tc**3)/(2*np.pi)
        b=1.0+0.146e-2*Tc*N-4.89e-2*N-2.97E-2*N**2+5.64e-3*(N**3)
        trelax=rel*b
    
        # Static dielectric constant of brine:
        # From Ulaby et al. 1986 vol. III E66
        eps=88.045-0.4147*Tc+6.295e-4*Tc**2+1.075e-5*Tc**3
        a=1.0-0.255*N+5.15e-2*N**2-6.89e-3*N**3
        eb0 = eps*a
    
        # The Debye relaxation law:
        # High-frequency limit of the dielectric constant:
        epsiwoo=4.9 
        # Relative permittivity:
        epsi_brine = epsiwoo+(eb0-epsiwoo)/(1.+(2*np.pi*trelax*freqHz)**2)
        # Loss factor:
        epsii_brine = 2*np.pi*trelax*freqHz*(eb0-epsiwoo)/(1.+(2*np.pi*trelax*freqHz)**2) + condbrine/(2*np.pi*freqHz*e0)
        
    elif method == 'StogrynDesargant':
        # Alternative calculation: From Stogryn and Desargant, 1985, Eq.7.
        # Derived from fit to data in temperature range: -2.8--25oC
        
        # Brine conductivity: Using here a different expression:
        condbrine = brine.conductivity(T,0,'StogrynDesargant1985')[0]
        # If using the conductivity equations from Stogryn 1971 instead:
        #sal = brine.salinity(T,'Assur')  
        #condbrine, N = brine.conductivity(T,sal,'Stogryn1971')
        
        # High-frequency limit of the dielectric constant:
        epsiwoo = (82.79 + 8.19*Tc**2)/(15.68 + Tc**2) 
        # Note: Values are somewhat different from 4.9 (the value used above)
        
        # Static dielectric constant of brine (Eq.10):
        epsib0 = (939.66 - 19.068*Tc)/(10.737 - Tc)
        
        # Relaxation time (in nano-seconds) (Eq.12): 
        trelax = (1./(2*np.pi)) * (0.10990 + 0.13603e-2*Tc + 0.20894e-3*Tc**2 + 0.28167e-5*Tc**3) 
        trelax = trelax/1e9 # in seconds
        
        # The Debye relaxation law:
        ebrine = epsiwoo + (epsib0 - epsiwoo)/(1.+2*np.pi*freqHz*trelax*1j) - (1j*condbrine)/(2*np.pi*e0*freqHz)
        # Splitting up in real and imaginary parts:
        epsi_brine = ebrine.real
        epsii_brine = -ebrine.imag
    
    # For temperatures close to OoC: 
    # Using the saline-water Double-Debye Dielectric Model, D3M.
    #mask = Tc>-1.8
    #epsi_brine[mask],epsii_brine[mask]=salinewater(T[mask],sal[mask],freq)
    # However, not calculated correctly. 

    # Always positive:
    epsi_brine = np.clip(epsi_brine,0,None)
    epsii_brine = np.clip(epsii_brine,0,None)
  
    return epsi_brine-1j*epsii_brine

def salinewater(T,sal,freq):
    """
    The Double-Debye Dielectric Model, D3M, developed for saline water.
    Applicability range: 0<T<30oC, 0<sal<40 psu, 0<freq<1000GHz.
    
    D3M was developed by Ellison, and reported in Matzler (2006, p. 431-455) 
    (Ulaby&Long, p. 125). 
        
    eps = salinewater(T,sal,freq)
        eps_brine:    dielectric constant
        T:            temperature [K]
        sal:          salinity [psu]
        freq:         frequency [GHz]
    """
    
    # Temperature in celcius:
    Tc = T-273.15

    # Parameters for the D3M model:
    a1 = 0.46606917e-2
    a2 = -0.26087876e-4
    a3 = -0.63926782e-5
    a4 = 0.63000075e1
    a5 = 0.26242021e-2
    a6 = -0.42984155e-2
    a7 = 0.34414691e-4
    a8 = 0.17667420e-3
    a9 = -0.20491560e-6
    a10 = 0.58366888e3
    a11 = 0.12684992e3
    a12 = 0.69227972e-4
    a13 = 0.38957681e-6
    a14 = 0.30742330e3
    a15 = 0.12634992e3
    a16 = 0.37245044e1
    a17 = 0.92609781e-2
    a18 = -0.26093754e-1
        
    ns = 1e-9 # Error in equation, p. 126?
    epsw0 = 87.85306 * np.exp(-0.00456992*Tc-a1*sal-a2*sal**2-a3*sal*Tc)
    epsw1 = a4 * np.exp(-a5*Tc-a6*sal-a7*sal*Tc)
    tau_w1 = (a8+a9*sal) * np.exp(a10/(Tc+a11)) * ns
    tau_w2 = (a12+a13*sal) * np.exp(a14/(Tc+a15)) * ns
    epswoo = a16 + a17*Tc + a18*sal
    e0 = 8.85419e-12 # Vaccuum electric permittivity
    
    # Conductivity: 
    sig = 2.903602+8.607e-2*Tc + 4.738817e-4*Tc**2 - 2.991e-6*Tc**3 + 4.3041e-9*Tc**4
    p = sal * (37.5109+5.45216*sal+0.014409*sal**2) / (1004.75+182.283*sal+sal**2)
    alpha0 = (6.9431+3.2841*sal-0.099486*sal**2) / (84.85+69.024*sal+sal**2)
    alpha1 = 49.843-0.2276*sal+0.00198*sal**2
    q = 1+alpha0*(Tc-15)/(Tc+alpha1)
    sigma = p*q*sig
    
    # Double Debye-model:
    epsi = epswoo + (epsw0-epsw1)/(1+(2*np.pi*freq*tau_w1)**2) + \
        (epsw1-epswoo)/(1+(2*np.pi*freq*tau_w2)**2)
        
    epsii = 2*np.pi*freq*tau_w1*(epsw0-epsw1)/(1+(2*np.pi*freq*tau_w1)**2) + \
        2*np.pi*freq*tau_w2*(epsw1-epswoo)/(1+(2*np.pi*freq*tau_w2)**2) + \
        sigma/(2*np.pi*e0*freq)
        
    return epsi-1j*epsii

def randomneedles(ehost, einc, volinc,method='PolderVanSanten'):
    """
    The effective dielectric constant of a mixture of clean sea ice 
    (host material) with inclusions of random brine needles. 
    The dielectric constant is calculated according to the Polder-Van Santen 
    equations (e.g. Shokr (1998) Arctic sea ice in the microwave c-band, IEEE TGRS 36(2), 463-78)
    Similar equation is also given in Ulaby and Long (2014), E4.36, p. 313.  

    [eps_mix] = randomneedles(epsi_host, epsii_host, epsi_inc, epsii_inc, volinc)
        eps_mix:    mixture dielectric constant
        ehost:      dielectric constant of host material, i.e. ice
        einc:       dielectric constant of brine inclusions
        volinc:     volume of brine inclusions in mixture
        method:     'PolderVanSanten' (default), 'TingaVossBlossey'
    """
    
    if method == 'PolderVanSanten':
        # Initialize:
        for i in range(len(ehost)):
            # Convergence measures:
            measurea=1.2 
            measureb=1.2 
            # Iterate to solve for emix:
            m=0
            while measurea>0.001 and measureb>0.001:
                m = m+1
                if m>30: break
        
                # Initial estimate of mixture permittivity: 
                # Following Vant et al. (1978), J. Appl. Phys. 49(3), 1264-80.
                emix = 3.05+0.72*volinc[i]-(0.024+0.33*volinc[i])*1j
    
                f1=(einc[i]-ehost[i])/(einc[i]+emix)
                f2=5.0*emix+einc[i]
                emix_new = ehost[i]+volinc[i]/3.*f1*f2
        
                # Check for convergence:
                measurea=abs(emix.real-emix_new.real)
                measureb=abs(emix.imag-emix_new.imag)
                emix = emix_new
    elif method=='TingaVossBlossey':
        emix = ehost + volinc/3*(einc-ehost)*(ehost*(5+volinc)+(1-volinc)*einc)/(ehost*(1+volinc)+einc*(1-volinc))          
    return emix

def sphericalinclusions(ehost,einc,volinc,method='ImprovedBorn'):
    """
    It can e.g. be computed using the Polder-Van Santen/de Loor mixing formulae for 
    spherical inclusions, or the improved born approximation as given by 
    C. Mätzler (1998). J. Appl. Phys. 83(11),6111-7. 
    
    There is significant difference, and likely the Polder Van Santen 
    implementation is incorrect. 
    
    [epsi_mix,epsii_mix] = spherical inclusions(epsi_host, epsii_host, epsi_inc, epsii_inc, volinc)
        epsi_mix:   relative permittivity of mixture (real part of dielectric constant) 
        epsii_mix:  dielectric loss factor of mixture (imaginary part of dielectric constant)
        epsi_host:  relative permittivity of host material, i.e. ice
        epsii_host: dielectric loss factor of host material, i.e. ice
        epsi_inc:   relative permittivity of brine inclusions
        epsii_inc:  dielectric loss factor of brine inclusions
        volinc:     volume of brine inclusions in mixture
        method:     'PolderVanSanten' (not recommended), 'ImprovedBorn' (default), 'TingaVossBlossey'
    """

    if method == 'PolderVanSanten':
        for i in range(len(ehost)):
            # Convergence measures:
            measurea=1.2 
            measureb=1.2 
            # Iterate to solve for emix:
            m=0
            while measurea>0.001 and measureb>0.001:
                m = m+1
                if m>30: break
            
                # Initial estimate of mixture permittivity: 
                # Following Vant et al. (1978), J. Appl. Phys. 49(3), 1264-80.
                emix = 3.05+0.72*volinc[i]-(0.024+0.33*volinc[i])*1j
                # New estimate:
                emix_new = ehost[i]+3*volinc[i]*emix*(einc[i]-ehost[i])/(einc[i]+2*emix)
        
                # Check for convergence:
                measurea=abs(emix.real-emix_new.real)
                measureb=abs(emix.imag-emix_new.imag)
                emix = emix_new       
                     
    elif method == 'TingaVossBlossey':
        emix = ehost + 3*volinc*ehost*(einc-ehost)/(2*ehost+einc-volinc*(einc-ehost))
        
    elif method == 'ImprovedBorn':
        # Effective permittivity of Polder and Van Santen (according to Matzler's notes)
        emix = 0.25 * (2*ehost-einc+3*volinc*(einc-ehost) + 
                       np.sqrt((2*ehost-einc+3*volinc*(einc-ehost))**2 + 8*ehost*einc))

        
    elif method == 'ImprovedBorn_shells':
        # according to Matzler, 1998
        emix = ehost + volinc*(einc-ehost)*(2+ehost/einc)/(3-volinc*(1-ehost/einc))

        
    return emix
   
def salineicewithairbubbles(sal,T,rho,freq):
    """
    The dielectric constant of saline sea ice with air bubbles.
    If no air, use rho = []
    Background information: Ulaby et al. 1986 vol. III.

    eps_mix = salineicewithairbubbles(sal,T,freq)
        eps_mix:    dielectric constant
        T:          temperature [K]
        rho:        density [g/cm³]
        freq:       frequency [GHz]
        
    Uses: epureice, ebrine, brine.volume, sphericalinclusions, (randomneedles)    
    """
    
    #%% Dielectrical constant of saline sea ice:
    # Dielectric constant of pure ice:
    eps_ice = epureice(T,freq,'Hufford_Matzler')
    # Dielectric constant of brine:
    eps_brine = ebrine(T,freq,'StogrynDesargant')   
    # Volume of brine in ice:
    volbrine = brine.volume(T,sal,'original') 
    #volbrine = brine.volume(T,sal,'Frankenstein') 
    
    # Assuming spherical brine inclusions:
    eps_spheres = sphericalinclusions(eps_ice,eps_brine,volbrine,'ImprovedBorn')

    #%% Dielectric constant of saline sea ice with air bubbles:
    if len(rho) > 0:
        # Air volume: 
        rho_ice = 0.926 #??? why this number??
        vol_air=(rho_ice-rho)/rho_ice
        vol_air[vol_air<0]=0
    
        # Permittivity of air: 
        eps_air = 1.0
        
        # Permittivity of saline sea ice with air bubbles as spherical inclusions:
        eps_mix = sphericalinclusions(eps_spheres,eps_air,vol_air)
        
    return eps_mix
