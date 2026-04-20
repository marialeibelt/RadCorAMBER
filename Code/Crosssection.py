import numpy as np

# =========================
# Konstanten
# =========================
alpha = 1/137.036  # Feinstrukturkonstante
M = 0.938272081    # Protonmasse [GeV]
m_l = 0.1056583745 # Leptonmasse (Myon) [GeV]

# =========================
# Sachs-Formfaktoren (Dipol-Ansatz) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Approximation
# =========================
def GE(Q2):
    # elektrischer Formfaktor
    return 1 / (1 + Q2 / 0.71)**2

def GM(Q2):
    # magnetischer Formfaktor (mit Protonmagnetmoment μp ≈ 2.793)
    mu_p = 2.792847
    return mu_p * GE(Q2)

# =========================
# Kinematik
# =========================
def Q2_lab(E, theta):
    """
    Q^2 im Laborsystem:
    E = Eingangsenergie des Leptons [GeV]
    theta = Streuwinkel [rad]
    """
    return 4 * E**2 * np.sin(theta/2)**2 / (1 + (2*E/M)*np.sin(theta/2)**2)

def epsilon(Q2, E, theta): #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Approximation
    """
    Photon-Polarisation ε
    """
    tau = Q2 / (4 * M**2)
    return 1 / (1 + 2*(1 + tau)*np.tan(theta/2)**2)

# =========================
# Mott-Cross-Section
# =========================
def sigma_mott(E, theta):
    """
    Mott-Wirkungsquerschnitt [GeV^-2 sr^-1]
    """
    return (alpha**2 * np.cos(theta/2)**2) / (4 * E**2 * np.sin(theta/2)**4)

# =========================
# Born Cross Section
# =========================
def born_cross_section(E, theta):
    """
    dσ/dΩ für elastische ep-Streuung im Born-Niveau
    E: Eingangsenergie [GeV]
    theta: Streuwinkel [rad]
    """
    Q2 = Q2_lab(E, theta)
    tau = Q2 / (4 * M**2)
    eps = epsilon(Q2, E, theta)

    sigma_R = GM(Q2)**2 + (eps / tau) * GE(Q2)**2

    return sigma_mott(E, theta) * (tau / (eps * (1 + tau))) * sigma_R


# =========================
# Beispiel
# =========================
if __name__ == "__main__":
    E = 1.0  # GeV
    theta = np.deg2rad(60)

    dsigma = born_cross_section(E, theta)

    print(f"dσ/dΩ = {dsigma:.5e} GeV^-2 sr^-1")