import numpy as np

# =========================
# Konstanten
# =========================
alpha = 1/137.035999084
M = 0.938272      # Protonmass
m_mu = 0.105658   # Muonmass

# =========================
# Minkowski metric (+,-,-,-)
# =========================
def dot(p, q):
    return p[0]*q[0] - np.dot(p[1:], q[1:])

# =========================
# Mandelstam variables
# =========================
def mandelstam(p1, p2, k1, k2):
    s = dot(p1 + k1, p1 + k1)
    t = dot(k1 - k2, k1 - k2)   # = q^2 (negative)
    u = dot(p1 - k2, p1 - k2)
    return s, t, u

# =========================
# Sachs Formfaktoren (dipole placeholder)
# =========================
def GE(Q2):
    return 1.0 / (1.0 + Q2 / 0.71)**2

def GM(Q2):
    mu_p = 2.793
    return mu_p * GE(Q2)

# =========================
# ε fully invariant form
# (from ν = (s-u)/4)
# =========================
def epsilon_invariant(s, u, t):
    Q2 = -t
    tau = Q2 / (4 * M**2)
    nu = (s - u) / 4.0

    eps_num = nu**2 - (M**4 * tau * (1 + tau))
    eps_den = nu**2 + (M**4 * tau * (1 + tau)) * (1 - 2 * m_mu**2 / Q2)

    return eps_num / eps_den

# =========================
# Born cross section (fully invariant)
# =========================
def born_ep(p1, k1, p2, k2):
    """
    p1: incoming proton
    k1: incoming lepton
    p2: outgoing proton
    k2: outgoing lepton
    """

    s, t, u = mandelstam(p1, p2, k1, k2)
    Q2 = -t

    tau = Q2 / (4 * M**2)

    # ν invariant
    nu = (s - u) / 4.0

    # ε (fully invariant)
    eps = epsilon_invariant(s, u, t)

    # reduced Born structure function
    sigma_R = GM(Q2)**2 + (eps / tau) * GE(Q2)**2

    # full Born prefactor (covariant form)
    pref = alpha**2 / (Q2**2)

    # flux factor structure (McMule-like normalization)
    flux = 1.0 / (16.0 * np.pi * (s - M**2 - m_mu**2)**2)

    return flux * pref * (nu**2 + M**4 * tau * (1 + tau)) * sigma_R