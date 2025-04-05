import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

# User-defined parameters
kT = 0.05  # Temperature in eV (example value)
F_As = 0.1  # As flux (user-defined)

# Fixed parameters from the document
a0 = 1.0            # cell
Omega_Ga = 1      # cell^3
C0 = 1.0            # cell^-2
kr = 0.1            # cell^2 / time
D_Ga = 1.0          # diffusion coefficient Ga
E_As_minus_E_Ga = 0.0  # eV (from table)
E_a_minus_E_Ga = 0.1   # eV (from table)
D_As = np.exp(-E_As_minus_E_Ga / kT)
w = 2.0             # width of Gaussian flux (assumed)
R_inf = 400.0       # domain size (cells)
theta = 60.0        # contact angle in degrees
rho0 = 100/R_inf         # initial rho = Rd/R_inf
h0 = 1.0            # cell height per layer

# Derived parameters
theta_rad = np.radians(theta)
cos_t = np.cos(theta_rad)
sin_t = np.sin(theta_rad)
cos_3t = np.cos(3 * theta_rad)
sin_3t = np.sin(3 * theta_rad)
numerator_B = 8 - 9 * cos_t + cos_3t
denominator_B = 3 * sin_t - sin_3t
B_theta = numerator_B / denominator_B

# Compute tau_As based on kT
tau_As = np.exp(E_a_minus_E_Ga / kT)

# Compute sigma
sigma = (2 * Omega_Ga * D_Ga * C0) / (B_theta * R_inf**3)

# Determine total simulation time T based on droplet depletion
rho_end = 0.01  # target final rho (user can adjust)
term_initial = rho0**3 * (3 * np.log(rho0) - 1)
term_final = rho_end**3 * (3 * np.log(rho_end) - 1)
T = (term_final - term_initial) / (9 * sigma)

# Grid setup
Nr = 400
dr = R_inf / Nr
r = np.linspace(0, R_inf, Nr + 1)

dt = dr**2/10  # Adjust dt based on stability
T_total = T - dt  # Ensure termination before full depletion
Nt = int(T_total / dt)

# Initialize arrays
c_ga = np.zeros(Nr + 1)
c_as = np.zeros(Nr + 1)
h = np.zeros(Nr + 1)

# Initial droplet radius
Rd = rho0 * R_inf
rho = rho0
time = 0.0

# For storing results
Rd_history = []
time_history = []

def update_rho(rho_old, time_current, sigma):
    RHS = rho0**3 * (3 * np.log(rho0) - 1) + 9 * sigma * time_current
    for _ in range(100):
        p = rho_old**3 * (3 * np.log(rho_old) - 1) - RHS
        if abs(p) < 1e-6:
            break
        p_prime = 3 * rho_old**2 * (3 * np.log(rho_old) - 2)
        if p_prime == 0:
            break
        rho_new = rho_old - p / p_prime
        rho_new = np.clip(rho_new, 1e-6, 0.99)
        rho_old = rho_new
    return rho_old

# Finite difference parameters
alpha = D_Ga * dt / dr**2
beta = D_As * dt / dr**2 
omega = C0 * kr * dt
kappa = dt * F_As / C0
gamma = dt / tau_As
epsilon = (F_As * tau_As) / C0
eta = dr**2/w**2
upsilon = 2*alpha*eta

def update(c_ga, c_as, flux_term):
    c_ga_next = np.zeros_like(c_ga)
    c_as_next = np.zeros_like(c_as)
    
    # Handle j=0
    j = 0
    c_ga_next[j] = c_ga[j] + 4 * alpha * (c_ga[j+1] - c_ga[j]) + upsilon * flux_term[j] - omega * c_ga[j] * c_as[j]
    c_as_next[j] = c_as[j] + 4 * beta * (c_as[j+1] - c_as[j]) - gamma*c_as[j] + kappa - omega * c_ga[j] * c_as[j]
    
    # Handle j from 1 to Nr-1
    for j in range(1, Nr):
        term_ga = alpha * ((1 + 1/(2*j)) * c_ga[j+1] - 2 * c_ga[j] + (1 - 1/(2*j)) * c_ga[j-1])
        c_ga_next[j] = c_ga[j] + term_ga + upsilon * flux_term[j] - omega * c_ga[j] * c_as[j]
        
        term_as = beta * ((1 + 1/(2*j)) * c_as[j+1] - 2 * c_as[j] + (1 - 1/(2*j)) * c_as[j-1])
        c_as_next[j] = c_as[j] + term_as - gamma*c_as[j] + kappa - omega * c_ga[j] * c_as[j]
    
    # Apply boundary condition at j=Nr
    c_ga_next[Nr] = 0.0
    c_as_next[Nr] = epsilon
    
    return c_ga_next, c_as_next

# Time loop
for step in range(Nt):
    # Compute F_Ga parameters
    x = Rd / w
    p_val = R_inf / w
    a_p = 3.545
    b_p = 0.187 / (p_val - 3.156)
    denominator = (a_p * x + b_p) * np.log(R_inf / Rd)
    flux_term = np.exp(-(r/w - x)**2) / denominator

    # Update concentrations (using Euler scheme)
    c_ga, c_as = update(c_ga, c_as, flux_term)

    # Update droplet radius
    time += dt
    rho = update_rho(rho, time, sigma)
    Rd = rho * R_inf
    Rd_history.append(Rd)
    time_history.append(time)

    # Update height profile
    dh = omega * h0 * c_ga * c_as
    h += dh

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(time_history, Rd_history, 'b-')
plt.xlabel('Time')
plt.ylabel('Droplet Radius (cells)')
plt.title('Droplet Radius Evolution')
plt.grid(True)

plt.figure(figsize=(10, 6))
plt.plot(r, h, 'r-')
plt.xlabel('Radial Position (cells)')
plt.ylabel('Height (cells)')
plt.title('Quantum Ring Height Profile')
plt.grid(True)

plt.show()