mod file;

fn main() {
    // User-defined parameters
    let kt = 0.05;          // Temperature in eV
    let f_as = 0.1;         // As flux

    // Fixed parameters
    let omega_ga = 1.0;     // cell^3
    let c0 = 1.0;           // cell^-2
    let kr = 0.1;           // cell^2 / time
    let d_ga = 1.0;         // diffusion coefficient Ga
    let e_as_minus_e_ga: f64 = 0.0; // eV
    let e_a_minus_e_ga: f64 = 0.05;  // eV
    let w: f64 = 2.0;            // width of Gaussian flux
    let r_inf = 400.0;      // domain size (cells)
    let theta = 60.0;       // contact angle in degrees
    let rho0_initial: f64 = 100.0 / r_inf; // initial rho
    let h0 = 1.0;           // cell height per layer

    // Derived parameters
    let theta_rad: f64 = theta / 180.0 * 3.14;
    let numerator_b = 8.0 - 9.0 * theta_rad.cos() + (3.0 * theta_rad).cos();
    let denominator_b = 3.0 * theta_rad.sin() - (3.0 * theta_rad).sin();
    let b_theta = numerator_b / denominator_b;
    let tau_as = (e_a_minus_e_ga / kt).exp();
    let sigma = (2.0 * omega_ga * d_ga * c0) / (b_theta * r_inf.powi(3));

    // Simulation time calculation
    let rho_end: f64 = 0.01;
    let term_initial = rho0_initial.powi(3) * (3.0 * rho0_initial.ln() - 1.0);
    let term_final = rho_end.powi(3) * (3.0 * rho_end.ln() - 1.0);
    let t_total = ((term_final - term_initial) / (9.0 * sigma)) - 1e-6;

    // Grid setup
    let nr = 400;
    let dr = r_inf / nr as f64;
    let r: Vec<f64> = (0..=nr).map(|i| i as f64 * dr).collect();
    let dt = dr.powi(2) / 10.0;
    let nt = (t_total / dt) as usize;

    // Number of snapshots
    let ns: usize = 10;

    // Initialize arrays
    let mut c_ga = vec![0.0; nr + 1];
    let mut c_as = vec![0.0; nr + 1];
    let mut h = vec![0.0; nr + 1];
    let mut h_history = vec![vec![0.0; nr + 1]; ns+1];
    let mut ga_history = vec![vec![0.0; nr + 1]; ns+1];
    let mut as_history = vec![vec![0.0; nr + 1]; ns+1];
    let mut rho = rho0_initial;
    let mut time = 0.0;

    // Precompute coefficients
    let alpha = d_ga * dt / dr.powi(2);
    let beta = (-e_as_minus_e_ga / kt).exp() * dt / dr.powi(2);
    let omega = c0 * kr * dt;
    let kappa = dt * f_as / c0;
    let gamma = dt / tau_as;
    let epsilon = (f_as * tau_as) / c0;
    let upsilon = 2.0 * alpha * dr.powi(2) / w.powi(2);

    // Simulation loop
    let mut rd_history = Vec::new();
    let mut time_history = Vec::new();

    let mut js = 0;

    for jt in 0..nt {
        // Compute flux term
        let rd = rho * r_inf;
        let x = rd / w;
        let p_val = r_inf / w;
        let a_p = 3.545;
        let b_p = 0.187 / (p_val - 3.156);
        let denominator = (a_p * x + b_p) * (r_inf / rd).ln();
        let flux_term: Vec<f64> = r.iter()
            .map(|ri| (-(ri / w - x).powi(2)).exp() / denominator)
            .collect();

        // Update concentrations
        let (new_c_ga, new_c_as) = update_concentrations(
            &c_ga,
            &c_as,
            &flux_term,
            alpha,
            beta,
            omega,
            gamma,
            kappa,
            epsilon,
            upsilon,
            nr,
        );
        c_ga = new_c_ga;
        c_as = new_c_as;

        // Update droplet radius
        time += dt;
        rho = update_rho(rho, time, sigma, rho0_initial);
        rd_history.push(rho * r_inf);
        time_history.push(time);

        // Update height profile
        for j in 0..=nr {
            h[j] += omega * h0 * c_ga[j] * c_as[j];
        }

        if jt % (nt / ns) == 0 {
            h_history[js] = h.clone();
            ga_history[js] = c_ga.clone();
            as_history[js] = c_as.clone();
            js += 1;
        }
    }

    h_history.push(r.clone());
    ga_history.push(r.clone());
    as_history.push(r);

    file::save_columns_to_file(&h_history, "results", "height.dat");
    file::save_columns_to_file(&ga_history, "results", "c_ga.dat");
    file::save_columns_to_file(&as_history, "results", "c_as.dat");

}

fn update_concentrations(
    c_ga: &[f64],
    c_as: &[f64],
    flux_term: &[f64],
    alpha: f64,
    beta: f64,
    omega: f64,
    gamma: f64,
    kappa: f64,
    epsilon: f64,
    upsilon: f64,
    nr: usize,
) -> (Vec<f64>, Vec<f64>) {
    let mut c_ga_next = vec![0.0; c_ga.len()];
    let mut c_as_next = vec![0.0; c_as.len()];

    // Handle j=0
    c_ga_next[0] = c_ga[0] + 4.0 * alpha * (c_ga[1] - c_ga[0]) 
        + upsilon * flux_term[0] 
        - omega * c_ga[0] * c_as[0];
    
    c_as_next[0] = c_as[0] + 4.0 * beta * (c_as[1] - c_as[0]) 
        - gamma * c_as[0] 
        + kappa 
        - omega * c_ga[0] * c_as[0];

    // Handle interior points
    for j in 1..nr {
        let j_f64 = j as f64;
        let coeff = 1.0 / (2.0 * j_f64);
        
        let term_ga = alpha * (
            (1.0 + coeff) * c_ga[j + 1] 
            - 2.0 * c_ga[j] 
            + (1.0 - coeff) * c_ga[j - 1]
        );
        
        let term_as = beta * (
            (1.0 + coeff) * c_as[j + 1] 
            - 2.0 * c_as[j] 
            + (1.0 - coeff) * c_as[j - 1]
        );

        c_ga_next[j] = c_ga[j] + term_ga + upsilon * flux_term[j] 
            - omega * c_ga[j] * c_as[j];
        
        c_as_next[j] = c_as[j] + term_as - gamma * c_as[j] + kappa 
            - omega * c_ga[j] * c_as[j];
    }

    // Boundary conditions
    c_ga_next[nr] = 0.0;
    c_as_next[nr] = epsilon;

    (c_ga_next, c_as_next)
}

fn update_rho(mut rho: f64, time: f64, sigma: f64, rho0: f64) -> f64 {
    let rhs = rho0.powi(3) * (3.0 * rho0.ln() - 1.0) + 9.0 * sigma * time;
    
    for _ in 0..100 {
        let p = rho.powi(3) * (3.0 * rho.ln() - 1.0) - rhs;
        if p.abs() < 1e-6 {
            break;
        }
        
        let p_prime = 3.0 * rho.powi(2) * (3.0 * rho.ln() - 2.0);
        if p_prime.abs() < 1e-12 {
            break;
        }
        
        rho = (rho - p / p_prime).clamp(1e-6, 0.99);
    }
    
    rho
}