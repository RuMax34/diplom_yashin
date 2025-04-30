mod file;
mod erf;

fn main() {
    // Fixed parameters
    let a0: f64 = 0.565/2.0; // nm
    let omega_ga = a0.powi(3);     // nm^3
    let c0 = 1.0/a0.powi(2);           // nm^-2
    let e_ga = 0.5;  // eV
    let e_as: f64 = 0.5; // eV
    let e_a: f64 = 1.0;  // eV

    // User-defined parameters
    let temperature = 300.0; // temperature in C
    let kt = 8.62e-5*(273.0+temperature); // Temperature in eV
    let nu0 = 2e6*kt/4.136; // GHz or ns^-1

    let flux_relative = 1e-10; // relative flux, correct order, do not change the order
    let f_as = c0*nu0*flux_relative;         // As flux, nm^(-2)/ns

    let d_ga = nu0/c0*(-e_ga/kt).exp(); // diffusion coefficient
    let d_as = nu0/c0*(-e_as/kt).exp(); // nm^2 / ns
    let tau_as = (e_a / kt).exp()/nu0;
    let diff_length_as = (d_as*tau_as).sqrt(); // nm
    println!("As diffusion length is {} nm", diff_length_as);
    let kr = d_ga;           // nm^2 / ns ! just a guess, may need to change
    let w: f64 = 3.0;            // width of Gaussian flux, nm
    let r_inf = 300.0;      // domain size, nm
    let theta = 60.0;       // contact angle in degrees
    let rd0: f64 = 30.0; // nm
    let rho0_initial: f64 = rd0 / r_inf; // initial rho
    let h0 = a0;           // cell height per layer

    // Derived parameters
    let theta_rad: f64 = theta / 180.0 * std::f64::consts::PI;
    let numerator_b = 8.0 - 9.0 * theta_rad.cos() + (3.0 * theta_rad).cos();
    let denominator_b = 3.0 * theta_rad.sin() - (3.0 * theta_rad).sin();
    let b_theta = numerator_b / denominator_b;
    let sigma = (2.0 * omega_ga * d_ga * c0) / (b_theta * r_inf.powi(3));

    // Simulation time calculation
    let rho_end: f64 = 0.001;
    let term_initial = rho0_initial.powi(3) * (3.0 * rho0_initial.ln() - 1.0);
    let term_final = rho_end.powi(3) * (3.0 * rho_end.ln() - 1.0);
    let t_total = (term_final - term_initial) / (9.0 * sigma);

    // Grid setup
    let nr = 400;
    let dr = r_inf / nr as f64;
    let r: Vec<f64> = (0..=nr).map(|i| i as f64 * dr).collect();
    let dt = dr.powi(2) / d_ga / 10.0;
    let nt = (t_total / dt) as usize;

    // Number of snapshots
    let ns: usize = 10;

    // Initialize arrays
    let mut c_ga: Vec<f64> = r.iter().map(|ri| (-(ri / w - rd0 / w).powi(2)).exp()).collect();
    let mut c_as = vec![(f_as * tau_as) / c0; nr + 1];
    let mut h = vec![0.0; nr + 1];
    let mut h_history = vec![vec![0.0; nr + 1]; ns+1];
    let mut ga_history = vec![vec![0.0; nr + 1]; ns+1];
    let mut as_history = vec![vec![0.0; nr + 1]; ns+1];
    let mut time_history = vec![0.0; nt];
    let mut rd_history = vec![0.0; nt];
    let mut tau_ga_history = vec![0.0; nt];
    let mut rho = rho0_initial;
    let mut time = 0.0;

    // Precompute coefficients
    let alpha = d_ga * dt / dr.powi(2);
    let beta =  d_as * dt / dr.powi(2);
    let omega = c0 * kr * dt;
    let kappa = dt * f_as / c0;
    let gamma = dt / tau_as;
    let epsilon = (f_as * tau_as) / c0;
    let upsilon = 2.0 * d_ga * dt / w.powi(2);

    let mut js = 0;

    for jt in 0..nt {
        // Compute flux term
        let rd = rho * r_inf;
        time_history[jt] = jt as f64 * dt;
        rd_history[jt] = rd / rd0;
        let x = rd / w;
        let p_val = r_inf / w;
        let a_p = 3.545;
        let b_p = 0.187 / (p_val - 3.156);
        let denominator = (a_p * x + b_p) * (r_inf / rd).ln();
        // let denominator = ((-x*x).exp()-(-(p_val-x).powi(2)).exp() + 2.0/std::f64::consts::FRAC_2_SQRT_PI*x*(erf(x)+erf(p_val-x))) * (r_inf / rd).ln();
        let flux_term: Vec<f64> = r.iter().map(|ri| (-(ri / w - x).powi(2)).exp() / denominator).collect();
        tau_ga_history[jt] = denominator;
        // Update concentrations
        (c_ga, c_as) = update_concentrations(
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

        // Update droplet radius
        time += dt;
        rho = update_rho(time, sigma, rho, term_initial);

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

    file::save_columns_to_file(&vec![time_history, rd_history, tau_ga_history], "results", "droplet.dat");
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

    //c_ga_next[0] = f64::max(ga_droplet_edge[0], c_ga_next[0]);
    
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
        
        //c_ga_next[j] = f64::max(ga_droplet_edge[j], c_ga_next[j]);

        c_as_next[j] = c_as[j] + term_as - gamma * c_as[j] + kappa 
            - omega * c_ga[j] * c_as[j];
    }

    // Boundary conditions
    c_ga_next[nr] = 0.0;
    c_as_next[nr] = epsilon;

    (c_ga_next, c_as_next)
}

fn update_rho(time: f64, sigma: f64, rho0: f64, term_initial: f64) -> f64 {
    let rhs = term_initial + 9.0 * sigma * time;
    let mut rho = rho0;
    for _ in 0..100 {
        let p = rho.powi(3) * (3.0 * rho.ln() - 1.0) - rhs;
        if p.abs() < 1e-6 {
            break;
        }
        
        let p_prime = 3.0 * rho.powi(2) * (3.0 * rho.ln() - 2.0);
        if p_prime.abs() < 1e-12 {
            break;
        }
        
        rho = rho - p / p_prime;
    }

    rho
}