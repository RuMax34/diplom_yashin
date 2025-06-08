use indicatif::ProgressBar;

mod file;

fn main() {
    // Fixed parameters
    let a0: f64 = 0.565/2.0; // nm
    let omega_ga = a0.powi(3);     // nm^3
    let c0 = 1.0/a0.powi(2);           // nm^-2
    let e_ga = 0.5;  // eV
    let e_as: f64 = 0.5; // eV
    let e_a: f64 = 1.0;  // eV

    // User-defined parameters
    let temperature = 250.0; // temperature in C
    let kt = 8.62e-5*(273.0+temperature); // Temperature in eV
    let nu0 = 2e6*kt/4.136; // GHz or ns^-1

    let flux_relative = 4.0*1e-10; // relative flux, correct order, do not change the order
    let f_as = c0*nu0*flux_relative;         // As flux, nm^(-2)/ns

    let d_ga = nu0/c0*(-e_ga/kt).exp(); // diffusion coefficient
    let d_as = nu0/c0*(-e_as/kt).exp(); // nm^2 / ns
    let tau_as = (e_a / kt).exp()/nu0;
    let diff_length_as = (d_as*tau_as).sqrt(); // nm
    println!("As diffusion length is {} nm", diff_length_as);
    let kr = d_ga;           // nm^2 / ns ! just a guess, may need to change
    let w: f64 = 3.0;            // width of Gaussian flux, nm
    let r_inf = 200.0;      // domain size, nm
    let theta = 60.0;       // contact angle in degrees
    let rd_initial: f64 = 30.0; // nm
    let h0 = a0;           // cell height per layer

    // Derived parameters
    let theta_rad: f64 = theta / 180.0 * std::f64::consts::PI;
    let numerator_b = 8.0 - 9.0 * theta_rad.cos() + (3.0 * theta_rad).cos();
    let denominator_b = 3.0 * theta_rad.sin() - (3.0 * theta_rad).sin();
    let b_theta = numerator_b / denominator_b;

    // Grid setup
    let nr = 1000;
    let dr = r_inf / nr as f64;
    let dr2 = dr*dr;
    let r: Vec<f64> = (0..=nr).map(|i| i as f64 * dr).collect();
    let dt = dr2 / d_ga / 10.0;

    // Precompute coefficients
    let alpha = d_ga * dt / dr2;
    let beta =  d_as * dt / dr2;
    let omega = c0 * kr * dt;
    let kappa = dt * f_as / c0;
    let gamma = dt / tau_as;
    let epsilon = (f_as * tau_as) / c0;
    // See Ga flux calculation below, same as before
    let flux_factor = dt * 2.0*d_ga/w.powi(2);
    let flux_radius2 = 4.0*d_ga*c0/f_as;
    let r03_droplet_coefficient = 2.0*omega_ga*omega*c0*dr2/b_theta;

    // Initial time limit
    let nt_max = 1_200_000;

    // Progress bar
    let pb = ProgressBar::new(nt_max as u64);

    // Minimal allowed droplet radius
    let rd_end: f64 = dr;

    // PRELIMINARY CALCULATION TO FIND FINAL TIME
    
    // Initialize arrays
    //let rj: Vec<f64> = (0..=nr).map(|j| j as f64 * dr).collect();
    let mut c_ga: Vec<f64> = vec![0.0; nr + 1];
    let mut c_as = vec![(f_as * tau_as) / c0; nr + 1];
    let mut c_ga_next = c_ga.clone();
    let mut c_as_next = c_as.clone();
    let mut rd = rd_initial;
    
    let mut jt_final = nt_max;

    for jt in 0..nt_max {
        pb.inc(1);

        // Update droplet radius
        rd = update_rho(rd, r03_droplet_coefficient, &c_ga, &c_as);

        if rd > rd_end {
            // Update concentrations
            update_concentrations(
                &c_ga,
                &c_as,
                &mut c_ga_next,
                &mut c_as_next,
                rd,
                w,
                dr,
                alpha,
                beta,
                omega,
                gamma,
                kappa,
                epsilon,
                flux_factor,
                flux_radius2,
                r_inf,
                nr,
            );
            c_ga.copy_from_slice(&c_ga_next);
            c_as.copy_from_slice(&c_as_next);
        } else {
            jt_final = jt;
            break;
        }
    }

    let nt = jt_final;

    pb.finish();

    println!("Completed first loop, nt = {}", nt);

    // Progress bar
    let pb = ProgressBar::new(nt as u64);

    // Number of snapshots
    let ns: usize = 10;

    // Initialize arrays again
    let mut c_ga: Vec<f64> = vec![0.0; nr + 1];
    let mut c_as = vec![(f_as * tau_as) / c0; nr + 1];
    let mut c_ga_next = c_ga.clone();
    let mut c_as_next = c_as.clone();
    let mut h = vec![0.0; nr + 1];
    let mut h_history = vec![vec![0.0; nr + 1]; ns+1];
    let mut ga_history = vec![c_ga.clone(); ns+1];
    let mut as_history = vec![c_as.clone(); ns+1];
    let mut time_history = vec![0.0; nt];
    let mut rd_history = vec![0.0; nt];
    let mut rd = rd_initial;


    let mut js = 0;

    for jt in 0..nt {
        pb.inc(1);

        // Update droplet radius
        rd = update_rho(rd, r03_droplet_coefficient, &c_ga, &c_as);

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

        if rd > rd_end {
            time_history[jt] = jt as f64 * dt;
            rd_history[jt] = rd;

            // Update concentrations
            update_concentrations(
                &c_ga,
                &c_as,
                &mut c_ga_next,
                &mut c_as_next,
                rd,
                w,
                dr,
                alpha,
                beta,
                omega,
                gamma,
                kappa,
                epsilon,
                flux_factor,
                flux_radius2,
                r_inf,
                nr,
            );
            c_ga.copy_from_slice(&c_ga_next);
            c_as.copy_from_slice(&c_as_next);
        } else {
            break;
        }
    }

    pb.finish();
    println!("Completed second loop, now saving files!");

    h_history.push(r.clone());
    ga_history.push(r.clone());
    as_history.push(r);

    file::save_columns_to_file(&vec![time_history, rd_history], "results", "droplet.dat");
    file::save_columns_to_file(&h_history, "results", "height.dat");
    file::save_columns_to_file(&ga_history, "results", "c_ga.dat");
    file::save_columns_to_file(&as_history, "results", "c_as.dat");

}

fn update_concentrations(
    c_ga: &[f64],
    c_as: &[f64],
    c_ga_next: &mut [f64],
    c_as_next: &mut [f64],
    rd: f64,
    w: f64,
    dr: f64,
    alpha: f64,
    beta: f64,
    omega: f64,
    gamma: f64,
    kappa: f64,
    epsilon: f64,
    flux_factor: f64,
    flux_radius2: f64,
    r_inf: f64,
    nr: usize,
) {

    // Handle j=0
    // Ga flux calculation using the formulas, derived before, now gives a correct result
    c_ga_next[0] = c_ga[0] + 4.0 * alpha * (c_ga[1] - c_ga[0])
        - omega * c_ga[0] * c_as[0] + flux_factor*wp_j(0.0, rd, w)/(3.545*rd/w+0.187*w/(r_inf-3.156*w))/(r_inf/rd).ln();
    
    c_as_next[0] = c_as[0] + 4.0 * beta * (c_as[1] - c_as[0]) 
        - gamma * c_as[0] 
        + kappa 
        - omega * c_ga[0] * c_as[0];

    // Handle interior points
    for j in 1..nr {
        let j_f64 = j as f64;
        let coeff = 1.0 / (2.0 * j_f64);
        
        // gallium atoms
        let term_ga = alpha * (
            (1.0 + coeff) * c_ga[j + 1] 
            - 2.0 * c_ga[j] 
            + (1.0 - coeff) * c_ga[j - 1]
        );
        // Ga flux calculation 3rd version from the notes
        let conc_ratio = flux_radius2/rd.powi(2);
        let rd_rinf_ratio2 = r_inf.powi(2)/rd.powi(2);
        let additional_factor = 1.0 + (rd_rinf_ratio2 - 1.0 + rd_rinf_ratio2.ln())/conc_ratio;
        c_ga_next[j] = c_ga[j] + term_ga
        - omega * c_ga[j] * c_as[j] + flux_factor*additional_factor*wp_j(j_f64*dr, rd, w)/(3.545*rd/w+0.187*w/(r_inf-3.156*w))/(r_inf/rd).ln();
        
        // as atoms
        let term_as = beta * (
            (1.0 + coeff) * c_as[j + 1] 
            - 2.0 * c_as[j] 
            + (1.0 - coeff) * c_as[j - 1]
        );
        c_as_next[j] = c_as[j] + term_as - gamma * c_as[j] + kappa 
            - omega * c_ga[j] * c_as[j];
    }

    // Boundary conditions
    c_ga_next[nr] = 0.0;
    c_as_next[nr] = epsilon;
}

fn update_rho(rd: f64, r03: f64, c_ga: &[f64], c_as: &[f64]) -> f64 {
    let sum = multiply_sum_with_index(c_ga, c_as);
    rd - r03/rd*sum
}

fn multiply_sum_with_index(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "Arrays must have the same length");
    a.iter()
        .zip(b.iter())
        .enumerate()
        .map(|(i, (&a_val, &b_val))| a_val * b_val * i as f64)
        .sum()
}

fn wp_j(x: f64, x0: f64, w: f64) -> f64 {
    let x2 = -((x-x0)/w).powi(2);
    x2.exp()
}
