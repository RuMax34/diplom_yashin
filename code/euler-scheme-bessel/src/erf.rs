/// Approximates the error function (erf) using Abramowitz and Stegun's approximation.
/// This approximation has a maximum error of `1.5e-7` for all real `x`.
pub fn erf_approx(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();

    // Constants for the approximation
    const P: f64 = 0.3275911;
    const A1: f64 = 0.254829592;
    const A2: f64 = -0.284496736;
    const A3: f64 = 1.421413741;
    const A4: f64 = -1.453152027;
    const A5: f64 = 1.061405429;

    let t = 1.0 / (1.0 + P * x);
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;

    let inner = A1 * t + A2 * t2 + A3 * t3 + A4 * t4 + A5 * t5;
    let erf = 1.0 - inner * (-x * x).exp();

    sign * erf
}

/// Approximates the error function (erf) using the exponential form with polynomial correction.
/// This approximation is based on: 
/// erf(x) ≈ sqrt(1 - exp(-(4x²/π) * P(x))), where P(x) is a rational polynomial.
pub fn _erf_approx_exp(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();

    if x == 0.0 {
        return 0.0;
    }

    // Compute the polynomial correction term P(x)
    let numerator = 1.0 + (10.0 - std::f64::consts::PI.powi(2)) / (5.0 * (std::f64::consts::PI - 3.0) * std::f64::consts::PI) * x.powi(2);
    let denominator = 1.0 + ((120.0 - 60.0 * std::f64::consts::PI + 7.0 * std::f64::consts::PI.powi(2)) 
                          / (15.0 * (std::f64::consts::PI - 3.0) * std::f64::consts::PI)) * x.powi(2);
    let p_x = numerator / denominator;

    // Compute the exponential term
    let exponent = -(4.0 * x.powi(2)) / std::f64::consts::PI * p_x;
    let erf = (1.0 - exponent.exp()).sqrt();

    sign * erf
}