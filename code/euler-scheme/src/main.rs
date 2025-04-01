mod file;
mod plot;

const BOLTZMANN: f64 = 8.62e-5;

struct GrowthConditions {
    temperature_celcius: f64,
    temperature_ev: f64, // kT, eV
    arsenic_flux: f64, // per nm^2 ?
}

impl GrowthConditions {
    fn new(temperature_celcius: f64, arsenic_flux: f64) -> GrowthConditions {
        GrowthConditions {
            temperature_celcius, 
            temperature_ev: BOLTZMANN*(temperature_celcius + 273.4), 
            arsenic_flux,
        }
    }
}

struct Crystal {
    lattice_parameter_a0: f64,
    lattice_parameter_h0: f64,
    max_concentration_c0: f64,
    arsenic_diffusion_energy: f64, // EAs-EGa, eV
    arsenic_desorption_energy: f64, // Ea-EGa, eV
    reaction_coefficient_kr: f64, // for k1*CGa*CAs
}

impl Crystal {
    fn default() -> Crystal {
        Crystal {
            lattice_parameter_a0: 1.0,
            lattice_parameter_h0: 1.0,
            max_concentration_c0: 1.0,
            arsenic_diffusion_energy: 0.0, // EAs-EGa, eV
            arsenic_desorption_energy: 0.5, // Ea-EGa, eV
            reaction_coefficient_kr: 0.1, // for k1*CGa*CAs
        }
    }
}

struct SurfaceAtoms {
    previous_concentration: Vec<f64>, // distance dependence
    current_concentration: Vec<f64>, // distance dependence
    diffusion_coefficient: f64,
    current_source_flux: Vec<f64>,
    desorption_time: f64,
    infinity_concentration: f64,
}

struct Droplet {
    gallium_molecular_volume: f64,
    contact_angle_theta: f64,
    form_factor_b: f64,
    radius_initial: f64, // in nm ?
    ga_flux_width_w: f64,
    radius_infinity: f64,
    current_radius: Vec<f64>, // time dependence
    number_of_iterations: Vec<usize>, // time dependence
    sigma_parameter: f64, // for Newton iterations
    c_parameter: f64, // for Newton iterations
    final_time: f64, // time need for total droplet depletion
}

impl Droplet {
    fn set_form_factor(&mut self){
        self.form_factor_b = (8.0-9.0*self.contact_angle_theta.cos()+(3.0*self.contact_angle_theta).cos())/(3.0*self.contact_angle_theta.sin()-(3.0*self.contact_angle_theta).sin())
    }
    fn get_flux(&self, distance: f64, time_index: usize) -> f64 {
        let rd = self.current_radius[time_index];
        (-(distance-rd).powi(2)/self.ga_flux_width_w.powi(2)).exp()/(3.545*rd/self.ga_flux_width_w+0.187*self.ga_flux_width_w/(self.radius_infinity-3.156*self.ga_flux_width_w))/(self.radius_infinity/rd).ln()
    }
    fn newton_function(radius_scaled: f64) -> f64 {
        radius_scaled.powi(3)*(1.0-3.0*radius_scaled.ln())
    }
    fn set_c_parameter(&mut self){
        self.c_parameter = Self::newton_function(self.radius_initial/self.radius_infinity);
    }
    fn set_sigma_parameter(&mut self, numerator: f64){
        self.sigma_parameter = 18.0*numerator/self.form_factor_b/self.radius_infinity.powi(3)
    }
    fn set_final_time(&mut self){
        self.final_time = self.c_parameter/self.sigma_parameter
    }
}

struct Newton {
    max_iterations: usize,
    tolerance: f64,
    current_value: f64,
    next_value: f64,
    answer: f64,
    actual_iterations: usize,
    actual_error: f64,
}

impl Newton {
    fn new(max_iterations: usize, tolerance: f64) -> Newton {
        Newton {
            max_iterations,
            tolerance,
            current_value: 0.0,
            next_value: 0.0,
            answer: 0.0,
            actual_iterations: 0,
            actual_error: 0.0,
        }
    }
    fn run(&mut self, function_ratio: &dyn Fn(f64) -> f64, initial_guess: f64){
        self.current_value = initial_guess;
        for it in 0..self.max_iterations {
            self.actual_error = function_ratio(self.current_value);
            self.next_value = self.current_value - self.actual_error;
            self.actual_error = self.actual_error.abs();
            self.actual_iterations = it+1;
            self.answer = self.current_value;
            if self.actual_error < self.tolerance*self.current_value.abs() {
                break;
            }
        }
    }
}

struct ExtraParameters {
    alpha: f64,
    beta: f64,
    gamma: f64,
    epsilon: f64,
    eta: f64,
    kappa: f64,
    omega: f64,
    upsilon: f64,
}

struct Grid {
    distance_step_dr: f64,
    time_step_dt: f64,
    number_of_points_r: usize,
    number_of_time_steps: usize,
    r_points: Vec<f64>,
    t_moments: Vec<f64>,
}



fn main() {
    println!("Hello, world!");
}
