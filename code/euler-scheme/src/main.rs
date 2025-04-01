mod file;
mod plot;

struct Variables {
    temperature_celsius: f64,
    temperature_ev: f64, // kT, eV
    arsenic_flux: f64, // per nm^2 ?
    initial_radius: f64, // in nm ?
}

struct InitialParameters {
    lattice_parameter_a0: f64,
    lattice_parameter_h0: f64,
    gallium_molecular_volume: f64,
    arsenic_diffusion_energy: f64, // EAs-EGa, eV
    arsenic_desorption_energy: f64, // Ea-EGa, eV
    reaction_coefficient_kr: f64, // for k1*CGa*CAs
    ga_flux_width_w: f64,
}

fn main() {
    println!("Hello, world!");
}
