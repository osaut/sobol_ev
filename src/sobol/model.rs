use params::Params;
use std::num::ln;

mod params;


pub trait Model {
  fn setup(params: &Params) -> Self;
  fn run(&self) -> f64;
}


 pub struct Gompertz {
  alpha: f64,
  K : f64,
  C0: f64
 }

impl Model for Gompertz {
  #[allow(dead_code)]
  fn setup(params: &Params) -> Gompertz {
    Gompertz{ alpha: params.get_float("alpha"), K: params.get_float("K"), C0 : params.get_float("C0")}
  }

  #[allow(dead_code)]
  fn run(&self) -> f64 {
    let tmax=10f64;
    let dt=0.001;
    let mut t=0.0;
    let mut C=self.C0;
    for _ in range(0, (tmax/dt) as uint) {
      C=C+self.alpha*dt*ln(self.K/C)*C;
      t=t+dt;
    }
    C
  }
}
