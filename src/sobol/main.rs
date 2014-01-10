#[crate_id = "sobol#0.1"];
#[comment = "Sobol indices evaluator"];
#[license = "MIT"];

use params::{Params, Value, Float, Range};
use model::{Model,Gompertz};
use std::vec;
use std::hashmap::HashMap;
pub mod model;
pub mod params;

fn sample(pr: &Params, num_real : uint) -> ~[Params] {
    let mut samples : ~[Params] = ~[];
    for _ in range(0, num_real) {
        samples.push(pr.realize());
    }
    samples
}
#[test]
fn test_sample() {
    let mut ranges=Params::new();
    ranges.insert("alpha", Range(0.0, 1.0));
    let samples=sample(&ranges, 10);

    assert_eq!(samples.len(), 10);
    for pp in samples.iter() {
        let value=pp.get_float("alpha");
        assert!((value>=0.0)&&(value<=1.0));
    }
}

#[allow(dead_code)]
fn eval_model(pr: &Params) -> f64 {
    let model : Gompertz = Model::setup(pr);
    model.run()
}

fn get_couples(n : uint) -> ~[(uint, uint)] {
    let mut res : ~[(uint, uint)]=~[];
    for i in range(0, n) {
        for j in range(i, n) {
            if (i!=j) {
                res.push( (i,j));
            }
        }
    }
    res
}

#[test]
fn test_get_couples() {
    let res=get_couples(3);

    assert_eq!(res[0], (0,1));
    assert_eq!(res[1], (0,2));
    assert_eq!(res[2], (1,2));
}

// First order Sobol indices
#[allow(dead_code)]
fn calc_sobol_1(pr: &Params, nsamp : uint) -> HashMap<~str, f64> {
    let samples_1=sample(pr, nsamp);
    let samples_2=sample(pr, nsamp);

    let big_one : ~[Params] = vec::append(samples_1.clone(), samples_2);

    // Average and deviation computaion
    let resvec=big_one.map(eval_model);
    let avg=resvec.iter().fold(0.0, |acc, &item| acc+item)/(resvec.len() as f64);
    let var=resvec.iter().fold(0.0, |acc, &item| acc+item*item)/(resvec.len() as f64)-avg*avg;

    let mut Sp: HashMap<~str, f64> = HashMap::new();
    for (k,_v) in pr.pp.iter() {
        let mut sump=0.0f64;
        for snum in range(0, nsamp) {
            let mut new_p=samples_2[snum].clone();
            new_p.set(k.to_owned(),Float(samples_1[snum].get_float(k.to_owned())));

            sump=sump+eval_model(&new_p)*eval_model(&samples_1[snum]);
        }
        let Up=sump/(nsamp as f64);
        Sp.insert(k.to_owned(), (Up-avg*avg)/var);
    }

    Sp
}

// Second order Sobol indices
#[allow(dead_code)]
fn calc_sobol_2(pr: &Params, sob1 : &HashMap<~str, f64>, nsamp : uint) -> ~[(~str, ~str, f64)] {
    let samples_1=sample(pr, nsamp);
    let samples_2=sample(pr, nsamp);
    let couples=get_couples(pr.len());

    let big_one : ~[Params] = vec::append(samples_1.clone(), samples_2);

    // Average and deviation computaion
    let resvec=big_one.map(eval_model);
    let avg=resvec.iter().fold(0.0, |acc, &item| acc+item)/(resvec.len() as f64);
    let var=resvec.iter().fold(0.0, |acc, &item| acc+item*item)/(resvec.len() as f64)-avg*avg;

    let mut Spr : ~[(~str, ~str, f64)]=~[];
    let keys=pr.keys();

    for &cpl in couples.iter() {
        let (i,j) = cpl;
        let ref key1=keys[i]; let ref key2=keys[j];
        let mut sump=0.0f64;
        for snum in range(0, nsamp) {
            let mut new_p=samples_2[snum].clone();
            new_p.set(key1.to_owned(),Float(samples_1[snum].get_float(key1.to_owned())));
            new_p.set(key2.to_owned(),Float(samples_1[snum].get_float(key2.to_owned())));

            sump=sump+eval_model(&new_p)*eval_model(&samples_1[snum]);
        }
        let Upr=sump/(nsamp as f64);
        let Vp=*sob1.find_equiv(key1).unwrap();
        let Vr=*sob1.find_equiv(key2).unwrap();
        Spr.push((key1.to_owned(), key2.to_owned(), (Upr-avg*avg-Vr-Vp)/var));
    }

    Spr
}


//
// Entry point
//
#[allow(dead_code)]
fn main() {
  let mut ranges=Params::new();
  ranges.insert("alpha", Range(0.01, 1.0)); ranges.insert("K", Range(0.01, 1.0)); ranges.insert("C0", Range(0.001, 0.01));

  let samples=sample(&ranges,5000);
  let resvec=samples.map(eval_model);

  let avg=resvec.iter().fold(0.0, |acc, &item| acc+item)/(resvec.len() as f64);
  let var=resvec.iter().fold(0.0, |acc, &item| acc+item*item)/(resvec.len() as f64)-avg*avg;

  println!("Mean value = {}, Standard Deviation = {}\n", avg, var);
  let sob_1=calc_sobol_1(&ranges,15000);
  println!("First order Sobol indices");
  for (k, v) in sob_1.iter() {
    println!("\t{:s} = {:f}", *k, *v);
  }
  let sob_2=calc_sobol_2(&ranges, &sob_1, 7500);
  println!("\nSecond order Sobol indices");
  for vv in sob_2.iter() {
    let (ref k1,ref k2,v) = *vv;
    println!("\t{:s} x {:s} = {:f}", *k1, *k2, v);
  }
}
