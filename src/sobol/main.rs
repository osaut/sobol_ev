#[crate_id = "sobol#0.1"];
#[comment = "Sobol indices evaluator"];
#[license = "MIT"];

extern crate collections;

use params::{Params, Float, Range};
use model::{Model,Gompertz};
use std::vec;
use std::vec::Vec;
mod model;
mod params;

fn sample(pr: &Params, num_real : uint) -> Vec<Params> {
    let mut samples : Vec<Params> = Vec::new();
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
    model.run(pr.get_float("tmax"))
}

fn get_couples(n : uint) -> Vec< (uint, uint) > {
    let mut res : Vec< (uint, uint) > = Vec::new();
    for i in range(0, n) {
        for j in range(i, n) {
            if i!=j {
                res.push( (i,j));
            }
        }
    }
    res
}

#[test]
fn test_get_couples() {
    let res=get_couples(3);

    assert_eq!(res.get(0), &(0,1));
    assert_eq!(res.get(1), &(0,2));
    assert_eq!(res.get(2), &(1,2));
}

// First order Sobol indices
#[allow(dead_code)]
fn calc_sobol_1(pr: &Params, nsamp : uint) -> collections::HashMap<~str, f64> {
    let samples_1=sample(pr, nsamp);
    let samples_2=sample(pr, nsamp);

    let big_one : Vec<Params> = vec::append(samples_1.clone(), samples_2.as_slice());

    // Average and deviation computaion
    let mut resvec=big_one.iter().map(eval_model);
    let avg=resvec.fold(0f64, |acc, item| acc+item)/(resvec.len() as f64);
    let var=resvec.fold(0f64, |acc, item| acc+item*item)/(resvec.len() as f64)-avg*avg;

    let mut sp: collections::HashMap<~str, f64> = collections::HashMap::new();
    let vkeys=pr.varying_keys();
    for k in vkeys.iter() {
        let mut sump=0.0f64;
        for snum in range(0, nsamp) {
            let mut new_p=samples_2.get(snum).clone();
            new_p.set(k.to_owned(),Float(samples_1.get(snum).get_float(k.to_owned())));

            sump=sump+eval_model(&new_p)*eval_model(samples_2.get(snum));
        }
        let up=sump/(nsamp as f64);
        sp.insert(k.to_owned(), (up-avg*avg)/var);
    }

    sp
}

// Second order Sobol indices
#[allow(dead_code)]
fn calc_sobol_2(pr: &Params, sob1 : &collections::HashMap<~str, f64>, nsamp : uint) -> Vec<(~str, ~str, f64)> {
    let samples_1=sample(pr, nsamp);
    let samples_2=sample(pr, nsamp);
    let keys=pr.varying_keys();
    let couples=get_couples(keys.len());

    let big_one : Vec<Params> = vec::append(samples_1.clone(), samples_2.as_slice());

    // Average and deviation computaion
    let mut resvec=big_one.iter().map(eval_model);
    let avg=resvec.fold(0f64, |acc, item| acc+item)/(resvec.len() as f64);
    let var=resvec.fold(0f64, |acc, item| acc+item*item)/(resvec.len() as f64)-avg*avg;

    let mut spr : Vec<(~str, ~str, f64)> = Vec::new();

    for &cpl in couples.iter() {
        let (i,j) = cpl;
        let ref key1=keys.get(i); let ref key2=keys.get(j);
        let mut sump=0.0f64;
        for snum in range(0, nsamp) {
            let mut new_p=samples_2.get(snum).clone();
            new_p.set(key1.to_owned(),Float(samples_1.get(snum).get_float(key1.to_owned())));
            new_p.set(key2.to_owned(),Float(samples_1.get(snum).get_float(key2.to_owned())));

            sump=sump+eval_model(&new_p)*eval_model(samples_1.get(snum));
        }
        let upr=sump/(nsamp as f64);
        let vp=*sob1.find_equiv(*key1).unwrap();
        let vr=*sob1.find_equiv(*key2).unwrap();
        spr.push((key1.to_owned(), key2.to_owned(), (upr-avg*avg-vr-vp)/var));
    }

    spr
}

// Total Sobol indices
fn calc_sobol_total(pr: &Params, nsamp: uint) -> collections::HashMap<~str, f64> {
    let samples_1=sample(pr, nsamp);
    let samples_2=sample(pr, nsamp);

    let big_one : Vec<Params> = vec::append(samples_1.clone(), samples_2.as_slice());

    // Average and deviation computaion
    let mut resvec=big_one.iter().map(eval_model);
    let avg=resvec.fold(0.0, |acc, item| acc+item)/(resvec.len() as f64);
    let var=resvec.fold(0.0, |acc, item| acc+item*item)/(resvec.len() as f64)-avg*avg;

    let mut sp: collections::HashMap<~str, f64> = collections::HashMap::new();
    let vkeys=pr.varying_keys();
    for k in vkeys.iter() {
        let mut sump=0.0f64;
        for snum in range(0, nsamp) {
            let mut new_p=samples_1.get(snum).clone();
            new_p.set(k.to_owned(),Float(samples_2.get(snum).get_float(k.slice_from(0))));

            sump=sump+eval_model(&new_p)*eval_model(samples_1.get(snum));
        }
        let up=sump/(nsamp as f64);
        sp.insert(k.to_owned(), 1.0-(up-avg*avg)/var);
    }

    sp
}

//
// Entry point
//
#[allow(dead_code)]
fn main() {
  let mut ranges=Params::new();
  ranges.insert("alpha", Range(0.01, 1.0)); ranges.insert("K", Range(0.01, 1.0));
  ranges.insert("C0", Range(0.001, 0.01)); ranges.insert("tmax", Float(10f64));

  let samples=sample(&ranges,5000);
  let mut resvec=samples.iter().map(eval_model);

  let avg=resvec.fold(0f64, |acc, item| acc+item)/(resvec.len() as f64);
  let var=resvec.fold(0f64, |acc, item| acc+item*item)/(resvec.len() as f64)-avg*avg;

  println!("Mean value = {}, Standard Deviation = {}\n", avg, var);

  println!("First order Sobol indices");
  let sob_1=calc_sobol_1(&ranges,15000);
  for (k, v) in sob_1.iter() {
    println!("\t{:s} = {:f}", *k, *v);
  }

  println!("\nSecond order Sobol indices");
  let sob_2=calc_sobol_2(&ranges, &sob_1, 7500);
  for vv in sob_2.iter() {
    let (ref k1,ref k2,v) = *vv;
    println!("\t{:s} x {:s} = {:f}", *k1, *k2, v);
  }

  println!("\nTotal Sobol indices");
  let sob_t=calc_sobol_total(&ranges, 15000);
  for (k, v) in sob_t.iter() {
    println!("\t{:s} = {:f}", *k, *v);
  }


}

