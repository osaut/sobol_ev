extern crate collections;
extern crate rand;

use std::vec::Vec;

//use collections::HashMap;

#[deriving(Clone)]
pub enum Value {
  Float(f64),
  Range(f64,f64)
}

// Set of parameters
#[deriving(Clone)]
pub struct Params {
  pp : collections::HashMap<~str, Value>
}

impl Params {
  #[allow(dead_code)]
  pub fn new() -> Params {
    Params{ pp: collections::HashMap::new() }
  }

  #[allow(dead_code)]
  pub fn insert(&mut self, tag: &str, value: Value) {
    self.pp.insert(tag.to_owned(), value);
  }

  #[allow(dead_code)]
  pub fn get(&self, tag: &str) -> Value {
    *self.pp.find_equiv(&tag).unwrap()
  }

  #[allow(dead_code)]
  pub fn get_float(&self, tag: &str) -> f64 {
    let value=self.pp.find_equiv(&tag);
    match value.unwrap() {
      &Float(f) => { f },
      _ => { fail!(); }
    }
  }

  #[allow(dead_code)]
  pub fn set(&mut self, tag: &str, val: Value) {
    if self.pp.contains_key_equiv(&tag) {
      self.pp.insert(tag.to_owned(), val);
    }
    else {
      fail!();
    }
  }

  #[allow(dead_code)]
  pub fn realize(&self) -> Params {
    let mut hash : collections::HashMap<~str, Value> = collections::HashMap::new();
    for (key, value) in self.pp.iter() {
      match *value {
        Float(f) => { hash.insert(key.to_owned(), Float(f)); },
        Range(f1, f2) => {
          let t=rand::random::<f64>();
          hash.insert(key.to_owned(), Float(t*f1+(1.0-t)*f2));
        }
      }
    }
    Params{ pp: hash }
  }


  #[allow(dead_code)]
  pub fn print(&self) {
    for (key, value) in self.pp.iter() {
      match *value {
        Float(f) => {println!("{} = {:f}", *key, f);}
        Range(f1,f2) => {println!("{} in [{:f}, {:f}]", *key, f1, f2);}
      }
    }
  }

  #[allow(dead_code)]
  pub fn len(&self) -> uint {
    self.pp.len()
  }

  pub fn keys<'s>(&'s self) -> Vec< &'s str > {
    let mut res : Vec < &str> = Vec::new();
    for (key, _) in self.pp.iter() {
      res.push(key.slice_from(0));
    }
    res
  }

  pub fn varying_keys<'s>(&'s self) -> Vec<&'s str> {
   let mut res : Vec<&str> = Vec::new();
   for (key, value) in self.pp.iter() {
      match *value {
        Float(_) => {},
        Range(_,_) => { res.push(key.slice_from(0)); }
      }
   }
   res
  }
}

#[test]
fn test_keys() {
  let mut pp=Params::new();
  pp.insert("alpha", Float(1.0)); pp.insert("beta", Float(1.0));

  let lst=pp.keys();
  assert!((lst==vec!(&"alpha", &"beta"))||(lst==vec!(&"beta", &"alpha")));
}

#[test]
fn test_varying_keys() {
 let mut pp=Params::new();
 pp.insert("alpha", Range(0.0, 1.0)); pp.insert("beta", Float(1.0));

 let lst=pp.varying_keys();
 assert_eq!(lst.len(), 1u);
 assert_eq!(lst, vec!(&"alpha"));
}

#[test]
fn test_realize() {
  let mut pp=Params::new();
  pp.insert("alpha", Float(1.0)); pp.insert("beta", Range(0.0,1.0));
  let ppr=pp.realize();

  assert_eq!(ppr.len(), 2u);
  assert_eq!(ppr.get_float("alpha"), 1.0);
  assert!( (ppr.get_float("beta")>=0.0) && (ppr.get_float("beta")<=1.0));
}
