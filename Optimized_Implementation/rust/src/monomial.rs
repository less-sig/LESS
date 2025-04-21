use std::ops::Mul;
use crate::fq::Fq;

#[derive(Debug)]
pub struct Monomial<const N: usize> {
    coeffs: [Fq; N],
    perms: [u16; N],
}

#[derive(Debug)]
pub struct MonomialIS<const N: usize> {
    perms: [u16; N],
}

impl<const N: usize> Monomial<N> {
    ///
    pub fn new() -> Self {
        let mut t: [u16; N] = [0u16; N];
        for i in 0..N {
            t[i] = i as u16;
        }

        Monomial {
            coeffs: [Fq(0); N],
            perms: t,
        }
    }

    /// 
    pub fn random(&mut self) {

    }

    /// l = l*r
    pub fn mul(l: &mut Self, r: &Self) {
        for i in 0..N {
            let t: usize = r.perms[i] as usize;
            l.perms[i] = l.perms[t] ;
            l.coeffs[i] = Fq::mul(l.coeffs[t], r.coeffs[i]);
        }
    }

    /// l = r**-1
    pub fn inv(l: &mut Self, r: &Self) {
        for i in 0..N {
            let t = r.perms[i] as usize;
            l.perms[t] = i as u16;
            l.coeffs[t] = Fq::inv(r.coeffs[i])
        } 
    }
}


impl<const N: usize> MonomialIS<N> {
    // const BYTES: usize = (N + 7) / 8;


    ///
    pub fn new() -> Self {
        let mut t: [u16; N] = [0u16; N];
        MonomialIS {
            perms: t,
        }
    }

    /// translate form `CheckCanonicakAction`
    pub fn from_bytes(b: [bool; N]) -> Self {
        let mut t: [u16; N] = [0u16; N];

        for i in 0..N {
            if b[i] == true {
                t[i] = i as u16;
            }
        }

        MonomialIS {
            perms: t,
        }
    }

    ///
    pub fn compress(&self, k: u32) -> [bool; N] {
        let mut t = [false; N];
        for i in 0..k as usize {
            let p: usize = self.perms[i] as usize;
            t[p] = true;
        }

        t
    }
}
