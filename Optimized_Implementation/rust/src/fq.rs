use std::ops::{Add, Mul, Sub};


#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Fq(pub u8);

const FQ127_INV_TABLE: [u8; 128] = [0, 1, 64, 85, 32, 51, 106, 109, 16, 113, 89, 104, 53, 88, 118, 17, 8, 15, 120, 107, 108, 121, 52, 116, 90, 61, 44, 80, 59, 92, 72, 41, 4, 77, 71, 98, 60, 103, 117, 114, 54, 31, 124, 65, 26, 48, 58, 100, 45, 70, 94, 5, 22, 12, 40, 97, 93, 78, 46, 28, 36, 25, 84, 125, 2, 43, 102, 91, 99, 81, 49, 34, 30, 87, 115, 105, 122, 33, 57, 82, 27, 69, 79, 101, 62, 3, 96, 73, 13, 10, 24, 67, 29, 56, 50, 123, 86, 55, 35, 68, 47, 83, 66, 37, 11, 75, 6, 19, 20, 7, 112, 119, 110, 9, 39, 74, 23, 38, 14, 111, 18, 21, 76, 95, 42, 63, 126, 0];


impl Fq {
    pub const Q: u8 = 127;
    pub const Q_BITS:u8 = 7;

    /// a mod q, a < 256
    fn conditional_sub_raw(a: u8) -> u8 {
        let sub_q = a - Fq::Q;
        let mask = (-((sub_q >> Fq::Q_BITS) as i8)) as u8;
        let ret = (mask & Fq::Q) + sub_q;
        return ret;
    }

    /// a mod q, a < 256
    //pub fn conditional_sub(a: Fq) -> Fq {
    //    return Fq(Fq::conditional_sub_raw(a.0));
    //}
    
    /// a mod q,  a < 256**2
    pub fn red(a: u16) -> Fq {
        return Fq(Fq::conditional_sub_raw(((a >> Fq::Q_BITS) as u8) + ((a as u8) & Fq::Q)));
    }

    pub fn add(a: Fq, b: Fq) -> Fq {
        return Fq(Fq::conditional_sub_raw(a.0 + b.0));
    }
    
    pub fn sub(a: Fq, b: Fq) -> Fq {
        return Fq(Fq::conditional_sub_raw(a.0 + Fq::Q - b.0));
    }

    pub fn mul(a: Fq, b: Fq) -> Fq {
        return Fq::red((a.0 as u16) * (b.0 as u16));
    }

    pub fn inv(a: Fq) -> Fq {
        Fq(FQ127_INV_TABLE[a.0 as usize])
    }
}

impl Add for Fq {
    type Output = Fq;
    fn add(self, r: Fq) -> Fq {
        Fq::add(self, r)
    }
}

impl Sub for Fq {
    type Output = Fq;
    fn sub(self, r: Fq) -> Fq {
        Fq::sub(self, r)
    }
}

impl Mul for Fq {
    type Output = Fq;
    fn mul(self, r: Fq) -> Fq {
        Fq::mul(self, r)
    }
}
