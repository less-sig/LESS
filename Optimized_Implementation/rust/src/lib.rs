#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Fq(u8);

impl Fq {
    pub const Q: u8 = 127;
    pub const Q_BITS:u8 = 7;

    fn conditional_sub_raw(a: u8) -> u8 {
        let sub_q = a - Fq::Q;
        let mask = (-((sub_q >> Fq::Q_BITS) as i8)) as u8;
        let ret = (mask & Fq::Q) + sub_q;
        return ret;
    }

    /// du kek
    pub fn conditional_sub(a: Fq) -> Fq {
        return Fq(Fq::conditional_sub_raw(a.0));
    }

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

}

pub struct Matrix<const N: usize, const M: usize> {

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
    }
}
