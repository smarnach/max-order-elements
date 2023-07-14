use std::collections::BTreeMap;

use either::Either;
use num_integer::{gcd, lcm};
use num_modular::{ModularInteger, MontgomeryInt, VanillaInt};
use num_prime::nt_funcs::factorize64;
use num_traits::Pow;

const M_MIN: u64 = 100_000_001;
const M_MAX: u64 = 110_000_000;
const M_STEP: u64 = 2;
const M_COUNT: u64 = (M_MAX - M_MIN + M_STEP - 1) / M_STEP;

fn main() {
    let phi: f64 = (5.0f64.sqrt() + 1.0) / 2.0;
    let mut sum = 0;
    let mut max = 0;
    for m in (M_MIN..M_MAX).step_by(M_STEP as _) {
        let ring = ModularRing::new(m);
        let (_, count) = ring.max_order_element_close_to(m as f64 / phi);
        sum += count;
        max = max.max(count);
    }
    println!("Average numbers checked: {}", sum as f64 / M_COUNT as f64);
    println!("Maximum numbers checked: {}", max);
}

#[derive(Debug)]
pub struct ModularRing {
    m: u64,
    exponent: u64,
    exp_factors: BTreeMap<u64, usize>,
    one: Either<VanillaInt<u64>, MontgomeryInt<u64>>,
}

impl ModularRing {
    pub fn new(m: u64) -> Self {
        let factors = factorize64(m);
        let exponent = factors
            .iter()
            .map(|(&p, &k)| match p {
                2 if k > 2 => 2u64.pow((k - 2) as _),
                _ => (p - 1) * p.pow((k - 1) as _),
            })
            .fold(1, lcm);
        let exp_factors = factorize64(exponent);
        let one = match m & 1 {
            0 => Either::Left(VanillaInt::new(1, &m)),
            1 => Either::Right(MontgomeryInt::new(1, &m)),
            _ => unreachable!(),
        };
        Self {
            m,
            exponent,
            exp_factors,
            one,
        }
    }

    pub fn has_max_order(&self, a: u64) -> bool {
        if gcd(a, self.m) != 1 {
            return false;
        }

        #[inline]
        fn helper<T>(this: &ModularRing, a: u64, one: &T) -> bool
        where
            T: ModularInteger<Base = u64> + Pow<u64, Output = T> + Copy,
        {
            let a = one.convert(a);
            this.exp_factors
                .keys()
                .all(|&p| a.pow(this.exponent / p) != *one)
        }

        match self.one {
            Either::Left(ref one) => helper(self, a, one),
            Either::Right(ref one) => helper(self, a, one),
        }
    }

    pub fn max_order_element_close_to(&self, x: f64) -> (u64, usize) {
        let mut a = x.round() as i64;
        let mut delta = if (a as f64) < x { 1 } else { -1 };
        let mut count = 1;
        loop {
            if self.has_max_order(a as _) {
                return (a as _, count);
            }
            a += delta;
            delta = -delta - delta.signum();
            count += 1;
        }
    }
}
