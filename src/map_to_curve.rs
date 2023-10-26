use ff::{Field, PrimeField};
use crypto_bigint::{U256, Encoding, CheckedSub};

#[allow(non_snake_case)]
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct CurveParams<F: Field> {
    q: U256,
    A: F,
    B: F, 
    Z: F,
    Zu: F,
    Zv: F,
}

// pub fn mod2<F: Field + PrimeField<Repr = [u8;32]>>(n: F) -> F {
//     let in_int = U256::from_le_bytes(n.to_repr());
//     let out_int = in_int.checked_rem(&U256::from(2u64));
//     assert!(out_int.is_some().unwrap_u8() == 1u8);
//     let out_int = out_int.unwrap();
//     let out = F::from_repr(U256::to_le_bytes(&out_int));
//     assert!(out.is_some().unwrap_u8() == 1u8);
//     out.unwrap()
// }

// pub fn direct_map<F: Field + PrimeField<Repr = [u8;32]>>(r: F, params: CurveParams<F>) -> (F, F) {
//     let mut u = r.square();
//     let mut t1 = u * params.Z;
//     let mut v = t1 + F::ONE;
//     let t2 = v.square();
//     let mut t3 = params.A.square();
//     t3 = t3 * t1;
//     t3 = t3 - t2;
//     t3 = t3 * params.A;
//     t1 = t2 * v;

//     println!("v is {:?}", v);
//     println!("t1 is {:?}", t1);

//     let (s, mut t1) = F::sqrt_ratio(&F::ONE, &(t3 * t1));
//     let s: bool = s.into();
//     u = u * params.Zu;
//     v = v * params.Zv;

//     if s {
//         u = F::ONE;
//         v = F::ONE;
//     }

//     v = v * t3;
//     v = v * t1;
//     t1 = t1.square();

//     u  = u * (params.A).neg();
//     u  = u * t3;
//     u  = u * t2;
//     u  = u * t1;
//     t1 = -v;

    
//     let is_neg = {
//         let dbl = v.double();
//         let tmp = mod2(dbl);
//         tmp == F::ONE
//     };

//     if s ^ is_neg {
//         v = t1;
//     }

//     (u, v)

// }

pub fn legendre<F: Field + PrimeField>(f: F) -> i8 {
    if f == F::ZERO {
        0
    } 
    else {
        let (c, _) = F::sqrt_ratio(&f, &F::ONE);
        if c.unwrap_u8() == 1u8 {
            1
        } else {
            -1
        }
    }
}

pub fn direct_map<F: Field + PrimeField<Repr = [u8;32]>>(r: F, params: CurveParams<F>) -> (F, F) {
    let num = params.A.neg();
    let den = (F::ONE + params.Z * r.square()).invert();
    assert!(den.is_some().unwrap_u8() == 1u8);
    let den = den.unwrap();
    let w = num * den;

    let f = w.square() * w + params.A * w.square() + params.B * w;
    let e = legendre(f);
    if e == 0 {
        let inv_2 = F::from(2u64).invert();
        assert!(inv_2.is_some().unwrap_u8() == 1u8);
        let inv_2 = inv_2.unwrap();
        let u = params.A.neg() * inv_2;
        let v = F::ZERO;
        (u, v)
    } else if e == 1 {
        let u = w;
        let sq = u.square() * u + params.A * u.square() + params.B * u;
        let num = sq.sqrt();
        assert!(num.is_some().unwrap_u8() == 1u8);
        let num = num.unwrap();
        let v = num.neg();
        (u, v)
    } else {
        let u = w.neg() + params.A.neg();
        let sq = u.square() * u + params.A * u.square() + params.B * u;
        let v = sq.sqrt();
        assert!(v.is_some().unwrap_u8() == 1u8);
        let v = v.unwrap();
        (u, v)
    }
}

pub fn legendre2<F: Field + PrimeField<Repr = [u8; 32]>>(f: F, params: CurveParams<F>) -> F {
    let char_min_1 = params.q.checked_sub(&U256::from(1u64)).unwrap();
    let exp = char_min_1.checked_div(&U256::from(2u64)).unwrap().to_words();
    f.pow(exp)
}


pub fn direct_map2<F: Field + PrimeField<Repr = [u8;32]>>(r: F, params: CurveParams<F>) -> (F, F) {
    let num = params.A.neg();
    let den = (F::ONE + params.Z * r.square()).invert();
    assert!(den.is_some().unwrap_u8() == 1u8);
    let den = den.unwrap();
    let w = num * den;

    let f = w.square() * w + params.A * w.square() + params.B * w;
    let e = legendre2(f, params);
    
    let u = e * w - (F::ONE - e) * params.A * F::TWO_INV;

    let sqrt = (u.square() * u + params.A * u.square() + params.B * u).sqrt();
    assert!(sqrt.is_some().unwrap_u8() == 1u8);
    let sqrt = sqrt.unwrap();
    let v = e.neg() * sqrt;

    (u, v)
}



#[cfg(test)]
mod tests {
    use std::ops::Neg;

    use super::*;
    use bp_ed25519::field::{Fe25519 as Fp, FIELD_MODULUS};

    #[allow(non_snake_case)]
    fn gen_params() -> CurveParams<Fp> {
        let sqrt_m1 = (FIELD_MODULUS - Fp::ONE).sqrt();
        assert!(sqrt_m1.is_some().unwrap_u8() == 1u8);
        let sqrt_m1 = sqrt_m1.unwrap();
        let Z = Fp::from(2u64);
        let Zu = Z.neg() * sqrt_m1;
        let Zv = Zu.sqrt();
        assert!(Zv.is_some().unwrap_u8() == 1u8);
        let Zv = Zv.unwrap();
        CurveParams { 
            q: U256::from_be_hex("7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed"),
            A: Fp::from(486662u64),
            B: Fp::ONE, 
            Z: Fp::from(2u64), 
            Zu: Zu, 
            Zv: Zv 
        }
    }

    #[test]
    fn test_direct_map() {
        let params = gen_params();
        // let r0 = Fp::ZERO;

        let r1_bytes: [u8; 32] = hex::decode("990b30e04e1c3620b4162b91a33429bddb9f1b70f1da6e5f76385ed3f98ab131")
            .expect("fail to parse")
            .try_into()
            .unwrap()
        ;
        let r1 = Fp::from_repr(r1_bytes).unwrap();
        println!("r1 is {:?}", r1);

        // let r1_int = U256::from_le_hex("673a505e107189ee54ca93310ac42e4545e9e59050aaac6f8b5f64295c8ec02f");
        // println!("r1 int is {:?}", r1_int);
        let (u, v) = direct_map(r1, params);

        println!("u is {:?}", u);
        println!("v is {:?}", v);

        println!("==============================");
        let (u, v) = direct_map2(r1, params);

        println!("u is {:?}", u);
        println!("v is {:?}", v);
    }

    #[test]
    fn get_leg() {
        let a = Fp::from(486662u64).neg();
        let (e, _) = Fp::sqrt_ratio(&a, &Fp::ONE);
        println!("e is {:?}", e);
    }

}