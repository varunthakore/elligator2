// Reference : https://elligator.org/map

use ff::{Field, PrimeField};
use crypto_bigint::{U256, CheckedSub};

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

pub fn legendre<F: Field + PrimeField<Repr = [u8; 32]>>(f: F, params: CurveParams<F>) -> F {
    let char_min_1 = params.q.checked_sub(&U256::from(1u64)).unwrap();
    let exp = char_min_1.checked_div(&U256::from(2u64)).unwrap().to_words();
    f.pow(exp)
}

pub fn direct_map<F: Field + PrimeField<Repr = [u8;32]>>(r: F, params: CurveParams<F>) -> (F, F) {
    let num = params.A.neg();
    let den = (F::ONE + params.Z * r.square()).invert();
    assert!(den.is_some().unwrap_u8() == 1u8);
    let den = den.unwrap();
    let w = num * den;

    let f = w.square() * w + params.A * w.square() + params.B * w;
    let e = legendre(f, params);
    
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

    // Generate curve parameters for Curve25519. Reference : https://elligator.org/formulas
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

    // Test Vectors copied from https://elligator.org/vectors/
    #[test]
    fn test_direct_map() {
        let params = gen_params();

        let r_bytes: [u8; 32] = hex::decode("66665895c5bc6e44ba8d65fd9307092e3244bf2c18877832bd568cb3a2d38a12")
            .expect("fail to parse")
            .try_into()
            .unwrap()
        ;
        let r = Fp::from_repr(r_bytes).unwrap();
        
        // Expected Output
        let u_bytes: [u8; 32] = hex::decode("04d44290d13100b2c25290c9343d70c12ed4813487a07ac1176daa5925e7975e")
            .expect("fail to parse")
            .try_into()
            .unwrap()
        ;
        let v_bytes: [u8; 32] = hex::decode("c35aa4226513c49a3c12b48d47f7e176e64c122345cad87c4e3ec9a72c828900")
            .expect("fail to parse")
            .try_into()
            .unwrap()
        ;
        let u_exp = Fp::from_repr(u_bytes).unwrap();
        let v_exp = Fp::from_repr(v_bytes).unwrap();
        

        let (u, v) = direct_map(r, params);
        assert_eq!(u, u_exp);
        assert_eq!(v, v_exp);
    }

}