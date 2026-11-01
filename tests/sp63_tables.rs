//! Табличные значения материалов, сверенные с текстом СП 63.13330.2018 (с Изменением № 1)
//! 2026-11-01: табл. 6.7, 6.8, 6.11, 6.13, 6.14, 6.15.

use femgeo::sp63::materials::{CalcType, ConcreteClass, RebarClass};

#[test]
fn concrete_tables() {
    for (c, rb, rbt, rbn, rbtn, eb) in [
        (ConcreteClass::B15, 8.5, 0.75, 11.0, 1.10, 24_000.0),
        (ConcreteClass::B25, 14.5, 1.05, 18.5, 1.55, 30_000.0),
        (ConcreteClass::B30, 17.0, 1.15, 22.0, 1.75, 32_500.0),
        (ConcreteClass::B40, 22.0, 1.40, 29.0, 2.10, 36_000.0),
        (ConcreteClass::B60, 33.0, 1.80, 43.0, 2.75, 39_500.0),
    ] {
        let (u, s) = (c.chars(CalcType::Uls), c.chars(CalcType::Sls));
        assert_eq!((u.rb, u.rbt, u.eb), (rb * 1e3, rbt * 1e3, eb * 1e3), "{c:?}");
        assert_eq!((s.rb, s.rbt), (rbn * 1e3, rbtn * 1e3), "{c:?}");
    }
}

#[test]
fn rebar_tables() {
    for (r, rs, rsc, rsn, rsw) in [
        (RebarClass::A240, 210.0, 210.0, 240.0, 170.0),
        (RebarClass::A400, 340.0, 340.0, 390.0, 280.0),
        (RebarClass::A500, 435.0, 400.0, 500.0, 300.0),
    ] {
        let (u, s) = (r.chars(CalcType::Uls), r.chars(CalcType::Sls));
        assert_eq!((u.rs, u.rsc, s.rs, r.rsw()), (rs * 1e3, rsc * 1e3, rsn * 1e3, rsw * 1e3), "{r:?}");
        assert_eq!(u.es, 2.0e8);
    }
}
