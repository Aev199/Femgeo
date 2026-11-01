//! Побитовый отпечаток подбора Капра–Мори (ПС1 и ПС2) на 60 строках усилий.
//!
//! Защищает от незаметного изменения результата при оптимизации. Меняется только намеренно:
//! при переносе изменений OpenCS (обновление `upstream/opencs.lock`) — с объяснением в коммите.
//! Значение зависит от реализации cos/sin платформы, поэтому проверяется на Linux x86_64.

use femgeo::sp63::capra_maury::{self, PlateForces, PlateSection, SlsContext, SlsOptions, UlsContext};
use femgeo::sp63::materials::{CalcType, ConcreteClass, ConcreteDiagramType, RebarClass};

fn fingerprint() -> u64 {
    let (cu, ru) = (ConcreteClass::B25.chars(CalcType::Uls), RebarClass::A500.chars(CalcType::Uls));
    let (cs, rs) = (ConcreteClass::B25.chars(CalcType::Sls), RebarClass::A500.chars(CalcType::Sls));
    let sec = PlateSection { h: 0.2, a_top: 0.03, a_bot: 0.035, ds_x: 0.012, ds_y: 0.016, concrete_diagram: ConcreteDiagramType::L3 };
    let uls = UlsContext::new(&sec, &cu, &ru);
    let sls = SlsContext::new(&sec, &cs, &rs, Some((&cu, &ru)), SlsOptions::default());
    let mut h: u64 = 0xcbf29ce484222325;
    let mut mix = |x: f64| {
        h ^= x.to_bits();
        h = h.wrapping_mul(0x100000001b3);
    };
    for i in 0..60 {
        let i = i as f64;
        let f = PlateForces { nx: -400.0 + 13.0 * i, ny: 120.0 - 7.0 * i, nxy: 30.0 - i, mx: -40.0 + 1.7 * i, my: 25.0 - 0.9 * i, mxy: 6.0 - 0.3 * i };
        for d in [capra_maury::design_uls(&uls, &f, 10.0), capra_maury::design_sls(&sls, &f, 10.0)] {
            for x in [d.top.asx, d.top.asy, d.bottom.asx, d.bottom.asy] {
                mix(x);
            }
            for fc in &d.facets {
                mix(fc.as_top);
                mix(fc.as_bot);
                mix(fc.eta);
                mix(fc.acrc_mm);
            }
        }
    }
    h
}

#[test]
#[cfg(all(target_os = "linux", target_arch = "x86_64"))]
fn design_fingerprint_is_unchanged() {
    assert_eq!(format!("{:016x}", fingerprint()), "45291358643b36ac");
}

/// Параллельный пакетный подбор совпадает с последовательным побитово.
#[test]
fn batch_equals_sequential() {
    let (c, r) = (ConcreteClass::B30.chars(CalcType::Uls), RebarClass::A500.chars(CalcType::Uls));
    let sec = PlateSection { h: 0.3, a_top: 0.04, a_bot: 0.05, ds_x: 0.016, ds_y: 0.016, concrete_diagram: ConcreteDiagramType::L3 };
    let ctx = UlsContext::new(&sec, &c, &r);
    let rows: Vec<PlateForces> = (0..97)
        .map(|i| {
            let i = i as f64;
            PlateForces { nx: -300.0 + 7.0 * i, ny: 40.0, nxy: 10.0 - 0.2 * i, mx: 60.0 - 1.3 * i, my: -20.0 + 0.5 * i, mxy: 5.0 }
        })
        .collect();
    let par = capra_maury::design_uls_batch(&ctx, &rows, 10.0);
    for (f, d) in rows.iter().zip(&par) {
        let s = capra_maury::design_uls(&ctx, f, 10.0);
        assert_eq!(s.total().to_bits(), d.total().to_bits());
    }
}
