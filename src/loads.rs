//! Сочетания нагрузок по СП 20.13330.2016, п. 6 (основные сочетания).
//!
//! Коэффициенты сочетаний (п. 6.4): длительные ψl1 = 1,0, остальные 0,95; кратковременные
//! ψt1 = 1,0, ψt2 = 0,9, остальные 0,7 — перебираются все варианты ведущей длительной и
//! ведущей и второй кратковременной нагрузки. ПС1 — расчётные значения (×γf), ПС2 —
//! нормативные; длительное сочетание ПС2 — постоянные, длительные и длительные части
//! кратковременных (ψдл). Особые сочетания не формируются.

/// Вид нагрузки.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Permanent,
    Long,
    Short,
}

/// Загружение: вид, γf, длительная доля (для кратковременных).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Case {
    pub kind: Kind,
    pub gamma_f: f64,
    pub psi_long: f64,
}

/// Сочетание: коэффициенты к нормативным значениям загружений (по номеру загружения).
pub type Combination = Vec<f64>;

/// Сочетания: ПС1, ПС2 полные, ПС2 длительное.
#[derive(Clone, Debug, PartialEq)]
pub struct Combinations {
    pub uls: Vec<Combination>,
    pub sls: Vec<Combination>,
    pub sls_long: Combination,
}

fn main_combinations(cases: &[Case], design: bool) -> Vec<Combination> {
    let idx = |k: Kind| -> Vec<usize> { cases.iter().enumerate().filter(|(_, c)| c.kind == k).map(|(i, _)| i).collect() };
    let (perm, long, short) = (idx(Kind::Permanent), idx(Kind::Long), idx(Kind::Short));
    let g = |i: usize| if design { cases[i].gamma_f } else { 1.0 };
    let leads_l: Vec<Option<usize>> = if long.is_empty() { vec![None] } else { long.iter().map(|&i| Some(i)).collect() };
    let mut pairs: Vec<(Option<usize>, Option<usize>)> = vec![];
    match short.len() {
        0 => pairs.push((None, None)),
        1 => pairs.push((Some(short[0]), None)),
        _ => {
            for &a in &short {
                for &b in &short {
                    if a != b {
                        pairs.push((Some(a), Some(b)));
                    }
                }
            }
        }
    }
    let mut out = vec![];
    for lead_l in &leads_l {
        for &(t1, t2) in &pairs {
            let mut c = vec![0.0; cases.len()];
            for &i in &perm {
                c[i] = g(i);
            }
            for &i in &long {
                c[i] = g(i) * if Some(i) == *lead_l { 1.0 } else { 0.95 };
            }
            for &i in &short {
                let psi = if Some(i) == t1 { 1.0 } else if Some(i) == t2 { 0.9 } else { 0.7 };
                c[i] = g(i) * psi;
            }
            if !out.contains(&c) {
                out.push(c);
            }
        }
    }
    out
}

pub fn combinations(cases: &[Case]) -> Combinations {
    let sls_long = cases
        .iter()
        .map(|c| match c.kind {
            Kind::Permanent | Kind::Long => 1.0,
            Kind::Short => c.psi_long,
        })
        .collect();
    Combinations { uls: main_combinations(cases, true), sls: main_combinations(cases, false), sls_long }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn case(kind: Kind, g: f64, psi: f64) -> Case {
        Case { kind, gamma_f: g, psi_long: psi }
    }

    #[test]
    fn coefficients_of_main_combinations() {
        let cs = [case(Kind::Permanent, 1.1, 0.0), case(Kind::Long, 1.2, 0.0), case(Kind::Long, 1.3, 0.0), case(Kind::Short, 1.2, 0.35), case(Kind::Short, 1.4, 0.0), case(Kind::Short, 1.3, 0.5)];
        let c = combinations(&cs);
        // 2 ведущих длительных × 3·2 пар кратковременных.
        assert_eq!(c.uls.len(), 12);
        assert!(c.uls.iter().all(|x| x[0] == 1.1));
        assert!(c.uls.contains(&vec![1.1, 1.2, 1.3 * 0.95, 1.2, 1.4 * 0.9, 1.3 * 0.7]));
        assert!(c.sls.contains(&vec![1.0, 0.95, 1.0, 0.7, 1.0, 0.9]));
        assert_eq!(c.sls_long, vec![1.0, 1.0, 1.0, 0.35, 0.0, 0.5]);
    }

    #[test]
    fn single_short_and_no_variable() {
        let c = combinations(&[case(Kind::Permanent, 1.1, 0.0), case(Kind::Short, 1.3, 0.3)]);
        assert_eq!(c.uls, vec![vec![1.1, 1.3]]);
        let c = combinations(&[case(Kind::Permanent, 1.05, 0.0)]);
        assert_eq!((c.uls.len(), c.sls.len()), (1, 1));
    }
}
