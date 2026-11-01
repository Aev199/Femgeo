//! Диаграммы σ(ε) материалов: кусочно-линейные ветви сжатия и растяжения.
//!
//! Источник: OpenCS `CScore/Diagramm.cs` (метод `Sig`) и `CSmath/Lspline.cs`
//! (интерполяция и экстраполяция), коммит `387e157`.

/// Вид материала: от него зависят правила ветвей в [`Diagram::sig`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaterialKind {
    Concrete,
    Rebar,
}

/// Кусочно-линейная ветвь σ(ε): узлы по возрастанию ε, без повторов.
///
/// Вне диапазона узлов — продолжение крайнего отрезка (как `LSpline`): для диаграмм СП 63
/// крайние отрезки — площадки, поэтому за ними напряжение постоянно.
#[derive(Clone, Debug)]
pub struct Branch {
    x: Vec<f64>,
    y: Vec<f64>,
}

impl Branch {
    /// Ветвь по точкам (ε, σ). Точки сортируются по ε; повторяющиеся ε — ошибка.
    pub fn new(points: &[(f64, f64)]) -> Result<Self, String> {
        if points.len() < 2 {
            return Err("ветвь диаграммы требует не менее двух точек".into());
        }
        let mut p = points.to_vec();
        p.sort_by(|a, b| a.0.total_cmp(&b.0));
        for w in p.windows(2) {
            if w[0].0 == w[1].0 {
                return Err("ветвь диаграммы содержит повторяющиеся деформации".into());
            }
        }
        Ok(Self {
            x: p.iter().map(|q| q.0).collect(),
            y: p.iter().map(|q| q.1).collect(),
        })
    }

    /// Узлы ветви по деформации.
    pub fn xs(&self) -> &[f64] {
        &self.x
    }

    /// Последний (наибольший) узел ветви.
    pub fn last_x(&self) -> f64 {
        *self.x.last().expect("ветвь не пустая")
    }

    fn slope(&self, i: usize) -> f64 {
        let n = self.x.len();
        let i = i.min(n - 2);
        (self.y[i + 1] - self.y[i]) / (self.x[i + 1] - self.x[i])
    }

    /// (σ, dσ/dε) в точке `e`.
    pub fn eval(&self, e: f64) -> (f64, f64) {
        let n = self.x.len();
        if e < self.x[0] {
            let b = self.slope(0);
            return (self.y[0] + b * (e - self.x[0]), b);
        }
        if e > self.x[n - 1] {
            let b = self.slope(n - 2);
            return (self.y[n - 1] + b * (e - self.x[n - 1]), b);
        }
        let mut i = 0;
        while i < n - 1 && e >= self.x[i + 1] {
            i += 1;
        }
        let b = self.slope(i);
        (self.y[i] + b * (e - self.x[i]), b)
    }
}

/// Диаграмма материала: ветвь сжатия (ε < 0) и растяжения (ε > 0).
#[derive(Clone, Debug)]
pub struct Diagram {
    pub kind: MaterialKind,
    pub compression: Branch,
    pub tension: Branch,
}

impl Diagram {
    /// Напряжение и касательный модуль при деформации `eps`.
    ///
    /// Бетон: сжатие учитывается всегда, растяжение — по `tension_concrete`; за последним
    /// узлом растяжения (εbt,ult) бетон с трещиной не работает. Арматура: растяжение всегда,
    /// сжатие — по `rebar_compression`.
    pub fn sig(&self, eps: f64, tension_concrete: bool, rebar_compression: bool) -> (f64, f64) {
        match self.kind {
            MaterialKind::Concrete => {
                if eps > 0.0 && tension_concrete {
                    if eps > self.tension.last_x() {
                        (0.0, 0.0)
                    } else {
                        self.tension.eval(eps)
                    }
                } else if eps < 0.0 {
                    self.compression.eval(eps)
                } else {
                    (0.0, 0.0)
                }
            }
            MaterialKind::Rebar => {
                if eps < 0.0 && rebar_compression {
                    self.compression.eval(eps)
                } else if eps > 0.0 {
                    self.tension.eval(eps)
                } else {
                    (0.0, 0.0)
                }
            }
        }
    }

    /// Только напряжение.
    pub fn stress(&self, eps: f64, tension_concrete: bool, rebar_compression: bool) -> f64 {
        self.sig(eps, tension_concrete, rebar_compression).0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn branch_interpolates_and_extrapolates_by_end_segments() {
        let b = Branch::new(&[(0.0, 0.0), (1.0, 2.0), (3.0, 2.0)]).unwrap();
        assert_eq!(b.eval(0.5), (1.0, 2.0));
        assert_eq!(b.eval(2.0), (2.0, 0.0));
        assert_eq!(b.eval(5.0), (2.0, 0.0));
        assert_eq!(b.eval(-1.0), (-2.0, 2.0));
        assert_eq!(b.eval(1.0), (2.0, 0.0));
    }

    #[test]
    fn branch_sorts_and_rejects_duplicates() {
        let b = Branch::new(&[(3.0, 2.0), (0.0, 0.0), (1.0, 2.0)]).unwrap();
        assert_eq!(b.xs(), &[0.0, 1.0, 3.0]);
        assert!(Branch::new(&[(0.0, 0.0), (0.0, 1.0)]).is_err());
    }
}
