//! Расчёты железобетонных сечений по СП 63.13330.
//!
//! Перенесено на Rust из OpenCS / CScore (<https://github.com/rumata-ap/OpenCS>, GPL-3.0),
//! коммит `387e15766278d30a01c4f501becd8e8273dd6ac8`. Источник указан в шапке каждого модуля.
//!
//! Единицы — кПа, м, кН/м, кН·м/м. Знаки пластины: «+» N — растяжение, «+» M растягивает
//! грань +Z (верх).

pub mod capra_maury;
pub mod diagram;
pub mod materials;
pub mod mesh_lp;
pub mod strip;
pub mod strip_design;
