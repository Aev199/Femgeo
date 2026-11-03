//! Femgeo: подбор расстановки свай и армирования конструкций.
//!
//! Этап 1 — расчёты сечений по СП 63.13330 ([`sp63`]); этап 2 — линейный МКЭ ([`fem`]).

pub mod fem;
pub mod aster;
pub mod building;
pub mod geobase;
pub mod loads;
pub mod optimize;
pub mod piles;
pub mod reinforce;
pub mod report;
pub mod service;
pub mod sp63;
