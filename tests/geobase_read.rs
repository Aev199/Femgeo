//! Чтение объекта geobase на синтетической базе той же схемы (`create_db.py`).

use femgeo::geobase::Object;
use rusqlite::Connection;

fn db(dir: &std::path::Path) -> std::path::PathBuf {
    let p = dir.join("obj.db");
    let c = Connection::open(&p).unwrap();
    c.execute_batch(
        "CREATE TABLE Soils (id INTEGER PRIMARY KEY, name TEXT NOT NULL, number TEXT, genesis TEXT, description TEXT, soil_type TEXT, notes TEXT, soil_category TEXT);
         CREATE TABLE ConstitutiveModels (id INTEGER PRIMARY KEY, name TEXT NOT NULL UNIQUE, software_compatibility TEXT);
         CREATE TABLE SoilProperties (id INTEGER PRIMARY KEY, soil_id INTEGER, model_id INTEGER, property_name TEXT, property_value REAL, unit TEXT, notes TEXT);
         CREATE TABLE boreholes (id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT NOT NULL, x_coord REAL, y_coord REAL, z_coord REAL);
         CREATE TABLE borehole_layers (id INTEGER PRIMARY KEY AUTOINCREMENT, borehole_id INTEGER, layer_number INTEGER, soil_id INTEGER, depth REAL, layer_name TEXT);
         INSERT INTO Soils (id, name) VALUES (1, 'Суглинок'), (2, 'Песок');
         INSERT INTO ConstitutiveModels (id, name) VALUES (1, 'Mohr-Coulomb'), (2, 'Hardening Soil');
         INSERT INTO SoilProperties (soil_id, model_id, property_name, property_value) VALUES
           (1, 1, 'Е', 12000), (1, 1, 'ν', 0.35), (2, 1, 'E', 30000), (2, 1, 'ν', 0.3), (1, 2, 'E50', 9000);
         INSERT INTO boreholes (id, name, x_coord, y_coord, z_coord) VALUES (1, 'Скв-1', 0, 0, 150.0), (2, 'Скв-2', 100, 0, 148.0);
         INSERT INTO borehole_layers (borehole_id, layer_number, soil_id, depth) VALUES
           (1, 1, 1, 0.0), (1, 2, 1, 6.0), (1, 3, 2, 30.0), (2, 1, 1, 4.0), (2, 2, 2, 25.0);",
    )
    .unwrap();
    p
}

#[test]
fn reads_columns_by_nearest_borehole() {
    let dir = std::env::temp_dir().join(format!("femgeo_gb_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = db(&dir);
    let o = Object::open(&path, None).unwrap();
    assert_eq!(o.model, "Mohr-Coulomb");
    let c = o.column_at(10.0, 5.0).unwrap();
    assert_eq!(c.top, 150.0);
    // Слой нулевой толщины пропущен; подошвы — отметки устья минус глубина.
    assert_eq!(c.layers.len(), 2);
    assert_eq!((c.layers[0].bottom, c.layers[0].e, c.layers[0].nu), (144.0, 12000.0, 0.35));
    assert_eq!((c.layers[1].bottom, c.layers[1].e), (120.0, 30000.0));
    let c2 = o.column_at(80.0, 0.0).unwrap();
    assert_eq!(c2.top, 148.0);
    // Модель без ν — понятная ошибка.
    let hs = Object::open(&path, Some("Hardening Soil")).unwrap();
    assert!(hs.column_at(0.0, 0.0).unwrap_err().contains("не заданы"));
    std::fs::remove_dir_all(&dir).ok();
}
