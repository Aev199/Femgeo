//! Согласованность `upstream/opencs.lock`: каждый перенесённый файл указан в карте, существует и
//! ссылается в шапке на закреплённый коммит OpenCS.

use std::collections::BTreeSet;
use std::path::Path;

fn lock() -> (String, Vec<(String, Vec<String>)>) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(root.join("upstream/opencs.lock")).expect("upstream/opencs.lock");
    let mut commit = String::new();
    let mut map = Vec::new();
    for line in text.lines().map(str::trim).filter(|l| !l.is_empty() && !l.starts_with('#')) {
        let (k, v) = line.split_once(" = ").expect("строка вида «ключ = значение»");
        match k {
            "repository" => {}
            "commit" => commit = v.to_string(),
            _ => map.push((k.to_string(), v.split(',').map(|s| s.trim().to_string()).collect())),
        }
    }
    (commit, map)
}

#[test]
fn every_ported_file_is_mapped_and_cites_pinned_commit() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let (commit, map) = lock();
    assert_eq!(commit.len(), 40, "полный SHA");
    let short = &commit[..7];
    let mapped: BTreeSet<_> = map.iter().map(|(k, _)| k.clone()).collect();
    for (ours, theirs) in &map {
        assert!(!theirs.is_empty(), "{ours}: нет файлов OpenCS");
        let text = std::fs::read_to_string(root.join(ours)).unwrap_or_else(|_| panic!("{ours} не найден"));
        let head: String = text.lines().take(12).collect::<Vec<_>>().join("\n");
        assert!(head.contains(short), "{ours}: в шапке нет коммита {short}");
    }
    // Все файлы расчётов сечений перенесены из OpenCS; в `fem` — только указанные в карте.
    for entry in std::fs::read_dir(root.join("src/sp63")).unwrap() {
        let p = entry.unwrap().path();
        let rel = format!("src/sp63/{}", p.file_name().unwrap().to_string_lossy());
        if rel.ends_with("mod.rs") {
            continue;
        }
        assert!(mapped.contains(&rel), "{rel} не указан в upstream/opencs.lock");
    }
}
