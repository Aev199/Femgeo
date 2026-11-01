// Оболочка приложения: одна команда `call` передаётся в `Service::dispatch`, как и в
// отладочном мосте. Долгие операции (расчёт) идут вне потока окна и сообщают этапы событием
// `progress`.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use femgeo::service::Service;
use std::sync::Mutex;
use tauri::{Emitter, Manager};

struct App(Mutex<Service>);

#[tauri::command]
async fn call(app: tauri::AppHandle, command: String, args: serde_json::Value) -> Result<serde_json::Value, String> {
    tauri::async_runtime::spawn_blocking(move || {
        let state = app.state::<App>();
        let mut service = state.0.lock().map_err(|_| "состояние приложения недоступно".to_string())?;
        service.dispatch(&command, args, &mut |stage| {
            let _ = app.emit("progress", stage);
        })
    })
    .await
    .map_err(|e| e.to_string())?
}

/// Переносная программа: данные WebView2 — в папке рядом с exe, если в неё можно писать,
/// иначе в локальных данных пользователя.
fn data_dir(app: &tauri::App) -> std::path::PathBuf {
    let beside = std::env::current_exe().ok().and_then(|exe| exe.parent().map(|d| d.join("femgeo-data")));
    if let Some(dir) = beside {
        if std::fs::create_dir_all(&dir).is_ok() && std::fs::write(dir.join(".write-test"), b"").is_ok() {
            let _ = std::fs::remove_file(dir.join(".write-test"));
            return dir;
        }
    }
    app.path().app_local_data_dir().unwrap_or_else(|_| std::env::temp_dir().join("femgeo"))
}

fn main() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .setup(|app| {
            let data = data_dir(app);
            app.manage(App(Mutex::new(Service::new())));
            tauri::WebviewWindowBuilder::new(app, "main", tauri::WebviewUrl::default())
                .title("Femgeo")
                .inner_size(1400., 900.)
                .min_inner_size(1000., 640.)
                .theme(Some(tauri::Theme::Light))
                .data_directory(data.join("webview"))
                .build()?;
            Ok(())
        })
        .invoke_handler(tauri::generate_handler![call])
        .run(tauri::generate_context!())
        .expect("ошибка запуска приложения");
}
