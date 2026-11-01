// Доступ к ядру: команда Tauri `call` в приложении, отладочный мост
// (`examples/ui_server.rs`) по HTTP при проверке интерфейса. Обе ведут в `Service::dispatch`.
import { invoke } from '@tauri-apps/api/core';
import { listen } from '@tauri-apps/api/event';
import { open, save } from '@tauri-apps/plugin-dialog';

export const desktop = typeof window !== 'undefined' && '__TAURI_INTERNALS__' in window;

export async function call(command, args = {}) {
  if (desktop) return invoke('call', { command, args });
  const response = await fetch(`./api/${command}`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(args),
  });
  const body = await response.json();
  if ('error' in body) throw new Error(body.error);
  return body.ok;
}

export async function onProgress(handler) {
  if (desktop) await listen('progress', (event) => handler(event.payload));
}

const FILTERS = {
  geobase: [{ name: 'Объект GeoBase', extensions: ['db'] }],
  project: [{ name: 'Проект Femgeo', extensions: ['femgeo.json', 'json'] }],
  csv: [{ name: 'Таблица CSV', extensions: ['csv'] }],
  dxf: [{ name: 'Чертёж DXF', extensions: ['dxf'] }],
  html: [{ name: 'Отчёт HTML', extensions: ['html'] }],
};

/// Путь к файлу: диалог в приложении, запрос пути при проверке через мост.
export async function pickFile(kind, saving = false) {
  if (desktop) {
    const options = { filters: FILTERS[kind] };
    const path = saving ? await save(options) : await open({ ...options, multiple: false });
    return path || null;
  }
  return window.prompt(saving ? 'Путь для сохранения' : 'Путь к файлу') || null;
}
