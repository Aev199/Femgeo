# CLAUDE.md — Femgeo

Начало работы в новом чате: прочитать этот файл, затем `docs/DEV_STATE.md`
(состояние, допущения, следующие шаги), при необходимости `docs/SPEC.md`
(постановка) и `AGENTS.md` (расчётные правила, перенос OpenCS, цикл разработки).
Состояние восстанавливать по репозиторию, а не по истории чата.

## Пользователь

- Инженер-конструктор/геотехник, пишет по-русски, работает с телефона: ответы
  краткие, итог пакета — что сделано, тесты, метрики, блокер, SHA.
- Бюджет ограничен — экономить токены: не перечитывать большие файлы целиком,
  править точечно, не выводить логи.
- CI запускать только при необходимости (разрешено). Exe не собирать без
  просьбы. Автору OpenCS ничего не сообщать.
- Решения, которые за пользователем: влить `stage1-sp63-plates` в `main`
  (нужно для ручной сборки exe), ветка `license-gpl3` в `topo_reconstruct_rs`.

## Устройство

```
src/sp63/       расчёты сечений СП 63 (перенос OpenCS, карта upstream/opencs.lock)
src/fem/        МКЭ: shell.rs (MITC4/T3), mod.rs (модель, сборка, faer), condense.rs
src/piles/      СП 24: осадка и взаимовлияние свай (mod.rs), правила расстановки (layout.rs)
src/loads.rs    сочетания СП 20
src/optimize/   подбор свай: Problem::regular (ряды), Problem::frontier (свободно), Вудбери
src/reinforce.rs армирование здания, изополя DXF
src/mxt.rs      здание из MIDAS .mxt (экспорт topo_reconstruct_rs)
src/geobase.rs  чтение объекта geobase (SQLite)
src/service.rs  команды окна: defaults, run, open_geobase, save/load_project, экспорты
app/            окно Tauri 2 (src-tauri) и интерфейс (ui: index.html, src/main.js, plan.js, app.css; style.css — дизайн-код GeoBase из topo)
examples/       ui_server.rs (отладочный мост), bench_*.rs, profile_sp63.rs
tests/          интеграционные тесты; common/ — фикстура OpenCS
```

Единицы: кН, м, кПа. Плиты: «+» M растягивает верх; CScore: N < 0 — сжатие.

## Команды

```sh
cargo test --release                       # все тесты (~2 мин)
cargo test --release --test optimize_raft  # выборочно
cargo run --release --example bench_service -- 48 36   # замер расчёта
cd app/ui && npm ci && npm run build       # интерфейс → app/ui/dist
cargo build --release --example ui_server  # мост: target/release/examples/ui_server 8787 app/ui/dist
```

Проверка окна — Playwright (Chromium предустановлен, `require('playwright')`),
скрипты держать вне репозитория (например `/tmp/uitest`), скриншоты смотреть.

## Особенности среды

- `pkill` в одной команде с другими обрывает всю команду (код 144) — запускать
  `pkill -f "ui_server 8787" || true` отдельным вызовом; мост перезапускать
  после каждой пересборки.
- .NET недоступен (установка запрещена) — эталоны C# только из тестов OpenCS
  или через CI.
- crates.io и npm доступны; Gmsh в контейнере нет.
- Связанные репозитории (только чтение): `Aev199/topo_reconstruct_rs`
  (геометрия ЛИРА, актуальная ветка `claude/building-geometry-recovery-vditix`),
  `Aev199/geobase`, `Aev199/cscore-pile-poc` (эталоны АРБАТ),
  `rumata-ap/OpenCS` (источник переноса).

## Правила, которые легко нарушить

- Несходимость — `indeterminate`, не успех; допуски тестов не ослаблять под
  результат — сначала понять причину (в истории проекта причиной не раз были
  исходные данные теста, а не код).
- Перенос из OpenCS — по алгоритму источника; оптимизации не меняют результат
  побитово (`tests/sp63_fingerprint.rs`); новый перенесённый файл — в
  `upstream/opencs.lock`.
- Нормы цитировать с пунктом; значения по памяти помечать и сверять
  (тексты СП 24 и СП 63 найдены: meganorm, копия КонсультантПлюс).
- Приватные модели ЛИРА и объекты geobase не коммитить.
- Коммиты — завершёнными пакетами, документация — с `[skip ci]`; в сообщениях
  коммитов не указывать модель.
