# pbrt-r3 Performance Optimization Notes (room-path-low)

このメモは、`room-path-low.pbrt` を使った `pbrt-r3` の性能調査・最適化の記録です。
主に `perf` + flamegraph でボトルネックを追跡し、局所最適化を段階的に適用しました。

## 計測条件

- Scene: `pbrt-v3-scenes/simple/room-path-low.pbrt`
- Threads: `--nthreads 24`
- `--quick` は未使用
- 基本計測: `perf record -F 99 -m 1024 --call-graph fp`（必要に応じて `dwarf`）
- 可視化: `FlameGraph` (`perf script` -> `stackcollapse-perf.pl` -> `flamegraph.pl`)

## 主要な観測

- `pbrt-r3` は `pbrt-v3` 比較で遅く、特にメモリランタイム（`memmove`/`free`）が目立つ。
- メモリランタイム寄与は概ね以下の傾向:
- `memcpy/memmove` が最大
- `free` が次点
- `malloc` は比較的小さい

## 実施した最適化

1. QBVH traversal の `nodes_to_visit` を固定長スタック化
- `Vec` 由来の push/pop と再配置コストを削減。

2. Light sampling surface fast path を追加
- `sample_li_surface` / `pdf_li_surface` を経由し、不要な変換を減らす。

3. `FilmTile::add_sample_filter` の一時確保削減
- 毎サンプルの scratch 生成を抑制。

4. `VisibilityTester` の軽量化
- `BaseInteraction` 丸ごとではなく、`unoccluded` に必要な endpoint データ中心に整理。

5. `material_bump` の `SurfaceInteraction::clone()` 廃止
- bump 評価用に必要フィールドのみを持つ軽量 `SurfaceInteraction` を構築。
- `compute_scattering_functions` 経路の copy/drop を削減。

6. `MixMaterial` の `si.clone()` 廃止
- `m1` 計算結果を保持しつつ、`m2` 実行前に必要フィールドのみリセットして再利用。

7. `SeparableBSSRDF::sample_sp` の clone 削減
- `base = si.clone()` を必要フィールド更新に置換。
- `chain[selected].clone()` を `swap_remove` へ変更。

8. `PathIntegrator::li` の所有移動削減
- `scene.intersect` 結果を `as_ref/as_mut` で扱い、`unwrap` move を減らす。
- 次 ray を `next_ray` で保持し、分岐後に確定代入。

9. `BSDF::num_components` の O(1) 化
- `BSDF` にフラグ組み合わせカウンタを追加し、線形走査を削減。
- `li` の非スペキュラ判定を O(1) メソッドへ置換。

## 結果サマリ（代表）

- `material_bump` clone 削減後:
- 実行時間で改善（約数%）
- `SurfaceInteraction::clone` / `drop_in_place<SurfaceInteraction>` の寄与低下

- `PathIntegrator::li` 所有移動削減後:
- `memmove` 寄与が有意に減少
- 全体時間は同等〜微改善（試行差あり）

- `BSDF::num_components` O(1) 化後（3回測定）:
- 実行時間中央値: 約 `11.66s`
- `memcpy/memmove` self 中央値: 約 `4.66%`
- `free` self 中央値: 約 `3.03%`
- `malloc` self 中央値: 約 `0.38%`

## 学び

- `MemoryArena` は `malloc/free` 抑制には効くが、`memmove` には直接効きにくい。
- `memmove` を減らすには「確保」より「コピー/再配置/所有移動」の削減が重要。
- Rust では全面的な out-parameter 化より、hot path に限定した in-place 化が実用的。

## 次の候補

- `SpotLight::sample_li_surface` と `VisibilityTester` 間のデータ受け渡しをさらに簡素化。
- `SurfaceInteraction` の重いフィールド分離（幾何本体と可変付帯情報の分離）を検討。
- 変更評価は 1 回ではなく 3 回以上の中央値で判断する。

