# C++化・工程3: 固定T ProblemとNewtonSessionの分離

2026-09-28。固定T Newtonの反復制御をC++へ移し、物理的な入力と可変な反復状態を分離した。
driver、Model、solverは独立しており、driverがModelやsolverを所有する設計にはしていない。

## 配置と責務

| ファイル | 実装 |
|---|---|
| `model/fixed_rate_problem.hpp` / `.cc` | `chemistry::FixedRateProblem`。初期状態、固定反応係数、T、dt、環境値、物理的な境界値を保持。RHS・gross生成消滅・残差・Jacobianの評価を提供 |
| `solver/fixed_rate_solver.hpp` | statelessな `chemistry::FixedRateSolver`、`FixedRateOptions`、caller所有の `NewtonScratch` |
| `solver/fixed_rate_solver.cc` | solveごとに生成する非公開の `NewtonSession` と反復ループ。conservation chart、未知数一覧、LU有効性、clipping、収束確認、deterministic/random restartを管理 |
| `solver/chemistry_newton_internal.h` | 不変データ用 `ChemNewtonReference` と、可変データ用 `ChemNewtonSystem` を分離 |
| `solver/fixed_rate_c.h` | 移行中のC実装との内部bridge。public APIではない |
| `solver/chemistry_solver.c` | energy・explicit solver、物理・保存則のhelperを保持。既存の固定T C入口はC++の実装へ委譲 |
| `cooling/chemistry_cooling.h` / `.c` | 化学加熱用の初期組成引数を `const double *` に変更。式は不変 |

`FixedRateProblem`は環境のrun/global値もコピーするが、network/tableの実体はModelから借用する。
ModelはProblemより長生きする必要がある。Problemのcopy/moveは禁止し、内部の借用参照を安定させる。
`valid()`は既存の初期化処理の成否を表す。内部APIとして有効なrun/globalポインタと評価済みrateを要求する。

初期組成・化学加熱用初期組成・floor・上下限はProblem側で保持し、Newton側にはconst pointerで渡す。
初期推定、保存則の切り替え、再試行によって物理的なold stateを書き換えることはできない。
SessionはProblemが用意した初期推定をコピーし、solveごとに反復履歴・LU有効性を初期化する。
保存則の再選択とLU等の無効化はSessionのメソッドで結び付けた。

`NewtonScratch`は行列・ベクトル・pivot配列のみを保持し、再利用できる。solve中のheap確保は追加していない。
出力に使う`ChemWorkspace`は既存のC評価・診断・化学加熱のために残している。
solverはProblemの初期物理状態・温度・rateから開始する。反復上限で失敗した際は物理的な初期状態を返す。
初期化・反復中の既存failure処理は維持し、export時のあらゆる失敗を新たにtransactional化したわけではない。

## 維持した数値処理

反応・Jacobian・保存則・初期推定の物理計算には既存のC helperを使う。
Newtonの計算順、LAPACK、LU再利用条件、収束確認回数、residual gate、clipping hold、
再試行順序、random seedと更新順を維持した。C++17で必要なenumの明示変換を追加している。
現在のsolverは具体的なchemistry型とC helperに依存する。
汎用 `NewtonSolver<Problem, LinearBackend>`、独立したlinear backend、統一したSolveReportは未実装。
energy solverと時間進行driverも引き続きCである。

## 検証

- 既存energy-step replay 4件が成功。
  C版 `782d710` のログと比較し、実行時間 `elapsed_s` のみを除いてログ全体が完全一致。
  [比較記録](chemistry_cpp_stage3_validation_2026-09-28.json)参照。
- 新しい `check-chemistry-fixed-rate` が成功。
  入力組成・rate・環境をcaller側で書き換えてもProblemが不変であること、
  同じProblemを独立したSessionで繰り返し解いた際の一致、
  反復上限0での状態復元、その後の正常な再利用、不正な初期組成の拒否を確認。
- 上記とC/C++ API・Model所有権テストをUBSanで実行して成功。
- profileとrandom-restart診断を有効にした構成でも固定T契約テストが成功。
- 新規production C++翻訳単位は `-Wall -Wextra -Werror -pedantic` で構文確認済み。
- 既存chemistryチェック、各ionization、CO解離、元素正規化のテストが成功。
- AREPOの公開API・依存検査、既存明示/陰解法回帰の2構成が成功。
  LU再利用・解析Jacobianを使う1ZONE構成に加え、これらを指定しないAREPO回帰構成でもビルド・実行した。
- AREPO側へ同期し、共有chemistry/radiation 73ファイルがbyte単位で一致。

既知の初期化回帰は引き続き3件失敗。AREPO全実行ファイルは未ビルド。
ASanのランタイム初期化停止は工程2の記録を参照。今回はUBSanを使用した。
ARCHEの接続実装・検証は行っていない。変更は未commit。

```sh
make -j4 BUILDDIR=build-cpp-stage3 EXEC=build-cpp-stage3/1ZONE \
  check-chemistry-fixed-rate check-chemistry-cpp build-cpp-stage3/1ZONE
make -j4 BUILDDIR=build-cpp-stage3-ubsan \
  OPTIMIZATION='-O1 -fsanitize=undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer' \
  check-chemistry-fixed-rate check-chemistry-cpp check-chemistry-model
make check-energy-replay
```
