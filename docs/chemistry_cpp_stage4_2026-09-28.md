# C++化 工程4: BE行列とlinear backendの分離

固定T Newtonから数値処理を分離した。数値アルゴリズム・演算順・再利用条件は変更していない。

## ファイルと責務

- `src/chemistry/numerics/backward_euler.hpp`: 次元とaccessorを受け、
  `(I - dt J) delta = dt R - total_step` の行列・右辺を作る。
  column-major配置、ループ順、右辺のfmaを保持。種名・N_SP・chemistry型への依存なし。
  accessorはtemplate引数で、種から縮約座標への変換はchemistry呼出し側に閉じる。
- `src/chemistry/numerics/lapack_dense.hpp/.cc`: 状態を持たない `numerics::LapackDense`。
  `solve` (dgesv)、`factorize` (dgetrf)、`solve_factored` (dgetrs) を提供。
  LP64のintとdoubleを使用する既存LAPACK ABIを保持。戻り値はLAPACK INFO。
  バッファを所有せず、呼出し側が正の次元と十分な容量を保証する内部API。
  分解済み行列のsolveには、同じ次元で成功したfactorizeが必要。
- `src/chemistry/solver/fixed_rate_solver.cc`: 種選択・gather/scatter・profiling・
  LU再利用条件を維持し、上の数値部品を利用。
  LU有効フラグはNewtonSession、行列・pivot・右辺はNewtonScratchが引き続き所有。
  chart変更・restartでの無効化もSessionが担う。backend側に隠れたcacheはない。
- `src/chemistry/solver/chemistry_solver.c`, `fixed_rate_c.h`:
  固定T向けCのlinear solve関数と宣言を削除。他のC経路のdgesv呼出しは維持。
- `src/chemistry/tests/numerics_dense.cc`: chemistry設定なしでもビルドできる数値テスト。
- ビルドsource/header依存・core依存検査・AREPO共有一致検査にnumericsを追加。

この段階は独立した数値部品の抽出であり、FixedRateSolver自体の汎用化は未完了。
保存則・clipping・chart再選択と反復制御の境界、統一SolveReport、energyとdriverのC++化が残る。
Model・solver・driverの所有関係は変更していない。ARCHEへの接続作業は行っていない。

## 検証

- chemistryのinclude/configなしで、数値テストをC++17・厳格警告設定でビルド・実行。
  非対称2次元行列の配置と解、同じLUで異なる右辺を解く処理、特異行列の失敗、1次元の解を確認。
- `check-chemistry-numerics`, `check-chemistry-fixed-rate`, C/C++ API、Model寿命テスト成功。
- LU再利用を無効にした構成でも固定T・数値テスト成功。
- 固定T・数値テストをUBSanで実行して成功。
- energy replay全4件がC版782d710とelapsed_s以外のログ全体で一致。
- AREPOモジュールのAPI・core依存検査・明示/陰解法の2構成の回帰成功。
- AREPO共有chemistry/radiationは76ファイル、差分0。

詳細: `chemistry_cpp_stage4_validation_2026-09-28.json`。
既存の初期化回帰3件の失敗は未解決（今回は再実行していない）。
AREPO全体の実行ファイルとARCHE接続は検証対象外。変更は未commit。
