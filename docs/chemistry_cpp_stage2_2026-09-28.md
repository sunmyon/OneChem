# C++化・工程2: Modelの所有権と読み取り専用の参照

2026-09-28。工程1に続き、run全体の資源所有をC++へ移した。

## 変更

- `chemistry::Model` が反応network、反応係数・分子冷却・dust・continuum各tableを所有する。
  constructorは非公開で、factoryはstatusと `unique_ptr<const Model>` を返す。
  本体のcopy/moveは禁止。所有ハンドルは移動でき、借用中のModelのアドレスは変わらない。
- テーブルの生成順と生成処理は既存の `InitChem` を利用し、destructorから `DestroyChem` を呼ぶ。
  初期化途中の失敗でも、公開contextの追加確保に失敗した場合でも、所有ハンドルが資源を解放する。
  `ChemRunContext` はC初期化処理の内部storageとして残る。数値処理の書き換えは行っていない。
- 計算側の `ChemEnvironment::run` は `ChemRunView` を参照する。
  viewはコピー可能な非所有の値型で、network/tableはconst pointerで参照する。
  Modelを破棄する前に借用を終える。constだけでthread safetyを保証する設計ではない。
- 公開API実装を `.cc` に移し、opaqueな `ChemistryContext` がModelを所有する。
  公開C ABIは維持。caller所有のworkspaceにはplacement newでC++の寿命を開始し、
  destructor不要であることをstatic assertionで確認する。
- standaloneも同じ公開context生成・破棄へ移行した。既存C adapterは内部bridge経由でviewを借用する。
  replay等の早期returnもcleanupへ集約し、借用を解除してからcontextを破棄する。
  adapterで変更するtimescale factorはviewの値コピーに反映し、Model内の値は変更しない。
- AREPOの共有ソース・build fragment・関連する既存回帰テストを追従させた。
  Cの翻訳単位とC++17の翻訳単位を個別にコンパイルする。

既存のC値型をそのまま利用できる部分に、重複するC++ wrapperは追加していない。
固定T Newton、energy solver、driverの数値制御自体は次工程以降。

## 検証

以下が成功した。

- 1ZONE混在ビルド、`check-chemistry-cpp`。
- `check-chemistry-model`: 所有ハンドル移動後の参照安定性、複数Modelの独立性、
  table_rootの寿命、invalid config、Model確保失敗、公開context確保失敗、
  後段のdust table読込失敗を繰り返した後の正常生成・破棄。
- 上記2テストをUBSanで実行。`-fno-sanitize-recover=undefined` を指定して成功。
- 既存のchemistryチェック、high/metal-ionization、photoionization、CO解離、元素正規化。
- energy-step replay 4件。C版基準 `782d710` のログに対し、実行時間 `elapsed_s` 以外は
  反復診断・最終温度・出力組成を含めて完全一致。
  [比較記録](chemistry_cpp_stage2_validation_2026-09-28.json)参照。
- AREPO側の実際のmake fragmentで共有core一式とC adapterをコンパイル。
  C callerをC linkerと既存のC++ runtime library指定でリンクし、context生成・破棄を実行。
- AREPOの公開ヘッダ・依存チェック、および既存 `explicit-implicit-check` の2構成。
- 共有chemistry/radiation 68ファイルのbyte一致。

AREPOの非energy構成のテストでは、既存の `chemistry_print_rerun_line` が
`CHEMISTRY_ENERGY_IMPLICIT` 内で定義されている一方、化学solverからは常に呼ばれるビルド不具合が見つかった。
関数本体を変えず条件コンパイルの外へ移した。energy構成の計算と出力は変わらない。

## 未解決・検証範囲

- 既知の初期化回帰テストは引き続き `short-step solve failed` が3件。
- ASan付きバイナリはランタイム初期化で停止し、時間制限内にテストへ進まなかった。
  サンドボックス外、MPIをリンクしない構成でも解消せず、ASanは未検証。
  UBSanと通常テストの成功を、ASanやリーク検証の代替結果としては扱わない。
- AREPO全実行ファイルは未ビルド。検証対象は共有chemistry一式と既存C adapter・関連テスト。
- ARCHEへの接続実装・接続検証は行っていない。
- 両作業ツリーの変更は未commit。

再実行例（1ZONE root）:

```sh
make -j4 BUILDDIR=build-cpp-stage2 EXEC=build-cpp-stage2/1ZONE \
  check-chemistry-cpp check-chemistry-model build-cpp-stage2/1ZONE
make -j4 BUILDDIR=build-cpp-stage2-ubsan \
  OPTIMIZATION='-O1 -fsanitize=undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer' \
  check-chemistry-cpp check-chemistry-model
make check-energy-replay
```
