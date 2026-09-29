# C++化・工程1: ビルドとC API境界

2026-09-28。基準コミットは1ZONE `782d710`、AREPO `f7d751c3e4`。
設計は [chemistry_cpp_arche_design_2026-09-28.md](chemistry_cpp_arche_design_2026-09-28.md) を参照。

## 実装

- 関数を宣言する公開ヘッダ3個に、それぞれ独立した `extern "C"` ガードを追加。
  Cの型レイアウト、関数名、引数、ABI versionは変更していない。
- 非数値処理の `chemistry_status.c` を内容不変で `.cc` に移行。
  C側からC++実装を呼ぶ経路を実際の実行ファイルに導入した。
- portable構成では `.c` をCC、`.cc` をCXXでコンパイルし、CXXでリンクする。
  C++17を指定し、C用の言語標準指定をC++コンパイルへ渡さない。
  legacy構成のリンクはCCのまま。初期化テストもCでコンパイルしてからリンクする。
- AREPOのchemistry専用make fragmentにも混在コンパイルを追加。
  `CHEMISTRY_CXX` は既存 `CPPC`、未設定なら `mpicxx`。必要に応じてホストと同じ
  MPI/toolchainのwrapperを指定できる。C用の `-std=` はC++コンパイル時だけ除去する。
  pattern ruleを用い、AREPO全体のコンパイラ、リンク設定、既定targetは変えない。
- parityチェックとソース走査の拡張子フィルタにC++を追加。

数値処理、solver構造、閾値、反応・冷却式は変更していない。
Modelなどの所有権classの実装は次工程。

## 検証

```sh
make -C src/chemistry check high-ionization-option-check metal-ionization-check \
  photoionization-check co-photodissociation-check carbon-oxygen-normalization-check
make BUILDDIR=build-cpp-stage1 EXEC=build-cpp-stage1/1ZONE check-chemistry-cpp
make check-energy-replay
```

上記はすべて成功。
`check-chemistry-cpp` はC++から実際のC APIをリンクしてcontext生成・破棄、workspace、
species照会、step準備、opacityを実行する。Cの別翻訳単位からC++のstatus実装を呼び、
C++で埋めた公開値型をCで読む。既存AREPO adapterもCでコンパイルしてリンクし、
context lifecycleと不正入力のstatus返却を確認する。
残る公開関数のリンク確認には不正入力を使い、数値経路は既存replayで検証する。

変更前後のenergy-step replay 4件は、実行時間 `elapsed_s` のみを除くと、
反復診断・最終温度・全出力組成を含むログ全体が完全一致した。
compiler、feature flagsのhash、比較したログのhashと最終値は
[検証記録](chemistry_cpp_stage1_validation_2026-09-28.json)に保存した。
生ログはローカルの `build-cpp-baseline/` と `build-replay-*.log` にある。

既知の `check-chemistry-initialization` の `short-step solve failed` 3件は引き続き失敗。
移行前にも確認済みの問題であり、成功扱いにはしていない。
同targetの絶対パスBUILDDIRの実行不具合は修正した。

## AREPO同期と検証範囲

`feature/1zone-chemistry` 作業ツリーの公開ヘッダ、status実装、make fragmentを同期。
共有chemistry/radiation 65ファイルはbyte単位で一致する。
AREPO側の実際のmake fragmentを最小のビルド用Makefileに読み込み、C11 callerと
C++17 status実装を個別コンパイル・リンク・実行する検証も成功した。

AREPOの全実行ファイルは未ビルド。この作業ツリーには実行用 `Config.sh` がなく、
今回確認したのは既存adapterのリンク・実行とchemistry用ビルド規則まで。
ARCHEへの接続実装・接続検証は今回の対象外。
変更は両作業ツリーで未commit。
