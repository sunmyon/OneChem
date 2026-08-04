# chemistry_new 整理・拡張準備のまとめ

## 目的

この作業の目的は、現在の H/He/D chemistry の計算結果と速度を保ちながら、C/O chemistry を追加できる形へ `chemistry_new` を整理することである。

全面的な書き直しは行わず、現在の solver 選択と時間積分の流れを維持したまま、次の部分を段階的に一般化した。

- species 定義
- reaction network
- explicit/implicit chemistry
- 中間生成物の平衡近似
- 保存則による implicit 系の縮約
- opacity table の探索
- dust 温度計算

ここまでの主要な変更は commit `6e5f0b3` (`Generalize chemistry network and optimize dust cooling`) として `origin/main` に反映済みである。

## 設計上の結論

### solver と rate/cooling の関係

reaction rate や cooling rate を solver から完全に独立させる設計は採用しない。現在の計算は solver の選択によって、次のように評価順序そのものが変わる。

1. timestep を推定し、solver を選択する。
2. chemistry explicit の場合は `chem_exp + radiative cooling` を評価する。
3. chemistry implicit / energy explicit の場合は `chem_imp + radiative cooling` を評価する。
4. energy implicit の場合は、温度反復の各 trial で `chem_imp + radiative cooling` を再評価する。

したがって、solver、reaction rate、chemical heating/cooling、radiative cooling の間には本質的な依存がある。今回の整理では、この依存を消すのではなく、データと責任範囲を明確にした。

### adapter は導入しない

外部状態と chemistry の間に互換層を積み重ねる方式は採用していない。`chemistry_new` の移行を直接進め、データのコピーや変換は必要な入口・出口に限定する。

## reaction network と species

### species 定義の一元化

species の列挙、名前、電荷、元素組成を `CHEM_SPECIES_LIST` から生成する形にした。

- `chemistry_new/species_definitions.h`: species の宣言元
- `chemistry_new/species.h`: species metadata の公開interface
- `chemistry_new/species.c`: metadata table と名前解決
- `chemistry_new/chemistry.h`: `enum SPECIES`、`N_SP`、`TGCHEM_NUM_ABUNDANCES`

現状では次が常に一致する。

```text
N_SP == TGCHEM_NUM_ABUNDANCES == enum SPECIES の要素数
```

これは C/O species を追加する際の二重管理を減らすための暫定的な唯一の定義源である。最終的に AREPO など外部コード側を唯一の定義源とする場合は、この一覧自体を外部定義から生成する必要がある。

### reaction table

純粋な forward reaction の化学量論を `chemistry_new/tables/reaction_network.tsv` から読み込むようにした。

初期化時に次を検査する。

- species 名が定義済みであること
- 行の形式が正しいこと
- 元素数と電荷が反応前後で保存されること
- reaction が一件以上読み込まれること

inverse reaction と photo-reaction は forward table に混在させない。critical density、shielding、chemical heating/cooling と密接に結びつく処理も、単純な reaction table へ無理に押し込めない。

### species ごとの reaction index

network 初期化時に、各 species について次の reaction index を構築する。

- その species を reactant として読む反応
- その species の生成・破壊に影響する反応

explicit 更新、production/destruction、chemical timestep、Jacobian column は、network 全体を毎回走査せず、この index を利用する。

## explicit chemistry

production/destruction と chemical timestep は `ChemReaction` から計算する。forward と inverse の両方を explicit 評価へ含める。

species 更新方法は compile-time option で選択できる。

```text
EXPLICIT_CHEMISTRY_UPDATE=SEQUENTIAL    # default
EXPLICIT_CHEMISTRY_UPDATE=SIMULTANEOUS
```

`SEQUENTIAL` は従来の順序依存更新を維持する。`SIMULTANEOUS` は同じ旧状態から全 species の production/destruction を評価する。収束性と従来挙動を優先し、default は `SEQUENTIAL` のままとした。

## 中間生成物の平衡近似

現在は HM を対象としているが、実装は複数中間体を扱える形にしている。

```text
INTERMEDIATE_EQUILIBRIUM=0    # 全speciesを非平衡値として保持
INTERMEDIATE_EQUILIBRIUM=1    # 指定した中間生成物を代数的平衡として扱う
```

主な実装内容は次の通りである。

- 中間体に関係する reaction index だけから生成率・破壊率を計算する。
- 単一中間体は直接平衡値を求める。
- 複数中間体は収束判定付きの反復で同時に更新する。
- explicit chemistry でも中間体の平衡値を更新する。
- implicit chemistry では中間体を Newton の独立変数から外す。
- 中間体平衡値の連鎖微分を解析 Jacobian に組み込み、一般の複数中間体行列として解く。
- 解析的縮約が失敗した場合、または20反復以降も収束しない場合は numerical Jacobian へ切り替える。

これにより、C/O network で複数の短寿命中間体を追加しても、同じ枠組みを利用できる。

## 保存則による implicit 系の縮約

全speciesを厳密な保存行列へ入れるのではなく、現在の abundance hierarchy に合わせ、支配的な項だけから次の関係を課している。

```text
H  + H+  + 2 H2       = H total
He + He+ + 2価He      = He total
D  + D+  + HD         = D total
e                     = H+ + He+ + 2 (2価He)
```

実際の helium 保存式では各 ion の原子数係数は1で、最後の式だけ二価Heの電荷係数が2である。

HI、HeI、DI、electron を従属変数として Newton 系から外し、独立変数の更新後に保存則から再構成する。中間体平衡と保存則を同じ algebraic species 系として扱い、解析 Jacobian を縮約する。

positivity と上限を満たすため、Newton step 全体に減衰係数を適用する。ただし常に強く減衰するのではなく、trial state が制約を満たす最大限の step を使う。

### 保存則へ全speciesを含めた場合の速度低下

保存則の一般化を試した段階では、元素を含む全speciesと中間生成物を保存式へ入れた。しかし、この方法では計算速度が大幅に低下した。

原因は保存式の加算回数ではない。保存則へ中間生成物を含めると、平衡中間体と保存則の従属speciesが代数的に結合する。その結果、独立変数を一つ更新するたびに、次の連鎖をJacobianへ反映する必要が生じる。

```text
independent species
  -> conservation-dependent species
  -> intermediate-equilibrium abundance
  -> reaction production/destruction
```

これにより、本来は関係するreactionだけを読めばよかったJacobian columnの結合範囲が広がり、縮約行列が密になる。また、Newtonのtrial step、line search、numerical Jacobianへのfallbackでも、保存則と中間体平衡を繰り返し再構成する必要がある。したがって、保存式自体が数個の加算であっても、solver全体の評価回数と一回当たりの処理量が増える。

一方、現在のnetworkでは、次の項は主要な保存量に対して十分小さい。

- H保存に対するHDは約 `1e-5` 程度
- H保存に対するHMなどの中間生成物はさらに小さい
- charge neutralityに対するD+はH+より十分小さい
- charge neutralityに対する負イオンなどの中間生成物も小さい

これらを保存則へ含めても物理精度の改善はほとんど得られず、計算コストだけが増える。そのため現在は、H保存からHDと中間生成物を外し、charge neutralityからD+と中間生成物を外している。D保存ではD、D+、HDを保持する。

この経験から、C/O chemistry追加時にも「元素を含む全speciesを自動的に保存式へ入れる」方針は採用しない。主要reservoirだけで従属speciesを再構成し、trace speciesと平衡中間体はreaction network側で扱う。C/O abundanceが支配的になり、無視した項が必要な精度を超える場合に限って保存式へ追加する。

## opacity table の改善

### continuum opacity

temperature と density の格子探索を線形走査から二分探索へ変更した。一般の不等間隔tableを保ったまま、探索コストを減らしている。

### dust opacity

dust opacity table は log temperature で等間隔なので、`log(Tdust)` から格子番号を直接計算する。

- Rossland opacity と Planck opacity を一度の探索で取得する。
- Planck opacity の温度微分も同時に計算する。
- 推定格子の丸め誤差に備え、隣接格子だけ補正する。

旧二分探索との比較では、1-zone の進化結果は bitwise に一致した。dust 関連処理の全体時間比は、二重table探索を行っていた時点の約12.5%から、おおむね6.4--7.0%まで低下した。

## dust 温度計算

### 平衡式

dust 温度は、dust-gas collision と放射収支の平衡

```text
F(Tdust) = Lambda_dust(Tdust) - Lambda_coll(Tdust) = 0
```

から求める。

gas の energy exchange には、低密度側では `Lambda_coll`、高密度側では `Lambda_dust` を使用する。低密度では `Tdust` が放射温度に、高密度では `Tdust` がgas温度に極めて近づくため、全密度域で片方だけを使うと、近い温度同士の差による丸め誤差が反対側の極限で顕在化する。現在は二つの温度差を比較し、数値的に条件のよい側を選ぶ。

### 低密度で発生していた false heating

従来実装では、低密度領域でdust由来の非物理的なgas加熱が発生し、温度進化に約 `1.4e4 K` の一時的なpulseが現れる場合があった。これはdustの物理過程ではなく、根探索の有限精度をgasのenergy exchangeとして扱ったことが原因だった。

低密度ではgas-dust collisionが非常に弱いため、dust平衡温度はほぼ放射温度に等しい。

```text
Tdust = Trad + delta T
delta T ~ 1e-11 K
```

一方、従来の初期bracketは `[0, 2000 K]` で、二分法を常に30回行っていた。30回後にも残る区間幅は

```text
2000 K / 2^30 = 1.86e-6 K
```

であり、必要な `delta T ~ 1e-11 K` より約5桁粗い。したがって、返された `Tdust` では

```text
Lambda_dust - Lambda_coll != 0
```

という有限のresidualが残る。特に `Tdust` と `TCMB` がほぼ等しいときは、`Tdust^4 - TCMB^4` の直接差も桁落ちしやすい。

旧実装は、この未収束な `Lambda_dust` を

```text
edot_dust = -Lambda_dust * rho * metallicity * dust_to_gas_ratio
```

としてgas energyへ戻していた。本来gasとdustの間で交換されるenergyは `Lambda_coll` だけであるため、放射平衡式に残った数値誤差が、そのまま実在するgas heatingとして解釈されていた。

現在は次のように修正している。

1. `Tdust^4 - TCMB^4` を

   ```text
   (Tdust - TCMB)(Tdust + TCMB)(Tdust^2 + TCMB^2)
   ```

   と因数分解し、近接温度での桁落ちを抑える。
2. residualだけでなくbracket幅でも収束を判定する。
3. 決定した最終 `Tdust` でbalanceを再評価する。
4. gasのenergy exchangeには、低密度側で `-Lambda_coll`、高密度側で `-Lambda_dust` を使用する。
5. bracket幅で終了するとき、得られたNewton点をbracket中点で上書きしない。

この修正により、根がmachine precisionに近すぎて放射residualを完全にはゼロにできない場合でも、低密度gasへfalse heatingが入らなくなった。

### 不連続性の除去

従来の collision factor は 159、275、425、680 K にhard stepを持っていた。このため residual が不連続になり、境界付近では数学的な根が存在しない場合があった。

各境界を相対幅 `1e-3` の `tanh` transition で連続化し、その微分も解析 residual derivative に含めた。324ケース比較では、この連続化によるgas温度進化の差は最大相対約 `9.4e-7`、最大絶対約 `9.3e-4 K` だった。

hard stepを残した試験では、二分法のbracketが275 Kなどの不連続点へ収縮し、`Tdust`がその境界に張り付く場合があった。これはbracket両端でresidualの符号が異なっていても、その間に `F(Tdust) = 0` を満たす温度が存在するとは限らないためである。

```text
F(T_boundary - 0) < 0
F(T_boundary + 0) > 0
```

というjumpがあれば、二分法は通常の連続関数と同様に区間を狭めるが、得られるのはzero crossingではなく不連続点である。局所微分も定義できないため、Newton法へ移っても安定した改善は期待できない。今回のcollision factorについては、hard stepを物理的に厳密な相転移として扱う必要性が低く、`tanh`で連続化するのが適切と判断した。

ただし、境界に張り付く解が常に数値的な失敗というわけではない。次の場合は境界解がphysical、または少なくとも採用している物理モデルの範囲内で正しい場合がある。

- collisionが消える極限で `Tdust = Trad` になる。
- gas-dust couplingが非常に強い極限で `Tdust` が `Tgas` に近づく。
- dustの存在温度やopacity tableの有効範囲に上限を課し、非制約平衡解がその範囲外にある。
- grain組成の相転移やsublimationを真の不連続modelとして与え、境界の両側からのenergy収支が境界へ戻す向きになる。

最後の場合は通常の `F = 0` ではなく、境界の片側極限を使って安定性を判定する必要がある。例えば

```text
F(T_boundary - 0) <= 0
F(T_boundary + 0) >= 0
```

で、温度変化が両側から境界へ向かうなら、その境界は一般化された平衡解として扱える。一方、単にtableが2000 Kで終わっているだけなら、`Tdust = 2000 K` はtable上の境界解であり、それだけではsublimationを物理的に解いたことにはならない。

### bracket と Newton 法

初期区間は物理的に根を挟む

```text
[min(Tgas, Trad), max(Tgas, Trad)]
```

を基本とし、opacity table の上限に合わせて `Tdust <= 2000 K` とする。

```text
Trad^4 = TCMB^4 + stellar_flux / (4 sigma)
```

初期値は区間中央ではなく `Trad` とし、最初から解析微分によるNewton候補を計算する。候補がbracket外に出る場合、または直前stepに対して大きすぎる場合は区間中央へ戻す safeguarded Newton とした。

停止条件は次の二つである。

```text
|F| <= 1e-5 (|Lambda_dust| + |Lambda_coll|)
bracket width <= max(1e-4 K, 1e-5 Tdust)
```

二つ目は `Tdust` と `TCMB` がほぼ等しく、両energy termが同時にゼロへ近づく場合に、相対residualだけで無意味な反復を続けることを防ぐ。

### 二分法とNewton法の将来的な切り替え

現在のsafeguarded Newtonは十分高速で安定しているが、C/O coolingやgrainの相転移を追加すると、residualの滑らかさや許容温度域が変わる可能性がある。将来的には、単なる固定反復回数ではなく、局所的な関数形に応じて次のように切り替えるのがよい。

1. `Trad`から開始し、物理的な許容区間でbracketを作る。
2. residualが連続で解析微分が有限なら、最初からNewton候補を計算する。
3. Newton候補がbracket内にあり、stepが十分小さく、residualまたはbracket幅を実際に減らす場合だけ採用する。
4. 候補が区間外へ出る、微分が小さすぎる、または改善しない場合は二分法へ戻る。
5. table cell境界や急なtransition付近では、区間幅がtransition幅またはtable分解能程度になるまで二分法を優先する。
6. bracketが既知の不連続境界へ収縮した場合は、通常の収束とみなす前に片側residualを評価する。
7. 片側residualが境界へ戻す向きならphysical boundary solutionとして受理し、そうでなければmodel不連続または探索失敗として記録する。

これはNewton法を主として使い、二分法でbracketを保証する現在の方法を、Brent法に近い受理判定へ発展させる方針である。重要なのは、区間が狭くなったという理由だけで不連続点を根とみなさないことと、physical boundary solutionを無理に内部の `F = 0` へ置き換えないことである。

将来の性能評価では、平均反復回数だけでなく、次も個別に数えるべきである。

- Newton候補の採用・棄却回数
- 二分法へfallbackした回数
- smooth root、physical boundary、table上限の各終了回数
- 不連続点を挟んだbracketの検出回数
- 各終了理由ごとのenergy residual

### 反復数と速度

物理bracketと `Trad` 初期値を導入した比較結果は次の通りである。

| 指標 | 変更前 | 変更後 |
|---|---:|---:|
| 平均反復回数 | 13.3777 | 6.4073 |
| dust balance評価回数 | 約396万回 | 約213万回 |
| balance評価回数の削減 | - | 46.3% |

診断runでは全体実行時間が約5.23秒から4.93--4.99秒、平均4.96秒へ短縮した。ただしadaptive evolutionの経路とdust関数の呼び出し総数も変わるため、この約5%は厳密な同一軌道benchmarkではない。

### HM平衡・非平衡の全体時間配分

現行commit `59f6f33` を独立したbuild directoryでclean buildし、同じ36条件surveyについてHM非平衡とHM平衡を改めてprofileした。

```text
HM非平衡: INTERMEDIATE_EQUILIBRIUM=0
HM平衡:   INTERMEDIATE_EQUILIBRIUM=1
compiler optimization: 現行Makefileと同じ -O0
```

wall timeはprofilerを付けず各modeを3回測定した。HM非平衡は `4.36, 4.37, 4.33 s`、HM平衡は `4.94, 4.92, 4.94 s` で、平均はそれぞれ `4.35 s` と `4.93 s` だった。HM平衡による全体時間の増加は約 `13.3%` である。

時間配分はmacOS sampling profilerを1 ms間隔でprocess終了まで実行し、各mode 2回分を合算した。HM非平衡は7,405、HM平衡は8,427のmain-function sampleを得た。表の秒数はsampling比率をprofilerなしの平均wall timeへ掛けた概算である。

| 処理 | HM非平衡 | 概算時間 | HM平衡 | 概算時間 |
|---|---:|---:|---:|---:|
| reaction network Jacobian・RHS・化学量論演算 | 57.35% | 2.497 s | 52.56% | 2.593 s |
| 保存則・algebraic Jacobian縮約 | 11.63% | 0.506 s | 9.66% | 0.477 s |
| HM平衡値評価・中間体判定 | 0.92% | 0.040 s | 10.51% | 0.519 s |
| reaction rate評価（dustを除く） | 2.39% | 0.104 s | 2.44% | 0.121 s |
| dust温度・dust opacity | 4.62% | 0.201 s | 4.12% | 0.203 s |
| radiative・chemical cooling | 1.55% | 0.068 s | 1.31% | 0.064 s |
| LAPACK/BLAS線形方程式 | 0.95% | 0.041 s | 1.90% | 0.094 s |
| implicit solverのその他の処理 | 7.62% | 0.332 s | 6.35% | 0.313 s |
| energy solver・thermodynamics | 0.68% | 0.029 s | 0.58% | 0.029 s |
| math library (`pow`, `exp`, `log`など) | 6.41% | 0.279 s | 5.67% | 0.280 s |
| 出力I/O | 3.97% | 0.173 s | 3.42% | 0.169 s |
| その他 | 1.92% | 0.083 s | 1.48% | 0.073 s |
| **合計** | **100%** | **4.353 s** | **100%** | **4.933 s** |

samplingはtop-of-stackを重複しない分類へまとめた。`pow`、`exp`、`log`などは複数の物理kernelから呼ばれるため、呼び出し元へ推定配分せずmath libraryとして独立させた。したがって、例えばreaction networkの実効コストにはmath library行の一部も含まれる。

HM平衡化による約 `0.58 s` の増加のうち、約 `0.48 s` はHM平衡値評価・中間体判定の増加で直接説明できる。線形方程式の概算時間も約 `0.05 s` 増えており、中間体応答を含むalgebraic Jacobian処理の影響が見える。一方、dust温度計算は両modeとも約 `0.20 s` で、全体時間に占める割合の差は分母となる総時間の違いによる。radiative cooling、thermodynamics、I/Oも絶対時間はほぼ変化していない。

現時点で最大の計算コストは、両modeともreaction networkのJacobian・RHS評価である。HM平衡modeをさらに高速化する場合は、dustやcoolingではなく、HM平衡値をNewton iteration中に何度再評価しているか、同じabundance/rateに対する再計算を省けるかを先に調べるべきである。

## 検証

commit前に、独立したbuild directoryから全objectを作り直し、次を確認した。

- `INTERMEDIATE_EQUILIBRIUM=0`: compile成功、36条件完走
- `INTERMEDIATE_EQUILIBRIUM=1`: compile成功、36条件完走
- 両modeの主要な `.dat` 出力に NaN/Inf なし
- `git diff --check` 成功

補助的な `_HM` diagnostic出力には、未使用の列に既存のNaNが残っている。主要な進化状態には現れていないが、diagnostic出力を今後利用する場合は別途整理が必要である。

また、出力先を開けなかった場合に `FILE *` が `NULL` のまま `fprintf()` され、SIGSEGVになる問題を修正した。現在は対象パスと `strerror(errno)` を表示し、終了コード1で停止する。

## build時の注意

現在のMakefileは、compile-time optionの値そのものをobjectの依存関係として追跡しない。例えば `INTERMEDIATE_EQUILIBRIUM=0` で作ったobjectを残したまま `=1` へ変更しても、source timestampが変わらなければ全objectは再構築されない。

optionを切り替えて比較する場合は、clean buildまたはmodeごとに異なる `BUILDDIR` を使う必要がある。

## C/O chemistry追加前の残課題

### 必須

1. `species_definitions.h` に C/O species を追加する。
2. 純forward reactionを `reaction_network.tsv` に追加する。
3. forward rate tableへ対応するrate式を追加する。
4. inverse reaction、photo-reaction、critical density依存をそれぞれ既存の専用経路へ追加する。
5. C/Oを含む保存則とcharge constraintで、どのspeciesを従属変数にするか決める。
6. 平衡近似するC/O中間体を `intermediate_equilibrium_species` に追加する。
7. elemental conservation、Jacobian、explicit/implicit一致を小networkから段階的に検証する。

### 後回しでよいもの

- dust温度solverの追加最適化
- 一般の不等間隔dust opacity table対応
- 2000 K以上のdust sublimationを含む物理モデル
- 外部species catalogからの自動生成

`Tdust = 2000 K` は現在のopacity table上限による境界解であり、2000 K以上の物理的平衡解を表現しているわけではない。dust sublimationを扱う段階では再検討が必要である。

## 現時点の判断

H/He/D networkの一般化、species metadata、reaction index、explicit/implicit更新、中間体平衡、保存則縮約までが一つの枠組みに入った。dust温度計算も、不連続性、bracket、収束判定、table探索、energy exchangeの主要な問題を解消した。

したがって、C/O chemistry追加前のdust部分は現時点で十分であり、次の主作業はC/O speciesとreactionを小さなsubsetから追加し、現在の一般化が実際に機能するか検証することである。

---

## 追記: C/O chemistry導入後の収束・energy・cooling診断（2026-08-04）

この節は、C/O chemistryとCO/OH/H2O line coolingを導入した後に現れた温度のぎざぎざ、chemical failure、dust cooling spikeを調べた記録である。上の「C/O chemistry追加前の残課題」より新しく、dust energy exchangeとC/O networkに関しては本節を現在の結論とする。

### 対象species

非平衡speciesとして次を追加した。

```text
C, C+, CO, CO2, O, OH, H2O, O2
```

短寿命中間体として次を扱う。

```text
HM, CO+, CH, CH2, O+, OH+, H2O+, H3O+, O2+
```

`INTERMEDIATE_EQUILIBRIUM=1` ではこれらをchemistry内部の代数的平衡speciesとし、Newtonの独立変数から外す。`INTERMEDIATE_EQUILIBRIUM=0` ではHMを含む全speciesを非平衡変数として時間発展させる。

中間体の結合構造は一様ではない。

- `CO+`, `O+`, `OH+`, `H2O+`, `H3O+`, `O2+` は、非平衡speciesを固定すれば逐次的または一意に決まる。
- `CH` と `CH2` は相互にcoupleするため、小さな連立平衡系になる。
- `HM` も平衡modeでは他の中間体と同じ代数的speciesとして扱う。

式の本数が少なくても、主要reservoirがほぼCOまたはO2へ移る領域では、解がpositivity boundや元素reservoirの境界に非常に近い。この条件の悪さが、単純な反応式に見えるにもかかわらずNewton収束が悪化した主因である。

### `delta_H2` とchemical heating

implicit chemistry後のH2結合energyは、基本的に

```text
delta_H2 = abundance_new[H2] - abundance_old[H2]
edot_H2 = 4.48 eV * delta_H2 * n_nuclei / dt
```

として評価する。これはH2の生成・破壊を経路ごとに完全分解するのではなく、step全体のH2存在量変化に対応する結合energyをgasへ戻す扱いである。

現在の判断は、最終的にH2が形成されるなら結合energyはdepositされるため、このstate-difference型の評価は基本的に維持する、というものである。低密度側ではcritical densityを使い、HM経路、dust形成、FUV解離、three-body形成のうち放射へ逃げる部分を補正している。高密度側では熱化が支配的として `4.48 eV * delta_H2` を直接使う。

C/O反応を経由した場合に、同じH2変化へ対応する余剰energyの一部がline photonまたは内部励起へ逃げる可能性は残る。ただし、今回観測した温度ぎざぎざについては、chemical solve失敗とabundance jumpが明確に同期していた。C/O経路のenergy partitionを原因と仮定して補正する前に、まず非線形solveを連続にする方針を採った。現時点ではC/O経路のenergy partitionに対する定量的補正は入れていない。

### Newton boundは物理timestepの短縮ではない

この点は実装思想として明確に区別する。

```text
physical timestep dt: 変更しない
Newton correction delta: abundance bound内へ成分ごとに補正する
initial guess: 同じdtのまま変更して再探索する
```

一時期試したNewton vector全体への `step_scale *= 0.5` は採用しない。全成分を一括縮小すると、支配的reservoirがboundへ近づく問題を解消せず、単にbound付近へゆっくり張り付かせる可能性がある。また、これは化学stepそのものの半減と誤解されやすい。

現在の `chemistry_apply_implicit_newton_step()` は、各非代数speciesについて次の旧来型処理を行う。

```text
candidate < 0.1 * current
    -> delta = -0.9 * current

candidate > upper_bound
    -> そのspeciesのdeltaだけを上限側へclip
```

bound外のNewton候補は、それ自体をchemical failureとはみなさない。補正した非平衡speciesから保存則従属speciesと平衡中間体を再構成し、その代数的再構成が有限かつbound内に入らない場合に限って、現在のguessを棄却する。物理timestep `dt` は維持する。

診断中に「最初にboundを破るspeciesがCI」と表示されたことは、CIが物理的に負になることや、CI反応だけが誤っていることを直接意味しない。Newton correctionを現在の独立変数順で検査したとき、最初にCIのtrial値が設定上限を越えた、という意味である。C/O保存則を使うmodeでは、低密度側はCI/O、高密度側はCO/O2を従属speciesに切り替え、主要reservoirがほぼCO/O2になる領域で境界変数を直接Newton更新しない構造も用意している。今回の図を作ったbuildは `CARBON_OXYGEN_CONSERVATION=0` であり、この保存則切替の影響は含まない。

### 収束しないstepとguess reset

`Z=10^-3, z=0` の `n ~ 5e6--5e7 cm^-3` で失敗しやすいchemical stepを抜き出し、同じ `dt` に対して元素reservoirを保存するrandom guessを多数与えた。目的はstepを短くすることではなく、そのstepが収束する先を調べることだった。

代表的な収束先は次の通りで、異なるrandom guessからほぼ同じ枝へ到達した。

| n [cm^-3] | T [K] | CO | CO2 | O2 | OH | H2O |
|---:|---:|---:|---:|---:|---:|---:|
| 6.9e6 | 119--122 | 3.118e-7 | 4.20e-9 | 1.138e-7 | 2.15e-8 | 5.5e-9 |
| 1.53e7 | 218 | 3.102e-7 | 5.79e-9 | 1.140e-7 | 3.69e-10 | 2.45e-8 |
| 2.53e7 | 296--304 | 3.102e-7 | 5.79e-9 | 1.140e-7 | ~2e-14 | 2.52e-8 |

abundanceは `XH` 規格化である。結論は「多数の全く異なるrootがある」のではなく、炭素のほぼ全量がCO、残った酸素の主要部分がO2へ入る境界近傍の一つの枝へ収束する、というものだった。

この結果をdeterministicな `iter_reset` guessへ反映した。

- 炭素をほぼCOへ置くguess。
- COに使われなかった酸素を主にO2へ置くguess。
- samplingで得た少量のCO2/H2Oを含むguess。
- それでも診断対象stepが失敗する場合だけ、compile-time diagnosticとして元素保存random restartを使える。

production runの思想は、同一のphysical timestepで複数guessを試すことである。random restartはroot探索と診断用であり、physical timestep halvingの代替ではない。

### 解析Jacobianと数値Jacobianへのfallback

平衡中間体をNewton変数から外した場合、その応答を無視するとJacobianが不完全になる。現在は複数中間体に一般化した連鎖微分を用い、代数的中間体と保存則従属speciesを縮約した解析Jacobianを作る。

解析Jacobianで20反復以降も `err_max > 1e-4` の場合は、全代数speciesを再構成しながら数値Jacobianへfallbackする。これは「一回でもboundを越えたら数値Jacobianへ移る」という意味ではない。解析近似で停滞した場合にだけ、より完全な微分を評価する保険である。

保存則へtrace speciesと全中間体を入れると、保存式自体は数個の加算でも、代数的連鎖微分が密になり速度が大幅に低下した。そのため現在は主要reservoirだけを保存則へ入れ、C/O保存則を外した診断modeも用意して、network solveと保存則縮約の問題を分離している。

### chemistry収束とenergy residualの結合

energy implicit modeでは、各trial温度について

```text
chemistry implicit -> chemical heating -> radiative cooling -> energy residual
```

を評価する。したがって、chemistryが異なるguess枝へ飛ぶ、bound処理後にrollbackする、または中間体平衡が未収束のまま返ると、temperatureだけを変えた一次元energy residualも不連続に見える。今回 `Z=10^-3, z=0` で見えたO/OHと温度の同時spikeは、この結合を示していた。

step halvingでは、同じ不適切なguessからboundへ近づくだけで、目的のchemical rootを保証しない。component-wise boundとCO/O2寄りのreset guessを戻した後、同じphysical timestepのままchemical warningと細かな温度ぎざぎざが消えた。したがって、この領域の主問題はenergy solver単体ではなく、energy residual内部で呼ばれるchemical nonlinear solveのbranch選択だった。

### molecular cooling on/offで変わる進化path

`MOLECULAR_LINE_COOLING=1/0` の二つをclean buildし、同じ36条件を完走させた。ここでoffにするのはCO/OH/H2O line coolingだけであり、C/O chemical networkそのものは両方で有効である。したがって、二つのrunのspecies差はnetworkの有無ではなく、coolingによる温度feedbackから生じる。

![molecular cooling on/off temperature paths](figures/molecular_cooling_on_off_paths.png)

`z=0` の代表値は次の通りだった。

| Z/Zsun | 最大 `|Ton-Toff|/Ton` | 最大差のn [cm^-3] | Ton [K] | Toff [K] | final `(Ton-Toff)/Ton` |
|---:|---:|---:|---:|---:|---:|
| 1e-3 | 1.557 | 2.53e7 | 112.8 | 288.4 | -8.38e-4 |
| 1e-2 | 0.240 | 1.62e5 | 51.8 | 64.2 | -7.68e-2 |
| 1e-1 | 0.714 | 8.08e3 | 25.8 | 44.3 | -1.57e-1 |
| 1 | 0.101 | 5.23e15 | 1403 | 1545 | -7.81e-2 |

`Z=10^-3` では `n ~ 1e7--1e8 cm^-3` の温度peakが最も大きく変わる。CO/OH/H2O coolingを入れると、offで約288 Kへ上がる枝が約113 Kへ抑えられる。一方、十分高密度ではdust/continuumと断熱加熱が支配的になり、二つのpathは再び近づく。

### C/O chemistry off + 固定CII/OI coolingとの比較

第三の比較として、C/O chemical networkを完全に外し、CIIとOIのfine-structure coolingだけを固定solar abundanceで残した。

```text
CARBON_OXYGEN_CHEMISTRY=0
MOLECULAR_LINE_COOLING=0
nCII = n_nuclei * XH * Y_C_sun * Z
nOI  = n_nuclei * XH * Y_O_sun * Z
nCI  = 0
```

現行sourceの `CARBON_OXYGEN_CHEMISTRY=0` 分岐は、CIとCIIの両方へ炭素全量を同時に与える。このままでは炭素coolingを二重に数えるため、今回の比較用一時buildだけでCIを0にした。本体sourceは変更していない。

上の温度軌道図には、この第三経路を緑点線で追加した。

完全C/O chemistry + molecular cooling onを基準にした比較は次の通りである。

| Z/Zsun | full C/O + mol. on 対 fixed CII/OI: 最大相対差 | n [cm^-3] | Tfull [K] | Tfixed [K] | final相対差 |
|---:|---:|---:|---:|---:|---:|
| 1e-3 | 0.361 | 6.63e2 | 74.7 | 47.7 | -8.73e-4 |
| 1e-2 | 1.755 | 3.29e2 | 7.39 | 20.4 | -7.68e-2 |
| 1e-1 | 14.50 | 1.00e-1 | 36.4 | 564 | -1.57e-1 |
| 1 | 4.007 | 2.03e-3 | 26.3 | 132 | -7.77e-2 |

相対差は `|Tfull-Tfixed|/Tfull`、final相対差は `(Tfull-Tfixed)/Tfull` である。

molecular coolingを切ったC/O chemistry runと固定CII/OI modelを比べると、`Z=10^-3` では `n=2.07e7 cm^-3` で、前者が277 K、後者が102 Kとなる。固定modelは炭素と酸素をCO/O2/H2Oへ移さないため、CII/OI coolingを中密度まで過大に残し、C/O chemistry runに現れる温度peakを強く抑える。

cooling成分の違いを次に示す。

![C/O chemistry cooling comparison](figures/carbon_chemistry_cooling_comparison_Zm3_z0.png)

左はC/O chemistryとmolecular coolingを両方有効、中は同じchemistryでmolecular coolingだけを無効、右はC/O chemistryを外した固定CII/OI modelである。固定modelではCII/OI abundanceが化学的に減らないため、そのcoolingがfull chemistryより広い密度域に残る。

ただし、この第三経路はcooling項だけを差し替えた完全なcontrolled experimentではない。C/O chemistryを切ると次も同時に変わる。

- CIIがcharge neutralityへ寄与しなくなる。
- metal electronがH/He chemistryとelectron collision coolingへfeedbackしなくなる。
- C/O reactionがH、H2、electron abundanceへ与える小さなfeedbackも消える。
- 全炭素・酸素を固定atomic reservoirとして扱い、CO/O2/H2Oへのdepletionを無視する。

特に `Z >= 10^-1` の低密度側で第三経路が大きく高温になるのは、単純にCO coolingがないことだけでは説明できない。固定CIIを仮定しながら、そのionizationに対応する電子をcharge neutralityへ加えていないlegacy近似の影響も含む。したがって第三経路は「旧来の固定metal cooling modelとの比較」と解釈し、C/O reactionだけの純粋な効果とはみなさない。

36条件surveyのwall timeは、同じ `-O0` buildで概ね次の通りだった。

| model | wall time |
|---|---:|
| C/O chemistry + molecular cooling | 34.3 s |
| C/O chemistry + molecular cooling off | 32.1 s |
| C/O chemistry off + fixed CII/OI | 5.0 s |

C/O networkと中間体平衡が計算時間の大部分を占め、molecular cooling table評価そのものの追加costは比較的小さい。

### 3モデルの計算時間配分と収束反復数

上の3モデルを `CHEMISTRY_PROFILE` 付きでclean buildし、同じ36条件についてcounterを集計した。絶対wall timeはprofileなしのrun、時間割合とcall数はprofile付きrunを使う。

![three-model performance comparison](figures/three_model_performance_comparison.png)

#### 全体時間とstep数

| 指標 | C/O + molecular | C/O + molecular off | C/O off + fixed CII/OI |
|---|---:|---:|---:|
| profileなしwall time [s] | 34.32 | 32.08 | 4.97 |
| fixed CII/OIに対する時間比 | 6.90 | 6.45 | 1.00 |
| collapse/chemistry step総数 | 54,972 | 54,972 | 54,972 |
| 1 case当たりcollapse step | 1,527 | 1,527 | 1,527 |
| energy residual評価総数 | 294,558 | 275,383 | 275,995 |
| energy trial / collapse step | 5.358 | 5.010 | 5.021 |
| Newton RHS評価総数 | 1,189,316 | 1,117,803 | 884,543 |
| Newton反復 / implicit trial（成功・失敗込み） | 4.038 | 4.059 | 3.205 |
| 全implicit trialのNewton反復中央値 | 4 | 4 | 4 |
| 全implicit trialの90 percentile | 5 | 5 | 4 |
| 全implicit trialの95 percentile | 9 | 9 | 4 |
| 全implicit trialの99 percentile | 24 | 25 | 4 |
| 全implicit trialの最大値 | 104 | 96 | 18 |
| 4反復以下で終了したtrial | 89.46% | 88.35% | 99.43% |
| resetを1回以上通過したtrial | 18,238 (6.19%) | 18,272 (6.64%) | 3 (0.001%) |
| reset trialが占める全Newton反復 | 23.73% | 25.18% | 0.005% |
| 70反復後のreset event | 0 | 0 | 0 |
| `chemistry_status=0` trial | 17,333 (5.88%) | 17,402 (6.32%) | 1 (<0.001%) |
| status=0 trialが消費したNewton反復 / trial | 14.88 | 14.82 | 13.00 |
| Newton反復 / collapse step | 21.635 | 20.334 | 16.091 |
| Newton反復 / case | 33,037 | 31,050 | 24,571 |
| 36条件の未収束warning | 0 | 0 | 0 |

collapse step数が全モデルで同じなのは、1-zoneのdensity evolution gridが同じだからである。モデル差は主に、一つのcollapse step内でenergy residualを何回評価するかと、一つのchemical trialを何回のNewton反復で解くかに現れる。

ここで「1回のchemistry」をどの階層で数えるかに注意が必要である。energy root探索中のcandidate temperatureを含む全implicit chemistry呼出しを直接計数しても、中央値は4反復だった。前stepのabundanceが初期値になり、1-zoneのcollapse stepが細かいため、C/Oモデルでも約89%のtrialは4反復以下で終了する。

ただし分布には長いtailがある。C/Oモデルの95 percentileは9、99 percentileは24、最大は96--104反復である。約6%のtrialがboundまたはalgebraic reconstructionを契機とするresetを通り、call数としては少数ながら全Newton反復の約24--25%を消費する。このsurveyでは `itr - 70*(flag_reset+1) > 0` による遅延resetは一度も発生せず、観測されたresetはすべてbound/algebraic failure直後の即時経路だった。したがって「resetがなければrunが止まる」ことと「全stepの中央値が4」であることは矛盾しない。runの成否は少数の難しいstepで決まるが、中央値は多数の容易なstepで決まるためである。

外側の一つのcollapse/chemistry stepではenergy root探索が温度を変えてimplicit solveを約5回呼ぶため、合計では20--22 Newton反復となる。またC/Oモデルではcandidate temperatureに対するtrialの約6%が `chemistry_status=0` で戻り、これらは平均約15反復を消費する。

#### 排他的な上位時間配分

次の4分類は重複しないように、chemistry totalから親子関係を差し引いて作った。

| 処理 | C/O + molecular | C/O + molecular off | fixed CII/OI |
|---|---:|---:|---:|
| implicit chemistry | 97.96% | 98.93% | 93.79% |
| radiative cooling | 1.69% | 0.72% | 4.56% |
| energy residual内のその他 | 0.21% | 0.21% | 1.04% |
| energy solver + driver overhead | 0.14% | 0.15% | 0.61% |

fixed CII/OIモデルでradiative coolingの割合が増えるのは、radiative coolingの絶対時間が増えたからではない。絶対時間は約0.228 sでmolecular-off C/O modelの約0.229 sとほぼ同じだが、implicit chemistryが大幅に軽くなったため割合が上がっている。

#### implicit chemistry内部のhot spot

以下は親処理を含むinclusive timeであり、例えばintermediate solveはJacobianやNewton stepの内部からも呼ばれる。そのため列を足して100%にはならない。

| 処理 | C/O + molecular | C/O + molecular off | fixed CII/OI |
|---|---:|---:|---:|
| Newton Jacobian | 54.34% | 55.06% | 41.34% |
| intermediate equilibrium | 27.07% | 27.57% | 10.82% |
| Newton step・algebraic再構成 | 13.47% | 13.50% | 8.95% |
| Newton RHS | 5.96% | 6.00% | 12.92% |
| reaction rate評価 | 2.11% | 2.04% | 10.47% |
| linear solve | 1.91% | 1.91% | 2.74% |
| molecular cooling table | 0.95% | 0% | 0% |
| dust温度 | 1.00% | 0.92% | 5.73% |

C/O chemistryを切るとNewton反復数は一trial当たり `4.06 -> 3.20`、約21%減る。しかしwall timeは `32.08 -> 4.97 s`、約6.45倍短くなる。したがって高速化の主因は反復数だけではなく、各反復で扱うspecies数、C/O reactionのJacobian、複数平衡中間体の代数的応答が消えることである。

molecular coolingを有効にすると、table評価自体は0.326 s、chemistry totalの約0.95%に過ぎない。一方でenergy residual評価が19,175回、Newton RHSが71,513回増え、wall timeは約6.96%増えた。追加costの大半はcooling関数そのものではなく、変化したthermal path上でenergy/chemistry solveのtrial数が増えたことによる。

再生成用ファイル:

- `figures/three_model_performance_comparison.dat`
- `figures/three_model_performance_comparison.gp`

追加の再生成用ファイル:

- `figures/carbon_chemistry_cooling_comparison_Zm3_z0.gp`
- `figures/no_co_chem_CII_OI_Z*_z0.dat`
- `figures/no_co_chem_CII_OI_Zm3_z0_cool.dat`

各cooling成分は次の図のようになる。off側でCO/OH/H2O曲線が消える一方、chemical networkは残り、dustと他のradiative coolingがenergy evolutionを担う。

![molecular cooling rate decomposition](figures/molecular_cooling_rates_Zm3_z0.png)

再生成用ファイル:

- `figures/molecular_cooling_on_off_paths.gp`
- `figures/molecular_cooling_rates_Zm3_z0.gp`
- `figures/molecular_on_Z*_z0.dat`
- `figures/molecular_off_Z*_z0.dat`

### C/O species進化

`Z=10^-3, z=0` のC/O species進化を、molecular cooling on/offで比較した。

![C/O species evolution](figures/carbon_oxygen_species_on_off_Zm3_z0.png)

共通する大きな流れは次の通りである。

1. 低密度ではC+とOが主要reservoirである。
2. 密度上昇とshielding・分子形成により、炭素はC/COへ、酸素はOH/H2O/O2へ再配分される。
3. `n ~ 1e4--1e7 cm^-3` でCOが主要な炭素reservoirになる。
4. 高密度では炭素のほぼ全量がCOへ入り、少量のCO2が残る。
5. 酸素はCOに使われた分を除き、主にO2とH2Oへ入る。

cooling on/offで最終的なreservoir構造は大きくは変わらないが、温度依存rateを通じて遷移密度とOH/H2Oの一時的な分配が変わる。このfeedbackが、温度peak付近で二つのthermal pathを分けている。

再生成用ファイル:

- `figures/carbon_oxygen_species_on_off_Zm3_z0.gp`
- `figures/molecular_on_Zm3_z0.dat`
- `figures/molecular_off_Zm3_z0.dat`

### CMB z=20低密度dust spike

`output_MT_HD_ME0_EIMP_CMBz20.dat` の低密度側では、gas温度のぎざぎざと同じ位置にdust coolingの孤立spikeがあった。

代表例:

| n [cm^-3] | 修正前dust cooling | 安定側の値 |
|---:|---:|---:|
| 2.60e-3 | 1.86e-4 | 3.33e-10 |
| 6.38e-2 | 2.02e-4 | 5.29e-9 |

低密度では真の平衡が

```text
Tdust - Tradiation ~ 1e-11 K
```

にあり、最初のNewton更新は正しい点へ到達していた。しかしdouble精度で達成できる相対residualが約 `1e-4` なのに、residual条件は `1e-5` だった。その直後にbracket幅条件が働き、Newton点を中点で上書きして `Tdust-Tradiation ~ 1e-7 K` の偽値を返す場合があった。`Lambda_dust` をgas energyへ使うと、この微小温度誤差が大きな放射cooling spikeへ増幅された。

一方、高密度では `Tdust ~ Tgas` となるため、`Lambda_coll` の `Tgas-Tdust` が桁落ちし、`Lambda_dust` の方が安定していた。したがって片方を全密度域で使うのではなく、次の条件数による切替を採用した。

```text
|Tgas - Tdust| >= |Tdust - Tradiation|
    -> Lambda_coll   （低密度側）
otherwise
    -> Lambda_dust   （高密度側）
```

正しく収束していれば両者は等しいため、切替自体は物理modelの変更ではない。有限精度でより大きく分離された温度差を使う側を選んでいる。

追加修正:

- `Tdust^4 - TCMB^4 - stellar_flux/(4 sigma)` の差を直接作らず、`Tradiation`を定義して `Tdust^4-Tradiation^4` を因数分解する。
- bracket幅で終了するとき、現在のNewton点を中点で上書きしない。
- 一度試した「全反復中の最小相対residual点を採用」は、high-density opacity transitionで別枝を選んだため撤回した。
- 温度幅許容値は以前の `max(1e-4 K, 1e-5 Tdust)` に戻した。

単一条件で再計算すると、低密度pathは次のように滑らかになった。

```text
n=1.00e-3  T=277.24 K
n=1.28e-3  T=175.40 K
n=1.66e-3  T=105.44 K
n=2.03e-3  T= 80.40 K
n=2.60e-3  T= 68.73 K
n=3.18e-3  T= 65.81 K
n=4.08e-3  T= 64.70 K
```

高密度側では `Tdust ~ Tgas` の枝を維持した。dust opacity tableでPlanck opacityが0になる高温sublimation領域では、`Lambda_dust`側がdust coolingを0へ落とす。このため、単純に常時 `Lambda_coll` を使うと、opacity上dustが消えている領域にもgas-dust collision coolingを残す不整合が起きる。

### 今回の診断から得た実装原則

1. chemical failureをphysical timestep halvingで隠さない。同じstepでguessを変えてrootを探す。
2. Newton candidateのbound越えを即failureとしない。各成分をbound内へ戻し、代数speciesを再構成して続行する。
3. 中間体平衡は解析的に縮約するが、停滞時だけ数値Jacobianを使う。
4. chemical heatingのpath依存補正を入れる前に、abundance solveとenergy residualの連続性を確保する。
5. 平衡式の左右の項が理論上等しくても、低密度・高密度の両極限で同じ側だけを数値的に使わない。
6. bracket幅が小さいことと、連続なrootへ収束したことを同一視しない。
7. molecular coolingのon/off比較ではchemical network自体を共通にし、temperature feedbackによるspecies pathの差を分離する。

今後、dust solverをさらに一般化する場合は

```text
Tdust = Tradiation + x (Tgas - Tradiation), 0 <= x <= 1
```

の `x` を解く形が候補になる。これなら低密度の `Tdust-Tradiation` と高密度の `Tgas-Tdust` を、それぞれ近い浮動小数点数の直接減算なしに構成できる。ただし現在のhybrid energy exchangeで今回の低密度・高密度の両問題は解消しているため、これは将来の整理項目とする。

## C/O chemistry追加後のJacobian・中間体処理の高速化

### 反応数が約2倍なのに計算時間が約6倍だった理由

C/O chemistryなしでは、平衡中間体はHMのみで、Newtonの独立変数は4、代数変数は4である。C/O chemistryありでは、平衡中間体がHMを含めて9種となり、Newtonの独立変数は12、代数変数も12となる。active reaction数はおよそ30から66への増加（約2.2倍）だが、次元増加により次の処理が同時に増えていた。

- 各Newton iterationにおけるfull reaction Jacobianの構築
- 9中間体の平衡解を求める9元連立方程式
- 9中間体の独立変数に対する応答を求める、多右辺の9元連立方程式
- 12独立変数と12代数変数の縮約（Schur update）
- 12元Newton行列の線形解法

したがって、Jacobian評価は「反応を足すだけ」では終わらない。反応由来の微分を構築した後に、平衡中間体と保存則による連鎖微分を独立変数空間へ縮約する費用が支配的だった。

最適化前の `-O0` profile（C/O chemistry + molecular cooling、36条件）では、全体35.014 sに対して以下を要していた。各項には包含関係があるので割合を単純加算してはならない。

| 処理 | 時間 [s] | 全体比 |
|---|---:|---:|
| implicit chemistry | 34.309 | 98.0% |
| Newton Jacobian全体 | 19.091 | 54.5% |
| solver/algebraic column生成 | 6.471 | 18.5% |
| 中間体平衡solve全体 | 9.533 | 27.2% |
| 中間体式assembly | 8.169 | 23.3% |
| 中間体linear solve | 1.007 | 2.9% |
| algebraic response | 0.476 | 1.4% |
| 中間体response solve | 0.552 | 1.6% |
| dense Schur update | 7.184 | 20.5% |
| main Newton linear solve | 0.664 | 1.9% |

反応率評価そのものはC/Oなしに対して約1.25倍に留まっていた。主因はrate coefficientではなく、平衡中間体を含むJacobian縮約である。

### 採用した変更

1. **full reaction Jacobianを反応ごとの1回の走査で構築**

   以前は微分するspeciesごとにreaction listを走査していた。各反応のrate積を一度評価し、その反応に現れる全reactant columnへ寄与を加える形へ変更した。

2. **中間体平衡系だけを直接assembly**

   full RHSと複数のfull Jacobian columnを構築せず、各中間体に影響するreaction indexから9本の式と9x9行列だけを作る。反応の加算順序は従来と同じに保った。

3. **Schur updateを構造的な非零項だけに限定**

   solver equationごとに、実際にその行へ寄与する代数speciesだけを抽出した。代数変数の応答自体は変えず、明らかなゼロとの積和だけを除いた。

4. **production buildの既定最適化を`-O2`へ変更**

   `OPTIMIZATION ?= -O2` とし、比較時は `OPTIMIZATION=-O0` で従来条件を再現できる。

`-O0`で構造変更だけを比較した結果は次の通りである。

| 項目 | 変更前 [s] | 変更後 [s] | 変化 |
|---|---:|---:|---:|
| profile total | 35.014 | 25.003 | -28.6% |
| Newton Jacobian | 19.091 | 11.624 | -39.1% |
| 中間体solve | 9.533 | 6.481 | -32.0% |
| 中間体assembly | 8.169 | 5.107 | -37.5% |
| Schur update | 7.184 | 2.366 | -67.1% |

36個の標準出力 `.dat` を比較し、`-O0`では出力精度10桁の範囲で温度進化が完全一致した。

### 採用しなかった変更

9中間体の構造依存関係からscalar/2x2 blockへ分解する一般的なblock solveも試した。しかし、小行列管理のoverheadにより中間体response solveが約0.56 sから約2.69 sへ悪化した。さらにhigh-metallicity, z=5の一部で温度が最大約26.7%変化し、有限精度で異なる平衡枝を選択した。この実装は完全に撤回し、現在は従来通りdense 9x9 LAPACK solveを用いている。

### 最終的な3モデル比較

36条件を連続実行した最終`-O2` buildの結果である。wall timeはprofile instrumentationなし、内訳はprofile buildで別に測定した。

| model | wall time [s] | profile total [s] | implicit chemistry [s] | intermediate solve [s] | Newton Jacobian [s] |
|---|---:|---:|---:|---:|---:|
| C/O chemistry + molecular cooling | 7.614 | 7.871 | 7.471 | 2.347 | 3.373 |
| C/O chemistry, molecular cooling off | 6.909 | 7.324 | 7.118 | 2.254 | 3.246 |
| CO chemistry off, CII/OI coolingのみ | 1.652 | 1.678 | 1.504 | 0.150 | 0.570 |

最適化前のprofileなし`-O0` wall timeは、それぞれ34.32 s、32.08 s、4.97 sだった。したがって最終buildではC/O + molecular modelが約4.51倍、molecular offが約4.64倍、CO chemistry off modelが約3.01倍高速になった。

最終状態でもC/O + molecular modelはCO chemistry off modelの約4.61倍、molecular offは約4.18倍である。molecular line cooling自体の追加費用は約10.2%であり、残りの大半はC/O networkと9中間体の平衡縮約に由来する。C/O + molecular profileでは、中間体平衡solveが全体の29.8%、Newton Jacobianが42.9%を占めた。呼び出し回数は中間体平衡solveが約284.5万回、中間体response solveが約116.9万回であり、小さな連立方程式でも累積費用は大きい。

最終`-O2`と構造変更後`-O0`の36条件比較では、温度の最大相対差は `3.1e-8` で、warningは発生しなかった。この差はcompiler optimizationによる浮動小数点演算順序の差の範囲であり、進化pathは一致している。

今後さらに高速化する場合の主対象は、reaction rate式ではなく、同一評価内での中間体行列factorizationの再利用、または現在の9x9系に特化した低overhead solverである。ただし平衡枝を変えないことが必須であり、一般的なblock分解は今回の検証では不採用となった。

### 中間体block分解の再検証

前回のblock分解で解が変わった原因を確認するため、reaction tableから中間体依存graphを構築し、dense解とblock解を全呼び出しで直接比較した。現在のnetworkの構造は次の通りだった。

```text
HM                          scalar
CO+                         scalar
CH <-> CH2                  2x2 coupled block
O+ -> OH+ -> H2O+ -> H3O+   directed chain
O2+                         scalar
```

この構造は中間体について線形な任意のreaction tableから機械的に抽出できる。ただし、directed chainをscalar SCCへ分解した前進代入は数値的に安全ではなかった。低電子存在量領域ではH3O+ equationの対角成分が実測で約 `4e-31` まで小さくなり、上流のH2O+項と右辺のほぼ完全な相殺からH3O+ responseが決まっていた。

代表的なJacobian response equationは概略

```text
5.50e-3 d(H2O+)/dx - 4.05e-31 d(H3O+)/dx = -2.11e-9
```

となっていた。この場合、scalar前進代入とdense LUはどちらも小さなresidualを持ち得るが、丸め誤差により `d(H3O+)/dx` がfactor 2程度異なる例があった。これはgraphの取り違えではなく、平衡中間体response系そのものの極端な悪条件性による。

そこで、無向に結合したreaction familyを一つのpivot可能blockとして

```text
{HM}, {CO+}, {CH,CH2}, {O+,OH+,H2O+,H3O+}, {O2+}
```

に分ける方法も検証した。しかし4x4 blockとfull 9x9 LAPACKでは悪条件系における丸め順序が異なり、多右辺responseの一致は得られなかった。

多右辺responseを従来のdense 9x9に残し、平衡abundanceの1右辺solveだけをblock化する構成では、36条件について最大相対温度差は `4.2e-7` まで抑えられた。ただし極微量中間体には数%の相対差が残った。また同一環境でのprofileなしwall timeは

| solve | wall time [s] |
|---|---:|
| dense 9x9 | 7.322 |
| abundanceのみblock | 7.414 |

となり、block版の方が約1.3%遅かった。小blockを複数回呼ぶoverheadが、Accelerate/LAPACKによる9x9 solve一回を上回ったためである。このためblock solver実装は撤回し、production codeはdense 9x9へ戻した。

今後の安全な候補はsparsity分解ではなく、同じ9x9 factorizationを平衡abundance solveと直後のresponse solveで再利用することである。これならglobal pivotingと現在の数値解を維持できる。ただし、両地点で構築される行列がbitwiseに一致することを先に確認する必要がある。もう一つの候補は中間体変数・平衡式のscalingだが、H3O+ responseの数値解自体を変える可能性があるため、性能変更ではなく数値安定性の問題として別に検証する。
