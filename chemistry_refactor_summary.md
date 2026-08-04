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

gas の energy exchange には、平衡探索終了時の `Lambda_coll` を使用する。低密度で `Tdust` が CMB 温度に極めて近い場合、有限の温度誤差を持つ `Lambda_dust` をgas heatingへ流用すると、非物理的な加熱が生じるためである。

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
4. gasのenergy exchangeには `-Lambda_coll` を使用し、放射residualをgasへ流さない。

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

最終のclean buildによる36条件surveyでは、通常species modeが約4.67秒、中間体平衡modeが約5.23秒だった。この二値は異なるsolver modeなので、互いのdust最適化率を表すものではない。

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
