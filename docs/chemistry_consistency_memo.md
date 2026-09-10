# chemistry / chemistry_old consistency memo

## 目的

`chemistry_old/` と新しい `chemistry/` の 1ZONE 計算を、まずは

- 外部 LW/J21 なし
- 0 metal
- HD あり
- energy implicit / chemistry implicit

の条件で比較した。

この条件では本来、主要な熱・化学進化はほぼ一致するはずなので、差の原因を切り分けた。

## 現状の主な変更

### ビルド・分岐

- `Makefile` に `USE_CHEMISTRY_OLD` と `ENABLE_LW` を追加した。
- `main.c` は `CHEMISTRY_OLD` の有無で old/new chemistry を切り替える。
- new 側では `chemistry_adapter.c/h` を通して、既存の `SphP` scalar fields と `chemdata_in/out` を接続する。

### adapter

- `sphp_p.elec` などの scalar abundance と `chemdata_out.Abund[]` の相互変換は `chemistry/` の外、`chemistry_adapter.c` に置いた。
- `ThermalEnergy` は chemistry 呼び出し前に entropy から再構成する。
- `molecular_weight` は `Abund[]` から計算する。
- `chemdata_in` に外側で更新済みの `SphP[i].Jeans_Length` を渡すようにした。
- 出力用の `tau_cont` 再計算も、new 側で `get_core_radius()` を呼び直さず、外側由来の `SphP[i].Jeans_Length` を使うようにした。

### D/HD abundance sum

- `timescales.c` の分子量計算で、`XHD` 定数ではなく、
  `DI + DII + HDI (+ HDII / DM)` を足すように old/new の方針を揃えた。
- ただし、この変更の影響は非常に小さく、今回の主な差の原因ではなかった。

### chemistry 内部の整合

- `chemical_cooling()` で `nH = XH * rho_gas / PROTONMASS` に修正した。
  これは初期の数 % 差の主因だった。
- energy implicit / Brent solver の trial 後に、`chemdata_out` だけでなく `chemdata_state` も restore するようにした。
  これにより、trial 温度で更新された `tau_cont` や `Jeans_Length` が次の residual 評価へ漏れなくなった。
- `Gamma` 評価は old の legacy な状態更新順に寄せた。
  具体的には、residual 内で `T_gas_new` そのものではなく、現在の `Utherm/Gamma/Abund` から作る温度で `compute_gamma_from_abund()` を呼ぶようにした。
- HD explicit branch で、old 側にある `k[50]` の DI destruction / DII creation 項を new 側にも入れた。
  今回の implicit run では実質的な差は出なかったが、整合性修正として必要。

## Jeans_Length について

今回いちばん議論が必要なのは `Jeans_Length` の扱い。

old 側では、1 timestep ごとに外側の collapse/update 処理で `SphP[i].Jeans_Length` が更新される。その後、chemistry に入る時点では、その外側で確定済みの値が基準 state になる。

一方で new 側は、chemistry 内部や adapter 出力側で `get_core_radius(rho, T)` を呼び直していた。energy implicit の trial 温度ごとに length が変わり得るため、old の基準 state とずれた。

この差は高密度で `tau_cont` に直接効く。

```c
tau_cont = opacity(T, rho, Jeans_Length);
f_reduce = 1. / (1. + 3./4. * tau_cont * tau_cont);
edot_H2   *= exp(-tau_cont);
edot_HD   *= exp(-tau_cont);
edot_cont *= f_reduce;
```

trace では、修正前の new は高密度で内部 `tau_cont` が old より約 9% 小さかった。その結果、冷却項がずれて、最終温度で約 1% の差まで積み上がった。

old/new consistency check としては、外側由来の `Jeans_Length` を使うのが正しい。

ただし、物理モデルとしては new 側の「implicit residual の候補温度に応じて `Jeans_Length` を再評価する」方が自然に見える面もある。したがってこれは「new が間違っていた」というより、

- old と同じ計算を再現する mode
- new の物理的に self-consistent な mode

を分けて考えるべき箇所。

現状の変更は、まず old との consistency を優先して old-compatible な基準に揃えている。

## 比較結果

最終確認は current source から old/new を再ビルドして実行した。

```text
old: USE_CHEMISTRY_OLD=1 ENABLE_LW=0
new: USE_CHEMISTRY_OLD=0 ENABLE_LW=0
```

`output_MT_HD_EIMP_CMBz0.dat` の比較:

![rho-T comparison](figures/rho_T_comparison_nolw_0metal.png)

図は `figures/rho_T_comparison_nolw_0metal.gp` で生成した。
プロット用の比較データは `figures/rho_T_comparison_nolw_0metal.dat` に置いた。

```text
T      max rel 9.295196e-05, final rel 9.295196e-05
Gamma  max rel 7.379382e-06, final rel 6.095914e-06
mu     max rel 3.878409e-05, final rel 3.878409e-05
tau    max rel 2.902560e-05, final rel 9.702490e-06
H2I    max rel 4.914973e-05, final rel 4.876332e-05
HDI    max rel 3.952224e-05, final rel 3.891933e-05
```

`elec/HM/DII` は相対差が大きく見えるが、絶対値は非常に小さい。

```text
elec ~ 1e-12
HM   ~ 1e-17
DII  ~ 1e-17
```

主要な熱進化、H2/HD abundance、`mu`、`tau_cont` にはほぼ効いていない。

## 今後の方針メモ

1. old との consistency check では、現在の old-compatible `Jeans_Length` 扱いを使う。
2. 物理的に new の self-consistent residual を採用するかは別途検討する。
3. その場合は `Jeans_Length` 再評価の有無を compile-time flag または runtime option として分けるのがよい。
4. 今回の 0 metal / no LW が十分揃ったので、次は LW あり、metal ありなど、差が物理的に出やすい設定で比較する。

## LW あり、0 metal の比較

`ENABLE_LW=1` で old/new を再ビルドし、`SURVEY_LYMANWERNER` の出力から
`J21 = 0.1, 1, 10, 100, 1000` を比較した。

この比較中に、`set_J21()` の呼び出し順の問題を見つけた。
`set_J21()` は `All.J21` から background rate を作るが、もともとは
`All.J21 = All.LymanWerner * 1.e-21` が `initialize_cells()` 内で設定されていた。
そのため survey loop では、出力ファイル名の J21 と実際の background rate が 1 ケース分ずれていた。

このため、survey loop 内で `set_J21()` を呼ぶ前に `All.J21` を更新するようにした。

```c
All.LymanWerner = LW[il];
All.J21 = All.LymanWerner * 1.e-21;
set_J21();
```

修正後の比較図:

![rho-T LW comparison](figures/rho_T_comparison_lw_0metal.png)

再生成用:

- `figures/rho_T_comparison_lw_0metal.gp`
- `figures/rho_T_comparison_J21_*_0metal.dat`

比較 summary:

| J21 | max dT/T | final dT/T | final T new | final T old | final dH2I/H2I | final H2I new | final H2I old | final dtau/tau | final tau new | final tau old |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 9.659046e-03 | 9.659046e-03 | 2.119271e+03 | 2.098997e+03 | -4.724193e-03 | 4.794384e-01 | 4.817141e-01 | 9.626194e-02 | 2.033508e+00 | 1.854947e+00 |
| 1 | -1.334370e-02 | 9.640831e-03 | 2.119271e+03 | 2.099035e+03 | -4.715772e-03 | 4.794384e-01 | 4.817100e-01 | 9.626003e-02 | 2.033508e+00 | 1.854951e+00 |
| 10 | 2.037161e-02 | 2.037161e-02 | 3.557831e+03 | 3.486799e+03 | -2.458111e-01 | 1.236266e-02 | 1.639200e-02 | 5.505023e-01 | 1.396763e+01 | 9.008453e+00 |
| 100 | 2.037160e-02 | 2.037160e-02 | 3.557831e+03 | 3.486799e+03 | -2.458110e-01 | 1.236266e-02 | 1.639200e-02 | 5.505020e-01 | 1.396763e+01 | 9.008455e+00 |
| 1000 | 2.037158e-02 | 2.037158e-02 | 3.557831e+03 | 3.486799e+03 | -2.458108e-01 | 1.236267e-02 | 1.639199e-02 | 5.505013e-01 | 1.396762e+01 | 9.008456e+00 |

現時点の観察:

- `J21=0.1, 1` では high-density 側で最終温度差が約 1%。
- `J21>=10` では強 LW 側に入り、最終温度差は約 2%、H2/HD abundance 差は約 25%。
- `J21=10, 100, 1000` はほぼ同じ軌道になっており、この設定では強 LW 側で飽和しているように見える。
- 差は主に high-density 側の `tau_cont` 差と結びついている。`J21>=10` では final `tau_cont` が new で約 14、old で約 9 になっている。

次に調べるべき候補:

1. LW shielding / photodissociation rate の old/new mapping。
2. `All.rate_Hminus_pd`, `All.rate_H2_diss`, `All.J21` が chemistry 内部で同じ値として入っているか。
3. 強 LW 側で H2 が減った後の `tau_cont` 差が、abundance 差由来なのか、opacity/state 更新順由来なのか。

## J21=1000 のずれの原因

主因は LW shielding column に使う `Jeans_Length` の扱いだった。

old 側では `SHIELD_LW` 有効時、`compute_abundances_d()` の入口で
`set_core_radius(ithis, T_gas)` を呼び、その時点の `T_gas` から
`SphP[ithis].Jeans_Length` を更新する。その後、同じ値が H2/HI/HD の
column density と `set_reaction_rate_coeff()` 内の continuum opacity に使われる。

new 側では adapter が外側由来の `SphP[i].Jeans_Length` を
`chemdata_in.Jeans_Length` として渡していた。この値は chemistry 入口の
current `T_gas` ではなく、外側の更新順に依存した値なので、LW shielding
column が old とずれる。`J21=1000` では H2 が非常に小さい期間に
self-shielding が効き始めるため、この lag が H2/HD の履歴に積み上がり、
最終的に温度で約 2%、H2/HD で約 25% の差になっていた。

対応として、`chemistry/` は触らず、adapter 側で `SHIELD_LW` のときだけ
old と同じ挙動にした。具体的には `chemdata_in.Jeans_Length = 0` を渡し、
new chemistry 内部の `get_core_radius(rho_cgs, T_gas)` で current `T_gas`
から再計算させる。さらに出力用の `tau_cont` 再計算でも同じ current-T
Jeans length を使い、`SphP[i].Jeans_Length` へ反映する。

`SHIELD_LW` なしでは old はこの入口再計算をしないので、no-LW では従来通り
外側由来の `SphP[i].Jeans_Length` を渡す。これにより no-LW の一致は保たれる。

修正後の `J21=0.1, 1, 10, 100, 1000` 比較:

| J21 | max dT/T | final dT/T | max dH2/H2 | final dH2/H2 | max dtau/tau | final dtau/tau |
|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 1.045704e-04 | 1.045704e-04 | -5.045119e-05 | -4.904363e-05 | -2.382958e-02 | -8.374243e-03 |
| 1 | 8.653152e-05 | 8.653152e-05 | -4.330200e-05 | -4.057927e-05 | -4.839600e-02 | -8.375971e-03 |
| 10 | 1.282787e-04 | 1.639326e-07 | 5.744958e-03 | -2.311188e-06 | -1.834471e-02 | -1.287514e-02 |
| 100 | 1.280482e-04 | 1.307216e-07 | 5.743714e-03 | -1.843522e-06 | -1.834444e-02 | -1.287570e-02 |
| 1000 | 1.268874e-04 | 1.335035e-07 | 5.723324e-03 | -1.883419e-06 | -1.834197e-02 | -1.287565e-02 |

`J21=1000` の最終値では、温度・gamma・mu・H2/HD はほぼ一致した。
`tau_cont` は最大で約 1.8%、最終で約 1.3% 残るが、これは熱化学進化には
大きく残っていない。必要なら次に `calc_continuum_opacities()` への入力
(`rho_cgs`, `T_gas`, `l_core`, abundance) を old/new で trace して詰める。

## J21=0、metallicity sweep

`J21=0`、`DUST_COOLING`、`H2_FORMATION_ON_DUST`、`METAL_COOLING_LW` を有効にし、
`Z/Zsun = 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 1` を old/new で比較した。

出力と再生成用ファイル:

- `figures/rho_T_comparison_metal_J0.png`
- `figures/rho_T_comparison_metal_J0.gp`
- `figures/rho_T_comparison_Z_*_J0.dat`
- `figures/metal_j0_comparison_summary.md`

![rho-T metal J0 comparison](figures/rho_T_comparison_metal_J0.png)

比較 summary:

| log10 Z | max dT/T | n at max | final dT/T | final T new | final T old | max dH2/H2 | final dH2/H2 | max dTdust/Tdust | final dTdust/Tdust | final dtau/tau |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -5 | 5.463226e-08 | 5.233e+15 | 4.930502e-08 | 1.133759e+03 | 1.133759e+03 | -1.611705e-09 | 0.000000e+00 | 5.369407e-08 | 4.931829e-08 | 3.443874e-09 |
| -4 | 2.429434e-07 | 8.207e+15 | 2.027559e-07 | 1.438183e+03 | 1.438183e+03 | 9.256893e-10 | -4.400514e-10 | 2.431646e-07 | 2.020605e-07 | 3.965419e-07 |
| -3 | 9.911646e-08 | 6.538e+07 | -4.237962e-09 | 1.368582e+03 | 1.368582e+03 | -7.181263e-09 | 0.000000e+00 | 5.618502e-08 | -4.091826e-09 | 3.157659e-08 |
| -2 | -1.315500e-07 | 6.392e+15 | 7.492909e-08 | 1.480066e+03 | 1.480066e+03 | -2.800550e-10 | -2.800550e-10 | -1.310370e-07 | 7.553717e-08 | 1.740013e-07 |
| -1 | 2.442720e-07 | 8.207e+15 | 2.095956e-07 | 1.567304e+03 | 1.567303e+03 | -1.921009e-09 | -1.921009e-09 | -3.594709e-07 | 2.091490e-07 | 4.096579e-07 |
| 0 | 2.487662e-07 | 8.207e+15 | 2.251407e-07 | 1.632757e+03 | 1.632757e+03 | -3.843794e-09 | -3.843794e-09 | 2.486397e-07 | 2.247732e-07 | 3.828686e-07 |

注: `dT/T`, `dH2/H2`, `dTdust/Tdust`, `dtau/tau` は `(new - old_fixed) / new`。

修正内容:

- old 側の `dust_temperature()` で `get_dust_opacities(T_dust, &kappa_R, &kappa_P)` と呼ぶように修正した。関数の引数順は `Rossland, Planck` なので、修正前は Planck/Rossland opacity が入れ替わっていた。
- old 側で continuum/dust opacity 後に `sph_temp.f_reduce` だけを更新していた箇所を、出力が読む `SphP[ithis].f_reduce` にも戻すようにした。
- new 側では chemistry 本体内で dust 込みの `tau_cont` が計算されるが、adapter の `out2particle_chemdata_1zone()` 後に `calc_continuum_opacities()` だけで `SphP[i].tau_cont` を再計算しており、dust tau を落としていた。`chemistry_adapter.c` 側で `DUST_COOLING` 時は `dust_temperature()` まで呼ぶようにして、出力用 `SphP[i].tau_cont` も chemistry 内部と同じ形にそろえた。
- new 側の `chemistry/chemistry.h` に `T_CMB0=2.728` の fallback があり、`allvars.h` の `T_CMB0=2.73` より先に効いていた。`allvars.h` の共通 `T_CMB0` 定義を chemistry include より前に移して、old/new とも `T_CMB0=2.73` にそろえた。

結果:

- 以前の 25--35% 程度の温度差は dust opacity 順の修正で消えた。
- その後に残っていた `Z/Zsun=1` 低密度側の `1.25e-3` 程度の差は `T_CMB0=2.728` vs `2.73` の不一致が原因だった。
- `T_CMB0` を共通化した後は、全 metallicity で最大 `|dT/T|` は `2.5e-7` 程度まで落ちた。
- `tau_cont`, `f_reduce`, `Tdust`, `H2I` も final relative difference はおおむね `1e-7` 以下から数 `1e-7` 程度。

## metallicity + LW sweep

`DUST_COOLING`、`H2_FORMATION_ON_DUST`、`METAL_COOLING_LW` を有効にし、
`Z/Zsun = 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 1` と
`J21 = 0.1, 1, 10, 100, 1000` の 30 ケースを old/new で比較した。

最初の比較では、metallicity が入った LW 有りケースで大きな差が出た。
原因は photoelectric heating の `G0` の受け渡しだった。

- `set_J21()` で `All.G0` は J21 で規格化済みなのに、初期化で
  `SphP[i].G0 = All.J21 * All.G0` としていた。old 側の
  `dust_PEheating()` は `SphP[ithis].G0` を読むため、old の PE heating が
  実質ゼロになっていた。
- new 側の `chemistry/react_rate_FUV.c` は `calc_G0_rate()` の中で既に
  `/1.6e-3` している。一方 old 側の `chemistry_old/stellar_spectrum.c` は
  raw flux を返し、caller 側で `/1.6e-3` していた。共通の `main.c:set_J21()`
  で常に `/1.6e-3` すると、new 側だけ `All.G0` がさらに 625 倍になる。

対応:

- `initialize_cells()` では `SphP[i].G0 = All.G0` とする。
- `main.c:set_J21()` では `CHEMISTRY_OLD` のときだけ
  `All.G0 = Norm * g0 / 1.6e-3`、new 側では `All.G0 = Norm * g0` とする。

修正後の heatmap:

![metallicity and LW heatmap](figures/metal_lw_heatmap_g0unit.png)

出力と再生成用ファイル:

- `figures/metal_lw_heatmap_g0unit.png`
- `figures/metal_lw_heatmap_g0unit.gp`
- `figures/metal_lw_heatmap_g0unit.dat`

J21 ごとの最悪ケース:

| J21 | worst log10 Z | max \|dT/T\| | n at max | final dT/T | final T new | final T old | max \|dH2/H2\| | max \|dtau/tau\| |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | -2 | 2.732312e-06 | 1.277e+01 | 8.595571e-08 | 1.455401e+03 | 1.455401e+03 | 5.974527e-06 | 1.603181e-01 |
| 1 | -3 | 5.961682e-06 | 6.325e+01 | 6.390374e-09 | 1.345774e+03 | 1.345774e+03 | 1.395821e-05 | 2.946482e-01 |
| 10 | -3 | 1.980585e-04 | 5.164e+02 | 6.166906e-09 | 1.345894e+03 | 1.345894e+03 | 3.182410e-04 | 2.630368e-01 |
| 100 | -3 | 4.126456e-04 | 3.275e+04 | 5.568202e-09 | 1.346934e+03 | 1.346934e+03 | 1.173324e-03 | 3.545509e-01 |
| 1000 | -3 | 3.414378e-04 | 3.999e+04 | 3.911899e-09 | 1.354841e+03 | 1.354841e+03 | 4.504686e-03 | 2.432786e-01 |

`max |dT/T|` の全 grid:

| J21 \ log10 Z | -5 | -4 | -3 | -2 | -1 | 0 |
|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 5.13e-08 | 1.75e-07 | 1.25e-07 | 2.73e-06 | 2.40e-07 | 2.70e-07 |
| 1 | 2.40e-07 | 9.88e-07 | 5.96e-06 | 4.78e-06 | 1.53e-06 | 4.91e-06 |
| 10 | 1.29e-04 | 1.40e-04 | 1.98e-04 | 3.35e-06 | 4.90e-07 | 5.57e-05 |
| 100 | 1.29e-04 | 1.40e-04 | 4.13e-04 | 2.87e-06 | 3.96e-06 | 4.98e-05 |
| 1000 | 1.28e-04 | 1.39e-04 | 3.41e-04 | 1.59e-06 | 2.34e-05 | 6.46e-05 |

まとめ:

- `G0` の受け渡しをそろえる前の大きな差は PE heating の不一致が原因。
- 修正後、全 30 ケースで `max |dT/T|` は最大 `4.13e-4`。
- final の温度差は、ほとんどのケースで `1e-7` 以下。最大でも
  `6.46e-5` 程度。
- `tau_cont` は relative difference としては数十% 出る場所があるが、
  温度進化への影響は最終的に小さい。次に詰めるなら、`tau_cont` の残差を
  `calc_continuum_opacities()` と dust opacity の入力値ごとに trace する。

## compute_abundances_d / implicit solver cleanup

`chemistry/chemistry.c` と `chemistry_old/chemistry_200601.c` の内部契約を整理した。

変更点:

- `radiative_cooling_d_nometal()` の入口で `set_reaction_rate_coeff()` を呼ぶようにした。これにより cooling 単独呼び出しでも `Jeans_Length`, self-shielding, `tau_cont`, `f_reduce`, dust 関連 state が古いまま使われない。
- old 側にも `update_radiation_state_old()` を追加し、`compute_abundances_d()` 入口にあった LW shielding 更新を `set_reaction_rate_coeff()` 側へ寄せた。
- subcycle の上限 `time3` を毎回 `dt_tot` に戻さず、残り時間 `a_end - a_now` から計算するようにした。
- chemistry の失敗 return を cooling の return で上書きしないようにした。
- implicit chemistry が非収束のときは abundance を書き戻さず、旧 state に戻して失敗 return するようにした。
- `switch_func[]` の analytic Jacobian dispatch table は implicit solve ごとではなく static one-time init にした。
- `solve_brent()` は実際には無条件 bisection だったため、死んだ secant / inverse quadratic branch を消して `solve_bisection()` に改名した。アルゴリズム自体は変えていない。
- `give_g()` は trial evaluation ごとに外側 state を破壊しないよう、local copy 上で chemistry/cooling を評価し、収束後の最後だけ `commit_g()` で反映する形にした。
- `CHEMISTRY_IMPLICIT` で `CHEM_DIFF_ANALYTIC` が無効な build は、明示的に `ALLOW_CHEM_DIFF_NUMERIC` を付けない限り compile error にした。accidental slow build を防ぐため。
- `ENERGY_IMPLICIT` の inner loop で、計算済みの `dT` を使わず `T_new *= 1.001` していた分岐を `T_new_temp + dT` に修正した。

検証:

- build:
  - new: `make -C /Users/j/Desktop/codex/1ZONE USE_CHEMISTRY_OLD=0 ENABLE_LW=1 EXTRA_OPT='-DSURVEY_METALLICITY -DMETALLICITY=0 -DMETAL_COOLING_LW -DDUST_COOLING -DH2_FORMATION_ON_DUST -DHAVE_MPI'`
  - old: `make -C /Users/j/Desktop/codex/1ZONE USE_CHEMISTRY_OLD=1 ENABLE_LW=1 EXTRA_OPT='-DSURVEY_METALLICITY -DMETALLICITY=0 -DMETAL_COOLING_LW -DDUST_COOLING -DH2_FORMATION_ON_DUST -DHAVE_MPI'`
- run dirs:
  - new: `/private/tmp/1zone-run-new-cleanup`
  - old: `/private/tmp/1zone-run-old-cleanup`
- 30 ケース (`Z/Zsun = 10^-5 ... 1`, `J21 = 0.1 ... 1000`) は new/old とも完走。

修正後の worst:

| quantity | worst max relative diff | final relative diff | case |
|---|---:|---:|---|
| T | 3.916437e-04 | 6.615609e-09 | `output_MT_HD_ME-3_EIMP_B_100_CMBz0.dat` |
| tau | 1.981535e-03 | 9.289538e-09 | `output_MT_HD_ME-3_EIMP_B_1000_CMBz0.dat` |
| H2 | 5.833666e-03 | -2.000148e-10 | `output_MT_HD_ME-4_EIMP_B_1000_CMBz0.dat` |
| Tdust | 8.212853e-05 | 6.912939e-09 | `output_MT_HD_ME-3_EIMP_B_100_CMBz0.dat` |

今回の cleanup 後も、温度の最大差は `4e-4` 未満、final 温度差は worst case でも `1e-7` より十分小さい。

### 追記: Jeans_Length の受け渡し廃止

`chemdata_in` で `Jeans_Length` を渡す設計はやめた。new chemistry 側では
`update_radiation_state()` が常に
`get_core_radius(rho_cgs, T_gas)` で `sph_state->Jeans_Length` を更新する。

対応:

- `chemdata_in.Jeans_Length` を削除した。
- `chemistry_adapter.c` の `USE_OUTER_JEANS_LENGTH_IN_CHEM` 分岐を削除した。
- adapter 側の出力用 opacity 再計算でも、常に現在の `rho_cgs, T_gas` から
  `get_core_radius()` を呼び、その値を `SphP[i].Jeans_Length` に反映する。

これで new 側の関数契約は、「Jeans length は外側から渡される入力ではなく、
chemistry/cooling の評価温度から内部で決まる state」として一本化された。
