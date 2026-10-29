# isotope_fw_metric_analysis(日本語版)

[English](README.md)

IsoFresh データベースの淡水食物網(同位体が測定された分類群が 10 以上ある 606 件)を使い、緯度、気温、気温の季節性が食物網の構造をどう変えるかを、安定同位体の類似度ネットワークで解析したリポジトリです。

## 主な結果

- 極に近いほど、結合度、入れ子度、リンク強度が低下する(標準化効果 −0.2〜−0.3)。止水・流水に共通で、サンプリング努力(測定サンプル数)を調整しても変わらない。
- 結合度の低下の約半分は気候の直接効果で、残りは種数と餌ニッチ面積を介した効果である(piecewise SEM)。入れ子度とリンク強度は、主に種数を介して低下する。
- 結合度は 1 ℃あたり約 0.004 上昇する。この向きと大きさは、淡水魚の活性化エネルギー(Rall et al. 2012)を用いた Petchey, Brose & Rall (2010) のモデルの予測と一致する。
- モジュラリティ(Q)は極ほど高い。標準化モジュラリティ(z-Q)に全体の傾向がないのは、種数を介した低下と、気候の直接の上昇が打ち消し合うためである。
- 気温の季節性(月平均気温の標準偏差;夏と冬の気温差の大きさ)は、結合度への効果が最も大きく、年平均気温とは独立に効いている。

詳しい内容は `reports/report_concise.md` と `reports/IsoFresh_Report_concise.docx` を参照してください。

## フォルダ構成

| パス | 内容 |
|---|---|
| `R/` | 解析スクリプト(実行順に番号付き) |
| `results/` | モデルの出力(CSV、テキスト)。`server_rerun/` はサーバーでの再実行結果で、Mac の結果と一致したもの |
| `results/theory_results/` | 理論にもとづく検証(Petchey 2010、ニッチの充填と拡大、魚の雑食性) |
| `figures/` | 本文図と補足図(PNG と PDF) |
| `reports/` | レポート(Markdown と Word) |
| `data/` | IsoFresh の CSV を置く場所(リポジトリには含まない) |

## スクリプト

| スクリプト | 内容 |
|---|---|
| `01_network_analysis_full_pipeline.R` | 同位体類似度ネットワークの作成、指標とヌルモデルによる標準化、WorldClim・HydroLAKES の共変量の取得、混合モデル、ブートストラップ媒介、piecewise SEM までの一括処理 |
| `02_make_manuscript_figures.R` | 本文図と補足図の作成 |
| `03_sampling_effort_sensitivity.R` | サンプリング努力(`N_Sample` の合計、log10)を共変量に加えた感度解析と、サンプル数の多い網に絞った解析 |
| `04_sem_sampling_effort.R` | サンプリング努力を全方程式に入れた piecewise SEM(緯度、気温、季節性) |
| `05_mediation_bootstrap_sampling_effort.R` | 研究単位のクラスターブートストラップ(999 回)による間接効果の推定 |
| `06_sem_modularity_zQ.R` | Q と z-Q の効果の分解(SEM) |
| `07_sem_path_diagrams.R` | SEM の経路図(図8、図S10、図S11) |
| `08_theory_tests.R` | Petchey et al. 2010 の検証、温度×体サイズ、リンク–種数スケーリング×緯度、ニッチの充填と拡大、魚の栄養段階 |
| `09_theory_sem.R` | 魚の栄養段階を経路に加えた SEM |
| `10_theory_figures.R` | 理論検証の図(図5〜7) |

## 実行方法

IsoFresh の CSV 4 つ(`IsoFresh_Ecosystem.csv`、`IsoFresh_FoodWeb.csv`、`IsoFresh_Isotope_Fish.csv`、`IsoFresh_Isotope_Other.csv`)を `data/` に置いてから、次のように実行します。

```bash
mkdir work && cp data/IsoFresh_*.csv work/ && cd work
Rscript ../R/01_network_analysis_full_pipeline.R          # 重い処理(ネットワーク指標とヌルモデル)
for s in 02 03 04 05 06 07 08 09 10; do ISOFRESH_RAW=. Rscript ../R/${s}_*.R; done
```

主な環境変数:`ISOFRESH_RAW`(IsoFresh の CSV の場所)、`ISOFRESH_DATA`(`IsoFresh_Metrics_v2.RData` の場所)、`ISOFRESH_N_SIMUL`(ヌルモデルの反復数、既定 199)、`ISOFRESH_N_BOOT`(ブートストラップの反復数、既定 999)、`ISOFRESH_CORES`(並列コア数)。

## 方法の概要

- **ノード**:食物網ごとの分類群(δ13C と δ15N の平均値)。
- **リンクの重み**:1 −(ユークリッド距離 ÷ その食物網内の最大距離)。全食物網に共通の閾値(0.66;全食物網の類似度の平均)以上のペアをリンクとした。
- **指標**:結合度 L / [S(S−1)/2]、重み付き NODF とモジュラリティ(rnetcarto)をヌルモデル(重みの対称な並べ替え)で標準化した値、平均リンク強度、平均近接中心性。
- **気候**:緯度の絶対値、WorldClim 2.1 の年平均気温(BIO1)と気温の季節性(BIO4/100;月平均気温の標準偏差、℃)。
- **統計**:地域と出典論文を入れ子のランダム効果 `(1 | Region/Source_DOI)` とする混合モデル、FDR 補正(Benjamini–Hochberg)、piecewise SEM、研究単位のクラスターブートストラップ。

## 実行環境

主な解析は R 4.5.1(Ubuntu 24.04)で行った:lme4 2.0.6、lmerTest 3.2.1、piecewiseSEM 2.3.1、rnetcarto 0.2.6、igraph 2.2.1、vegan 2.7.2、terra 1.9.46、rfishbase 5.0.1。理論検証(08〜10)は R 4.4.2(macOS)で実行した。
