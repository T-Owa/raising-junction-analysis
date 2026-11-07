# RAISING downstream解析

バージョン：**1.0.0**。同梱の800 bpドナーと3′／5′側の最終プライマーを対象に、全リードペアの4分類、ドナー–ヒト接合の抽出・再評価、MH／ドナー末端欠失／介在配列のFigure用集計を行います。

ソースコード、ドナー配列・最終プライマー、合成データによる実行例を同梱しています。研究で取得したFASTQやサンプルごとの解析結果は含みません。対象のドナーとプライマーに合わせた解析であり、別の配列への適用にはコードの変更が必要です。

編集者・査読者向けの実行手順と解析・出力の対応は、英語の[Reviewer guide](docs/REVIEWER_GUIDE.md)にまとめています。

## まず確認するもの

- [README.md](README.md)：英語の導入・実行手順。
- [environment.yml](environment.yml)：必要なソフトウェアとバージョン。
- [Methods](docs/METHODS.md)：判定条件、分母、解釈の限界。
- [Validation](docs/VALIDATION.md)：実施したテストと検証範囲。
- [Reference provenance](docs/REFERENCE_PROVENANCE.md)：参照配列の出典と同定情報。

## インストールと合成データでの確認

Ubuntu／WSLのターミナルで、展開したフォルダーに移動して実行します。

```bash
micromamba create -y -f environment.yml --strict-channel-priority
micromamba activate raising-public
python tools/check_environment.py
python tests/run_all.py
```

実際のBWA／BLAST／samtoolsを使用する、合成FASTQからFigure集計までのテスト：

```bash
python tests/run_all.py --include-external
```

合成データのホスト配列はランダムに生成したもので、ヒトの実データではありません。実データへの感度・特異度を検証するベンチマークとは区別してください。

## 完了済みの接合解析から集計する

既に作成済みの接合解析の完了フォルダーに、次の4つが揃っていれば、全FASTQの再解析は不要です。

- `batch_summary.json`
- `all_identical_sequence_groups.tsv.gz` または `.tsv`
- `all_junction_reassessment.tsv.gz` または `.tsv`
- `all_template_comparisons.tsv.gz` または `.tsv`

```bash
python tools/prepare_figure_data.py \
  --input-dir /path/to/completed_junction_batch \
  --output-dir /path/to/new_figure_data
```

`/path/to/...`は実際の場所に置き換えます。条件名と反復番号を付ける場合は、`examples/metadata.csv`を実サンプル名に書き換えたコピーを用意し、`--metadata /path/to/metadata.csv`を追加します。サンプル名から条件・生物学的反復は推定しません。

Figureの主集計は、`figure_values.csv`の **nuclear_all・min_support=1・all_supported_pairs** を選び、3p／5pを分けます。支持2ペア以上も感度解析として併記できます。

## 全量解析を最初から再実行する場合

英語READMEの順に、新しい解析ルート内へraw FASTQと同じホスト参照を用意し、4分類→接合解析→集計を実行します。最初は4分類バッチに`--allow-index-build`を付けます。全ゲノムの初回索引作成には時間・ディスク容量を要します。

ソースコード、入力、解析設定を変更する場合は、**新しい解析フォルダーを使ってください。** 再開時には入力とコードの同一性を確認するため、互換性のない完了済み処理は再利用できません。各実行のmanifestは保存してください。

解析単位、分母、MH・欠失の定義と検出範囲は[Methods](docs/METHODS.md)を参照してください。

---

[MITライセンス](LICENSE) · [引用情報](docs/CODE_AVAILABILITY.md)
