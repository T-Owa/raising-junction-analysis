# RAISING解析スクリプト 公開用パッケージ

公開パッケージ版：**1.0.0**。今回使用した800 bpドナーと3′／5′側の最終プライマーを対象に、全リードペアの4分類、ドナー–ヒト接合の抽出・再評価、MH／ドナー末端欠失／介在配列のFigure用集計を行います。

個人名、PC名、個人の絶対パス、元サンプルの結果・ログ・リードは同梱していません。科学的な再現に必要なドナー配列・最終プライマー・判定条件は残しています。別ドナーにもそのまま使える汎用解析ソフトではありません。

## まず確認するもの

- `README.md`：英語の導入・実行手順。
- `environment.yml`：必要なソフトウェア。既存のraising環境とは別環境で導入できます。
- `docs/METHODS.md`：判定条件、分母、解釈の限界、Methods文章のたたき台。
- `docs/VALIDATION.md`：実施したテストと検証範囲。
- `docs/REFERENCE_PROVENANCE.md`：リファレンス情報で確定している点・未確定の点。
- `docs/CODE_AVAILABILITY.md`：公開後にURL／DOIを入れて使う文章。
- `docs/CITATION.cff.example`：公開する著者名・URLを入れる引用情報のテンプレート。

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

## これまでの結果を集計し直す場合

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

公開版では説明文などの変更でもコードハッシュが変わるため、**元の解析フォルダーに上書きして途中再開しないでください。** 元の結果・スクリプトは保持し、公開版の全量再実行には別フォルダーを使います。入力やコードの変更を検出する仕組みは維持しています。

## 論文での解釈

4分類は全入力リードペアの組成です。コンカテマー候補にはヒトへのマッピング支持を併せ持つペアも含まれます。Short donorは遊離ドナーや非コンカテマーを確定する分類ではありません。

接合集計は支持された配列パターン単位です。細胞数、独立分子数、独立挿入イベント数とは一致しません。MHは配列の特徴であり、それだけで修復経路は確定できません。介在配列を含む接合のMHは未評価（NA）として扱います。

今回の欠失はドナー末端の欠失です。ヒトゲノム側の欠失や挿入物全体の構造を確定するものではありません。元の解析基準を保っているため、単純接合の名目上のドナー欠失上限39 bp（3p）／45 bp（5p）、介在配列上限20 bpも変更していません。これより大きな変化が生物学的に存在しないという意味ではありません。

## 公開するとき

このフォルダーをコード用リポジトリに配置できます。著者表記・リポジトリURL・DOIは未登録のため、公開時にテンプレートを埋めます。MITライセンス案を同梱しています。論文で使用した版を固定してアーカイブし、その版のDOIとバージョンを引用してください。

**これから実行して生成されるログ・manifest・詳細証拠表には、再現性のため実行環境のパスやread ID／配列が記録されます。** コードを匿名化したことと、生成結果を匿名化したことは別です。実データや生成物はコード公開用フォルダーに混ぜず、共有用コピーを別に確認してください。
