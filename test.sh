#!/bin/bash

OUTPUT_DIR="./test"  # ここで出力先ディレクトリを指定します。デフォルトはカレントディレクトリです。

for YEAR in {2024..2027}; do
  for MONTH in {01..12}; do
    # 2025年の5月以降は処理しない
    if [ "$YEAR" -eq 2027 ] && [ "$MONTH" -gt 12 ]; then
      break
    fi

    # 日付をyyyy/mm形式にフォーマット
    DATE="$YEAR/$MONTH"

    # Pythonスクリプトを実行
    python ./altitude_calculation.py --date $DATE --output $OUTPUT_DIR --moon 70 --jupiter 70 --venus 70 --mars 70 --sun 50
  done
done
