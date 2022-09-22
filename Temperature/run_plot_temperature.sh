#!/bin/bash

#:<<"#__COMMENT__"
python3 plot_temperature.py \
    --channel_columns temp1 temp2 temp3 temp4 \
    --channel_labels \
        'CH1 2nd stage' \
        'CH2 1st stage' \
        'CH3 Wave-guide @ 40K' \
        'CH4 2nd stage far' \
    --data_dir '/data/lakeshore218' \
    --data_starts \
        '2022/08/24 08:48' '2022/09/21 17:43' \
    --data_stops \
        '2022/08/25 12:48' '2022/09/30 00:00' \
    --data_labels \
        '5th cooling' '6th cooling (Wt SIS)' \
    --data_linestyles "solid" "dashed" \
    --data_linewidth 2 \
    --xtype 'sec' \
    --grid --logy \
    --figsizeW 16 --figsizeH 4 \
    --outname 'temperature_compare' \
#__COMMENT__

:<<"#__COMMENT__"
python3 plot_temperature.py \
    --channel_columns temp1 \
    --channel_labels \
        '4K stage' \
    --data_dir '/data/lakeshore218' \
    --data_starts \
        '2022/08/24 08:48' \
    --data_stops \
        '2022/08/25 06:48' \
    --data_labels \
        '' \
    --data_linestyles "solid" \
    --data_linewidth 2 \
    --xtype 'sec' \
    --ymin 1.0 \
    --grid --logy \
    --figsizeW 10 --figsizeH 4 \
    --outname 'temperature_single_5thCooling' \
#__COMMENT__
