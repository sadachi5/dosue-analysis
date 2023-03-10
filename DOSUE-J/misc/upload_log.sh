#!/bin/bash
date >> upload_log.out
rsync -ruavv dosuedaq:/data/ms2840a/dosue-j/2023-03.log 2023-03.csv >> upload_log.out
