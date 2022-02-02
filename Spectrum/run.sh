#!/bin/bash

nLoop=1000
gain='1.e+6'
T=100
time=60
freq='20.e+9'
phi=6

#python3 spectrum.py -o 'freq-bin-width-fine_vs_fit-result_2e-10-signal_whiterandom.pdf' --test 0 -n $nLoop -c 2.e-10 -T $T -g $gain -t $time -f $freq -A $phi
#python3 spectrum.py -o 'freq-bin-width-fine_vs_fit-result_no-signal_whiterandom.pdf' --test 0 -n $nLoop -c 0 -T $T -g $gain -t $time -f $freq -A $phi
#python3 spectrum.py -o 'freq-bin-width_vs_fit-result_2e-10-signal_whiterandom.pdf' --test 0 -n $nLoop -c 2.e-10 -T $T -g $gain -t $time -f $freq -A $phi
#python3 spectrum.py -o 'freq-bin-width_vs_fit-result_no-signal_whiterandom.pdf' --test 0 -n $nLoop -c 0 -T $T -g $gain -t $time -f $freq -A $phi


#python3 spectrum.py -o 'dosue-k_c1e-9.pdf' --test 0 -b '0.5e+3,1.e+3,5.e+3,10.e+3' -n 1 -c '1e-9' -T 130 --gainDB 65 -t 60 -f $freq -A $phi
#python3 spectrum.py -o 'aho.pdf' --test 1 -n 1 -c 0 -T 130 --gainDB 65 -t 60 -f $freq -A $phi

# KUMODES
#python3 spectrum.py -o 'kumodes_c0_n1.pdf' --test 0 -b '1e+3' -n 1 -c 0 -T 190 --gainDB 65 -t 120 -f '28.e+9' -A $phi


# DOSUE-K
# chi = 10e-10
#python3 spectrum.py -o 'aho.pdf' --test 0 -b '1.e+3' -n 1 -c '1e-9' -T 130 --gainDB 65 -t 60 -f $freq -A $phi
#python3 spectrum.py -o 'dosue-k_c1e-9_n1000.pdf' --test 0 -b '0.5e+3,1.e+3,5.e+3,10.e+3,50.e+3,100.e+3' -n 1000 -c '1e-9' -T 130 --gainDB 65 -t 60 -f $freq -A $phi
#python3 spectrum.py -o 'dosue-k_c1e-9_n1000.pdf' --test 0 -b '0.5e+3,1.e+3,5.e+3,10.e+3,50.e+3,100.e+3' -n 1000 -c '1e-9' -T 130 --gainDB 65 -t 60 -f $freq -A $phi
#python3 spectrum.py -o 'dosue-k_c1e-9_n1000_fineBW.pdf' --test 0 -b '0.1e+3,0.2e+3,0.5e+3,0.75e+3,1.e+3' -n 1000 -c '1e-9' -T 130 --gainDB 65 -t 60 -f $freq -A $phi
#python3 spectrum.py -o 'dosue-k_c1e-9_n1000_fineBW2.pdf' --test 0 -b '1.e+3,3.e+3,5.e+3,7.e+3,9.e+3' -n 1000 -c '1e-9' -T 130 --gainDB 65 -t 60 -f $freq -A $phi
#python3 spectrum.py -o 'dosue-k_c1e-9_n1000_many.pdf' --test 0 -b '0.1e+3,0.5e+3,0.75e+3,1.e+3,3.e+3,5.e+3,7.e+3,10.e+3,50.e+3,100.e+3' -n 1000 -c '1e-9' -T 130 --gainDB 65 -t 60 -f $freq -A $phi
# chi = 0
#python3 spectrum.py -o 'dosue-k_c0_n1000.pdf' --test 0 -b '0.5e+3,1.e+3,5.e+3,10.e+3,50.e+3,100.e+3' -n 1000 -c '0.' -T 130 --gainDB 65 -t 60 -f $freq -A $phi
#python3 spectrum.py -o 'dosue-k_c0_n1000_fineBW.pdf' --test 0 -b '0.1e+3,0.2e+3,0.5e+3,0.75e+3,1.e+3' -n 1000 -c '0.' -T 130 --gainDB 65 -t 60 -f $freq -A $phi
#python3 spectrum.py -o 'dosue-k_c0_n1000_many.pdf' --test 0 -b '0.1e+3,0.5e+3,0.75e+3,1.e+3,3.e+3,5.e+3,7.e+3,10.e+3,50.e+3,100.e+3' -n 1000 -c '1e-9' -T 130 --gainDB 65 -t 60 -f $freq -A $phi


#python3 spectrum.py -o 'aho.pdf' --test 0 -b '0.3e+3' -n 1 -c '0.' -T 130 --gainDB 65 -t 60 -f $freq -A $phi -v 1
python3 spectrum.py -o 'aho2.pdf' --test 0 -b '0.3e+3' -n 1 -c '0.' -T 130 --gainDB 65 -t 60 -f $freq -A $phi --rebin 10 -v 1
