#!/bin/bash

# experiment number comes from first argument
EXPNO=$1

if [ -z "$EXPNO" ]; then
  echo "Usage: $0 <experiment_number>"
  exit 1
fi

echo "Running for experiment number: $EXPNO"

grep -i 'Number of plastic and viscous grid cells' \
  ../output_post_files/output_mu_test_"$EXPNO" \
  > outputs_VP/num_visc_plas_"$EXPNO".out

grep -i 'Number of viscous eta and viscous zeta grid cells' \
  ../output_post_files/output_mu_test_"$EXPNO" \
  > outputs_VP/num_mixed_visc_plas_"$EXPNO".out

grep 'min, max u' ../output_post_files/output_mu_test_"$EXPNO" >\
   ../output_post_files/min_max_vel_"$EXPNO".out