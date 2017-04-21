#~ python mmatch.py \
    #~ --sources "17:29:21.4 -30:56:02" "17:29:21.4 -30:56:02" \
    #~ --pawprints fits/v20100418_00895_st_cat.fits fits \
    #~ -o matched.csv \
    #~ -fcl bin/fitsio_cat_list \
    #~ --procs 0

python mmatch.py \
    -srcf sources.csv \
    --pawprints fits/v20100418_00895_st_cat.fits fits \
    -o matched.csv \
    -fcl bin/fitsio_cat_list \
    --procs 0

#~ python mmatch.py  -o matched.csv --pawprints fits/v20100418_00895_st_cat.fits fits
