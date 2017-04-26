#~ python mmatch.py \
    #~ --sources "17:29:21.4 -30:56:02" "17:29:21.4 -30:56:02" \
    #~ --pawprints fits/v20100418_00895_st_cat.fits fits \
    #~ -o matched.csv \
    #~ -f2m bin/vvv_flx2mag \
    #~ --procs 0

python mmatch.py \
    -srcf sources.csv \
    --pawprints fits/v20100418_00895_st_cat.fits fits \
    -o matched.csv \
    -f2m bin/vvv_flx2mag \
    -wd _temp \
    --procs 3

#~ python mmatch.py  -o matched.csv --pawprints fits/v20100418_00895_st_cat.fits fits
