tag=$1;
data_dir="/home/lhcb/rich-pd/pmt-test-facility/data/afterglow-data/"

mv wave_1.dat TR_0_0.dat $data_dir

cd $data_dir
mv wave_1.dat sig_$tag.dat
mv TR_0_0.dat trg_$tag.dat

