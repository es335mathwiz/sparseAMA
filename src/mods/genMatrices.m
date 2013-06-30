SPSolvePath='c:/users/m1gsa00/sp_solve';
mExchangePath='./matlabFileExchange';
addpath(SPSolvePath);
addpath(mExchangePath);
'gen firmvalue matrices'
pushd './modelez/firmvalue'
genInsAndOuts('firmvalue')
popd

