thePwd=pwd;
mExchangePath=[thePwd '/matlabFileExchange'];
addpath(SPSolvePath);
addpath(mExchangePath);
'gen firmvalue matrices'
pushd './modelez/firmvalue'
genInsAndOuts('firmvalue')
popd

