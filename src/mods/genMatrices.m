thePwd=pwd;
mExchangePath=[thePwd '/matlabFileExchange'];
addpath(SPSolvePath);
addpath(mExchangePath);
addpath(thePwd);
'gen firmvalue matrices'
pushd './modelez/firmvalue'
genInsAndOuts('firmvalue')
popd
'gen chrisDavidMatt matrices'
pushd './modelez/chrisDavidMatt'
genInsAndOuts('modelv1')
popd

