function genInsAndOuts(modName)
[parserRetCode,...
param_,np,modname,neq,nlag,nlead,eqname_,...
eqtype_,endog_,delay_,vtype_]=SPParser('./',modName);
setParams;
eval([modName,'_aim_matrices']);
hmat=cofh;hmat(1:neq,1:neq*(nlag+1))=hmat(1:neq,1:neq*(nlag+1))+cofg;
[hmatRp,hmatCi,hmatAi,ncol]=sparse_to_csr(sparse(hmat));
condn=1.0e-8;uprbnd=1+condn;
[bmat,rts,ia,nexact,nnumeric,lgroots,aimcode]= ...
SPAmalg(hmat,neq,nlag,nlead,condn,uprbnd);
[bmatRp,bmatCi,bmatAi,ncol]=sparse_to_csr(sparse(bmat));
save('denseMats.mat','hmat','hmatRp','hmatCi','hmatAi',...
'bmat','bmatRp','bmatCi','bmatAi','-v4');     
end       