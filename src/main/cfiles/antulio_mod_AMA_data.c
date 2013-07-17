void antulio_mod_AMA_data(int param_,const char *eqname_,int eqtype_,int endog_,int delay_,int vtype_) {

  // antulio.mod_AMA_data()
  //     This function will return various information about the AMA model,
  //     but will not compute the G and H matrices.

const char *eqname[8];
const char *param[20];
const char *endog[8];
int delay[8][1];
int vtype[8][1];


// char modname;
//modname  = "antulio.mod";
int neq = 8;
int np = 20;
int nlag = 1;
int nlead = 1;
int eqtype[8][1];

  eqname[1] = "YA";
  eqname[2] = "KA";
  eqname[3] = "YB";
  eqname[4] = "KB";
  eqname[5] = "LAT";
  eqname[6] = "AT";
  eqname[7] = "SMALLA";
  eqname[8] = "ONE";
  //eqname_ = eqname;

  eqtype[1][1] = 1;     eqtype[2][1] = 1;     eqtype[3][1] = 1;   
  eqtype[4][1] = 1;     eqtype[5][1] = 1;     eqtype[6][1] = 1;   
  eqtype[7][1] = 1;     eqtype[8][1] = 1;   
  //eqtype_ = eqtype;

  param[1] = "YAK";
  param[2] = "YAE";
  param[3] = "YAB";
  param[4] = "YAL";
  param[5] = "MAK";
  param[6] = "QA11";
  param[7] = "QA12";
  param[8] = "MAL";
  param[9] = "YBK";
  param[10] = "YBE";
  param[11] = "YBA";
  param[12] = "MBK";
  param[13] = "MBE";
  param[14] = "MBA";
  param[15] = "MUA2";
  param[16] = "P11";
  param[17] = "P21";
  param[18] = "P12";
  param[19] = "P22";
  param[20] = "RHO";
  //param_ = char(param);

  endog[1] = "YA";
  endog[2] = "KA";
  endog[3] = "YB";
  endog[4] = "KB";
  endog[5] = "LAT";
  endog[6] = "AT";
  endog[7] = "SMALLA";
  endog[8] = "ONE";
  // endog_ = char(endog);

  delay[1][1] = 0;     delay[2][1] = 0;     delay[3][1] = 0;   
  delay[4][1] = 0;     delay[5][1] = 0;     delay[6][1] = 0;   
  delay[7][1] = 0;     delay[8][1] = 0;   
  //delay_ = delay;

  vtype[1][1] = 1;     vtype[2][1] = 1;     vtype[3][1] = 1;   
  vtype[4][1] = 1;     vtype[5][1] = 1;     vtype[6][1] = 1;   
  vtype[7][1] = 1;     vtype[8][1] = 2;   
  //vtype_ = vtype;



}
