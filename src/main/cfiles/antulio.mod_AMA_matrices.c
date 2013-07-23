int antulio.mod_AMA_matrices(char paramnames, int params)

{
//  This script will compute the G and H matrices.
// it reads in paramnames
// and params
// and then uses them to assign values to the matrix



int g[8][16] = { { 0 } };
int  h[8][24] = { { 0 } } ;

  g[65] = g[65] + 1;
  g[9] = g[9] - (YAK*1);
  g[105] = g[105] - (YAE*1);
  g[81] = g[81] - (YAB*1);
  g[97] = g[97] - (YAL*1);
  g[74] = g[74] + 1;
  g[10] = g[10] - (MAK*1);
  g[106] = g[106] - (QA11*1);
  g[82] = g[82] - (QA12*1);
  g[98] = g[98] - (MAL*1);
  g[83] = g[83] + 1;
  g[27] = g[27] - (YBK*1);
  g[107] = g[107] - (YBE*1);
  g[3] = g[3] - (YBA*1);
  g[92] = g[92] + 1;
  g[28] = g[28] - (MBK*1);
  g[108] = g[108] - (MBE*1);
  g[4] = g[4] - (MBA*1);
  g[101] = g[101] + 1;
  g[109] = g[109] - (P21*1);
  g[85] = g[85] - (P22*1);
  h[165] = h[165] + (-1.0*((1.0*(MUA2^-1.0))*1));
  h[173] = h[173] - (P11*1);
  h[149] = h[149] - (P12*1);
  g[110] = g[110] + 1;
  g[46] = g[46] - (RHO*1);
  g[118] = g[118] - 1;
  g[119] = g[119] + 1;
  g[127] = g[127] - (0.0*1);
  g[128] = g[128] + 1;
  g[64] = g[64] - 1;

  cofg = g;
  cofh = h;
  return (g,h);
}
