/*
STO-3G
---------------
****
H     0 
S   3   1.00
      3.42525091             0.15432897       
      0.62391373             0.53532814       
      0.16885540             0.44463454       
****
N     0 
S   3   1.00
     99.1061690              0.15432897       
     18.0523120              0.53532814       
      4.8856602              0.44463454       
SP   3   1.00
      3.7804559             -0.09996723             0.15591627       
      0.8784966              0.39951283             0.60768372       
      0.2857144              0.70011547             0.39195739       
****
*/

typedef struct gto {
    // gaussian = A * exp(-a * r^2)
    // A = A * (2a/pi)^(3/4)    仅对1S轨道
    double l, m, n;
    double alpha;
    double coeff;
}GTO;

typedef struct _b {
    GTO gaussian[3];   // 3 表示一个基函数由3个gaussian函数构成
}BASIS;

void* read_basis(const char * );
void basis_set_output(const BASIS*, int, char* );
