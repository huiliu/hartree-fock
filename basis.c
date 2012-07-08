#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "basis.h"

void* read_basis(const char * file_name)
{
// 每个原子有几个基函数需要指定说明
    FILE *f;
    char *sparate = "****";
    char *atom_orbit_type = "SP";   // SP型一行有3个参数
    char symbol[5];
    double param;
    int state = 0, i;
    int gauss_num; //row_num每一行参数的个数，gauss_num    每一块gauss函数的个数
    int basis_num = 2, basis_total = 9, basis_i = 0, ib = 0;
    double tmp_alpha, tmp_coeff_1, tmp_coeff_2;     //用以临时存储从文件中读取的基函数参数信息

    // basis用于保存所有基函数
    BASIS *basis = calloc(sizeof(BASIS), basis_total);

    f = fopen(file_name, "r");
    
    while (1) {
        switch (state) {
            case 0: // initial state
                fscanf(f, "%s", symbol);
                // symbol == "****", start read atom information
                if (strcmp(symbol, sparate) == 0)   state = 1;
                //basis_i++;
                break;
            case 1: // read an atom's information
                if (fscanf(f, "%s", symbol) == EOF) // "H" read the atom symbol, initialize some parameter which i donn't know now
                    return basis;
                    //exit(EXIT_FAILURE);
                fscanf(f, "%lf", &param);     // 0

                if (strncmp(symbol, "N", 1) == 0)
                    basis_num = 5;
                else if (strncmp(symbol, "H", 1) == 0)
                    basis_num = 1;
                ib = 0;
                state = 2;  // start read basis set information
                break;
            case 2:
                if (ib >= basis_num) {
                    state = 0;
                    break;
                }else{
                    ib++;
                }

                fscanf(f, "%s", symbol);        // " "S" 3   1.00  "
                fscanf(f, "%d", &gauss_num);     // "  S "3"  1.00  " 每一块高斯函数的数目
                fscanf(f, "%lf", &param);         // "  S  3  "1.00" "    不清楚什么用途

                if (strcmp(symbol, atom_orbit_type) == 0) {
                    state = 4;  // SP 有三列数据，第一列为gaussian函数的指数，2,3列分别为2S, 2P中的组合系数
                }else{
                    state = 3;  // start read basis set data of S orbital
                }

                break;
            case 3:
                for (i = 0; i < gauss_num; i++) {
                    fscanf(f, "%lf", &tmp_alpha);
                    fscanf(f, "%lf", &tmp_coeff_1);

                    basis[basis_i].gaussian[i].alpha = tmp_alpha;
                    // 注意，此处的归一化只是针对1S轨道的 (2a/pi)^(3/4) 《量子化学》中册，P50
                    //basis[basis_i].gaussian[i].A = tmp_coeff_1 * pow(2*tmp_alpha/M_PI, (double)3.0/4);
                    basis[basis_i].gaussian[i].coeff = tmp_coeff_1;
#ifdef BEBUG_OUTPUT_BASIS_SET
                    printf("%20.10lf%20.10lf\n", basis[basis_i].gaussian[i].alpha, basis[basis_i].gaussian[i].coeff);
#endif

                    //basis[basis_i].gaussian[i].l = 0;
                    //basis[basis_i].gaussian[i].m = 0;
                    //basis[basis_i].gaussian[i].n = 0;
                    //printf("%15.10lf%15.10lf%15.10lf%15.10lf\n", tmp_coeff_1, pow(2*tmp_alpha/M_PI, (double)3.0/4), tmp_alpha, 2.0);
                }
                basis_i ++;
                ib++;
                state = 2;  //读完一组基函数信息，状态返回，读取下一组
                break;
            case 4: // 一个gaussian函数由3个参数决定的情形
                for (i = 0; i < gauss_num; i++) {
                    fscanf(f, "%lf", &tmp_alpha);
                    fscanf(f, "%lf", &tmp_coeff_1);
                    fscanf(f, "%lf", &tmp_coeff_2);
                    // 2S
                    basis[basis_i].gaussian[i].alpha = tmp_alpha;
                    basis[basis_i].gaussian[i].coeff = tmp_coeff_1;
                    // 2P
                    // Px
                    basis[basis_i+1].gaussian[i].alpha = tmp_alpha;
                    basis[basis_i+1].gaussian[i].coeff = tmp_coeff_2;
                    basis[basis_i+1].gaussian[i].l = 1;
                    // Py
                    basis[basis_i+2].gaussian[i].alpha = tmp_alpha;
                    basis[basis_i+2].gaussian[i].coeff = tmp_coeff_2;
                    basis[basis_i+2].gaussian[i].m = 1;
                    // Pz
                    basis[basis_i+3].gaussian[i].alpha = tmp_alpha;
                    basis[basis_i+3].gaussian[i].coeff = tmp_coeff_2;
                    basis[basis_i+3].gaussian[i].n = 1;
                }
                basis_i += 4;
                ib += 4;
                state = 2;  //读完一组基函数信息，状态返回，读取下一组
                break;
        }
    }
    fclose(f);
    return basis;
}

// 参数 count 表示一个基函数由count个gaussian函数组成
void basis_set_output(const BASIS* b, int count, char* msg)
{
    int i;
    double alpha, coeff, l, m, n;
    printf("%s", msg);
    for (i = 0; i < count; i++) {
        alpha = b->gaussian[i].alpha;
        coeff = b->gaussian[i].coeff;
        l = b->gaussian[i].l;
        m = b->gaussian[i].m;
        n = b->gaussian[i].n;
        printf("%12.8lf%12.8lf\t%g\t%g\t%g\n", alpha, coeff, l, m, n);
    }
}

//int main(int argc, char** argv)
//{
//    BASIS *b = read_basis("basis_set");    
//    return 0;
//}
