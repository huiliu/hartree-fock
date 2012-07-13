#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include "basis.h"
#include "overlap.h"

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
    printf("%s\n", msg);
    for (i = 0; i < count; i++) {
        alpha = b->gaussian[i].alpha;
        coeff = b->gaussian[i].coeff;
        l = b->gaussian[i].l;
        m = b->gaussian[i].m;
        n = b->gaussian[i].n;
        printf("%12.8lf%12.8lf\t%g\t%g\t%g\n", alpha, coeff, l, m, n);
    }
}

void atom_output(const ATOM_INFO* atom, int n)
{
    ATOM_INFO *a=NULL;
    int i, j;
    for (i = 0; i < n; i++) {
        a = (ATOM_INFO*)(atom + i);
        printf("%s %d %d%12.7lf%12.7lf%12.7lf\n", a->symbol, a->n, a->basis_count, a->c->data[0], a->c->data[1], a->c->data[2]); 
        for (j = 0; j < a->basis_count; j++)
            basis_set_output(a->basis + j, 3, "Basis Function:");
    }
}

// -------------------------------------------------------------------------------
// 读取基组方法的重新实现
#define Item_count              3
#define Atom_count_max          4
#define Basis_set_count_max     2

INPUT_INFO* parse_input(const char* file_name)
{
// parse the input file, translate the input file into useable data
    FILE *f;
    INPUT_INFO *input_information;
    int Atom_index = 0, basis_count = 0;
    char Item[Item_count][8] = {"$COORD", "$BASIS", "$END"};
    char BasisSetName[Basis_set_count_max][8] = {"STO-3G", "6-31G"};
    ATOM_INFO *atom; // save the coordination of atoms
    int state = 0, i;
    char input_item[9];
    input_information = malloc(sizeof(INPUT_INFO));

    f = fopen(file_name, "r");

    while (1) {
        switch (state) {
            case 0: // initial state
                if (fscanf(f, "%s", input_item) == EOF) {
                    exit(EXIT_FAILURE);
                    return NULL;
                }
                for (i = 0; i < Item_count; i++) {
                    // setup the state will go according the Item value
                    if (strcmp(input_item, Item[i]) == 0) {
                        state = i + 1;
                        break;
                    }
                }
                break;
            case 1: // read the coordination information
                fscanf(f, "%s", input_item);
                // if read the item "$END", continue the other state
                if (strcmp(input_item, Item[2]) == 0) {
                    state = 0;
                    break;
                }
                // save coordination
                input_information->gxyz = (gsl_vector**) \
                                realloc(sizeof(gsl_vector*), (Atom_index + 1));
                input_information->gxyz[Atom_index] = gsl_vector_alloc(3);
                // save atom information
                input_information->c = (ATOM_INFO**) \
                                    realloc(sizeof(ATOM_INFO *), Atom_index+1);
                atom = input_information->c[Atom_index];
                // read element symbol of atom
                strcpy(atom[Atom_index].symbol, input_item);
                // read element core electronics of atom
                fscanf(f, "%d", &atom[Atom_index].n);
                for (i = 0; i < 3; i++)
                    // read coordination of atom
                    fscanf(f, "%lf", atom[Atom_index].c->data + i);

                Atom_index++;
                break;
            case 2: // read the basis set information
                fscanf(f, "%s", input_item);

                for (i = 0; i < Basis_set_count_max; i++) {
                    if (strcmp(input_item, BasisSetName[i]) == 0)
                        break;
                }

                int j;
                switch (i) {
                    case 0: // STO-3G
                        // Set some parameters.
                        for (j = 0; j < Atom_index; j++) {
                            switch (atom[j].n) {
                                case 1:     // the basis count of atom j
                                    input_information->basis_count += 1;
                                    input_information->basis_set = realloc(
                                                sizeof(BASIS*), 
                                                input_information->basis_count);
                                    break;
                                case 2:
                                    break;
                                case 7:
                                    input_information->basis_count += 5;
                                    input_information->basis_set = calloc(
                                                sizeof(BASIS),
                                                input_information->basis_count);
                                    break;
                            }
                        }
                        break;
                    case 1: // 6-31G
                        break;
                } //end setup basis set information
                // read basis set
                // TODO:
                //      读完基函数直接退出，要修改的更友善
                readbasis(f, atom, Atom_index);
                input_information->n = Atom_index;
                return input_information;
                break;
            case 3: // block end
                break;
        }
    }
    input_information->n = Atom_index;
    return input_information;
}

// 是函数read_basis针对新的数据结构的升级版
int readbasis(FILE * f, ATOM_INFO* atom_list, int atom_count)
{
// read the part contain basis set in the input file
// 每个原子有几个基函数需要指定说明
    char *sparate = "****";
    char *atom_orbit_type = "SP";   // SP型一行有3个参数
    char symbol[5];
    double param;
    int state = 0, i;
    int gauss_num; //gauss_num    每一块gauss函数的个数
    int basis_num, basis_i = 0, ib = 0;
    double tmp_alpha, tmp_coeff_1, tmp_coeff_2;     //用以临时存储从文件中读取的基函数参数信息
    int atom_index = -1;

    // basis用于保存所有基函数
    //BASIS *basis = calloc(sizeof(BASIS), basis_total);
    BASIS *basis;

    //atom_output(atom_list, atom_count);

    while (1) {
        switch (state) {
            case 0: // initial state
                fscanf(f, "%s", symbol);
                // symbol == "****", start read atom information
                if (strcmp(symbol, sparate) == 0)   state = 1;
                break;
            case 1: // read an atom's information, such as the line 11 21 27 33 of this file
                //if (fscanf(f, "%s", symbol) == EOF) // 
                fscanf(f, "%s", symbol); // "H" read the atom symbol, initialize some parameter which i donn't know now
                if (strcmp(symbol, "$END") == 0) {
                    //return basis;
                    return 0;
                }

                fscanf(f, "%lf", &param);     // 0
                basis = atom_list[++atom_index].basis;
                basis_num = atom_list[atom_index].basis_count;
                basis_i = 0;
                // read the basis set count of every atom
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
                basis[basis_i].gauss_count = gauss_num;

                if (strcmp(symbol, atom_orbit_type) == 0) {
                    state = 4;  // SP 有三列数据，第一列为gaussian函数的指数，2,3列分别为2S, 2P中的组合系数
                }else{
                    state = 3;  // start read basis set data of S orbital
                }

                break;
            case 3: // 当前为读取S
                for (i = 0; i < gauss_num; i++) {
                    fscanf(f, "%lf", &tmp_alpha);
                    fscanf(f, "%lf", &tmp_coeff_1);

                    basis[basis_i].gaussian[i].alpha = tmp_alpha;
                    basis[basis_i].gaussian[i].coeff = tmp_coeff_1;
                    basis[basis_i].gaussian[i].norm = normalize_coeff(&basis[basis_i].gaussian[i]);
#ifdef DEBUG_OUTPUT_BASIS_SET
                    printf("%20.10lf%20.10lf\n", basis[basis_i].gaussian[i].alpha, basis[basis_i].gaussian[i].coeff);
#endif

                }
                basis_i ++;
                ib++;
                state = 2;  //读完一组基函数信息，状态返回，读取下一组
                break;
            case 4: // 读取SP标签：S与P轨道
                for (i = 0; i < gauss_num; i++) {
                    fscanf(f, "%lf", &tmp_alpha);
                    fscanf(f, "%lf", &tmp_coeff_1);
                    fscanf(f, "%lf", &tmp_coeff_2);
                    // 2S
                    basis[basis_i].gaussian[i].alpha = tmp_alpha;
                    basis[basis_i].gaussian[i].coeff = tmp_coeff_1;
                    basis[basis_i].gaussian[i].norm = normalize_coeff(&basis[basis_i].gaussian[i]);

                    basis[basis_i+1].gauss_count = basis[basis_i+2].gauss_count \
                    = basis[basis_i+3].gauss_count = basis[basis_i].gauss_count;
                    // 2P
                    // Px
                    basis[basis_i+1].gaussian[i].alpha = tmp_alpha;
                    basis[basis_i+1].gaussian[i].coeff = tmp_coeff_2;
                    basis[basis_i+1].gaussian[i].l = 1;
                    basis[basis_i+1].gaussian[i].norm = normalize_coeff(&basis[basis_i+1].gaussian[i]);
                    // Py
                    basis[basis_i+2].gaussian[i].alpha = tmp_alpha;
                    basis[basis_i+2].gaussian[i].coeff = tmp_coeff_2;
                    basis[basis_i+2].gaussian[i].m = 1;
                    basis[basis_i+2].gaussian[i].norm = normalize_coeff(&basis[basis_i+2].gaussian[i]);
                    // Pz
                    basis[basis_i+3].gaussian[i].alpha = tmp_alpha;
                    basis[basis_i+3].gaussian[i].coeff = tmp_coeff_2;
                    basis[basis_i+3].gaussian[i].n = 1;
                    basis[basis_i+3].gaussian[i].norm = normalize_coeff(&basis[basis_i+3].gaussian[i]);
                }
                basis_i += 4;
                ib += 4;
                state = 2;  //读完一组基函数信息，状态返回，读取下一组
                break;
        }
    }
    return 0;
}

double normalize_coeff(const GTO *g)
{
    double alpha = g->alpha;
    double l = g->l;
    double m = g->m;
    double n = g->n;

    return pow(2 * alpha / M_PI, 0.75) * sqrt(pow(4*alpha, l + m + n) / \
        (factorial_2(2*l-1) * factorial_2(2*m-1) * factorial_2(2*n-1)));
}
