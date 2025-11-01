#include <iostream>
#include <cmath>
#include <limits>

using namespace std;

int main() {
    int k = 3;
    int K = 4;

    // 使用double直接表示大数，避免移位溢出
    double q_bs = 2048;
    double q_ks = 262144;          // 1ULL << 18
    double B_ks = 32;             // 1ULL << 9
    double d_ks = 4;
    double n = 800;

    // 使用科学计数法表示非常大的数
    double Q = 1ULL<<60;
    //double Q = 1.4757395258967641e+20;  // 2^67
    //double B_bs = 8388608;               // 1ULL << 23
    double B_bs = 512;
    double d_bs = 7;
    double N = 4096;                     // 1ULL << 13

    double sigma_e_squared = 9;
    double sigma_k_squared = 0.5;

    //计算sigma_br^2
    double sigma_br_squared = (2 * n + 1) * d_bs * N * B_bs * B_bs * sigma_e_squared / 12;
    cout<<"sigma_br_squared = "<<sigma_br_squared<<" 6*sqrt = "<<6*sqrt(sigma_br_squared)<<" Q/16 = "<<Q / 16;
    if (6 * sqrt(sigma_br_squared) < Q / 16) {
        cout << " 符合要求" << endl;
    } else {
        cout << " 不符合要求," << endl;
    }

    //计算sigma_hp^2
    double term1 = k * K * (K + 1) * d_bs * N * N * B_bs * B_bs * sigma_br_squared / 48;
    double term2 = k * k * K * (K + 1) * d_bs * N * N * B_bs * B_bs * sigma_e_squared / 24;
    double term3 = K * d_bs * N * B_bs * B_bs * (2 * k * sigma_e_squared + sigma_br_squared) / 12;
    
    double sigma_hp_squared = term1 + term2 + term3;
    cout<<"sigma_hp_squared = "<<sigma_hp_squared<<" 6*sqrt = "<<6*sqrt(sigma_hp_squared)<<" Q/16 = "<<Q / 16;
    if (6 * sqrt(sigma_hp_squared) < Q / 16) {
        cout << " 符合要求" << endl;
    } else {
        cout << " 不符合要求" << endl;
    }

    //计算第二次模切换
    double sigma_MS_2_squared = ((q_ks * q_ks) / (Q * Q)) * sigma_hp_squared + (1 + k * K * N * sigma_k_squared) / 12;
    cout<<"sigma_MS_2_squared = "<<sigma_MS_2_squared<<" 6*sqrt = "<<6*sqrt(sigma_MS_2_squared)<<" q_ks/16 = "<<q_ks / 16;
    if (6 * sqrt(sigma_MS_2_squared) < q_ks / 16) {
        cout << " 符合要求" << endl;
    } else {
        cout << " 不符合要求" << endl;
    }

    //计算密钥切换
    double sigma_KS_squared = sigma_MS_2_squared + (k * K * d_ks * B_ks * B_ks * sigma_e_squared) / 12;
    cout<<"sigma_KS_squared = "<<sigma_KS_squared<<" 6*sqrt = "<<6*sqrt(sigma_KS_squared)<<" q_ks/16 = "<<q_ks / 16;
    if (6 * sqrt(sigma_KS_squared) < q_ks / 16) {
        cout << " 符合要求" << endl;
    } else {
        cout << " 不符合要求" << endl;
    }

    //第三次模切换
    double sigma_MS_3_squared = ((q_bs * q_bs) / (q_ks * q_ks)) * sigma_KS_squared + (1 + k * K * n * sigma_k_squared) / 12;
    cout<<"sigma_MS_3_squared = "<<sigma_MS_3_squared<<" 6*sqrt = "<<6*sqrt(sigma_MS_3_squared)<<" q_bs/16 = "<<q_bs / 16;
    if (6 * sqrt(sigma_MS_3_squared) < q_bs / 16) {
        cout << " 符合要求" << endl;
    } else {
        cout << " 不符合要求" << endl;
    }

    // double std_dev = sqrt(sigma_MS_3_squared);
    // cout << "能参与下次计算的噪声阈值为Q/8 = " << Q / 8 << endl;
    // cout << " 6倍标准差为 = " << 6 * std_dev << endl;

    // if (6 * std_dev < Q / 8) {
    //     cout << " 符合要求" << endl;
    // } else {
    //     cout << " 不符合要求" << endl;
    // }

    return 0;
}