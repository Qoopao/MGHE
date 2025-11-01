#include <bits/stdint-uintn.h>
#include "binfhe-constants.h"
#include "binfhecontext.h"
#include "lattice/hal/lat-backend.h"
#include "lwe-ciphertext-fwd.h"
#include "lwe-privatekey-fwd.h"
#include "lwe-privatekey.h"
#include "math/hal/nativeintbackend.h"
#include "math/ternaryuniformgenerator.h"
#include "rlwe-ciphertext.h"
#include "utils/inttypes.h"
#include "secret-sharing.h"
#include <cstdint>
#include <iostream>
#include <vector>
#include <string>
#include "rlwe-publickey-fwd-lyh.h"
#include "rlwe-publickey-lyh.h"


using namespace lbcrypto;



int main(){
    auto cc = BinFHEContext();
    cc.GenerateBinFHEContext(BIGP, XZDDF);
    auto params = cc.GetParams()->GetVectorNTRUParams();
    uint32_t N = params->GetN();
    NativeInteger Q = params->GetQ();

    int k = 1;
    int ki = 4;
    cc.inituser(k, ki);

    cc.RLWECRPGen(cc.GetParams()->GetVectorNTRUParams());
    for(int i = 0; i < k; ++i){
        for(int j = 0; j < ki; ++j){
            auto rlwesk = cc.RLWESecretKeyGen(i, j);
            auto rlwepk = cc.RLWEPublicKeyGen(rlwesk, i, j);
        }
        auto jointrlwepk = cc.RLWEJointPublicKeyGen(cc.GetParams()->GetVectorNTRUParams(), i); 
    }



    std::vector<std::vector<NativePoly>> Triples(3,std::vector<NativePoly>(ki));
    TernaryUniformGeneratorImpl<NativeVector> tug;
    std::vector<RLWECiphertext> C(ki);
    for(int i=0;i<ki;++i){       //参与方P_i生成随机多项式a_i，加密
        NativePoly a_poly(tug,params->GetPolyParams(),Format::COEFFICIENT);
        a_poly.SetFormat(Format::EVALUATION);
        Triples[0][i] = a_poly;
        C[i] = cc.AHE_Encrypt(cc.GetParams(),a_poly,0);
    }

    //P_0计算密文和，得到a的密文
    RLWECiphertext a_sum = C[0];
    for(int i=1;i<ki;++i){
        a_sum = cc.AHE_Add(a_sum,C[i]);
    }



    //各方计算E_i=Enc(0)，随机选取b_i，计算C_i'=Add(ScalarMult(E_i,b_i),a_sum)
    std::vector<RLWECiphertext> C_prime(ki);
    for(int i=0;i<ki;++i){      
        NativePoly b_poly(tug,params->GetPolyParams(),Format::COEFFICIENT);
        b_poly.SetFormat(Format::EVALUATION);
        Triples[1][i] = b_poly;
        NativePoly zero = NativePoly(params->GetPolyParams(),Format::EVALUATION,true);
        RLWECiphertext E_i = cc.AHE_Encrypt(cc.GetParams(),zero,0);
        C_prime[i] = cc.AHE_Add(cc.AHE_ScalarMult(a_sum,b_poly),E_i);
    }

    //除P_0，其他方随机选取c_i，计算C_i''=Enc(-c_i)
    std::vector<RLWECiphertext> C_prime_prime(ki); 
    for(int i=1;i<ki;++i){
        NativePoly c_poly(tug,params->GetPolyParams(),Format::COEFFICIENT);
        c_poly.SetFormat(Format::EVALUATION);
        Triples[2][i] = c_poly;
        auto c_poly_negate = c_poly.Negate();
        C_prime_prime[i] = cc.AHE_Encrypt(cc.GetParams(),c_poly_negate,0);
    }
    //P_0计算所有C_i'的和再加上所有C_i''的和，得到ab-\sum c_i的密文，然后解密得到c_0
    RLWECiphertext ab_minus_sum_c = C_prime[0];
    for(int i=1;i<ki;++i){
        ab_minus_sum_c = cc.AHE_Add(ab_minus_sum_c,C_prime[i]);
        ab_minus_sum_c = cc.AHE_Add(ab_minus_sum_c,C_prime_prime[i]);
    }

    auto c_0 = cc.AHE_Decrypt(cc.GetParams(),ab_minus_sum_c,0,ki);
    c_0.SetFormat(Format::EVALUATION);
    Triples[2][0] = c_0;


    //验证所有a的和*所有b的和=所有c的和
    // NativePoly a_sum_poly(params->GetPolyParams(),Format::EVALUATION,true);
    // NativePoly b_sum_poly(params->GetPolyParams(),Format::EVALUATION,true);
    // NativePoly c_sum_poly(params->GetPolyParams(),Format::EVALUATION,true);
    // for(int i=0;i<ki;++i){
    //     a_sum_poly += Triples[0][i];
    //     b_sum_poly += Triples[1][i];
    //     c_sum_poly += Triples[2][i];
    // }

    // if(c_sum_poly == (a_sum_poly*b_sum_poly)){
    //     std::cout<<"验证成功"<<std::endl;
    // }else{
    //     std::cout<<"验证失败"<<std::endl;
    //     // 只打印前5个元素
    //     std::cout<<"c_sum_poly: [";
    //     auto c_values = c_sum_poly.GetValues();
    //     for(int i=0; i<std::min(5, static_cast<int>(c_values.GetLength())); i++){
    //         std::cout<<c_values[i];
    //         if(i < std::min(5, static_cast<int>(c_values.GetLength())) - 1) std::cout<<" ";
    //     }
    //     std::cout<<"]"<<std::endl;
        
    //     std::cout<<"a_sum_poly*b_sum_poly: [";
    //     auto ab_values = (a_sum_poly*b_sum_poly).GetValues();
    //     for(int i=0; i<std::min(5, static_cast<int>(ab_values.GetLength())); i++){
    //         std::cout<<ab_values[i];
    //         if(i < std::min(5, static_cast<int>(ab_values.GetLength())) - 1) std::cout<<" ";
    //     }
    //     std::cout<<"]"<<std::endl;
    // }

    return 0;
}