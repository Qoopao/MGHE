#include <bits/stdint-uintn.h>
#include "binfhe-constants.h"
#include "binfhecontext.h"
#include "rlwe-privatekey-lyh.h"
#include "lattice/hal/lat-backend.h"
#include "lwe-ciphertext-fwd.h"
#include "lwe-privatekey.h"
#include "math/hal/nativeintbackend.h"
#include "enroll-shares.h"
#include "utils/inttypes.h"
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <system_error>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace lbcrypto;
namespace fs = std::filesystem;

std::vector<std::vector<std::vector<NativePoly>>> TriplesGen();

int main() {
    // Sample Program: Step 1: Set CryptoContext，注意，大参数下有几率自举失败，原因未知，可能时hybrid product噪声上溢？或者说Q不够大，导致自举后的密文噪声>q/16
    
    auto cc = BinFHEContext();
    cc.GenerateBinFHEContext(L110T_4, XZDDF);

    //组数以及组内人数
    int k = 3;
    int ki = 3;
    //注意要先生成这个目录，否则会出问题
    cc.inituser(k, ki);

    std::vector<std::vector<LWEPrivateKey>> sk(k, std::vector<LWEPrivateKey>(ki));
    auto LWEParams = cc.GetParams()->GetLWEParams();
    

    cc.RLWECRPGen(cc.GetParams()->GetVectorNTRUParams());
    for(int i = 0; i < k; ++i){
        cc.LWECRMGen(LWEParams->Getn(), FreshEncModulus, i);
        for(int j = 0; j < ki; ++j){
            auto lwesk = cc.LWESecretKeyGen(i, j);
            auto lwepk = cc.LWEPublicKeyGen(lwesk, i, j);
            sk[i][j]=lwesk;
            auto rlwesk = cc.RLWESecretKeyGen(i, j);
            auto rlwepk = cc.RLWEPublicKeyGen(rlwesk, i, j);
        }
        auto jointpk = cc.LWEJointPublicKeyGen(cc.GetParams()->GetLWEParams(), i);
        auto jointrlwepk = cc.RLWEJointPublicKeyGen(cc.GetParams()->GetVectorNTRUParams(), i);  
    }

   
    int m0=1;
    int m1=0;
    int m2=0;
    LWEPlaintext result1;
    LWEPlaintext result2;

    // Generate the bootstrapping keys (refresh and switching keys)
    std::cout << "Generating the bootstrapping keys..." << std::endl;
    for(int i = 0; i < k; ++i){
        for(int j = 0; j < ki; ++j){
            cc.MGNBTKeyGen(sk, i, j);
        }
    }
    
    //Enroll，第i组，4个类型的分片，第j个参与方的分片
    std::vector<EnrollShares> Enroll_shares_all(k);
    for(int i=0;i<k;++i){
        auto es = EnrollShares(cc, i, ki, LWEParams->Getn());
        Enroll_shares_all[i] = es;
    }

    //重新组织顺序
    std::vector<std::vector<NativePoly>> F_shares_all(k, std::vector<NativePoly>(ki));
    std::vector<std::vector<NativePoly>> F_inv_shares_all(k, std::vector<NativePoly>(ki));
    std::vector<std::vector<std::vector<NativePoly>>> sk_sum_shares_all(k, std::vector<std::vector<NativePoly>>(cc.GetParams()->GetLWEParams()->Getn(), std::vector<NativePoly>(ki)));
    std::vector<std::vector<NativePoly>> sk_inv_all_sum_shares_all(k, std::vector<NativePoly>(ki));

    for(int i=0;i<k;++i){
        F_shares_all[i] = Enroll_shares_all[i].GetF_shares();
        F_inv_shares_all[i] = Enroll_shares_all[i].GetF_inv_shares();
        sk_sum_shares_all[i] = Enroll_shares_all[i].Getsk_inv_sum_shares();
        sk_inv_all_sum_shares_all[i] = Enroll_shares_all[i].Getsk_all_sum_shares();
    }
        

    //i=0的时候恢复一下
    // NativePoly sum(cc.GetParams()->GetVectorNTRUParams()->GetPolyParams(), Format::EVALUATION, true);
    // for(int j = 0; j < ki; ++j){
    //     sum += sk_all_sum_shares_all[0][j];
    // }
    // sum.SetFormat(Format::COEFFICIENT);
    // cout<<"sum = "<<sum<<endl;

    
    cc.GetBinFHEScheme()->GroupLWEKSKGen(cc.GetParams(),k,ki);
    cc.GetBinFHEScheme()->GroupHPKGen(cc.GetParams(), F_shares_all,k,ki);
    cc.GetBinFHEScheme()->GroupEVKGen(cc, F_shares_all, F_inv_shares_all, sk_sum_shares_all, sk_inv_all_sum_shares_all,k,ki);
    
    //cc.GetBinFHEScheme()->GroupLWEKSKGen(cc.GetParams());

    std::cout << "Completed the key generation." << std::endl;
    std::cout<<std::endl;


    // Sample Program: Step 3: Encryption
    //自举代码中有些内容不按组号严格递增读取，故groupidxs需要使用set并仔细检查
    LWECiphertext ct1 = cc.MGHE_Encrypt(m0, 0);//加密是q/4*m
    LWECiphertext ct2 = cc.MGHE_Encrypt(m1, 1);
    LWECiphertext ct3 = cc.MGHE_Encrypt(m2, 2);
    
    LWECiphertext ctNAND1;
    LWECiphertext ctNAND2;


    // Sample Program: Step 4: Evaluation
    // std::cout << "Start  the  gate bootstrapping " << std::endl;
    // ctNAND1 = cc.MGEvalBinGate(NAND, ct1, ct2, k);
    // cc.MGHE_Decrypt(sk, ctNAND1, &result1, k);
    // std::cout << "Result of encrypted computation of ( ct1:"<<m0<<" NAND ct2:"<<m1<<" ) = " << result1 << std::endl;

    ctNAND1 = cc.MGEvalBinGate(NAND, ct2, ct3, k);
    cout<<"---------第一个NAND门计算完成-------------"<<endl;
    cc.MGHE_Decrypt(sk, ctNAND1, &result1, k);
    std::cout << "Result of encrypted computation of ( ct2:"<<m1<<" NAND ct3:"<<m2<<" ) = " << result1 << std::endl;
    
    ctNAND2 = cc.MGEvalBinGate(NAND, ct1, ctNAND1, k);
    cout<<"---------第二个NAND门计算完成-------------"<<endl;
    cc.MGHE_Decrypt(sk, ctNAND2, &result2, k);
    std::cout << "Result of encrypted computation of ct1:"<<m0<<" NAND ( ct2:"<<m1<<" NAND ct3:"<<m2<<" ) = " << result2 << std::endl;
    
    return 0;
}
