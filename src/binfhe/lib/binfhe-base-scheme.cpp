//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================


/*
 * Custom Modifications:
 * - [This code is the implementation of the algorithm in the paper https://eprint.iacr.org/2023/1564]
 * 
 * This modified section follows the terms of the original BSD 2-Clause License.
 * Other modifications are provided under the terms of the BSD 2-Clause License.
 * See the BSD 2-Clause License text below:
 */


//==================================================================================
// Additional BSD License for Custom Modifications:
//
// Copyright (c) 2023 Binwu Xiang,Kaixing Wang and other contributors
//
// All rights reserved.
//
// Author TPOC: wangkaixing22@mails.ucas.ac.cn
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================


#include "binfhe-base-scheme.h"
#include <bits/stdint-uintn.h>
#include <sys/types.h>
#include "ZZ.h"
#include "ZZ_pX.h"
#include "binfhe-base-params.h"
#include "binfhe-constants.h"
#include "cereal/archives/json.hpp"
#include "lattice/hal/lat-backend.h"
#include "lwe-ciphertext-fwd.h"
#include "lwe-ciphertext.h"
#include "lwe-privatekey-fwd.h"
#include "math/discreteuniformgenerator.h"
#include "math/hal/nativeintbackend.h"
#include "math/ternaryuniformgenerator.h"
#include "ntru-ciphertext.h"
#include "rlwe-privatekey-fwd-lyh.h"
#include "rlwe-privatekey-lyh.h"
#include "rlwe-publickey-fwd-lyh.h"
#include "rlwe-publickey-lyh.h"
#include "utils/inttypes.h"
#include "utils/serializable.h"
#include "utils/serial.h"
#include "utils/sertype.h"
#include "vntru-acc.h"
#include <cereal/archives/binary.hpp>
//#include "hpk-custom-serial.h"
#include "hpk-fwd.h"
#include "hpk.h"
#include "rlwe-ciphertext.h"
#include "secret-sharing.h"


#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <iomanip>
#include <filesystem>
#include <set>
#include <vector>
#include <filesystem>

#include <chrono>
#include <iomanip>

namespace fs = std::filesystem;

namespace lbcrypto {

// wrapper for KeyGen methods
RingGSWBTKey BinFHEScheme::KeyGen(const std::shared_ptr<BinFHECryptoParams>& params, ConstLWEPrivateKey& LWEsk,
                                  KEYGEN_MODE keygenMode = SYM_ENCRYPT) const {
    const auto& LWEParams = params->GetLWEParams();

    RingGSWBTKey
        ek;  
    LWEPrivateKey skN;  
    if (keygenMode == SYM_ENCRYPT) {
        skN = LWEscheme->KeyGen(LWEParams->GetN(), LWEParams->GetQ());
    }
    else if (keygenMode == PUB_ENCRYPT) {
        ConstLWEKeyPair kpN = LWEscheme->KeyGenPair(LWEParams);
        skN                 = kpN->secretKey;
        ek.Pkey             = kpN->publicKey;
    }
    else {
        OPENFHE_THROW(config_error, "Invalid KeyGen mode");
    }
    /*-----------KSkey-------------*/
    ek.KSkey = LWEscheme->KeySwitchGen(LWEParams, LWEsk, skN);

    const auto& RGSWParams = params->GetRingGSWParams();
    const auto& polyParams = RGSWParams->GetPolyParams(); 
    NativePoly skNPoly(polyParams);
    skNPoly.SetValues(skN->GetElement(), Format::COEFFICIENT);
    skNPoly.SetFormat(Format::EVALUATION);
    /*-----------BSkey-------------*/
    ek.BSkey = ACCscheme->KeyGenAcc(RGSWParams, skNPoly, LWEsk);

    return ek;
}



VectorNTRUBTKey BinFHEScheme::NKeyGen(const std::shared_ptr<BinFHECryptoParams>& params, ConstLWEPrivateKey& LWEsk,
                                      KEYGEN_MODE keygenMode = SYM_ENCRYPT) const {
    const auto& LWEParams = params->GetLWEParams();
    VectorNTRUBTKey ek; 

    const auto& VNTRUParams = params->GetVectorNTRUParams();
    const auto& polyParams  = VNTRUParams->GetPolyParams();  
    uint64_t Q{VNTRUParams->GetQ().ConvertToInt<uint64_t>()};   //这里改成了64
    uint64_t N{VNTRUParams->GetN()};                            //这里也是

    NativeVector NatVec(N, Q);
    NativeVector NatVec_inv(N, Q);
    if(VNTRUParams->m_f!=nullptr){
        NativeVector ff=*VNTRUParams->m_f;
        NativeVector invff=*VNTRUParams->m_invf;
        for(uint64_t i=0;i<N;++i){
            NatVec[i]=ff[i];
            NatVec_inv[i]=invff[i];
        }
    }else{
        Get_invertible_NativeVector(NatVec, NatVec_inv, Q, N);
        VNTRUParams->m_f=std::make_shared<NativeVector>(NatVec);
        VNTRUParams->m_invf=std::make_shared<NativeVector>(NatVec_inv);
    }
    
    LWEPrivateKey LWEskN = std::make_shared<LWEPrivateKeyImpl>(LWEPrivateKeyImpl(NatVec));

    ek.KSkey = LWEscheme->KeySwitchGen(LWEParams, LWEsk, LWEskN);

    NativePoly skNPoly(polyParams);
    skNPoly.SetValues(NatVec, Format::COEFFICIENT);
    
    NativePoly invskNPoly(polyParams);
    invskNPoly.SetValues(NatVec_inv, Format::COEFFICIENT);

    
    skNPoly.SetFormat(Format::EVALUATION);
    invskNPoly.SetFormat(Format::EVALUATION);
    
    ek.BSkey = NACCscheme->KeyGenAcc(VNTRUParams, skNPoly, invskNPoly, LWEsk);  //ek的Pkey不赋值

    return ek;
}

void BinFHEScheme::MGNKeyGen(const std::shared_ptr<BinFHECryptoParams>& params, std::vector<std::vector<LWEPrivateKey>>& LWEsk,int groupidx, int partyidx,
                                      KEYGEN_MODE keygenMode = SYM_ENCRYPT) const {

    //cout<<"MGNKeyGen: groupidx="<<groupidx<<", partyidx="<<partyidx<<endl;
    //生成了NTRU密钥f,f->s的KSK，还有论文里面的{evk，ksk}                                    
    const auto& LWEParams = params->GetLWEParams();
    

    const auto& VNTRUParams = params->GetVectorNTRUParams();
    const auto& polyParams  = VNTRUParams->GetPolyParams();  
    NativeInteger Q{VNTRUParams->GetQ()};   //这里128位时有问题
    uint64_t N{VNTRUParams->GetN()};                            //这里也是

    NativeVector NatVec(N, Q);
    NativeVector NatVec_inv(N, Q);

    Get_invertible_NativeVector(NatVec, NatVec_inv, Q, N);
    VNTRUParams->m_f=std::make_shared<NativeVector>(NatVec);
    VNTRUParams->m_invf=std::make_shared<NativeVector>(NatVec_inv);

    // 存储m_f和m_invf到对应目录
    std::string user_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user";
    std::string save_dir = user_dir + "/group_" + std::to_string(groupidx) + "/party_" + std::to_string(partyidx);
    namespace fs = std::filesystem;
    fs::create_directories(save_dir);
    
    // 存储m_f到skN文件
    std::ofstream skN_ofs(save_dir + "/skN");
    for (size_t i = 0; i < VNTRUParams->m_f->GetLength(); ++i) {
        skN_ofs << (*VNTRUParams->m_f)[i];
        if (i + 1 != VNTRUParams->m_f->GetLength()) skN_ofs << " ";
    }
    skN_ofs << '\n';
    skN_ofs.close();
    
    // 存储m_invf到invskN文件
    std::ofstream invskN_ofs(save_dir + "/invskN");
    for (size_t i = 0; i < VNTRUParams->m_invf->GetLength(); ++i) {
        invskN_ofs << (*VNTRUParams->m_invf)[i];
        if (i + 1 != VNTRUParams->m_invf->GetLength()) invskN_ofs << " ";
    }
    invskN_ofs << '\n';
    invskN_ofs.close();
    
}


void BinFHEScheme::GroupEVKGen(BinFHEContext &cc, 
        std::vector<std::vector<NativePoly>>& F_shares_all,
        std::vector<std::vector<NativePoly>>& F_inv_shares_all,
        std::vector<std::vector<std::vector<NativePoly>>>& sk_inv_sum_shares_all, 
        std::vector<std::vector<NativePoly>>& sk_all_sum_shares_all,
        int numofgroups, int numofparties) const {

    auto params = cc.GetParams();
    const auto& VNTRUParams = params->GetVectorNTRUParams();
    const auto& polyParams = VNTRUParams->GetPolyParams();
    NativeInteger Q = VNTRUParams->GetQ();
    uint64_t N = VNTRUParams->GetN();
    uint64_t q = params->GetLWEParams()->Getq().ConvertToInt<uint64_t>();

    cout<<"    Start Group EVK Generation..."<<endl;

    //Step1: 对于每个组读取各方的LWE私钥并相加得到各组的联合LWE私钥jsk
    for (int groupIdx = 0; groupIdx < numofgroups; ++groupIdx) {
        std::string group_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user/group_" + std::to_string(groupIdx);
        if (!fs::exists(group_dir)) {
            break;  // 如果组目录不存在，则跳过
        }
        
        //Step3: 调用MGKeyGenAcc_new生成GroupEVK
        VectorNTRUACCKey GroupEVK = NACCscheme->MGKeyGenAcc_new(cc, F_shares_all[groupIdx], F_inv_shares_all[groupIdx], sk_inv_sum_shares_all[groupIdx], sk_all_sum_shares_all[groupIdx], groupIdx, numofparties);
        //Step4：将GroupEVK写入各组下的GroupEVK文件
        std::string GroupEVK_path = group_dir + "/GroupEVK";
        std::ofstream ofs(GroupEVK_path, std::ios::binary);
        if (ofs.is_open()) {
            try {
                // 检查GroupEVK是否有效
                if (!GroupEVK) {
                    std::cerr << "Error: GroupEVK is nullptr for group " << groupIdx << std::endl;
                } else {
                    // 检查内部结构是否有效
                    bool isValid = true;
                    try {
                        const auto& elements = GroupEVK->GetElements();
                        if (elements.empty() || elements[0].empty() || elements[0][0].empty()) {
                            std::cerr << "Warning: GroupEVK has empty internal structure" << std::endl;
                        } else {
                            // 检查第一个元素是否有效
                            if (!elements[0][0][0]) {
                                std::cerr << "Warning: First element of GroupEVK is nullptr" << std::endl;
                            }
                        }
                    } catch (const std::exception& e) {
                        std::cerr << "Error accessing GroupEVK elements: " << e.what() << std::endl;
                        isValid = false;
                    }
                    
                    if (isValid) {
                        // 直接传递shared_ptr给Serialize函数，不要解引用
                        Serial::Serialize(GroupEVK, ofs, SerType::BINARY);
                        //std::cout << "Successfully wrote GroupEVK to " << GroupEVK_path << std::endl;
                    }
                }
            } catch (const std::exception& e) {
                std::cerr << "Error serializing GroupEVK for group " << groupIdx << ": " << e.what() << std::endl;
                // 输出更多调试信息
                std::cerr << "Exception type: " << typeid(e).name() << std::endl;
            }
            ofs.close();
        } else {
            std::cerr << "Failed to open file for writing GroupEVK: " << GroupEVK_path << std::endl;
        }
        
    }
    
    std::cout << "    Group EVK Generation done" << std::endl;
    
}

void BinFHEScheme::GroupHPKGen(const std::shared_ptr<BinFHECryptoParams>& params, std::vector<std::vector<NativePoly>>& F_shares_all, int maxGroups, int maxParties) const {
    std::cout << "    Start Group HPK Generation..." << std::endl;



    // 获取VectorNTRU参数
    const auto vectorNTRUParams = params->GetVectorNTRUParams();
    uint32_t digitsG = vectorNTRUParams->GetDigitsG();
    auto polyParams = vectorNTRUParams->GetPolyParams();
    auto Gpow = vectorNTRUParams->GetGPower();
    NativeInteger Q = vectorNTRUParams->GetQ();

    // 从user目录下读取rlweCRP文件中的所有多项式到p数组
    std::string user_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user";
    std::string crp_path = user_dir + "/rlweCRP";
    std::vector<NativePoly> p;
    if (fs::exists(crp_path)) {
        std::ifstream ifs(crp_path);
        std::string line;
            
        while (std::getline(ifs, line)) {
            std::istringstream iss(line);
            NativeVector crp_vec(vectorNTRUParams->GetN(), Q);
            size_t j = 0;
                
            // 读取这一行的所有值
            std::string value_str;
            while (iss >> value_str && j < vectorNTRUParams->GetN()) {
                crp_vec[j] = NativeInteger(value_str);
                j++;
            }
                
            NativePoly poly(polyParams, Format::COEFFICIENT, true);
            poly.SetValues(crp_vec, Format::COEFFICIENT);
            p.push_back(poly);
        }
            
        ifs.close();
    }else {
        std::cerr << "Warning: rlweCRP文件不存在于路径" << crp_path << std::endl;
    }

    for (int groupIdx = 0; groupIdx < maxGroups; groupIdx++) {
        // 检查组目录是否存在
        std::string groupDir = "/vscode/myProgram/RRRREALL/src/binfhe/user/group_" + std::to_string(groupIdx);
        if (!fs::exists(groupDir)) {
            break; // 结束HPK生成，因为组序号是递增的，不存在还有后续组的情况
        }
        
        // 创建HPK结构：HPK={I_0,I_1,I_2, H}
        // 使用嵌套的vector结构存储，第一个维度对应I_0,I_1,I_2,H，第二个维度对应不同digit
        HPK hpk = std::make_shared<HPKImpl>(Q, digitsG);

        // 生成I_0，digitsG个随机选取的多项式
        for(int partyIdx = 0; partyIdx < maxParties; partyIdx++){
            std::string partyDir = groupDir + "/party_" + std::to_string(partyIdx);
            if(!fs::exists(partyDir)){
                break;
            }

            // 生成I_0，digitsG个随机选取的多项式
            DiscreteUniformGeneratorImpl<NativeVector> dug;
            dug.SetModulus(Q);
            for (uint32_t i = 0; i < digitsG; i++) {
                NativePoly randomPoly(dug, polyParams, Format::COEFFICIENT);
                randomPoly.SetFormat(Format::EVALUATION);
                if(partyIdx == 0){
                    hpk->SetHPK(0, i, randomPoly);
                }else{
                    NativePoly tmp = hpk->GetHPK(0, i);
                    NativePoly tmp2 = tmp + randomPoly;
                    hpk->SetHPK(0, i, tmp2);
                }
                
            }
        }  

        //每个参与方生成固定的R
        std::vector<NativePoly> R_poly_all(maxParties);
        TernaryUniformGeneratorImpl<NativeVector> tug;
        for(int partyIdx = 0; partyIdx < maxParties; partyIdx++){
            NativeVector Rvec = tug.GenerateVector(vectorNTRUParams->GetN(),32);    //模数不能是NativeInteger，所以用这种方式生成
            Rvec.SwitchModulus(vectorNTRUParams->GetQ());
            R_poly_all[partyIdx] = NativePoly(polyParams, Format::COEFFICIENT, true);
            R_poly_all[partyIdx].SetValues(Rvec, Format::COEFFICIENT);
            R_poly_all[partyIdx].SetFormat(Format::EVALUATION);
        }
        
        for(int partyIdx = 0; partyIdx < maxParties; partyIdx++){
            std::string partyDir = groupDir + "/party_" + std::to_string(partyIdx);
            if(!fs::exists(partyDir)){
                break;
            }

            // 生成I_1，I_1 = p[i]*R + Gpow[i]*[F] + e, 注意I_1和I_2中的R是同一个
            for (uint32_t d = 0; d < digitsG; d++) {
                NativePoly e(vectorNTRUParams->GetDgg(), polyParams, Format::COEFFICIENT);
                p[d].SetFormat(Format::EVALUATION);

                F_shares_all[groupIdx][partyIdx].SetFormat(Format::EVALUATION);
                e.SetFormat(Format::EVALUATION);
                NativePoly I1 = p[d] * R_poly_all[partyIdx] + F_shares_all[groupIdx][partyIdx] * Gpow[d] + e ; 
                if(partyIdx == 0){
                    hpk->SetHPK(1, d, I1);
                }else{
                    NativePoly tmp = hpk->GetHPK(1, d);
                    NativePoly tmp2 = tmp + I1;
                    hpk->SetHPK(1, d, tmp2);
                }
            }
        }

        // 生成I_2：I_2 = -I_0*Z + R*Gpow[i] + e
        for(int partyIdx = 0; partyIdx < maxParties; partyIdx++){
            std::string partyDir = groupDir + "/party_" + std::to_string(partyIdx);
            if(!fs::exists(partyDir)){
                break;
            }

            //读取rlwesk
            std::string rlwesk_path = partyDir + "/rlwesk";
            NativePoly partyrlwesk(polyParams, Format::COEFFICIENT, true);
            try{
                std::ifstream z_ifs(rlwesk_path);
                cereal::JSONInputArchive archive(z_ifs);
                RLWEPrivateKey rlwesk = std::make_shared<RLWEPrivateKeyImpl>(NativePoly(polyParams, Format::COEFFICIENT, true));
                archive(rlwesk);
                partyrlwesk = rlwesk->GetElement();
            }catch(const std::exception& e){
                std::cerr << "Error: 读取" << rlwesk_path << "时发生异常: " << e.what() << std::endl;
            }
            
            // 计算I_2 = -I_0*Z + R*Gpow[i] + e   
            std::vector<NativePoly> I0 = hpk->GetHPK(0);
            partyrlwesk.SetFormat(Format::EVALUATION);
            for(uint32_t d = 0; d < digitsG; d++){
                NativePoly I2 = I0[d] * partyrlwesk;
                I2 = I2.Negate();
                NativePoly e(vectorNTRUParams->GetDgg(), polyParams, Format::COEFFICIENT);
                e.SetFormat(Format::EVALUATION);
                I2 = I2 + R_poly_all[partyIdx] * Gpow[d] + e;               
                if(partyIdx == 0){
                    hpk->SetHPK(2, d, I2);
                }else{
                    NativePoly tmp = hpk->GetHPK(2, d);
                    NativePoly tmp2 = tmp + I2;
                    hpk->SetHPK(2, d, tmp2);
                }
            }     
        }
       
        //生成H，H= -p[i]*Z + e 
        for(int partyIdx = 0; partyIdx < maxParties; partyIdx++){
            std::string partyDir = groupDir + "/party_" + std::to_string(partyIdx);
            if(!fs::exists(partyDir)){
                break;
            }

            //读取rlwesk
            std::string rlwesk_path = partyDir + "/rlwesk";
            NativePoly partyrlwesk(polyParams, Format::COEFFICIENT, true);
            try{
                std::ifstream z_ifs(rlwesk_path);
                cereal::JSONInputArchive archive(z_ifs);
                RLWEPrivateKey rlwesk = std::make_shared<RLWEPrivateKeyImpl>(NativePoly(polyParams, Format::COEFFICIENT, true));
                archive(rlwesk);
                partyrlwesk = rlwesk->GetElement();
            }catch(const std::exception& e){
                std::cerr << "Error: 读取" << rlwesk_path << "时发生异常: " << e.what() << std::endl;
            }
            partyrlwesk.SetFormat(Format::EVALUATION);

            for (uint32_t d = 0; d < digitsG; d++){
                NativePoly e(vectorNTRUParams->GetDgg(), polyParams, Format::COEFFICIENT);
                e.SetFormat(Format::EVALUATION);
                p[d].SetFormat(Format::EVALUATION);
                NativePoly H = p[d] * partyrlwesk; 
                H = H.Negate();             
                H += e;                  
                if(partyIdx == 0){
                    hpk->SetHPK(3, d, H);
                }else{
                    NativePoly tmp = hpk->GetHPK(3, d);
                    NativePoly tmp2 = tmp + H;
                    hpk->SetHPK(3, d, tmp2);
                }
            }
        }
        
        // 保存HPK到文件 
        std::string hpkFile = groupDir + "/hpk";
        try{
            ofstream hpkofs(hpkFile);
            cereal::JSONOutputArchive archive(hpkofs);
            archive(hpk);
        }
        catch (const std::exception& e) {
            std::cerr << "Error: Failed to save HPK for group " << groupIdx << " - " << e.what() << std::endl;
        }

    }

    std::cout << "    Group HPK Generation done" << std::endl;
    
}


void BinFHEScheme::GroupLWEKSKGen(const std::shared_ptr<BinFHECryptoParams>& params, int numofgroups, int numofparties) const {
    std::cout << "    Start Group LWEKSK Generation..." << std::endl;
    
    //必要的参数
    const size_t n(params->GetLWEParams()->Getn());
    const size_t N(params->GetLWEParams()->GetN());
    NativeInteger qKS(params->GetLWEParams()->GetqKS());
    NativeInteger::Integer value{1};
    NativeInteger::Integer baseKS(params->GetLWEParams()->GetBaseKS());
    const auto digitCount = 
        static_cast<size_t>(std::ceil(log(qKS.ConvertToDouble()) / log(static_cast<double>(baseKS))));
    std::vector<NativeInteger> digitsKS(digitCount);
    for (size_t i = 0; i < digitCount; ++i) {
        digitsKS[i]=value;
        value *= baseKS;
    }
    
    // 遍历所有组
    for (int groupIdx = 0; groupIdx < numofgroups; ++groupIdx) {
        std::string group_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user/group_" +
        std::to_string(groupIdx);
        // 检查组目录是否存在
        if (!fs::exists(group_dir)) {
            cerr << "Error: 组目录不存在 " << group_dir << endl;
            break;  // 如果组目录不存在，则退出
        }
    
        // Step1: 生成N个随机的矩阵A ∈ Z_{qKS}^{d×n}
        std::vector<std::vector<std::vector<NativeInteger>>> A_all(N, std::vector<std::vector<NativeInteger>>(digitCount, std::vector<NativeInteger>(n)));
        std::vector<std::vector<NativeInteger>> B_all(N, std::vector<NativeInteger>(digitCount));
        for(size_t l = 0; l < N; ++l){
            for(size_t d = 0; d < digitCount; ++d){
                B_all[l][d] = 0;
            }
        }


        for (size_t l = 0; l < N; ++l) {
            std::vector<std::vector<NativeInteger>> A(digitCount, std::vector<NativeInteger>(n));
            DiscreteUniformGeneratorImpl<NativeVector> dug;
            dug.SetModulus(qKS);
            for (size_t i = 0; i < digitCount; ++i) {
                NativeVector randomVector = dug.GenerateVector(n);
                for (size_t j = 0; j < n; ++j) {
                    A[i][j] = randomVector[j];
                }
            }
            A_all[l] = A;
        }
        
        // 步骤2: 读取组内各个party的sk（读取到n维NativeVector变量中），rlwesk（读取到N维NativeVector变量中）文件
        for (int partyIdx = 0; partyIdx < numofparties; ++partyIdx) {
            std::string party_dir = group_dir + "/party_" + std::to_string(partyIdx);
            if(!filesystem::exists(party_dir)){
                break;
            }
            
            // 读取sk文件
            // 创建NativeVector存储sk值，要注意这里模数要先设成原来的模数，后面switchmodulus才有效，否则会截断出问题！！！！！！
            NativeVector sk_vec(n, params->GetLWEParams()->Getq());
            std::string sk_path = party_dir + "/sk";
            if (fs::exists(sk_path)) {
                ifstream sk_ifs(sk_path);
                if (sk_ifs.is_open()) {
                    std::vector<std::string> values;
                    std::string val_str;
                    while (sk_ifs >> val_str) {
                        values.push_back(val_str);
                    }
                    sk_ifs.close();
                    
                    if (!values.empty()) {
                        for (size_t i = 0; i < values.size(); ++i) {
                            sk_vec[i] = NativeInteger(values[i]);
                        }
                    }else{
                        cout<<"读取sk异常"<<endl;
                    }
                }
            }
            sk_vec.SwitchModulus(qKS);
            
            // 读取rlwesk文件
            NativeVector rlwesk_vec(N, params->GetLWEParams()->GetQ());
            std::string rlwesk_path = party_dir + "/rlwesk";
            if (fs::exists(rlwesk_path)) {
                try{
                    std::ifstream rlwesk_ifs(rlwesk_path);
                    cereal::JSONInputArchive archive(rlwesk_ifs);
                    RLWEPrivateKey rlwesk = std::make_shared<RLWEPrivateKeyImpl>();
                    archive(rlwesk);
                    rlwesk_vec = rlwesk->GetElement().GetValues();
                }catch (const std::exception& e) {
                    std::cerr << "Error: Failed to read rlwesk for party " << partyIdx << " - " << e.what() << std::endl;
                }
            }
            rlwesk_vec.SwitchModulus(qKS);

            // 步骤3：对于每个A，记为A_l，计算b_l=-A_l*sk（矩阵乘法，结果为长度digitCount的NativeVector）
            DiscreteGaussianGeneratorImpl<NativeVector> dgg=params->GetLWEParams()->GetDgg();
            for (size_t l = 0; l < N; ++l) {
                // 计算b_l = -A_l * sk
                NativeVector B(digitCount, qKS);
            
                for (size_t d = 0; d < digitCount; ++d) {
                    NativeInteger tmp = NativeInteger(0);
                    for(size_t j = 0; j < n; ++j){
                        NativeInteger tmp2 = A_all[l][d][j].ModMul(sk_vec[j], qKS);
                        tmp += tmp2;
                        tmp = tmp.Mod(qKS);
                    }
                    B[d] = qKS - tmp;
                }
                
    
                // 生成噪声e
                NativeVector e = dgg.GenerateVector(digitCount, qKS);
                // 加上rlwe_vec[l] * digitsKS + e
                for (size_t d = 0; d < digitCount; ++d) {
                    B[d] = B[d].ModAdd((rlwesk_vec[l].ModMul(digitsKS[d], qKS)).ModAdd(e[d], qKS), qKS);
                }
                
                for(size_t d = 0; d < digitCount; ++d){
                    B_all[l][d] = B_all[l][d].ModAdd(B[d], qKS);
                }
            }

        }

        // 步骤4：记LWEKSK={b_l,A_l}_{1\le l\le N}存储到各个组的LWEKSK文件中
        // 保存LWEKSK到文件
        std::string lweksk_path = group_dir + "/LWEKSK";
        std::ofstream lweksk_ofs(lweksk_path);
        if (lweksk_ofs.is_open()) {
            // 先写入N, digitCount和n
            lweksk_ofs << N << " " << digitCount << " " << n << std::endl;
            
            // 写入所有的矩阵A和对应的b
            for (size_t l = 0; l < N; ++l) {
                // 写入矩阵A
                for (size_t i = 0; i < digitCount; ++i) {
                    for (size_t k = 0; k < n; ++k) {
                        lweksk_ofs << A_all[l][i][k];
                        if (k + 1 < n) {
                            lweksk_ofs << " ";
                        }
                    }
                    lweksk_ofs << std::endl;
                }
                
                // 写入对应的B向量
                for (size_t i = 0; i < digitCount; ++i) {
                    lweksk_ofs << B_all[l].at(i);
                    if (i + 1 < digitCount) {
                        lweksk_ofs << " ";
                    }
                }
                lweksk_ofs << std::endl;
            }
            
            lweksk_ofs.close();
            //std::cout << "    Successfully wrote LWEKSK to " << lweksk_path << std::endl;
        } else {
            std::cerr << "Failed to open file for writing LWEKSK: " << lweksk_path << std::endl;
        }

    }
    
    std::cout << "    Group LWEKSK Generation done" << std::endl;

}

//这里改了int的类型
void Get_invertible_NativeVector(NativeVector& NatVec, NativeVector& NatVec_inv, NativeInteger q_boot, uint64_t N) {
    //cout<<"q_boot:"<<q_boot<<endl;
    // 三值集合的均匀分布
    // uniform_int_distribution<int> ternary_sampler(-1,1);
    //正态分布
    normal_distribution<double> gaussian_sampler(0.0, 1);
    // 随机引擎
    default_random_engine rand_engine(std::chrono::system_clock::now().time_since_epoch().count());

    std::vector<NativeInteger> vec     = std::vector<NativeInteger>(N, q_boot);
    std::vector<NativeInteger> vec_inv = std::vector<NativeInteger>(N, q_boot);

    NativeInteger half_q_boot = q_boot / 2;
    //polynomial with the coefficient vector vec (will be generated later)
    ZZ_pX poly;
    
    
    ZZ_p::init(conv<ZZ>(q_boot.ToString().c_str()));
    //cout<<ZZ_p::modulus()<<endl;
    //element of Z_(q_boot)
    ZZ_p coef;
    // 使用字符串中转初始化ZZ_p上下文
    coef.init(conv<ZZ>(q_boot.ToString().c_str()));
    //the inverse of poly modulo poly_mod (will be generated later)
    ZZ_pX inv_poly;
    //random sampling
    //int sum = 0;
    while (true) {
        //create the polynomial with the coefficient vector of the desired form
        SetCoeff(poly, 0, gaussian_sampler(rand_engine));
        for (uint32_t i = 1; i < N; i++) {
            coef = gaussian_sampler(rand_engine);
            // if(coef == 0)
            // {
            //     sum++;
            // }
            SetCoeff(poly, i, coef);
        }
        //cout<<double(sum)/N<<endl;
        //test invertibility
        try {
            // static ZZ_pX get_def_poly()
            ZZ_pX def_poly;
            ZZ_p coef;
            coef.init(conv<ZZ>(q_boot.ToString().c_str()));
            coef = 1;
            SetCoeff(def_poly, 0, coef);
            SetCoeff(def_poly, N, coef);

            InvMod(inv_poly, poly, def_poly);
            break;
        }
        catch (...) {
            cout << "Polynomial " << poly << " isn't a unit" << endl;
            continue;
        }
    }
    //cout<<"poly:"<<poly<<endl;
    NativeInteger tmp_coef;
    for (long int i = 0; i <= deg(poly); i++) {       
        // 获取ZZ_p元素的内部ZZ表示，然后通过字符串转换为NativeInteger
        ZZ_p elem = poly[i];
        ZZ elem_zz = rep(elem);
        std::stringstream ss;
        ss << elem_zz;
        tmp_coef = NativeInteger(ss.str());
        //cout<<"elem,elem_zz,tmp_coef:"<<elem<<","<<elem_zz<<","<<tmp_coef<<endl;
        // if (tmp_coef > half_q_boot)
        //     tmp_coef -= q_boot;
        vec[i] = tmp_coef;
    }

    for (long int i = 0; i <= deg(inv_poly); i++) {
        // 获取ZZ_p元素的内部ZZ表示，然后通过字符串转换为NativeInteger
        ZZ_p elem = inv_poly[i];
        ZZ elem_zz = rep(elem);
        std::stringstream ss;
        ss << elem_zz;
        tmp_coef = NativeInteger(ss.str());
        // if (tmp_coef > half_q_boot)
        //     tmp_coef -= q_boot;
        vec_inv[i] = tmp_coef;
    }
    
    for (uint32_t i = 0; i < N; i++) {
        NativeInteger v     = vec[i];
        NativeInteger v_inv = vec_inv[i];
        if (v < 0)
            NatVec[i] = q_boot - typename NativeVector::Integer(-v);
        else
            NatVec[i] = typename NativeVector::Integer(v);
        if (v_inv < 0)
            NatVec_inv[i] = q_boot - typename NativeVector::Integer(-v_inv);
        else
            NatVec_inv[i] = typename NativeVector::Integer(v_inv);
    }
}

//注意这里用的是vec-NTRU的EK，所以这是XZD+的bs
LWECiphertext BinFHEScheme::EvalBinGate(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                        const VectorNTRUBTKey& EK, ConstLWECiphertext& ct1,
                                        ConstLWECiphertext& ct2) const {
    if (ct1 == ct2)
        OPENFHE_THROW(config_error, "Input ciphertexts should be independant");

    // By default, we compute XOR/XNOR using a combination of AND, OR, and NOT gates 迭代计算
    if ((gate == XOR) || (gate == XNOR)) {
        const auto& ctAND1 = EvalBinGate(params, AND, EK, ct1, EvalNOT(params, ct2));
        const auto& ctAND2 = EvalBinGate(params, AND, EK, EvalNOT(params, ct1), ct2);
        const auto& ctOR   = EvalBinGate(params, OR, EK, ctAND1, ctAND2);

        // NOT is free so there is not cost to do it an extra time for XNOR
        return (gate == XOR) ? ctOR : EvalNOT(params, ctOR);
    }

    LWECiphertext ctprep = std::make_shared<LWECiphertextImpl>(*ct1);

    //构造(5q/8,0)
    auto n = params->GetLWEParams()->Getn();
    NativeVector zero(n,0);
    uint32_t q = params->GetLWEParams()->Getq().ConvertToInt<uint32_t>();
    zero.SetModulus(q);
    NativeInteger temp_b= 5*q/8;
    LWECiphertext ct_temp = std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(zero), temp_b.Mod(q)));

    // the additive homomorphic operation for XOR/NXOR is different from the other gates we compute
    // 2*(ct1 - ct2) mod 4 for XOR, me map 1,2 -> 1 and 3,0 -> 0
    if ((gate == XOR_FAST) || (gate == XNOR_FAST)) {
        LWEscheme->EvalSubEq(ctprep, ct2);
        LWEscheme->EvalAddEq(ctprep, ctprep);
    }
    else {
        LWEscheme->EvalAddEq(ctprep, ct2);
        LWEscheme->EvalSubEq(ct_temp,ctprep);

    }
      
    auto acc{BootstrapGateCore(params, gate, EK.BSkey, ct_temp)};  //一个NativePoly，多项式

    NativePoly& accVec{acc->GetElements()};
    accVec = accVec.Transpose();
    accVec.SetFormat(Format::COEFFICIENT);
   
    // we add Q/8 to "b" to to map back to Q/4 (i.e., mod 2) arithmetic.
    const auto& LWEParams = params->GetLWEParams();
    NativeInteger Q{LWEParams->GetQ()};
    NativeInteger b{(Q >> 3) + 1};
    auto ctExt = std::make_shared<LWECiphertextImpl>(std::move(accVec.GetValues()), std::move(b));

    // Modulus switching to a middle step Q'
    auto ctMS = LWEscheme->ModSwitch(LWEParams->GetqKS(), ctExt);
    // Key switching
    auto ctKS = LWEscheme->KeySwitch(LWEParams, EK.KSkey, ctMS);
    // Modulus switching
    return LWEscheme->ModSwitch(ct1->GetModulus(), ctKS);
}

LWECiphertext BinFHEScheme::MGEvalBinGate(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                        const VectorNTRUBTKey& EK, ConstLWECiphertext& ct1,
                                        ConstLWECiphertext& ct2,int numofgroups) const {

                                        //这里的numofgroups是初始化时的组数，不是密文ct1或ct2的组数
    if (ct1 == ct2)
        OPENFHE_THROW(config_error, "Input ciphertexts should be independant");

    // By default, we compute XOR/XNOR using a combination of AND, OR, and NOT gates 迭代计算
    if ((gate == XOR) || (gate == XNOR)) {
        const auto& ctAND1 = MGEvalBinGate(params, AND, EK, ct1, EvalNOT(params, ct2),numofgroups);
        const auto& ctAND2 = MGEvalBinGate(params, AND, EK, EvalNOT(params, ct1), ct2,numofgroups);
        const auto& ctOR   = MGEvalBinGate(params, OR, EK, ctAND1, ctAND2,numofgroups);

        // NOT is free so there is not cost to do it an extra time for XNOR
        return (gate == XOR) ? ctOR : EvalNOT(params, ctOR);
    }

    
    auto n = params->GetLWEParams()->Getn();
    uint32_t q = params->GetLWEParams()->Getq().ConvertToInt<uint32_t>();

    //根据组数对两个密文补0
    std::set<int> idxs1=ct1->groupidxs;
    std::set<int> idxs2=ct2->groupidxs;


    std::set<int> idxs3;
    idxs3.insert(idxs1.begin(),idxs1.end());
    idxs3.insert(idxs2.begin(),idxs2.end());


    NativeVector newa1(numofgroups*n,params->GetLWEParams()->Getq());
    NativeVector newa2(numofgroups*n,params->GetLWEParams()->Getq());
    //置零
    for(uint32_t i=0;i<numofgroups*n;++i){
        newa1[i]=0;
        newa2[i]=0;
    }
    
    //cout<<"numofgroups="<<numofgroups<<endl;
    //两种情况，如果groupidxs只有一个元素，说明是第一次加密，密文长度仅有n+1；
    //如果groupidxs不止有一个元素，说明密文是经过一次自举的，长度有kn+1
    //补齐ct1
    int acum1=0;
    if(ct1->groupidxs.size()==1){
        for(int id:ct1->groupidxs){
            if(id>=numofgroups){cout<<"补齐ct1 error"<<endl; break;}
            for(uint32_t i=0;i<n;++i,++acum1){
                newa1[id*n+i]=ct1->GetA()[acum1];
            }
        }
    }else{
        for(int i=0;i<ct1->GetA().GetLength();++i){
            newa1[i]=ct1->GetA()[i];
        }
    }


    //补齐ct2
    int acum2=0;
    if(ct2->groupidxs.size()==1){
        for(int id:ct2->groupidxs){
            if(id>=numofgroups){cout<<"补齐ct2 error"<<endl; break;}
            for(uint32_t i=0;i<n;++i,++acum2){
                newa2[id*n+i]=ct2->GetA()[acum2];
            }
        }
    }else{
        for(int i=0;i<ct2->GetA().GetLength();++i){
            newa2[i]=ct2->GetA()[i];
        }
    }


    newa1.SetModulus(q);
    newa2.SetModulus(q);
    LWECiphertext newct1=std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(newa1), ct1->GetB().Mod(q)));
    LWECiphertext newct2=std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(newa2), ct2->GetB().Mod(q)));
    newct1->SetModulus(q);
    newct2->SetModulus(q);

// {    
//     // 添加解密逻辑，对newct1和newct2进行解密
//     LWEPlaintext result1 = 0, result2 = 0;

//     vector<vector<LWEPrivateKey>> sk111(numofgroups, vector<LWEPrivateKey>(2));

//     for (int groupIdx = 0; groupIdx < 10; ++groupIdx) {
//         std::string group_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user/group_" + std::to_string(groupIdx);
//         if (!fs::exists(group_dir)) {
//             continue;  // 如果组目录不存在，则跳过
//         }
        
//         // 初始化联合私钥jsk
//         LWEPrivateKey jsk = nullptr;
//         bool first_sk = true;
        
//         // 读取该组中所有party的私钥
//         for (int partyIdx = 0; partyIdx < 4; ++partyIdx) {
//             std::string party_sk_path = group_dir + "/party_" + std::to_string(partyIdx) + "/sk";
//             if (!fs::exists(party_sk_path)) {
//                 continue;  // 如果party私钥文件不存在，则跳过
//             }
            
//             // 从文件读取LWE私钥
//             std::ifstream sk_ifs(party_sk_path);
//             if (!sk_ifs.is_open()) {
//                 continue;
//             }
            
//             // 读取私钥向量值
//             std::vector<std::string> values;
//             std::string val_str;
//             while (sk_ifs >> val_str) {
//                 values.push_back(val_str);
//             }
//             sk_ifs.close();
            
//             if (values.empty()) {
//                 continue;
//             }
            
//             // 创建LWE私钥
//             LWEPrivateKey sk = std::make_shared<LWEPrivateKeyImpl>();
//             NativeVector sk_vec(values.size(), NativeInteger(q));
//             for (size_t i = 0; i < values.size(); ++i) {
//                 sk_vec[i] = NativeInteger(values[i]);
//             }
            
//             // 设置私钥向量（使用const_cast来绕过const限制）
//             LWEPrivateKeyImpl* sk_impl = const_cast<LWEPrivateKeyImpl*>(sk.get());
//             sk_impl->~LWEPrivateKeyImpl(); // 析构旧对象
//             new (sk_impl) LWEPrivateKeyImpl(sk_vec); // 放置新对象
            
//             // 将当前party的私钥添加到联合私钥jsk中
//             sk111[groupIdx][partyIdx] = sk;

//         }
        

//     }
//     // 对newct1和newct2进行解密
//     LWEscheme->MGHE_Decrypt(params->GetLWEParams(), sk111, newct1, &result1, numofgroups);
//     LWEscheme->MGHE_Decrypt(params->GetLWEParams(), sk111, newct2, &result2, numofgroups);
   
//     // 可选：输出解密结果
//     std::cout << "Decrypted result for newct1: " << result1 << std::endl;
//     std::cout << "Decrypted result for newct2: " << result2 << std::endl;


// } 




    LWECiphertext ctprep = std::make_shared<LWECiphertextImpl>(*newct1);

    //构造(5q/8,0,...,0) numofgroups*n个0
    NativeVector zero(n*numofgroups,0);
    zero.SetModulus(q);
    NativeInteger temp_b= 5*q/8;
    LWECiphertext ct_temp = std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(zero), temp_b.Mod(q)));

    // the additive homomorphic operation for XOR/NXOR is different from the other gates we compute
    // 2*(ct1 - ct2) mod 4 for XOR, me map 1,2 -> 1 and 3,0 -> 0
    
    if ((gate == XOR_FAST) || (gate == XNOR_FAST)) {
        LWEscheme->EvalSubEq(ctprep, newct2);
        LWEscheme->EvalAddEq(ctprep, ctprep);
        ctprep->groupidxs=idxs3;
    }
    else {        
        LWEscheme->EvalAddEq(ctprep, newct2);//ctprep=newct1+newct2
        LWEscheme->EvalSubEq(ct_temp,ctprep);//ct_temp=ct_temp-(ctprep)，即FHEW中的同态Nand门
        ct_temp->groupidxs=idxs3;
    }
    //输出ct1，ct2，newct1,newct2,ct_temp的A
    // cout<<"ct1->GetA()="<<ct1->GetA()<<endl;
    // cout<<"ct2->GetA()="<<ct2->GetA()<<endl;
    // cout<<"newct1->GetA()="<<newct1->GetA()<<endl;
    // cout<<"newct2->GetA()="<<newct2->GetA()<<endl;
    // cout<<"ct_temp->GetA()="<<ct_temp->GetA()<<endl;
    
    
    //cout<<"ct_temp"<<ct_temp->GetA();
    // cout<<"ct_temp->groupidxs="<<endl;
    // for(int id:ct_temp->groupidxs){
    //     cout<<id<<endl;
    // }

    //读取各个GroupEVK
    std::vector<VectorNTRUACCKey> EVK_ALL(numofgroups);
    for(int index=0;index<numofgroups;++index){
    //读取第index组的GroupEVK
        std::string group_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user/group_" + std::to_string(index);
        std::string GroupEVK_path = group_dir + "/GroupEVK";
    
        VectorNTRUACCKey ek;
        if (!Serial::DeserializeFromFile(GroupEVK_path, ek, SerType::BINARY)) {
            std::cerr << "Failed to read GroupEVK from file: " << GroupEVK_path << std::endl;
            // 可以选择抛出异常或返回空指针，但根据上下文我们继续执行，已经重定向试过了，evk的写入读取没问题
        }
        EVK_ALL[index]=ek;
    }
    
    //读取各组的hpk
    std::vector<HPK> HPK_ALL(numofgroups, std::make_shared<HPKImpl>(params->GetVectorNTRUParams()->GetQ(), params->GetVectorNTRUParams()->GetDigitsG()));
    for (int id = 0; id < numofgroups; ++id) {
        std::string hpk_path = std::string("/vscode/myProgram/RRRREALL/src/binfhe/user/group_") + std::to_string(id) + "/hpk";
        // 读取hpk数据 
        try {
            HPK hpkgroup = std::make_shared<HPKImpl>(params->GetVectorNTRUParams()->GetQ(), params->GetVectorNTRUParams()->GetDigitsG());
            std::ifstream hpk_file(hpk_path, std::ios::in | std::ios::binary);
            cereal::JSONInputArchive archive(hpk_file);
            archive(hpkgroup);
            HPK_ALL[id]=hpkgroup;
        }
        catch (const std::exception& e) {
            OPENFHE_THROW(deserialize_error, "Failed to deserialize hpk for group " + std::to_string(id) + ": " + e.what());
        }
  
    }

    //读取NegativeCRP
    // 从rlweCRP文件中读取p的值并取负
    std::vector<NativePoly> negativeCRP(params->GetVectorNTRUParams()->GetDigitsG(), NativePoly(params->GetVectorNTRUParams()->GetPolyParams(), Format::EVALUATION, true));
    std::string user_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user";
    std::string crp_path = user_dir + "/rlweCRP";
    if (std::filesystem::exists(crp_path)) {
        std::ifstream ifs(crp_path);
        for (uint32_t d = 0; d < params->GetVectorNTRUParams()->GetDigitsG(); ++d) {
            NativeVector p_vec(params->GetVectorNTRUParams()->GetN(), params->GetVectorNTRUParams()->GetQ());
            for (size_t j = 0; j < params->GetVectorNTRUParams()->GetN(); ++j) {
                std::string value_str;
                ifs >> value_str;
                p_vec[j] = NativeInteger(value_str);
            }
            NativePoly pd(params->GetVectorNTRUParams()->GetPolyParams(), Format::COEFFICIENT, true);
            pd.SetValues(p_vec, Format::COEFFICIENT);
            pd.SetFormat(Format::EVALUATION);
            pd = pd.Negate();
            negativeCRP[d] = pd; //读取取负没问题111
        }
        ifs.close();
    }
    else {
        // 如果rlweCRP文件不存在
        std::cerr << "警告：rlweCRP文件不存在于路径" << crp_path << std::endl;
    }

    //读取LWEKSK
    auto lweparams = params->GetLWEParams();
    const size_t nn(numofgroups * lweparams->Getn());
    const size_t N(lweparams->GetN());
    NativeInteger qKS(lweparams->GetqKS());
    NativeInteger::Integer baseKS(lweparams->GetBaseKS());
    const auto digitCount = static_cast<size_t>(std::ceil(log(qKS.ConvertToDouble()) / log(static_cast<double>(baseKS))));
    //读取每个组的A_all和B_all
    //A_all的第一个维度指示组索引，第二个维度是N，第三维度指示矩阵的行；B_all的第一个维度指示组索引第二个维度是N；
    std::vector<std::vector<std::vector<std::vector<NativeInteger>>>> A_all(numofgroups, std::vector<std::vector<std::vector<NativeInteger>>>(N, std::vector<std::vector<NativeInteger>>(digitCount, std::vector<NativeInteger>(lweparams->Getn()))));
    std::vector<std::vector<std::vector<NativeInteger>>> B_all(numofgroups, std::vector<std::vector<NativeInteger>>(N, std::vector<NativeInteger>(digitCount)));
    //这里空间应该开numofgroups，LWEKSKGen函数是对user文件下所有组生成的
    A_all.resize(numofgroups);
    B_all.resize(numofgroups);
    //cout<<"KSK numofgroups="<<numofgroups<<endl;
    // 为每个组加载LWEKSK, 注意组号严格递增
    for (int groupIdx = 0; groupIdx < numofgroups; ++groupIdx) {
        std::string group_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user/group_" + std::to_string(groupIdx);
        if (!std::filesystem::exists(group_dir)) {
            continue;  // 如果组目录不存在，则跳过
        }
        std::string lweksk_path = group_dir + "/LWEKSK";
        
        std::ifstream lweksk_ifs(lweksk_path);
        if (!lweksk_ifs.is_open()) {
            throw std::runtime_error("Failed to open LWEKSK file for group " + std::to_string(groupIdx) + ": " + lweksk_path);
        }
        
        size_t file_N, file_digitCount, file_n;
        lweksk_ifs >> file_N >> file_digitCount >> file_n;
        
        // 确保文件中的参数与当前参数匹配
        if (file_N != N || file_digitCount != digitCount) {
            throw std::runtime_error("LWEKSK parameters do not match for group " + std::to_string(groupIdx));
        }
        
        // 为当前组创建A和B的存储空间
        std::vector<std::vector<std::vector<NativeInteger>>> A_group(N, std::vector<std::vector<NativeInteger>>(digitCount, std::vector<NativeInteger>(lweparams->Getn())));
        std::vector<std::vector<NativeInteger>> B_group(N, std::vector<NativeInteger>(digitCount));
        
        // 读取矩阵A和向量B - 严格按照GroupLWEKSKGen的写入格式
        for (size_t l = 0; l < N; ++l) {
            // 读取矩阵A_all[l]的digitCount行
            for (size_t d = 0; d < digitCount; ++d) {
                // 读取params->Getn()个元素
                for (size_t j = 0; j < lweparams->Getn(); ++j) {
                    string val;
                    lweksk_ifs >> val;
                    A_group[l][d][j] = NativeInteger(val);
                }
            }
            
            // 读取向量B_all[l]的digitCount个元素
            for (size_t d = 0; d < digitCount; ++d) {
                string val;
                lweksk_ifs >> val;
                B_group[l][d] = NativeInteger(val);
            }
        }
        
        // 将当前组的A和B添加到A_all和B_all中
        A_all[groupIdx] = A_group;
        B_all[groupIdx] = B_group;
    }



//----------------------计时开始-------------------------

auto start = std::chrono::high_resolution_clock::now();

    //盲旋转
    auto accRLWEvec=MGBootstrapGateCore_new(params, gate, ct_temp, EVK_ALL, HPK_ALL, negativeCRP, numofgroups);  
    
    //sample extraction
    auto accVec = accRLWEvec;
    accVec[0].SetFormat(Format::COEFFICIENT);
    for(int i=1;i<=numofgroups;++i){
        accVec[i]=accVec[i].Transpose();
        accVec[i].SetFormat(Format::COEFFICIENT);
    }

    // we add Q/8 to "b" to to map back to Q/4 (i.e., mod 2) arithmetic.
    // 由于r(X)的设置，m的值只能是1或7（-1），所以加上Q/8之后可以将delta变成Q/4
    const auto& LWEParams = params->GetLWEParams();
    NativeInteger Q{LWEParams->GetQ()};
    NativeInteger b{(Q >> 3) + 1};
    b.ModAddEq(accVec[0][0],Q);

    NativeVector aa(numofgroups*LWEParams->GetN(),Q);
    for(int i=0;i<numofgroups;++i){ 
        for(uint32_t j=0;j<LWEParams->GetN();++j){
            aa[i*LWEParams->GetN()+j]=accVec[i+1][j];
        }
    }
    
    auto ctExt = std::make_shared<LWECiphertextImpl>(std::move(aa), std::move(b));
    ctExt->groupidxs=idxs3;


    // Modulus switching to a middle step Q', 模切换函数按照ctExt的长度进行模切换，所以不需要写一个MG模切换函数
    auto ctMS = LWEscheme->ModSwitch(LWEParams->GetqKS(), ctExt);
    ctMS->groupidxs=idxs3;

    // Key switching
    auto ctKS = LWEscheme->MGKeySwitch_new(LWEParams, ctMS, A_all, B_all, numofgroups);
    ctKS->groupidxs=idxs3;
    

    // Modulus switching
    auto ctMS2 = LWEscheme->ModSwitch(LWEParams->Getq(), ctKS);
    ctMS2->groupidxs=idxs3;

//-----------------------------计时结束-----------------------------
auto end = std::chrono::high_resolution_clock::now();

    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    double elapsed_time = duration_ms.count() / 1000.0;
    
    // 输出结果
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "---------------------------执行一次自举用时：" << elapsed_time << " 秒" << std::endl;
    
    return ctMS2;
}

//修改的··
LWECiphertext BinFHEScheme::EvalBinGateMD(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                        std::vector<VectorNTRUBTKey>& EKvec, ConstLWECiphertext& ct1,
                                        ConstLWECiphertext& ct2,vector<NativeVector>& fvec) const{
    if (ct1 == ct2)
        OPENFHE_THROW(config_error, "Input ciphertexts should be independant");

    //构造(5q/8,0)
    int k=EKvec.size();
    auto n = k*params->GetLWEParams()->Getn();
    NativeVector zero(n,0);
    uint32_t q = params->GetLWEParams()->Getq().ConvertToInt<uint32_t>();
    zero.SetModulus(q);
    NativeInteger temp_b= 5*q/8;
    LWECiphertext ct_temp = std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(zero), temp_b.Mod(q)));

    LWECiphertext ctprep = std::make_shared<LWECiphertextImpl>(*ct1);
    // the additive homomorphic operation for XOR/NXOR is different from the other gates we compute
    // 2*(ct1 - ct2) mod 4 for XOR, me map 1,2 -> 1 and 3,0 -> 0
    if ((gate == XOR_FAST) || (gate == XNOR_FAST)) {
        LWEscheme->EvalSubEq(ctprep, ct2);
        LWEscheme->EvalAddEq(ctprep, ctprep);
    }
    else {
        //计算NAND gate, 5q/8-ct1-ct2
        LWEscheme->EvalAddEq(ctprep, ct2);
        LWEscheme->EvalSubEq(ct_temp,ctprep);
    }

    //上面已经完成了计算与非门，接下来是bootstrapping


    //bootstrapping
    std::vector<VectorNTRUACCKey> bskvec;
    for(int i=0;i<k;++i){
        bskvec.push_back(EKvec[i].BSkey);
    }

    auto acc{BootstrapGateCoreMD(params, gate, bskvec, ct_temp,fvec)};

    //用F解密acc；k=1时，首项为1，ext后再解密的结果是0，首项为7解密结果为1；k=2时，F*acc*1/三角 的首项都是0，导致ext后解密的值都是1
    auto polyParams  = params->GetVectorNTRUParams()->GetPolyParams();
    NativePoly F(polyParams,Format::COEFFICIENT);
    F.SetValues(fvec[0],Format::COEFFICIENT);
    F.SetFormat(Format::EVALUATION);

    // F.SetFormat(Format::COEFFICIENT);
    // cout<<"F"<<endl;
    // for(uint32_t j=0;j<params->GetLWEParams()->GetN();++j){
    //     cout<<F[j]<<" ";
    // }
    // cout<<endl;
    // F.SetFormat(Format::EVALUATION);

    NativePoly accPolynomial{acc->GetElements()};
    NativePoly rresult=F*accPolynomial;
    rresult.SetFormat(Format::COEFFICIENT);

    // std::cout<<"------------F*accPoly-------------"<<std::endl;      
    // for(uint32_t j=0;j<params->GetLWEParams()->GetN();++j){
    //     std::cout<<rresult[j]<<" ";
    // }
    // std::cout<<std::endl;

    // std::cout<<"------------F*accPoly-rounding----------"<<std::endl;
    // for(uint32_t j=0;j<params->GetLWEParams()->GetN();++j){
    //     NativeInteger r=rresult[j];
    //     rresult[j]=r.MultiplyAndRound(NativeInteger(8),params->GetLWEParams()->GetQ()).ConvertToInt();
    //     cout<<rresult[j]<<" ";
    // }
    // std::cout<<std::endl;


   
    //acc的值是NTT形式的数组，这里用Transpose只是将多项式系数在NTT形式下进行排序取反？
    NativePoly& accVec{acc->GetElements()};
    accVec = accVec.Transpose();
    accVec.SetFormat(Format::COEFFICIENT);
    
    // we add Q/8 to "b" to to map back to Q/4 (i.e., mod 2) arithmetic.
    const auto& LWEParams = params->GetLWEParams();
    NativeInteger Q{LWEParams->GetQ()};
    NativeInteger b{(Q >> 3) + 1};
    auto ctlwe = std::make_shared<LWECiphertextImpl>(std::move(accVec.GetValues()), std::move(b));

    //用F_i解密
    {
        const NativeInteger& mod = (const NativeInteger&) q;    //用于密文的模切换等
        auto ctdet = LWEscheme->ModSwitch(q,ctlwe);              //经过ext的密文进行模切换
        NativeVector adet   = ctdet->GetA();                      //获取模切换后的a向量 (a decrypt test)
        F.SetFormat(Format::COEFFICIENT);                          //获取用于解密的F多项式的coefficient态
        
        uint32_t N = params->GetLWEParams()->GetN();
        NativeVector Fvec(N,Q);
        for(uint32_t index=0;index<params->GetLWEParams()->GetN();++index){
            Fvec[index]=F[index];
        }                                                          //得到F多项式的系数，存放在Fvec中

        
        NativeInteger mu = mod.ComputeMu();
        Fvec.SwitchModulus(mod);
        
        NativeInteger inner(0);
        for (size_t j = 0; j < params->GetLWEParams()->GetN(); ++j) {
            inner += adet[j].ModMulFast(Fvec[j], mod, mu);
        }
        inner.ModEq(mod);                                       //=inner%mod
    
        NativeInteger r = ctdet->GetB();   
        r.ModSubFastEq(inner, mod);                             // r = B- \sum{ai*si} % q
        r.ModAddFastEq((mod / (4 * 2)), mod);                   //rounding
        int result = ((NativeInteger(4) * r) / mod).ConvertToInt();
        cout<<"用F向量当作sk解密经过ext的密文结果是: "<<result<<endl;
    }


    
    

    //auto ctMS = LWEscheme->ModSwitch(LWEParams->GetqKS(), ctlwe);
    //auto ctKS = LWEscheme->KeySwitch(LWEParams, EKvec[0].KSkey, ctMS);
    //return LWEscheme->ModSwitch(ct1->GetModulus(), ctKS);
    return LWEscheme->ModSwitch(q,ctlwe);
    // return ctExt;

}



// Full evaluation as described in https://eprint.iacr.org/2020/086
LWECiphertext BinFHEScheme::EvalBinGate(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                        const RingGSWBTKey& EK, ConstLWECiphertext& ct1,
                                        ConstLWECiphertext& ct2) const {
    if (ct1 == ct2)
        OPENFHE_THROW(config_error, "Input ciphertexts should be independant");

    // By default, we compute XOR/XNOR using a combination of AND, OR, and NOT gates 迭代计算
    if ((gate == XOR) || (gate == XNOR)) {
        const auto& ctAND1 = EvalBinGate(params, AND, EK, ct1, EvalNOT(params, ct2));
        const auto& ctAND2 = EvalBinGate(params, AND, EK, EvalNOT(params, ct1), ct2);
        const auto& ctOR   = EvalBinGate(params, OR, EK, ctAND1, ctAND2);

        // NOT is free so there is not cost to do it an extra time for XNOR
        return (gate == XOR) ? ctOR : EvalNOT(params, ctOR);
    }

    LWECiphertext ctprep = std::make_shared<LWECiphertextImpl>(*ct1);
    // the additive homomorphic operation for XOR/NXOR is different from the other gates we compute
    // 2*(ct1 - ct2) mod 4 for XOR, me map 1,2 -> 1 and 3,0 -> 0
    if ((gate == XOR_FAST) || (gate == XNOR_FAST)) {
        LWEscheme->EvalSubEq(ctprep, ct2);
        LWEscheme->EvalAddEq(ctprep, ctprep);
    }
    else {
        // for all other gates, we simply compute (ct1 + ct2) mod 4
        // for AND: 0,1 -> 0 and 2,3 -> 1
        // for OR: 1,2 -> 1 and 3,0 -> 0
        LWEscheme->EvalAddEq(ctprep, ct2);
    }

    auto acc{BootstrapGateCore(params, gate, EK.BSkey, ctprep)};

    // the accumulator result is encrypted w.r.t. the transposed secret key
    // we can transpose "a" to get an encryption under the original secret key
    std::vector<NativePoly>& accVec{acc->GetElements()};
    accVec[0] = accVec[0].Transpose();  //对a重排列
    accVec[0].SetFormat(Format::COEFFICIENT);
    accVec[1].SetFormat(Format::COEFFICIENT);

    // we add Q/8 to "b" to to map back to Q/4 (i.e., mod 2) arithmetic.
    const auto& LWEParams = params->GetLWEParams();
    NativeInteger Q{LWEParams->GetQ()};
    NativeInteger b{(Q >> 3) + 1};
    b.ModAddFastEq(accVec[1][0], Q);

    auto ctExt = std::make_shared<LWECiphertextImpl>(std::move(accVec[0].GetValues()), std::move(b));
    // Modulus switching to a middle step Q'
    auto ctMS = LWEscheme->ModSwitch(LWEParams->GetqKS(), ctExt);
    // Key switching
    auto ctKS = LWEscheme->KeySwitch(LWEParams, EK.KSkey, ctMS);
    // Modulus switching
    return LWEscheme->ModSwitch(ct1->GetModulus(), ctKS);
}

// Full evaluation as described in https://eprint.iacr.org/2020/086
LWECiphertext BinFHEScheme::EvalBinGate(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                        const RingGSWBTKey& EK, const std::vector<LWECiphertext>& ctvector) const {
    // check if the ciphertexts are all independent
    for (size_t i = 0; i < ctvector.size(); i++) {
        for (size_t j = i + 1; j < ctvector.size(); j++) {
            if (ctvector[j] == ctvector[i]) {
                OPENFHE_THROW(config_error, "Input ciphertexts should be independent");
            }
        }
    }

    NativeInteger p = ctvector[0]->GetptModulus();

    LWECiphertext ctprep = std::make_shared<LWECiphertextImpl>(*ctvector[0]);
    ctprep->SetptModulus(p);
    if ((gate == MAJORITY) || (gate == AND3) || (gate == OR3) || (gate == AND4) || (gate == OR4)) {
        // we simply compute sum(ctvector[i]) mod p
        for (size_t i = 1; i < ctvector.size(); i++) {
            LWEscheme->EvalAddEq(ctprep, ctvector[i]);
        }
        auto acc = BootstrapGateCore(params, gate, EK.BSkey, ctprep);

        std::vector<NativePoly>& accVec = acc->GetElements();
        // the accumulator result is encrypted w.r.t. the transposed secret key
        // we can transpose "a" to get an encryption under the original secret key
        accVec[0] = accVec[0].Transpose();
        accVec[0].SetFormat(Format::COEFFICIENT);
        accVec[1].SetFormat(Format::COEFFICIENT);

        // we add Q/8 to "b" to to map back to Q/4 (i.e., mod 2) arithmetic.
        auto& LWEParams = params->GetLWEParams();
        NativeInteger Q = LWEParams->GetQ();
        NativeInteger b = Q / NativeInteger(2 * p) + 1;
        b.ModAddFastEq(accVec[1][0], Q);

        auto ctExt = std::make_shared<LWECiphertextImpl>(std::move(accVec[0].GetValues()), std::move(b));
        // Modulus switching to a middle step Q'
        auto ctMS = LWEscheme->ModSwitch(LWEParams->GetqKS(), ctExt);
        // Key switching
        auto ctKS = LWEscheme->KeySwitch(LWEParams, EK.KSkey, ctMS);
        // Modulus switching
        return LWEscheme->ModSwitch(ctvector[0]->GetModulus(), ctKS);
    }
    else if (gate == CMUX) {
        if (ctvector.size() != 3)
            OPENFHE_THROW(not_implemented_error, "CMUX gate implemented for ciphertext vectors of size 3");

        auto ccNOT   = EvalNOT(params, ctvector[2]);
        auto ctNAND1 = EvalBinGate(params, NAND, EK, ctvector[0], ccNOT);
        auto ctNAND2 = EvalBinGate(params, NAND, EK, ctvector[1], ctvector[2]);
        auto ctCMUX  = EvalBinGate(params, NAND, EK, ctNAND1, ctNAND2);
        return ctCMUX;
    }
    else {
        OPENFHE_THROW(not_implemented_error, "This gate is not implemented for vector of ciphertexts at this time");
    }
}

// Full evaluation as described in https://eprint.iacr.org/2020/086
LWECiphertext BinFHEScheme::Bootstrap(const std::shared_ptr<BinFHECryptoParams>& params, const RingGSWBTKey& EK,
                                      ConstLWECiphertext& ct) const {
    NativeInteger p = ct->GetptModulus();
    LWECiphertext ctprep{std::make_shared<LWECiphertextImpl>(*ct)};
    // ctprep = ct + q/4
    LWEscheme->EvalAddConstEq(ctprep, (ct->GetModulus() >> 2));

    auto acc{BootstrapGateCore(params, AND, EK.BSkey, ctprep)};

    // the accumulator result is encrypted w.r.t. the transposed secret key
    // we can transpose "a" to get an encryption under the original secret key
    std::vector<NativePoly>& accVec{acc->GetElements()};
    accVec[0] = accVec[0].Transpose();
    accVec[0].SetFormat(Format::COEFFICIENT);
    accVec[1].SetFormat(Format::COEFFICIENT);

    // we add Q/8 to "b" to to map back to Q/4 (i.e., mod 2) arithmetic.
    const auto& LWEParams = params->GetLWEParams();
    NativeInteger Q{LWEParams->GetQ()};
    NativeInteger b = Q / NativeInteger(2 * p) + 1;
    b.ModAddFastEq(accVec[1][0], Q);

    auto ctExt = std::make_shared<LWECiphertextImpl>(std::move(accVec[0].GetValues()), std::move(b));
    // Modulus switching to a middle step Q'
    auto ctMS = LWEscheme->ModSwitch(LWEParams->GetqKS(), ctExt);
    // Key switching
    auto ctKS = LWEscheme->KeySwitch(LWEParams, EK.KSkey, ctMS);
    // Modulus switching
    return LWEscheme->ModSwitch(ct->GetModulus(), ctKS);
}

// Evaluation of the NOT operation; no key material is needed
LWECiphertext BinFHEScheme::EvalNOT(const std::shared_ptr<BinFHECryptoParams>& params, ConstLWECiphertext& ct) const {
    NativeInteger q{ct->GetModulus()};
    uint32_t n{ct->GetLength()};

    NativeVector a(n, q);
    for (uint32_t i = 0; i < n; ++i)
        a[i] = ct->GetA(i) == 0 ? 0 : q - ct->GetA(i);

    return std::make_shared<LWECiphertextImpl>(std::move(a), (q >> 2).ModSubFast(ct->GetB(), q));
}

// Evaluate Arbitrary Function homomorphically
// Modulus of ct is q | 2N
LWECiphertext BinFHEScheme::EvalFunc(const std::shared_ptr<BinFHECryptoParams>& params, const RingGSWBTKey& EK,
                                     ConstLWECiphertext& ct, const std::vector<NativeInteger>& LUT,
                                     const NativeInteger& beta) const {
    auto ct1 = std::make_shared<LWECiphertextImpl>(*ct);
    NativeInteger q{ct->GetModulus()};
    uint32_t functionProperty{this->checkInputFunction(LUT, q)};

    if (functionProperty == 0) {  // negacyclic function only needs one bootstrap
        auto fLUT = [LUT](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
            return LUT[x.ConvertToInt()];
        };
        LWEscheme->EvalAddConstEq(ct1, beta);
        return BootstrapFunc(params, EK, ct1, fLUT, q);
    }

    if (functionProperty == 2) {  // arbitary funciton
        const auto& LWEParams = params->GetLWEParams();
        uint32_t N{LWEParams->GetN()};
        if (q.ConvertToInt() > N) {  // need q to be at most = N for arbitary function
            std::string errMsg =
                "ERROR: ciphertext modulus q needs to be <= ring dimension for arbitrary function evaluation";
            OPENFHE_THROW(not_implemented_error, errMsg);
        }

        // TODO: figure out a way to not do this :(

        // repeat the LUT to make it periodic
        std::vector<NativeInteger> LUT2 = LUT;
        LUT2.insert(LUT2.end(), LUT.begin(), LUT.end());

        NativeInteger dq{q << 1};
        // raise the modulus of ct1 : q -> 2q
        ct1->GetA().SetModulus(dq);

        auto ct2 = std::make_shared<LWECiphertextImpl>(*ct1);
        LWEscheme->EvalAddConstEq(ct2, beta);
        // this is 1/4q_small or -1/4q_small mod q
        auto f0 = [](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
            if (x < (q >> 1))
                return Q - (q >> 2);
            else
                return (q >> 2);
        };
        auto ct3 = BootstrapFunc(params, EK, ct2, f0, dq);
        LWEscheme->EvalSubEq2(ct1, ct3);
        LWEscheme->EvalAddConstEq(ct3, beta);
        LWEscheme->EvalSubConstEq(ct3, q >> 1);

        // Now the input is within the range [0, q/2).
        // Note that for non-periodic function, the input q is boosted up to 2q
        auto fLUT2 = [LUT2](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
            if (x < (q >> 1))
                return LUT2[x.ConvertToInt()];
            else
                return Q - LUT2[x.ConvertToInt() - q.ConvertToInt() / 2];
        };
        auto ct4 = BootstrapFunc(params, EK, ct3, fLUT2, dq);
        ct4->SetModulus(q);
        return ct4;
    }

    // Else it's periodic function so we evaluate directly
    LWEscheme->EvalAddConstEq(ct1, beta);
    // this is 1/4q_small or -1/4q_small mod q
    auto f0 = [](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
        if (x < (q >> 1))
            return Q - (q >> 2);
        else
            return (q >> 2);
    };
    auto ct2 = BootstrapFunc(params, EK, ct1, f0, q);
    LWEscheme->EvalSubEq2(ct, ct2);
    LWEscheme->EvalAddConstEq(ct2, beta);
    LWEscheme->EvalSubConstEq(ct2, q >> 2);

    // Now the input is within the range [0, q/2).
    // Note that for non-periodic function, the input q is boosted up to 2q
    auto fLUT1 = [LUT](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
        if (x < (q >> 1))
            return LUT[x.ConvertToInt()];
        else
            return Q - LUT[x.ConvertToInt() - q.ConvertToInt() / 2];
    };
    return BootstrapFunc(params, EK, ct2, fLUT1, q);
}

// Evaluate Homomorphic Flooring
LWECiphertext BinFHEScheme::EvalFloor(const std::shared_ptr<BinFHECryptoParams>& params, const RingGSWBTKey& EK,
                                      ConstLWECiphertext& ct, const NativeInteger& beta, uint32_t roundbits) const {
    const auto& LWEParams = params->GetLWEParams();
    NativeInteger q{roundbits == 0 ? LWEParams->Getq() : beta * (1 << roundbits + 1)};
    NativeInteger mod{ct->GetModulus()};

    auto ct1 = std::make_shared<LWECiphertextImpl>(*ct);
    LWEscheme->EvalAddConstEq(ct1, beta);

    auto ct1Modq = std::make_shared<LWECiphertextImpl>(*ct1);
    ct1Modq->SetModulus(q);
    // this is 1/4q_small or -1/4q_small mod q
    auto f1 = [](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
        if (x < (q >> 1))
            return Q - (q >> 2);
        else
            return (q >> 2);
    };
    auto ct2 = BootstrapFunc(params, EK, ct1Modq, f1, mod);
    LWEscheme->EvalSubEq(ct1, ct2);

    auto ct2Modq = std::make_shared<LWECiphertextImpl>(*ct1);
    ct2Modq->SetModulus(q);

    // now the input is only within the range [0, q/2)
    auto f2 = [](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
        if (x < (q >> 2))
            return Q - (q >> 1) - x;
        else if (((q >> 2) <= x) && (x < 3 * (q >> 2)))
            return x;
        else
            return Q + (q >> 1) - x;
    };
    auto ct3 = BootstrapFunc(params, EK, ct2Modq, f2, mod);
    LWEscheme->EvalSubEq(ct1, ct3);

    return ct1;
}

// Evaluate large-precision sign
LWECiphertext BinFHEScheme::EvalSign(const std::shared_ptr<BinFHECryptoParams>& params,
                                     const std::map<uint32_t, RingGSWBTKey>& EKs, ConstLWECiphertext& ct,
                                     const NativeInteger& beta, bool schemeSwitch) const {
    auto mod{ct->GetModulus()};
    const auto& LWEParams = params->GetLWEParams();
    auto q{LWEParams->Getq()};
    if (mod <= q) {
        std::string errMsg =
            "ERROR: EvalSign is only for large precision. For small precision, please use bootstrapping directly";
        OPENFHE_THROW(not_implemented_error, errMsg);
    }

    const auto& RGSWParams = params->GetRingGSWParams();
    const auto curBase     = RGSWParams->GetBaseG();
    auto search            = EKs.find(curBase);
    if (search == EKs.end()) {
        std::string errMsg("ERROR: No key [" + std::to_string(curBase) + "] found in the map");
        OPENFHE_THROW(openfhe_error, errMsg);
    }
    RingGSWBTKey curEK(search->second);

    auto cttmp = std::make_shared<LWECiphertextImpl>(*ct);
    while (mod > q) {
        cttmp = EvalFloor(params, curEK, cttmp, beta);
        // round Q to 2betaQ/q
        //  mod   = mod / q * 2 * beta;
        mod   = (mod << 1) * beta / q;
        cttmp = LWEscheme->ModSwitch(mod, cttmp);

        // if dynamic
        if (EKs.size() == 3) {
            // TODO: use GetMSB()?
            uint32_t binLog = static_cast<uint32_t>(ceil(GetMSB(mod.ConvertToInt()) - 1));
            uint32_t base{0};
            if (binLog <= static_cast<uint32_t>(17))
                base = static_cast<uint32_t>(1) << 27;
            else if (binLog <= static_cast<uint32_t>(26))
                base = static_cast<uint32_t>(1) << 18;

            if (0 != base) {  // if base is to change ...
                RGSWParams->Change_BaseG(base);

                auto search = EKs.find(base);
                if (search == EKs.end()) {
                    std::string errMsg("ERROR: No key [" + std::to_string(curBase) + "] found in the map");
                    OPENFHE_THROW(openfhe_error, errMsg);
                }
                curEK = search->second;
            }
        }
    }
    LWEscheme->EvalAddConstEq(cttmp, beta);

    if (!schemeSwitch) {
        // if the ended q is smaller than q, we need to change the param for the final boostrapping
        auto f3 = [](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
            return (x < q / 2) ? (Q / 4) : (Q - Q / 4);
        };
        cttmp = BootstrapFunc(params, curEK, cttmp, f3, q);  // this is 1/4q_small or -1/4q_small mod q
        LWEscheme->EvalSubConstEq(cttmp, q >> 2);
    }
    else {  // return the negated f3 and do not subtract q/4 for a more natural encoding in scheme switching
        // if the ended q is smaller than q, we need to change the param for the final boostrapping
        auto f3 = [](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
            return (x < q / 2) ? (Q - Q / 4) : (Q / 4);
        };
        cttmp = BootstrapFunc(params, curEK, cttmp, f3, q);  // this is 1/4q_small or -1/4q_small mod q
    }
    RGSWParams->Change_BaseG(curBase);
    return cttmp;
}

// Evaluate Ciphertext Decomposition
std::vector<LWECiphertext> BinFHEScheme::EvalDecomp(const std::shared_ptr<BinFHECryptoParams>& params,
                                                    const std::map<uint32_t, RingGSWBTKey>& EKs, ConstLWECiphertext& ct,
                                                    const NativeInteger& beta) const {
    auto mod         = ct->GetModulus();
    auto& LWEParams  = params->GetLWEParams();
    auto& RGSWParams = params->GetRingGSWParams();

    NativeInteger q = LWEParams->Getq();
    if (mod <= q) {
        std::string errMsg =
            "ERROR: EvalDecomp is only for large precision. For small precision, please use bootstrapping directly";
        OPENFHE_THROW(not_implemented_error, errMsg);
    }

    const auto curBase = RGSWParams->GetBaseG();
    auto search        = EKs.find(curBase);
    if (search == EKs.end()) {
        std::string errMsg("ERROR: No key [" + std::to_string(curBase) + "] found in the map");
        OPENFHE_THROW(openfhe_error, errMsg);
    }
    RingGSWBTKey curEK(search->second);

    auto cttmp = std::make_shared<LWECiphertextImpl>(*ct);
    std::vector<LWECiphertext> ret;
    while (mod > q) {
        auto ctq = std::make_shared<LWECiphertextImpl>(*cttmp);
        ctq->SetModulus(q);
        ret.push_back(std::move(ctq));

        // Floor the input sequentially to obtain the most significant bit
        cttmp = EvalFloor(params, curEK, cttmp, beta);
        mod   = mod / q * 2 * beta;
        // round Q to 2betaQ/q
        cttmp = LWEscheme->ModSwitch(mod, cttmp);

        if (EKs.size() == 3) {  // if dynamic
            uint32_t binLog = static_cast<uint32_t>(ceil(log2(mod.ConvertToInt())));
            uint32_t base   = 0;
            if (binLog <= static_cast<uint32_t>(17))
                base = static_cast<uint32_t>(1) << 27;
            else if (binLog <= static_cast<uint32_t>(26))
                base = static_cast<uint32_t>(1) << 18;

            if (0 != base) {  // if base is to change ...
                RGSWParams->Change_BaseG(base);

                auto search = EKs.find(base);
                if (search == EKs.end()) {
                    std::string errMsg("ERROR: No key [" + std::to_string(curBase) + "] found in the map");
                    OPENFHE_THROW(openfhe_error, errMsg);
                }
                curEK = search->second;
            }
        }
    }
    RGSWParams->Change_BaseG(curBase);
    ret.push_back(std::move(cttmp));
    return ret;
}

// private:
NTRUCiphertext BinFHEScheme::BootstrapGateCore(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                               ConstVectorNTRUACCKey& ek, ConstLWECiphertext& ct) const {
    if (ek == nullptr) {
        std::string errMsg =
            "Bootstrapping keys have not been generated. Please call BTKeyGen "
            "before calling bootstrapping.";
        OPENFHE_THROW(config_error, errMsg);
    }

    auto& LWEParams  = params->GetLWEParams();
    auto& NTRUParams = params->GetVectorNTRUParams();
    auto polyParams  = NTRUParams->GetPolyParams();

    // Specifies the range [q1,q2) that will be used for mapping
    NativeInteger p  = ct->GetptModulus();                                     //4 明文模数
    NativeInteger q  = ct->GetModulus();                                       //1024
    NativeInteger Q      = LWEParams->GetQ();             //
    NativeInteger Q2p    = Q / NativeInteger(2 * p) + 1;  //Q/8+1
    NativeInteger Q2pNeg = Q - Q2p;                       //7/8Q-1

    uint32_t N = LWEParams->GetN();
                   
    NativeVector m(N, Q);
    NativeVector new_m(N, Q);
    // Since q | (2*N), we deal with a sparse embedding of Z_Q[x]/(X^{q/2}+1) to
    // Z_Q[x]/(X^N+1)
    uint32_t factor = (2 * N / q.ConvertToInt());
    const NativeInteger& b = ct->GetB()*(2*NativeInteger(N)/q);// b嵌入到Z_2N,    0~2N
    //const NativeInteger& b = ct->GetB();

    for (size_t j = 0; j < N; ++j) {
        m[j] = j<N/2 ?  Q2p:Q2pNeg;         //m: ... + + + + - - - - ... 分成两半 m=△*r(X)
    }
    for (size_t j = 0; j < N; ++j) {        //左移b个单位  new_m=△*r(X) * X^{b}
        auto k = b.ConvertToInt()+j;
        if (k>=N && k<2*N )
        {
            new_m[k%N]=Q- m[j];             //Q-Q2p=Q2pNeg
        }
        else
        {
             new_m[k%N]= m[j];
        }
    }
    NativeInteger azero = ct->GetA()[0];                    //a_0
    uint32_t wzero = factor * azero.ConvertToInt() + 1;     //(2N/q * a_0) +1 == \omega_0
    uint32_t invw = ModInverse(wzero, 2 * N) % (2 * N);     // \omega_0 '
    NativePoly polym(polyParams);
    polym.SetValues(new_m, Format::COEFFICIENT);
    polym.SetFormat(EVALUATION);
    auto polym2{polym.AutomorphismTransform(invw)};         //△*r(X^{w'0}) * X^{-b*w'0}
    auto acc = std::make_shared<NTRUCiphertextImpl>(std::move(polym2));
    NACCscheme->EvalAcc(NTRUParams, ek, acc, ct->GetA());

    //应该是这样的，1是0，7是1
    std::vector<std::vector<NTRUCiphertext>> accV(1,std::vector<NTRUCiphertext>(1));
    accV[0][0]=acc;
    std::vector<NativeVector> fvec(1);
    fvec.push_back(*(params->GetVectorNTRUParams()->m_f));
    showacc(params, accV, 1, fvec);

    return acc;
}


NTRUCiphertext BinFHEScheme::MGBootstrapGateCore(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                               ConstVectorNTRUACCKey& ek, ConstLWECiphertext& ct) const {
    if (ek == nullptr) {
        std::string errMsg =
            "Bootstrapping keys have not been generated. Please call BTKeyGen "
            "before calling bootstrapping.";
        OPENFHE_THROW(config_error, errMsg);
    }

    auto& LWEParams  = params->GetLWEParams();
    auto& NTRUParams = params->GetVectorNTRUParams();
    auto polyParams  = NTRUParams->GetPolyParams();

    // Specifies the range [q1,q2) that will be used for mapping
    NativeInteger p  = ct->GetptModulus();                                     //4 明文模数
    NativeInteger q  = ct->GetModulus();                                       //1024
    NativeInteger Q      = LWEParams->GetQ();             //
    NativeInteger Q2p    = Q / NativeInteger(2 * p) + 1;  //Q/8+1
    NativeInteger Q2pNeg = Q - Q2p;                       //7/8Q-1

    uint32_t N = LWEParams->GetN();
                   
    NativeVector m(N, Q);
    NativeVector new_m(N, Q);
    // Since q | (2*N), we deal with a sparse embedding of Z_Q[x]/(X^{q/2}+1) to
    // Z_Q[x]/(X^N+1)
    uint32_t factor = (2 * N / q.ConvertToInt());
    const NativeInteger& b = ct->GetB()*(2*NativeInteger(N)/q);// b嵌入到Z_2N,    0~2N

    for (size_t j = 0; j < N; ++j) {
        m[j] = j<N/2 ?  Q2p:Q2pNeg;         //m: ... + + + + - - - - ... 分成两半 m=△*r(X)
    }
    for (size_t j = 0; j < N; ++j) {        //左移b个单位  new_m=△*r(X) * X^{b}
        auto k = b.ConvertToInt()+j;
        if (k>=N && k<2*N )
        {
            new_m[k%N]=Q- m[j];             //Q-Q2p=Q2pNeg
        }
        else
        {
             new_m[k%N]= m[j];
        }
    }
    NativeInteger azero = ct->GetA()[0];                    //a_0
    uint32_t wzero = factor * azero.ConvertToInt() + 1;     //(2N/q * a_0) +1 == \omega_0
    uint32_t invw = ModInverse(wzero, 2 * N) % (2 * N);     // \omega_0 '
    NativePoly polym(polyParams);
    polym.SetValues(new_m, Format::COEFFICIENT);
    polym.SetFormat(EVALUATION);
    auto polym2{polym.AutomorphismTransform(invw)};         //△*r(X^{w'0}) * X^{-b*w'0}
    auto acc = std::make_shared<NTRUCiphertextImpl>(std::move(polym2));
    NACCscheme->MGEvalAcc(NTRUParams, ek, acc, ct->GetA());

    std::vector<std::vector<NTRUCiphertext>> accV(1,std::vector<NTRUCiphertext>(1));
    accV[0][0]=acc;
    //std::vector<NativeVector> fvec(1);
    //fvec.push_back(*(params->GetVectorNTRUParams()->m_f));
    //showacc(params, accV, 1, fvec);


    return acc;
}

// 计算f(X^t)的系数
// 输入：f的系数向量a，模数Q，多项式度数N，变换指数t
// 输出：f(X^t)的系数向量
NativeVector compute_f_xt(const std::shared_ptr<BinFHECryptoParams> params, const NativeVector& a, int t) {
    uint32_t N = params->GetLWEParams()->GetN();
    NativeInteger Q = params->GetLWEParams()->GetQ();
    
    NativeVector b(N, Q);
    
    uint32_t t_unsigned = static_cast<uint32_t>(std::abs(t));
    uint32_t d = gcd(t_unsigned, N);
    
    for (uint32_t k_prime = 0; k_prime < N/d; k_prime++) {
        uint32_t k = d * k_prime;
        
        for (uint32_t i = 0; i < d; i++) {
            uint32_t j = (k_prime + i * (N/d)) % N;
            
            uint64_t j_t = static_cast<uint64_t>(j) * static_cast<uint64_t>(t_unsigned);
            uint64_t exponent = j_t / N;
            bool is_positive = (exponent % 2 == 0);
            
            // 使用NativeInteger的模运算方法
            if (is_positive) {
                b[k] = b[k].ModAdd(a[j], Q);
            } else {
                b[k] = b[k].ModSub(a[j], Q);
            }
        }
    }
    
    return b;
}

std::vector<NativePoly> BinFHEScheme::MGBootstrapGateCore_new(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                               ConstLWECiphertext& ct, std::vector<VectorNTRUACCKey>& EK_ALL, 
                                               std::vector<HPK>& HPK_ALL, std::vector<NativePoly> negativeCRP, int numofgroups) const {


    std::string user_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user";
    

    auto& LWEParams  = params->GetLWEParams();
    auto& NTRUParams = params->GetVectorNTRUParams();
    auto polyParams  = NTRUParams->GetPolyParams();
    

    // Specifies the range [q1,q2) that will be used for mapping
    NativeInteger p  = ct->GetptModulus();                                     //4 明文模数
    NativeInteger q  = ct->GetModulus();                                       //1024
    NativeInteger Q      = LWEParams->GetQ();             //
    NativeInteger Q2p    = Q / NativeInteger(2 * p) + 1;  //Q/8+1
    NativeInteger Q2pNeg = Q - Q2p;                       //7/8Q-1

    uint32_t N = LWEParams->GetN();
                   
    NativeVector m(N, Q);
    NativeVector new_m(N, Q);
    // Since q | (2*N), we deal with a sparse embedding of Z_Q[x]/(X^{q/2}+1) to
    // Z_Q[x]/(X^N+1)
    [[maybe_unused]] uint32_t factor = (2 * N / q.ConvertToInt());
    const NativeInteger& b = ct->GetB()*(2*NativeInteger(N)/q);// b嵌入到Z_2N,    0~2N

    for (size_t j = 0; j < N; ++j) {
        m[j] = j<N/2 ?  Q2pNeg:Q2p;         //m: ... + + + + - - - - ... 分成两半 m=△*r(X)
    }
    for (size_t j = 0; j < N; ++j) {        //new_m=△*r(X) * X^{(2N/q)*b}
        auto k = b.ConvertToInt()+j;
        if (k>=N && k<2*N )
        {
            new_m[k%N]=Q- m[j];             //Q-Q2p=Q2pNeg
        }
        else
        {
             new_m[k%N]= m[j];
        }
    }
    

    NativePoly polym(polyParams);
    polym.SetValues(new_m, Format::COEFFICIENT);
    polym.SetFormat(EVALUATION);
    auto acc = std::make_shared<NTRUCiphertextImpl>(std::move(polym));

    //使用numofgroups而不是ct->groupidxs.size()来初始化accV的大小
    std::vector<std::vector<NTRUCiphertext>> accV(numofgroups, std::vector<NTRUCiphertext>(params->GetVectorNTRUParams()->GetDigitsG()));
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < numofgroups; ++i) {
        for(size_t j=0;j<params->GetVectorNTRUParams()->GetDigitsG();++j){
            //获得X^{as}的密文
            accV[i][j]=sub_MGBootstrapGateCore_new(params, NAND, ct, EK_ALL, i, params->GetVectorNTRUParams()->GetGPower()[j]);
        }
    }


    //将acc变成一个MGRLWE密文
    std::vector<NativePoly> accRLWE;
    if (numofgroups > 0 && acc) {
        accRLWE.resize(numofgroups+1);
        // 设置第一个元素为acc的elements
        accRLWE[0] = acc->GetElements();
        // 为其余元素创建值为0的多项式
        for(int i = 1; i < numofgroups+1; ++i){
            // 创建一个全0的多项式，使用第一个元素的参数
            accRLWE[i] = NativePoly(params->GetVectorNTRUParams()->GetPolyParams(), Format::COEFFICIENT, true);
            accRLWE[i].SetFormat(Format::EVALUATION);
        }
    }
    // 确保i不超过accRLWE的大小，避免越界访问
    for(int i = 0; i < numofgroups; ++i){
        MGHybridProd(params, accRLWE, HPK_ALL, negativeCRP, accV[i], i);
    }

    
    //测试一下hybrid product是否正确
    //Step1: 读取各组个party的rlwesk，并将其相加，记为RZ
    // std::vector<NativePoly> RZ;
    // RZ.resize(numofgroups);
    // for(int i = 0; i < numofgroups; ++i){
    //     int group_idx = i; 
    //     std::string group_dir = user_dir + "/group_" + std::to_string(group_idx);
    //     NativePoly Z(polyParams, Format::EVALUATION, true);
        
    //     // 遍历组目录下的所有party目录，读取rlwesk文件并累加到Z 
    //     // 遍历组目录下的所有party目录
    //     for (const auto& entry : fs::directory_iterator(group_dir)) {
    //         if (entry.is_directory()) {
    //             std::string dir_name = entry.path().filename().string();
                
    //             // 检查是否是party目录（格式为"party_数字"）
    //             if (dir_name.rfind("party_", 0) == 0) {
    //                 std::string rlwesk_path = entry.path().string() + "/rlwesk";
    //                 RLWEPrivateKey rsk=std::make_shared<RLWEPrivateKeyImpl>();
    //                 try{
    //                     std::ifstream z_ifs(rlwesk_path);
    //                     cereal::JSONInputArchive archive(z_ifs);
    //                     archive(rsk);
    //                     NativePoly tempZ = rsk->GetElement();
    //                     tempZ.SetFormat(Format::EVALUATION);
    //                     //Z.SetFormat(Format::EVALUATION);
    //                     Z += tempZ;
    //                     //Z.SetFormat(Format::COEFFICIENT);
    //                 }
    //                 catch (const std::exception& e) {
    //                     OPENFHE_THROW(deserialize_error, "Failed to deserialize rlwesk for group " + std::to_string(group_idx) + " party " + dir_name + ": " + e.what());
    //                 }
                    
    //                 // 检查rlwesk文件是否存在
    //                 // if (fs::exists(rlwesk_path)) {
    //                 //     std::ifstream z_ifs(rlwesk_path);
    //                 //     NativeVector z_vec(params->GetVectorNTRUParams()->GetN(), Q);
                        
    //                 //     for (size_t j = 0; j < params->GetVectorNTRUParams()->GetN(); ++j) {
    //                 //         std::string value_str;
    //                 //         z_ifs >> value_str;
    //                 //         z_vec[j] = NativeInteger(value_str);
    //                 //     }
                        
    //                 //     z_ifs.close();
                        
    //                 //     // 累加当前party的rlwesk到Z
    //                 //     NativePoly tempZ(polyParams, Format::COEFFICIENT, true);
    //                 //     tempZ.SetValues(z_vec, Format::COEFFICIENT);
    //                 //     tempZ.SetFormat(Format::EVALUATION);
    //                 //     Z.SetFormat(Format::EVALUATION);
    //                 //     Z += tempZ;
    //                 //     Z.SetFormat(Format::COEFFICIENT);
    //                 //}
    //             }
    //         }
    //     }
    //     Z.SetFormat(Format::COEFFICIENT);
    //     // 将当前组的Z赋值给RZ
    //     RZ[i] = Z;
    // }

    // //Step2: 计算accRLWE[0] + \sum_{i=1}^{numofgroups} accRLWE[i] * RZ[i]
    // NativePoly RLWEB = accRLWE[0];
    // RLWEB.SetFormat(Format::EVALUATION);

    // for(int i = 1; i < numofgroups+1; ++i){
    //     accRLWE[i].SetFormat(Format::EVALUATION);
    //     RZ[i-1].SetFormat(Format::EVALUATION);
    //     RLWEB += (accRLWE[i] * RZ[i-1]);
    // }
    // // cout<<"accRLWE[0]:"<<accRLWE[0].GetValues()<<endl;
    // // cout<<"accRLWE[1]:"<<accRLWE[1].GetValues()<<endl;
    // //accRLWE[2].SetFormat(Format::COEFFICIENT);
    // //cout<<"这里accRLWE[2]:"<<accRLWE[2]<<endl;
    // //accRLWE[2].SetFormat(Format::EVALUATION);

    // //Step3: 打印结果
    // RLWEB.SetFormat(Format::COEFFICIENT);
    // auto rrr = RLWEB.GetValues();
    // //cout<<"before rounding:"<<rrr<<endl;
    // for(uint32_t j = 0; j < params->GetVectorNTRUParams()->GetN(); ++j){
    //     auto tmp = rrr[j].MultiplyAndRound(NativeInteger(8), Q);
    //     rrr[j] = tmp;
    // }
    // //以追加模式向rrr.log中写入下面的内容
    // string rrrpath = "/vscode/myProgram/RRRREALL/rrr.log";
    // std::ofstream rrr_ifs(rrrpath, std::ios::app);
    // rrr_ifs<<"RLWE Ciphertext Decrypt: "<<rrr<<endl;
    // rrr_ifs.close();

    // cout<<"RLWE Ciphertext Decrypt: "<<rrr<<endl;

    return accRLWE;
}

//ct是MGRLWE密文，acc是冗余的自举密文
void BinFHEScheme::MGHybridProd(const std::shared_ptr<BinFHECryptoParams>& params, std::vector<NativePoly>& ct,
                                     std::vector<HPK>& HPK_ALL, std::vector<NativePoly> negativeCRP, std::vector<NTRUCiphertext>& acc, int index) const{

    size_t k = ct.size() - 1;
    //cout<<"k:"<<k<<endl;
    auto vectorNTRUParams = params->GetVectorNTRUParams();
    auto polyParams = vectorNTRUParams->GetPolyParams();
    uint32_t digitsG = vectorNTRUParams->GetDigitsG();
    [[maybe_unused]] auto Q = vectorNTRUParams->GetQ();
    
    // 定义⊡算子操作的辅助函数

    auto multiplyAccumulate = [&](const NativePoly& ciphertext, const std::vector<NativePoly>& accumulator) {
        std::vector<NativePoly> dct(digitsG, NativePoly(polyParams, Format::COEFFICIENT, true));
        
        NativePoly tempCiphertext = ciphertext; // 创建副本以避免修改原始密文
        tempCiphertext.SetFormat(Format::COEFFICIENT);
        VectorNTRUAccumulator va;
        va.CompleteSignedDigitDecompose(vectorNTRUParams, tempCiphertext, dct);
        
        NativePoly sum(polyParams, Format::EVALUATION, true);
        for(uint32_t d = 0; d < digitsG; ++d){
            dct[d].SetFormat(Format::EVALUATION);
        }
        for (uint32_t d = 0; d < digitsG; ++d){
            auto tmp = dct[d] * accumulator[d];
            sum += tmp;
        }
        return sum;
    };//外积计算没问题111

    // // 从rlweCRP文件中读取p的值并取负
    // std::vector<NativePoly> negativeP(digitsG, NativePoly(polyParams, Format::EVALUATION, true));
    // std::string user_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user";
    // std::string crp_path = user_dir + "/rlweCRP";
    // if (std::filesystem::exists(crp_path)) {
    //     std::ifstream ifs(crp_path);
    //     for (uint32_t d = 0; d < digitsG; ++d) {
    //         NativeVector p_vec(vectorNTRUParams->GetN(), vectorNTRUParams->GetQ());
    //         for (size_t j = 0; j < vectorNTRUParams->GetN(); ++j) {
    //             std::string value_str;
    //             ifs >> value_str;
    //             p_vec[j] = NativeInteger(value_str);
    //         }
    //         NativePoly pd(polyParams, Format::COEFFICIENT, true);
    //         pd.SetValues(p_vec, Format::COEFFICIENT);
    //         pd.SetFormat(Format::EVALUATION);
    //         pd = pd.Negate();
    //         negativeP[d] = pd; //读取取负没问题111
    //     }
    //     ifs.close();
    // }
    // else {
    //     // 如果rlweCRP文件不存在
    //     std::cerr << "警告：rlweCRP文件不存在于路径" << crp_path << std::endl;
    // }
    
    
    // 创建y和u的容器
    std::vector<NativePoly> y(k + 1);
    std::vector<NativePoly> u(k + 1);


    // 对于0≤l≤k，计算y_l = ct[l] ⊡ acc；u_l = y_l ⊡ hpk_all[index][1]；v+=∑(y_l ⊡ hpk_all[index][3])
    for(size_t l = 0; l <= k; ++l){
        // 计算y_l = ct[l] ⊡ acc
        // 转换acc为vector<NativePoly>类型
        std::vector<NativePoly> acc_native(digitsG);
        for (uint32_t d = 0; d < digitsG; ++d) {
            acc_native[d]=acc[d]->GetElements();
        }
        y[l] = multiplyAccumulate(ct[l], acc_native);
        
        // 计算u_l = y_l ⊡ hpk[1]
        // 注意：这里hpk[1]是一个向量，我们需要为每个维度单独计算
        std::vector<NativePoly> hpk1_as_acc = HPK_ALL[index]->GetHPK(1);
        u[l] = multiplyAccumulate(y[l], hpk1_as_acc);
    }    

    // v = Σ_{l=0}^k y_l ⊡ H̃_l
    // 注意：l=0时使用H̃_0 = -p，l≥1时使用hpk_all[index][3]中的H̃_l
    NativePoly v = multiplyAccumulate(y[0], negativeCRP);
    // 对于l=1到k，使用存储在hpk_all[index][3]中的H̃_l
    for(size_t l = 1; l <= k; ++l){
        std::vector<NativePoly> H_l = HPK_ALL[l-1]->GetHPK(3); // H部分
        NativePoly yl_Hl = multiplyAccumulate(y[l], H_l);
        v += yl_Hl;
    }
    
    
    
    // 根据l的值处理c_l'
    for(size_t l = 0; l <= k; ++l){
        if(l == 0){
            // c_l'=u_0 + v ⊡ hpk[2]
            std::vector<NativePoly> hpk2_as_acc = HPK_ALL[index]->GetHPK(2);
            NativePoly v_hpk2 = multiplyAccumulate(v, hpk2_as_acc);
            ct[l] = u[l] + v_hpk2;
        } else if (l == (static_cast<size_t>(index+1))) {
            // c_l'=u_l + v ⊡ hpk[0]
            std::vector<NativePoly> hpk0_as_acc = HPK_ALL[index]->GetHPK(0);
            NativePoly v_hpk0 = multiplyAccumulate(v, hpk0_as_acc);
            ct[l] = u[l] + v_hpk0;
        } else {
            // c_l'=u_l
            ct[l] = u[l];
        }
        
        // 确保结果在模Q下
        // ct[l].SetFormat(Format::COEFFICIENT);
        // ct[l] = ct[l].Mod(Q);
        // ct[l].SetFormat(Format::EVALUATION);
    }
}

NTRUCiphertext BinFHEScheme::sub_MGBootstrapGateCore_new(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                     ConstLWECiphertext& ct, std::vector<VectorNTRUACCKey>& EK_ALL, int index, NativeInteger Delta) const{
 


    auto& LWEParams  = params->GetLWEParams();
    auto& NTRUParams = params->GetVectorNTRUParams();
    auto polyParams  = NTRUParams->GetPolyParams();

    auto n = params->GetLWEParams()->Getn();
    //[[maybe_unused]] int k = ct->GetA().GetLength()/n;

    // Specifies the range [q1,q2) that will be used for mapping
    NativeInteger q  = ct->GetModulus();                                       //1024
    NativeInteger Q = LWEParams->GetQ();         
    uint32_t N = LWEParams->GetN();

    
    NativeVector m(N, Q);
    NativeVector new_m(N, Q);
    uint32_t factor = (2 * N / q.ConvertToInt());
    const NativeInteger& b = NativeInteger(0);

    m[0] = Delta;
    for (size_t j = 1; j < N; ++j) {
        m[j] = 0;         
    }
    for (size_t j = 0; j < N; ++j) {        //左移b个单位  new_m=△*r(X) * X^{b}
        auto k = b.ConvertToInt()+j;
        if (k>=N && k<2*N )
        {
            new_m[k%N]=Q- m[j];            
        }
        else
        {
             new_m[k%N]= m[j];
        }
    }

    NativeInteger azero = ct->GetA()[index * n];  //这里的n是参数的n
    uint32_t wzero = factor * azero.ConvertToInt() + 1;
    uint32_t invw = ModInverse(wzero, 2 * N) % (2 * N);
    NativePoly polym(polyParams);
    polym.SetValues(m, Format::COEFFICIENT);
    polym.SetFormat(EVALUATION);

    auto polym2{polym.AutomorphismTransform(invw)};
    auto acc = std::make_shared<NTRUCiphertextImpl>(std::move(polym2)); 
    NACCscheme->MGEvalAcc_new(NTRUParams, EK_ALL[index], acc, ct->GetA(), params->GetLWEParams()->Getn(), index);
    
    return acc;                                                


}



NTRUCiphertext BinFHEScheme::sub_BootstrapGateCoreMD(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                     ConstVectorNTRUACCKey& ek,ConstLWECiphertext& ct,int index, NativeInteger tri) const{
    if (ek == nullptr) {
        std::string errMsg =
            "Bootstrapping keys have not been generated. Please call BTKeyGen "
            "before calling bootstrapping.";
        OPENFHE_THROW(config_error, errMsg);
    }

    auto& LWEParams  = params->GetLWEParams();
    auto& NTRUParams = params->GetVectorNTRUParams();
    auto polyParams  = NTRUParams->GetPolyParams();

    auto n = params->GetLWEParams()->Getn();
    //int k = ct->GetA().GetLength()/n;

    // Specifies the range [q1,q2) that will be used for mapping
    NativeInteger q  = ct->GetModulus();                                       //1024
    NativeInteger Q      = LWEParams->GetQ();             //
    uint32_t N = LWEParams->GetN();

    
    NativeVector m(N, Q);
    uint32_t factor = (2 * N / q.ConvertToInt());
    //设置r(X^{2N/q})=1 且 b=0，外面还要乘一个三角形
    //m[0]=1*(Q / NativeInteger(8) + 1);

    NativeVector fv=*(params->GetVectorNTRUParams()->m_f);
    //m[0]=1*tri;
    //m[0]=1;
    const NativeInteger& mod = (const NativeInteger&) Q;
    //NativeInteger mu = mod.ComputeMu();   
    //cout<<"mod="<<mod<<endl;
    //cout<<"tri="<<tri<<endl;
    for (size_t j = 0; j < N; ++j) {
        //cout<<fv[j]<<"---";
        m[j] = fv[j].ModMulEq(tri,mod);
        //cout<<m[j]<<" ";
        
    }

    NativeInteger azero = ct->GetA()[index*n];  //这里的n是参数的n
    uint32_t wzero = factor * azero.ConvertToInt() + 1;
    uint32_t invw = ModInverse(wzero, 2 * N) % (2 * N);
    NativePoly polym(polyParams);
    polym.SetValues(m, Format::COEFFICIENT);
    polym.SetFormat(EVALUATION);

    auto polym2{polym.AutomorphismTransform(invw)};
    auto acc = std::make_shared<NTRUCiphertextImpl>(std::move(polym2)); 
    NACCscheme->EvalAccMD(NTRUParams, ek, acc, ct->GetA(), n, index);
    return acc;
}


//修改的bootstrap核心
NTRUCiphertext BinFHEScheme::BootstrapGateCoreMD(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                               std::vector<VectorNTRUACCKey>& ekvec, ConstLWECiphertext& ct,vector<NativeVector>& fvec) const {
    if (ekvec.size() == 0) {
        std::string errMsg =
            "Bootstrapping keys have not been generated. Please call BTKeyGen "
            "before calling bootstrapping.";
        OPENFHE_THROW(config_error, errMsg);
    }

    //获取参数
    auto& LWEParams  = params->GetLWEParams();
    auto& NTRUParams = params->GetVectorNTRUParams();
    auto polyParams  = NTRUParams->GetPolyParams();

    NativeInteger p  = ct->GetptModulus();                                     //4 明文模数
    NativeInteger q  = ct->GetModulus();                                       //1024
    NativeInteger Q      = LWEParams->GetQ();             //
    NativeInteger Q2p    = Q / NativeInteger(2 * p) + 1;  //Q/8+1   = 文中的三角形
    NativeInteger Q2pNeg = Q - Q2p;                       //7/8Q-1
    uint32_t N = LWEParams->GetN();



    /*--------------------调用多次sub_bootstrapcore---------------------------*/
    int k=ekvec.size();     //参与方个数
    // cout<<k<<endl;
    
    //设置初始ACC=r(X^{2N/q}) * X^{-(2N/q)*b}        
    NativeVector m(N, Q);
    NativeVector new_m(N, Q);

    uint32_t factor = (2 * N / q.ConvertToInt());
    const NativeInteger& b = ct->GetB()*(2*NativeInteger(N)/q);// b嵌入到Z_2N,    0~2N

    for (size_t j = 0; j < N; ++j) {
       m[j] = j<N/2 ?  Q2p:Q2pNeg;         //m: ... + + + + - - - - ... 分成两半 m=r(X^{2N/q})
    }
    for (size_t j = 0; j < N; ++j) {        // new_m=r(X) * X^b，注意mod的是X^N+1 与 Q
        auto k = b.ConvertToInt()+j;
        if (k>=N && k<2*N )
        {
            new_m[k%N]=Q- m[j];             //Q-Q2p=Q2pNeg
        }
        else
        {
             new_m[k%N]= m[j];
        }
    }

    NativeInteger azero = ct->GetA()[0];                    //a_0
    uint32_t wzero = factor * azero.ConvertToInt() + 1;     //(2N/q * a_0) +1 == \omega_0
    uint32_t invw = ModInverse(wzero, 2 * N) % (2 * N);     // \omega_0 '
    NativePoly polym(polyParams);
    polym.SetValues(new_m, Format::COEFFICIENT);
    polym.SetFormat(EVALUATION);
    auto polym2{polym.AutomorphismTransform(invw)};         //这里同构一下，就得到了论文里面的初始ACC
    auto initacc = std::make_shared<NTRUCiphertextImpl>(std::move(polym2));

    uint32_t digit = params->GetVectorNTRUParams()->GetDigitsG();
    std::vector<std::vector<NTRUCiphertext>> acc(k,std::vector<NTRUCiphertext>(digit));
    acc[0][0]=initacc;
    NACCscheme->EvalAccMD(NTRUParams,ekvec[0],acc[0][0],ct->GetA(),ct->GetA().GetLength()/k, 0);//k=1时通过测试，说明evalaccmd没问题
    //第0个acc返回r(X)·X^{(第一组a)}
    
    for(int i=1;i<k;++i){
        for(uint32_t j=0;j<digit;++j){
            //1，B，B^2，B^3 分别传入
            //cout<<"j="<<j<<"时传入的值："<<params->GetVectorNTRUParams()->GetGPower()[j]<<endl;
            acc[i][j]=sub_BootstrapGateCoreMD(params, NAND, ekvec[i], ct, i, params->GetVectorNTRUParams()->GetGPower()[j]);
        }
    }


    showacc(params,acc,k,fvec);
    



    //acc[i]之间的乘法,目前acc都是evaluation状态
    NativePoly resultPoly(acc[0][0]->GetElements());
    resultPoly.SetFormat(Format::EVALUATION);

    NativePoly ret(params->GetVectorNTRUParams()->GetPolyParams(), Format::COEFFICIENT, true);
    NativePoly toBedecomPoly(acc[0][0]->GetElements());
    for(int i=1;i<k;++i){   //k大于2的时候这里有问题
        // if(i!=1){
        //     toBedecomPoly.SetFormat(Format::EVALUATION);
        //     for(uint32_t a=0; a<N; ++a){
        //         toBedecomPoly[a]=ret.GetValues()[a];
        //     }
        // }
        std::vector<NativePoly> dct(digit,NativePoly(params->GetVectorNTRUParams()->GetPolyParams(), Format::COEFFICIENT, true));  // d-1维N长多项式
        VectorNTRUAccumulator va;


        toBedecomPoly.SetFormat(Format::COEFFICIENT);
        va.CompleteSignedDigitDecompose(params->GetVectorNTRUParams(),toBedecomPoly,dct);

        ret.SetFormat(Format::EVALUATION);
        for (uint32_t d = 0; d < digit; ++d){
            NativePoly accPolydecomposition(acc[i][d]->GetElements());  //1,0 1,1 1,2 1,3 
            dct[d].SetFormat(Format::EVALUATION);
            auto tmp = dct[d]*accPolydecomposition;  
            ret += tmp;
        }

    }
    ret.SetFormat(Format::EVALUATION);
    auto result=std::make_shared<NTRUCiphertextImpl>(std::move(ret));
    return result;
}






RLWECiphertext BinFHEScheme::BootstrapGateCore(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                               ConstRingGSWACCKey& ek, ConstLWECiphertext& ct) const {
    if (ek == nullptr) {
        std::string errMsg =
            "Bootstrapping keys have not been generated. Please call BTKeyGen "
            "before calling bootstrapping.";
        OPENFHE_THROW(config_error, errMsg);
    }

    auto& LWEParams  = params->GetLWEParams();
    auto& RGSWParams = params->GetRingGSWParams();
    auto polyParams  = RGSWParams->GetPolyParams();

    // Specifies the range [q1,q2) that will be used for mapping
    NativeInteger p  = ct->GetptModulus();  //4 明文模数
    NativeInteger q  = ct->GetModulus();
    uint32_t qHalf   = q.ConvertToInt() >> 1;
    NativeInteger q1 = RGSWParams->GetGateConst()[static_cast<size_t>(gate)];  //3/8q
    NativeInteger q2 = q1.ModAddFast(NativeInteger(qHalf), q);                 //7/8

    // depending on whether the value is the range, it will be set
    // to either Q/8 or -Q/8 to match binary arithmetic
    NativeInteger Q      = LWEParams->GetQ();
    NativeInteger Q2p    = Q / NativeInteger(2 * p) + 1;  //Q/8+1
    NativeInteger Q2pNeg = Q - Q2p;                       //7/8Q-1

    uint32_t N = LWEParams->GetN();
    NativeVector m(N, Q);
    // Since q | (2*N), we deal with a sparse embedding of Z_Q[x]/(X^{q/2}+1) to
    // Z_Q[x]/(X^N+1)
    uint32_t factor = (2 * N / q.ConvertToInt());

    const NativeInteger& b = ct->GetB();
    for (size_t j = 0; j < qHalf; ++j) {
        NativeInteger temp = b.ModSub(j, q);
        if (q1 < q2)
            m[j * factor] = ((temp >= q1) && (temp < q2)) ? Q2pNeg : Q2p;
        else
            m[j * factor] = ((temp >= q2) && (temp < q1)) ? Q2p : Q2pNeg;
    }
    //m(x)-m(x^w)
    std::vector<NativePoly> res(2);
    // no need to do NTT as all coefficients of this poly are zero
    res[0] = NativePoly(polyParams, Format::EVALUATION, true);
    res[1] = NativePoly(polyParams, Format::COEFFICIENT, false);
    res[1].SetValues(std::move(m), Format::COEFFICIENT);
    res[1].SetFormat(Format::EVALUATION);

    // main accumulation computation
    // the following loop is the bottleneck of bootstrapping/binary gate
    // evaluation
    //主累加计算
    //下面的循环是引导/二进制门的瓶颈
    //评估
    auto acc = std::make_shared<RLWECiphertextImpl>(std::move(res));
    ACCscheme->EvalAcc(RGSWParams, ek, acc, ct->GetA());
    return acc;
}

// Functions below are for large-precision sign evaluation, 大精度符号评估
// flooring, homomorphic digit decomposition, and arbitrary 向下取整、同态数字分解以及任意函数评估相关
// funciton evaluation, from https://eprint.iacr.org/2021/1337
template <typename Func>
RLWECiphertext BinFHEScheme::BootstrapFuncCore(const std::shared_ptr<BinFHECryptoParams>& params,
                                               ConstRingGSWACCKey& ek, ConstLWECiphertext& ct, const Func f,
                                               const NativeInteger& fmod) const {
    if (ek == nullptr) {
        std::string errMsg =
            "Bootstrapping keys have not been generated. Please call BTKeyGen before calling bootstrapping.";
        OPENFHE_THROW(config_error, errMsg);
    }

    auto& LWEParams  = params->GetLWEParams();
    auto& RGSWParams = params->GetRingGSWParams();
    auto polyParams  = RGSWParams->GetPolyParams();

    NativeInteger Q = LWEParams->GetQ();
    uint32_t N      = LWEParams->GetN();
    NativeVector m(N, Q);
    // For specific function evaluation instead of general bootstrapping
    NativeInteger ctMod    = ct->GetModulus();
    uint32_t factor        = (2 * N / ctMod.ConvertToInt());
    const NativeInteger& b = ct->GetB();
    for (size_t j = 0; j < (ctMod >> 1); ++j) {
        NativeInteger temp = b.ModSub(j, ctMod);
        m[j * factor]      = Q.ConvertToInt() / fmod.ConvertToInt() * f(temp, ctMod, fmod);
    }
    std::vector<NativePoly> res(2);
    // no need to do NTT as all coefficients of this poly are zero
    res[0] = NativePoly(polyParams, Format::EVALUATION, true);
    res[1] = NativePoly(polyParams, Format::COEFFICIENT, false);
    res[1].SetValues(std::move(m), Format::COEFFICIENT);
    res[1].SetFormat(Format::EVALUATION);

    // main accumulation computation
    // the following loop is the bottleneck of bootstrapping/binary gate
    // evaluation
    auto acc = std::make_shared<RLWECiphertextImpl>(std::move(res));
    ACCscheme->EvalAcc(RGSWParams, ek, acc, ct->GetA());
    return acc;
}



// Full evaluation as described in https://eprint.iacr.org/2020/086
template <typename Func>
LWECiphertext BinFHEScheme::BootstrapFunc(const std::shared_ptr<BinFHECryptoParams>& params, const RingGSWBTKey& EK,
                                          ConstLWECiphertext& ct, const Func f, const NativeInteger& fmod) const {
    auto acc = BootstrapFuncCore(params, EK.BSkey, ct, f, fmod);

    std::vector<NativePoly>& accVec = acc->GetElements();
    // the accumulator result is encrypted w.r.t. the transposed secret key
    // we can transpose "a" to get an encryption under the original secret key
    accVec[0] = accVec[0].Transpose();
    accVec[0].SetFormat(Format::COEFFICIENT);
    accVec[1].SetFormat(Format::COEFFICIENT);

    auto ctExt      = std::make_shared<LWECiphertextImpl>(std::move(accVec[0].GetValues()), std::move(accVec[1][0]));
    auto& LWEParams = params->GetLWEParams();
    // Modulus switching to a middle step Q'
    auto ctMS = LWEscheme->ModSwitch(LWEParams->GetqKS(), ctExt);
    // Key switching
    auto ctKS = LWEscheme->KeySwitch(LWEParams, EK.KSkey, ctMS);
    // Modulus switching
    return LWEscheme->ModSwitch(fmod, ctKS);
}









//./combinetest测试用函数
void showacc(const std::shared_ptr<BinFHECryptoParams>& params, std::vector<std::vector<NTRUCiphertext>>& acc,int k,vector<NativeVector>& fvec){

    NativePoly testPo(params->GetVectorNTRUParams()->GetPolyParams(),Format::COEFFICIENT,true);
    testPo[0]=1;
    testPo.SetFormat(Format::EVALUATION);  

    for(int i=0;i<k;++i){
        //检查每个acc[i]
        // for(unsigned long j=0;j<acc[i].size();++j){
        //     std::cout<<"-----------acc["<<i<<"]["<<j<<"]------------------"<<std::endl;
        //     cout<<"-------"<<endl;
        //     NativePoly accPoly(acc[i][j]->GetElements());
        //     accPoly.SetFormat(Format::COEFFICIENT);
        //     for(uint32_t kk=0;kk<params->GetLWEParams()->GetN();++kk){
        //         std::cout<<accPoly[kk]<<"  ";
        //     }
        //     std::cout<<std::endl;
        //     if(i==0){break;}
        // }

        //std::cout<<"-----------showaac 里面的f1 vec------------------"<<std::endl;        
        std::shared_ptr<NativeVector> f1vec=params->GetVectorNTRUParams()->m_f;
        // for(uint32_t j=0;j<params->GetLWEParams()->GetN();++j){
        //    std::cout<<(*f1vec)[j]<<" ";
        // }
        //std::cout<<std::endl;
        // accPoly.SetFormat(Format::EVALUATION);

        NativePoly fpoly(params->GetVectorNTRUParams()->GetPolyParams());
        fpoly.SetValues((*f1vec),Format::COEFFICIENT);    //注意不能直接设置EVALUTION
        cout<<"----"<<endl;
        fpoly.SetFormat(Format::EVALUATION);

        NativePoly accPoly(acc[0][0]->GetElements());
        auto result=accPoly*fpoly;      //c0*f  c1*f
        result.SetFormat(Format::COEFFICIENT);
        // std::cout<<"---------------accPoly["<<i<<"]*fpoly---------------------"<<std::endl;
        // for(uint32_t j=0;j<params->GetLWEParams()->GetN();++j){
        //     std::cout<<result[j]<<" ";
        // }
        // std::cout<<std::endl;

        result.SetFormat(Format::EVALUATION);
        testPo*=result;
        result.SetFormat(Format::COEFFICIENT);

        std::cout<<"---------------accPoly["<<i<<"]*fpoly-rounding---------------------"<<std::endl;
        for(uint32_t j=0;j<params->GetLWEParams()->GetN();++j){
            NativeInteger r=result[j];
            result[j]=r.MultiplyAndRound(NativeInteger(8),params->GetLWEParams()->GetQ()).ConvertToInt();
        }

        for(uint32_t j=0;j<params->GetLWEParams()->GetN();++j){
            std::cout<<result[j]<<" ";
        }
        std::cout<<std::endl;
    }

    // testPo.SetFormat(Format::COEFFICIENT);
    // cout<<"testPo"<<endl;
    // for(uint32_t i=0;i<params->GetLWEParams()->GetN();++i){
    //     cout<<testPo[i]<<" ";
    // }
    // cout<<endl;

}


LWECiphertext ctExtend(LWECiphertext& ori, NativeVector& F, vector<LWEPrivateKey> skvec, NativeInteger mod){
    //生成etk，得到F_i的多密钥加密密文
    int k=skvec.size();    
    std::vector<NativeVector> s;
    for(int i=0;i<k;++i){
        auto tmp = skvec[i]->GetElement();
        s.push_back(tmp);
    }
    
    int n = k*s[0].GetLength();

    vector<NativeVector> etk;
    for(size_t i=0;i<F.GetLength();++i){
        DiscreteUniformGeneratorImpl<NativeVector> dug;
        dug.SetModulus(mod);
        NativeVector a = dug.GenerateVector(n);

        NativeInteger mu = mod.ComputeMu();
        NativeInteger b = F[i];
        for(int j=0;j<k;++j){
            for(int jj=0;jj<n/k;++jj){
                b-=a[(n/k)*j+jj].ModMulFast(s[j][jj], mod, mu);     //b+as=m
            }
        }

        NativeVector tmp(n+1);
        for(int ii=0;ii<n;++ii){
            tmp[ii]=a[ii];
        }
        tmp[n]=b;
        etk.push_back(tmp);
    }

    //ori=(c1,...,cn,c0)
    NativeInteger b = ori->GetB();
    NativeVector a(n/k);
    for(size_t i=0;i<F.GetLength();++i){
        b-=ori->GetA()[i]*etk[i][n];

    }
    
    LWECiphertext result = LWECiphertext();
    return result;
}


}
;  // namespace lbcrypto
