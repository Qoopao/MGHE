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


#include "lwe-pke.h"
#include <bits/stdint-uintn.h>
#include <math.h>
#include <sys/types.h>
#include "binfhe-constants.h"
#include "binfhecontext.h"
#include "lwe-privatekey-fwd.h"
#include "lwe-publickey-fwd.h"
#include "rlwe-ciphertext.h"
#include "utils/serial.h"
#include "vntru-acc.h"

#include "math/binaryuniformgenerator.h"
#include "math/discreteuniformgenerator.h"
#include "math/hal/nativeintbackend.h"
#include "math/ternaryuniformgenerator.h"

#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <ostream>
#include <vector>
#include <filesystem>
#include <fstream>
#include <sstream>
#include "cereal/archives/json.hpp"

using namespace std;

namespace lbcrypto {

// the main rounding operation used in ModSwitch (as described in Section 3 of
// https://eprint.iacr.org/2014/816) The idea is that Round(x) = 0.5 + Floor(x)
NativeInteger LWEEncryptionScheme::RoundqQ(const NativeInteger& v, const NativeInteger& q,
                                           const NativeInteger& Q) const {
    return NativeInteger(static_cast<BasicInteger>(
                             std::floor(0.5 + v.ConvertToDouble() * q.ConvertToDouble() / Q.ConvertToDouble())))
        .Mod(q);
}

//三元分布的私钥
LWEPrivateKey LWEEncryptionScheme::KeyGen(usint size, const NativeInteger& modulus) const {
    TernaryUniformGeneratorImpl<NativeVector> tug;
    return std::make_shared<LWEPrivateKeyImpl>(LWEPrivateKeyImpl(tug.GenerateVector(size, modulus)));
}

void LWEEncryptionScheme::LWECRMGen(usint size, const NativeInteger& modulus, int groupidx) const{
    TernaryUniformGeneratorImpl<NativeVector> tug;

    // 生成CRM矩阵
    std::vector<NativeVector> crm;
    DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(modulus);
    for (usint i = 0; i < size; ++i) {
        crm.push_back(dug.GenerateVector(size, modulus));
    }
    // 遍历group目录
    namespace fs = std::filesystem;
    std::string user_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user";
    for (const auto& entry : fs::directory_iterator(user_dir)) {
        if (entry.is_directory()) {
            fs::path group_path = entry.path();
            fs::path crm_path = group_path / "lweCRM";
            if (!fs::exists(crm_path)) {//如果不存在则写入，存在则不写入
                std::ofstream ofs(crm_path);
                for (const auto& vec : crm) {
                    for (size_t j = 0; j < vec.GetLength(); ++j) {
                        ofs << vec[j];
                        if (j + 1 != vec.GetLength()) ofs << " ";
                    }
                    ofs << '\n';
                }
                ofs.close();
            }
        }
    }
}


void LWEEncryptionScheme::RLWECRPGen(const std::shared_ptr<VectorNTRUCryptoParams>& VNTRUParams) const{
    
    // 生成所有组共享的CRP多项式向量
    auto digitCount = VNTRUParams->GetDigitsG();
    std::vector<NativePoly> crp_all;
    crp_all.reserve(digitCount);
    
    NativePoly crp(VNTRUParams->GetPolyParams(), Format::COEFFICIENT, true);
    DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(VNTRUParams->GetQ());    
    
    // 生成digitCount个随机多项式
    for (size_t i = 0; i < digitCount; ++i) {
        // 正确初始化每个NativePoly实例
        NativePoly poly_i(VNTRUParams->GetPolyParams(), Format::COEFFICIENT, true);
        NativeVector crpvec_i = dug.GenerateVector(VNTRUParams->GetN(), VNTRUParams->GetQ());
        poly_i.SetValues(crpvec_i, Format::COEFFICIENT);
        crp_all.push_back(poly_i);
    }
    // 检查user目录下是否存在rlweCRP文件
    namespace fs = std::filesystem;
    std::string user_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user";
    fs::path crp_path = user_dir + "/rlweCRP";
    
    // 如果rlweCRP文件不存在则写入，存在则不写入
    if (!fs::exists(crp_path)) {
        std::ofstream ofs(crp_path);
        // 写入digitCount个随机多项式
        for (size_t i = 0; i < digitCount; ++i) {
            const NativeVector& crpvec_i = crp_all[i].GetValues();
            for (size_t j = 0; j < crpvec_i.GetLength(); ++j) {
                ofs << crpvec_i[j];
                if (j + 1 != crpvec_i.GetLength()) ofs << " ";
            }
            // 在每个多项式之后添加换行符，除了最后一个
            if (i + 1 != digitCount) ofs << '\n';
        }
        ofs.close();
    }
}


//组版本的三元分布私钥
LWEPrivateKey LWEEncryptionScheme::LWESecretKeyGen(usint size, const NativeInteger& modulus, int groupIdx, int partyIdx) const {
    TernaryUniformGeneratorImpl<NativeVector> tug;
    namespace fs = std::filesystem;
    std::string user_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user";

    // 生成私钥
    NativeVector sk_vec = tug.GenerateVector(size, modulus);
    // 存储私钥到对应目录
    std::string sk_dir = user_dir + "/group_" + std::to_string(groupIdx) + "/party_" + std::to_string(partyIdx);
    fs::create_directories(sk_dir);
    std::ofstream sk_ofs(sk_dir + "/sk");
    for (size_t i = 0; i < sk_vec.GetLength(); ++i) {
        sk_ofs << sk_vec[i];
        if (i + 1 != sk_vec.GetLength()) sk_ofs << " ";
    }
    sk_ofs << '\n';
    sk_ofs.close();
    return std::make_shared<LWEPrivateKeyImpl>(LWEPrivateKeyImpl(sk_vec));
}

//高斯分布的私钥
LWEPrivateKey LWEEncryptionScheme::KeyGenGaussian(usint size, const NativeInteger& modulus) const {
    DiscreteGaussianGeneratorImpl<NativeVector> dgg;
    double STD_DEV = 3.19 ;
    dgg.SetStd(STD_DEV);
    return std::make_shared<LWEPrivateKeyImpl>(LWEPrivateKeyImpl(dgg.GenerateVector(size, modulus)));
}

//组版本的高斯分布私钥
LWEPrivateKey LWEEncryptionScheme::LWESecretKeyGenGaussian(usint size, const NativeInteger& modulus, int groupIdx, int partyIdx) const {
    DiscreteGaussianGeneratorImpl<NativeVector> dgg;
    double STD_DEV = 3.0 ;
    dgg.SetStd(STD_DEV);
    namespace fs = std::filesystem;
    std::string user_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user";

    // 生成私钥
    NativeVector sk_vec = dgg.GenerateVector(size, modulus);
    // 存储私钥到对应目录
    std::string sk_dir = user_dir + "/group_" + std::to_string(groupIdx) + "/party_" + std::to_string(partyIdx);
    fs::create_directories(sk_dir);
    std::ofstream sk_ofs(sk_dir + "/sk");
    for (size_t i = 0; i < sk_vec.GetLength(); ++i) {
        sk_ofs << sk_vec[i];
        if (i + 1 != sk_vec.GetLength()) sk_ofs << " ";
    }
    sk_ofs << '\n';
    sk_ofs.close();
    return std::make_shared<LWEPrivateKeyImpl>(LWEPrivateKeyImpl(sk_vec));
}


LWEPublicKey LWEEncryptionScheme::LWEPublicKeyGen(const std::shared_ptr<LWECryptoParams>& params, ConstLWEPrivateKey& sk, int groupIdx, int partyIdx) const {
    size_t dim = params->Getn();
    NativeInteger modulus = FreshEncModulus;

    // Read M from file, 模数是FreshEncModulus
    std::vector<NativeVector> M;
    std::string crm_path = "/vscode/myProgram/RRRREALL/src/binfhe/user/group_" + std::to_string(groupIdx) + "/lweCRM";
    std::ifstream infile(crm_path);
    if (!infile.is_open()) {
        throw std::runtime_error("Cannot open lweCRM file: " + crm_path);
    }
    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::vector<NativeInteger> row;
        std::string sval;
        while (iss >> sval) {
            row.push_back(NativeInteger(sval));
        }
        NativeVector nv(row.size(), modulus);
        for (size_t i = 0; i < row.size(); ++i) nv[i] = row[i];
        M.push_back(nv);
    }
    infile.close();
    if (M.size() != dim) {
        throw std::runtime_error("lweCRM row count does not match dim");
    }

    //检查M是否正确读取
    // std::cout<<"M--PKGEN--------:"<<std::endl;
    // std::cout<<"PKGEN modulus="<<modulus<<std::endl;
    // std::cout<<"PKGEN dim="<<dim<<std::endl;
    // for (const auto& vec : M) {
    //     for (size_t j = 0; j < vec.GetLength(); ++j) {
    //         std::cout << vec[j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // Generate error vector e
    DiscreteGaussianGeneratorImpl<NativeVector> dgg;
    dgg.SetStd(3.0);
    NativeVector e = dgg.GenerateVector(dim, modulus);
    //std::cout<<"e="<<e<<std::endl;

    // Compute r = -Ms + e，这里要改
    NativeVector r(dim, modulus);
    NativeVector ske = sk->GetElement();
    //将s转换成FreshEncModulus再与LWECRM计算
    ske.SwitchModulus(modulus);
    NativeInteger mu = modulus.ComputeMu();
    for (size_t j = 0; j < dim; ++j) {
        NativeInteger sum = 0;
        for (size_t i = 0; i < dim; ++i) {
            sum.ModAddFastEq(M[j][i].ModMulFast(ske[i], modulus, mu), modulus);
        }
        // r[j] = (e[j] + (modulus - sum)) % modulus
        r[j] = (e[j] + (modulus - sum.Mod(modulus))).Mod(modulus);
    }
    // 保存r到对应的pk文件
    std::string pk_path = "/vscode/myProgram/RRRREALL/src/binfhe/user/group_" + std::to_string(groupIdx) + "/party_" + std::to_string(partyIdx) + "/pk";
    

    LWEPublicKey lwepublickey = std::make_shared<LWEPublicKeyImpl>(std::move(M), std::move(r));
    
    try{
        std::ofstream ofs(pk_path);
        cereal::JSONOutputArchive archive(ofs);
        archive(lwepublickey);
    }catch(const std::exception& e){
        std::cerr << "错误：无法将LWEPublicKey写入文件" << pk_path << "，错误信息：" << e.what() << std::endl;
    }

    return lwepublickey;
}

LWEPublicKey LWEEncryptionScheme::LWEJointPublicKeyGen(const std::shared_ptr<LWECryptoParams>& params, int groupIdx) const {
    namespace fs = std::filesystem;
    using std::string;
    using std::vector;
    
    // 获取模数和维度
    NativeInteger modulus = FreshEncModulus;
    size_t dim = params->Getn();

    // 构建组目录路径
    string user_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user";
    string group_dir = user_dir + "/group_" + std::to_string(groupIdx);

    // 查找组内所有party目录
    vector<string> party_dirs;
    for (const auto& entry : fs::directory_iterator(group_dir)) {
        if (entry.is_directory() && entry.path().filename().string().find("party_") == 0) {
            party_dirs.push_back(entry.path().string());
        }
    }
    if (party_dirs.empty()) {
        throw std::runtime_error("未找到任何party目录: " + group_dir);
    }

    vector<LWEPublicKey> all_pks; 
    for (const auto& party_dir : party_dirs) {
        string pk_path = party_dir + "/pk";
        LWEPublicKey pkj = std::make_shared<LWEPublicKeyImpl>();
        try{
            std::ifstream pk_ofs(pk_path);
            cereal::JSONInputArchive archive(pk_ofs);
            archive(pkj);
            all_pks.push_back(pkj);
        }catch(const std::exception& e){
            std::cerr << "错误：无法将文件内容读入变量" << pk_path << "，错误信息：" << e.what() << std::endl;
        }
    }

    // 按元素累加并取模
    NativeVector joint_v(dim, modulus);
    for (size_t i = 0; i < dim; ++i) {
        for (const auto& pk : all_pks) {
            joint_v[i] += pk->Getv()[i];
        }
        joint_v[i] = joint_v[i].Mod(modulus);
    }

    // 将结果写入jointpk文件
    string jointpk_path = group_dir + "/jointpk";
    std::vector<NativeVector> M = all_pks[0]->GetA();
    LWEPublicKey jpk = std::make_shared<LWEPublicKeyImpl>(std::move(M), std::move(joint_v));
    try{
        std::ofstream jpk_ofs(jointpk_path);
        cereal::JSONOutputArchive archive(jpk_ofs);
        archive(jpk);
    }catch(const std::exception& e){
        std::cerr << "错误：无法将LWEJointPublicKey写入文件" << jointpk_path << "，错误信息：" << e.what() << std::endl;
    }
    
    // 返回LWEPublicKeyImpl
    return jpk;
}


// size is the ring dimension N, modulus is the large Q used in RGSW encryption of bootstrapping.
LWEKeyPair LWEEncryptionScheme::KeyGenPair(const std::shared_ptr<LWECryptoParams>& params) const {
    int size              = params->GetN();
    NativeInteger modulus = params->GetQ();

    // generate secret vector skN of ring dimension N
    LWEPrivateKey skN = KeyGen(size, modulus);
    if (params->GetKeyDist() == GAUSSIAN) {
        skN = KeyGenGaussian(size, modulus);
    }
    else {
        skN = KeyGen(size, modulus);
    }
    // generate public key pkN corresponding to secret key skN
    auto pkN = PubKeyGen(params, skN);

    auto lweKeyPair = LWEKeyPairImpl(pkN, skN);

    // return the public key (A, v), private key sk pair
    return std::make_shared<LWEKeyPairImpl>(lweKeyPair);
}

// size is the ring dimension N, modulus is the large Q used in RGSW encryption of bootstrapping.
LWEPublicKey LWEEncryptionScheme::PubKeyGen(const std::shared_ptr<LWECryptoParams>& params,
                                            ConstLWEPrivateKey& skN) const {
    size_t dim            = params->GetN();
    NativeInteger modulus = params->GetQ();

    DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(modulus);
    std::vector<NativeVector> A(dim);

    // generate random matrix A of dimension N x N
    for (size_t i = 0; i < dim; i++) {
        NativeVector a = dug.GenerateVector(dim);
        A[i]           = std::move(a);
    }

    // generate error vector e
    DiscreteGaussianGeneratorImpl<NativeVector> dgg;
    NativeVector e = dgg.GenerateVector(dim, modulus);

    // compute v = As + e
    NativeVector v   = e;
    NativeVector ske = skN->GetElement();
    NativeInteger mu = modulus.ComputeMu();

    for (size_t j = 0; j < dim; ++j) {
        for (size_t i = 0; i < dim; ++i) {
            v[j].ModAddFastEq(A[j][i].ModMulFast(ske[i], modulus, mu), modulus);
        }
    }
    // public key A, v
    return std::make_shared<LWEPublicKeyImpl>(LWEPublicKeyImpl(std::move(A), std::move(v)));
}


// classical LWE encryption
// a is a randomly uniform vector of dimension n; with integers mod q
// b = a*s + e + m floor(q/4) is an integer mod q
LWECiphertext LWEEncryptionScheme::Encrypt(const std::shared_ptr<LWECryptoParams>& params, ConstLWEPrivateKey& sk,
                                           LWEPlaintext m, LWEPlaintextModulus p, NativeInteger mod) const {
    if (mod % p != 0 && mod.ConvertToInt() & (1 == 0)) {
        std::string errMsg = "ERROR: ciphertext modulus q needs to be divisible by plaintext modulus p.";
        OPENFHE_THROW(not_implemented_error, errMsg);
    }

    NativeVector s = sk->GetElement();      //获取私钥的值
    uint32_t n     = s.GetLength();         //获取私钥的长度=格维度n
    s.SwitchModulus(mod);                   //把s切换到当前的模数，如果s=1 mod q，那么s=s' mod Q
    
    //for(uint32_t i=0;i<n;++i){
    //	std::cout<<s[i]<<std::endl;
    //}
    NativeInteger b = (m % p) * (mod / p) + params->GetDgg().GenerateInteger(mod);//这一句等于 △m + e  Dgg意思是离散高斯变量生成器
    // NativeInteger b = (m % p) * (mod / p);

    DiscreteUniformGeneratorImpl<NativeVector> dug;     //dug离散均匀变量生成器
    dug.SetModulus(mod);                                //设置好范围
    NativeVector a = dug.GenerateVector(n);

    NativeInteger mu = mod.ComputeMu();                 //计算好巴雷特参数用于快速模乘

    for (size_t i = 0; i < n; ++i) {
        b += a[i].ModMulFast(s[i], mod, mu);
    }

    //输出的密文（a,b）=(a,△m + e + as)：noise(m)=-(-b+as)，论文里是（a,b）=(a,as-(△m + e))：noise(m)=-b+as
    auto ct = std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(a), b.Mod(mod)));
    ct->SetptModulus(p);
    return ct;
}




//MGHE.Encrypt
// LWECiphertext LWEEncryptionScheme::MGHE_Encrypt(const std::shared_ptr<LWECryptoParams>& params,
//     LWEPlaintext m, int groupidx, LWEPlaintextModulus p, NativeInteger mod) const {
//     if (mod % p != 0 && mod.ConvertToInt() & (1 == 0)) {
//         std::string errMsg = "ERROR: ciphertext modulus q needs to be divisible by plaintext modulus p.";
//         OPENFHE_THROW(not_implemented_error, errMsg);
//     }
//     // 1. 读取联合公钥jointpk
//     size_t n = params->Getn();
//     NativeInteger modulus = params->Getq();
//     std::string group_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user/group_" + std::to_string(groupidx);
//     std::string jointpk_path = group_dir + "/jointpk";
//     std::ifstream jointpk_ifs(jointpk_path);
//     if (!jointpk_ifs.is_open()) {
//         throw std::runtime_error("无法打开jointpk文件: " + jointpk_path);
//     }
//     std::vector<NativeInteger> jointpk_vec;
//     std::string sval;
//     while (jointpk_ifs >> sval) {
//         jointpk_vec.push_back(NativeInteger(sval));
//     }
//     jointpk_ifs.close();
//     if (jointpk_vec.size() != n) {
//         throw std::runtime_error("jointpk长度与n不符");
//     }
//     //std::cout<<"pk (jointpk) = ";
//     //for(size_t i = 0; i < jointpk_vec.size(); ++i) {
//     //     std::cout<<jointpk_vec[i]<<" ";
//     // }
//     // std::cout<<std::endl;

//     // 2. 读取CRM矩阵M
//     std::string crm_path = group_dir + "/lweCRM";
//     std::ifstream crm_ifs(crm_path);
//     if (!crm_ifs.is_open()) {
//         throw std::runtime_error("无法打开lweCRM文件: " + crm_path);
//     }
//     std::vector<NativeVector> M;
//     std::string line;
//     while (std::getline(crm_ifs, line)) {
//         std::istringstream iss(line);
//         std::vector<NativeInteger> row;
//         while (iss >> sval) {
//             row.push_back(NativeInteger(sval));
//         }
//         NativeVector nv(row.size(), modulus);
//         for (size_t i = 0; i < row.size(); ++i) nv[i] = row[i];
//         M.push_back(nv);
//     }
//     crm_ifs.close();
//     if (M.size() != n) {
//         throw std::runtime_error("lweCRM行数与n不符");
//     }

//     // std::cout<<"M--ENC--------:"<<std::endl;
//     // std::cout<<"ENC modulus="<<modulus<<std::endl;
//     // std::cout<<"ENC n="<<n<<std::endl;
//     // for (const auto& vec : M) {
//     //     for (size_t j = 0; j < vec.GetLength(); ++j) {
//     //         std::cout << vec[j] << " ";
//     //     }
//     //     std::cout << std::endl;
//     // }

//     // 3. 生成随机向量v ∈ Zq^n (从二元分布)
//     BinaryUniformGeneratorImpl<NativeVector> bug;
//     NativeVector v = bug.GenerateVector(n, modulus);
//     //std::cout<<"v="<<v<<std::endl;

//     // 4. 采样噪声e1, e2
//     DiscreteGaussianGeneratorImpl<NativeVector> dgg;
//     double STD_DEV = 0.9;
//     dgg.SetStd(STD_DEV);
//     NativeInteger e1 = dgg.GenerateInteger(modulus); // 标量噪声
//     NativeVector e2 = dgg.GenerateVector(n, modulus); // 向量噪声
//     //std::cout<<"e1="<<e1<<std::endl;
//     //std::cout<<"e2="<<e2<<std::endl;

//     // 5. 计算b = <jointpk, v> + (q/4)*m + e1
//     NativeInteger b = 0;
//     NativeInteger mu = modulus.ComputeMu();
//     for (size_t i = 0; i < n; ++i) {
//         b.ModAddFastEq(jointpk_vec[i].ModMulFast(v[i], modulus, mu), modulus);
//     }
//     //std::cout<<"b after <pk,v> = "<<b<<std::endl;
//     b.ModAddFastEq((m % p) * (mod / p), modulus);
//     //std::cout<<"b after adding (q/p)*m = "<<b<<std::endl;
//     b.ModAddFastEq(e1, modulus);
//     //std::cout<<"b after adding e1 = "<<b<<std::endl;
//     // 6. 计算a = v^T * M + e2
//     NativeVector a(n, modulus);
//     for (size_t j = 0; j < n; ++j) {
//         NativeInteger sum = 0;
//         for (size_t i = 0; i < n; ++i) {
//             sum.ModAddFastEq(v[i].ModMulFast(M[i][j], modulus, mu), modulus);
//         }
//         a[j] = sum.ModAdd(e2[j], modulus);
//     }
//     //std::cout<<"a = "<<a<<std::endl;
//     // 7. 输出密文(a, b)
//     auto ct = std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(a), b.Mod(modulus)));
//     ct->SetptModulus(p);
//     ct->groupidxs.insert(groupidx);   //自举代码中有些内容不按组号严格递增读取，故需要使用set并仔细检查
//     return ct;
// }

//第二版MGEnc
/*
    用一个大模数FreshEncModulus来加密，然后再模切换回上下文中定义的q
*/
LWECiphertext LWEEncryptionScheme::MGHE_Encrypt(const std::shared_ptr<LWECryptoParams>& params,
    LWEPlaintext m, int groupidx, LWEPlaintextModulus p, NativeInteger mod) const {

    NativeInteger freshmod(FreshEncModulus);  // 使用无符号长整型位移创建2^24
    
    if ( (mod % p != 0 && mod.ConvertToInt() & (1 == 0)) || (freshmod % p != 0 && freshmod.ConvertToInt() & (1 == 0))){
        std::string errMsg = "ERROR: ciphertext modulus q1, q2 needs to be divisible by plaintext modulus p.";
        OPENFHE_THROW(not_implemented_error, errMsg);
    }


    // 1. 读取联合公钥jointpk，生成算法需要改
    size_t n = params->Getn();
    std::string group_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user/group_" + std::to_string(groupidx);
    std::string jointpk_path = group_dir + "/jointpk";
    LWEPublicKey jointpk = std::make_shared<LWEPublicKeyImpl>();
    try{
        std::ifstream jointpk_ifs(jointpk_path);
        cereal::JSONInputArchive archive(jointpk_ifs);
        archive(jointpk);
    }catch(const std::exception& e){
        throw std::runtime_error("无法解析jointpk文件: " + jointpk_path + "，错误信息: " + e.what());
    }
    


    // std::cout<<"pk (jointpk) = ";
    // for(size_t i = 0; i < jointpk_vec.size(); ++i) {
    //     std::cout<<jointpk_vec[i]<<" ";
    // }
    // std::cout<<std::endl;

    // 2. 读取CRM矩阵M
    // std::string crm_path = group_dir + "/lweCRM";
    // std::ifstream crm_ifs(crm_path);
    // if (!crm_ifs.is_open()) {
    //     throw std::runtime_error("无法打开lweCRM文件: " + crm_path);
    // }
    // std::vector<NativeVector> M;
    // std::string line;
    // while (std::getline(crm_ifs, line)) {
    //     std::istringstream iss(line);
    //     std::vector<NativeInteger> row;
    //     while (iss >> sval) {
    //         row.push_back(NativeInteger(sval));
    //     }
    //     NativeVector nv(row.size(), freshmod);  
    //     for (size_t i = 0; i < row.size(); ++i) nv[i] = row[i];
    //     M.push_back(nv);
    // }
    // crm_ifs.close();
    // if (M.size() != n) {
    //     throw std::runtime_error("lweCRM行数与n不符");
    // }

    // std::cout<<"M--ENC--------:"<<std::endl;
    // std::cout<<"ENC modulus="<<freshmod<<std::endl;
    // std::cout<<"ENC n="<<n<<std::endl;
    // for (const auto& vec : M) {
    //     for (size_t j = 0; j < vec.GetLength(); ++j) {
    //         std::cout << vec[j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // 3. 生成随机向量v ∈ Zq^n (从二元分布)
    BinaryUniformGeneratorImpl<NativeVector> bug;
    NativeVector v = bug.GenerateVector(n, freshmod);
    //std::cout<<"v="<<v<<std::endl;

    // 4. 采样噪声e1, e2
    DiscreteGaussianGeneratorImpl<NativeVector> dgg;
    double STD_DEV = 3.0;
    dgg.SetStd(STD_DEV);
    NativeInteger e1 = dgg.GenerateInteger(freshmod); // 标量噪声
    NativeVector e2 = dgg.GenerateVector(n, freshmod); // 向量噪声
    //std::cout<<"e1="<<e1<<std::endl;
    //std::cout<<"e2="<<e2<<std::endl;

    // 5. 计算b = <r, v> + (q/4)*m + e1
    NativeInteger b = 0;
    NativeInteger mu = freshmod.ComputeMu();
    for (size_t i = 0; i < n; ++i) {
        b.ModAddFastEq(jointpk->Getv()[i].ModMulFast(v[i], freshmod, mu), freshmod);
    }
    //std::cout<<"b after <pk,v> = "<<b<<std::endl;
    b.ModAddFastEq((m % p) * (freshmod / p), freshmod);
    //std::cout<<"b after adding (q/p)*m = "<<b<<std::endl;
    b.ModAddFastEq(e1, freshmod);
    //std::cout<<"b after adding e1 = "<<b<<std::endl;
    // 6. 计算a = v^T * M + e2
    NativeVector a(n, freshmod);
    for (size_t j = 0; j < n; ++j) {
        NativeInteger sum = 0;
        for (size_t i = 0; i < n; ++i) {
            sum.ModAddFastEq(v[i].ModMulFast(jointpk->GetA()[i][j], freshmod, mu), freshmod);
        }
        a[j] = sum.ModAdd(e2[j], freshmod);
    }
    //std::cout<<"a = "<<a<<std::endl;
    // 7. 输出密文模切换后的密文(a, b)
    auto ct = std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(a), b.Mod(freshmod)));
    auto ctMs = ModSwitch(mod, ct);
    ctMs->SetptModulus(p);
    ctMs->groupidxs.insert(groupidx);   //自举代码中有些内容不按组号严格递增读取，故需要使用set并仔细检查


    return ctMs;
}


LWECiphertext LWEEncryptionScheme::MGHE_ALLEncrypt(const std::shared_ptr<LWECryptoParams>& params,int k,
    LWEPlaintext m, LWEPlaintextModulus p, NativeInteger mod) const {
    if (mod % p != 0 && mod.ConvertToInt() & (1 == 0)) {
        std::string errMsg = "ERROR: ciphertext modulus q needs to be divisible by plaintext modulus p.";
        OPENFHE_THROW(not_implemented_error, errMsg);
    }
    size_t n = params->Getn();
    NativeInteger modulus = params->Getq();
    std::string group_dir = "/vscode/myProgram/RRRREALL/src/binfhe/user/group_0";
    std::string jointpk_path = group_dir + "/jointpk";
    std::ifstream jointpk_ifs(jointpk_path);
    if (!jointpk_ifs.is_open()) {
        throw std::runtime_error("无法打开jointpk文件: " + jointpk_path);
    }
    std::vector<NativeInteger> jointpk_vec;
    std::string sval;
    while (jointpk_ifs >> sval) {
        jointpk_vec.push_back(NativeInteger(sval));
    }
    jointpk_ifs.close();
    if (jointpk_vec.size() != n) {
        throw std::runtime_error("jointpk长度与n不符");
    }
    std::string crm_path = group_dir + "/lweCRM";
    std::ifstream crm_ifs(crm_path);
    if (!crm_ifs.is_open()) {
        throw std::runtime_error("无法打开lweCRM文件: " + crm_path);
    }
    std::vector<NativeVector> M;
    std::string line;
    while (std::getline(crm_ifs, line)) {
        std::istringstream iss(line);
        std::vector<NativeInteger> row;
        while (iss >> sval) {
            row.push_back(NativeInteger(sval));
        }
        NativeVector nv(row.size(), modulus);
        for (size_t i = 0; i < row.size(); ++i) nv[i] = row[i];
        M.push_back(nv);
    }
    crm_ifs.close();
    if (M.size() != n) {
        throw std::runtime_error("lweCRM行数与n不符");
    }
    BinaryUniformGeneratorImpl<NativeVector> bug;
    NativeVector v = bug.GenerateVector(n, modulus);
    DiscreteGaussianGeneratorImpl<NativeVector> dgg;
    double STD_DEV = 3.19;
    dgg.SetStd(STD_DEV);
    NativeInteger e1 = dgg.GenerateInteger(modulus);
    NativeVector e2 = dgg.GenerateVector(n, modulus);
    NativeInteger b = 0;
    NativeInteger mu = modulus.ComputeMu();
    for (size_t i = 0; i < n; ++i) {
        b.ModAddFastEq(jointpk_vec[i].ModMulFast(v[i], modulus, mu), modulus);
    }
    b.ModAddFastEq((m % (2 * p)) * (modulus / (2 * p)), modulus);
    b.ModAddFastEq(e1, modulus);
    NativeVector a_head(n, modulus);
    for (size_t j = 0; j < n; ++j) {
        NativeInteger sum = 0;
        for (size_t i = 0; i < n; ++i) {
            sum.ModAddFastEq(v[i].ModMulFast(M[i][j], modulus, mu), modulus);
        }
        a_head[j] = sum.ModAdd(e2[j], modulus);
    }
    NativeVector a(k * n, modulus);
    for (size_t i = 0; i < n; ++i) {
        a[i] = a_head[i];
    }
    for (size_t i = n; i < k * n; ++i) {
        a[i] = 0;
    }
    auto ct = std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(a), b.Mod(modulus)));
    ct->SetptModulus(p);
    return ct;
}

void LWEEncryptionScheme::MGHE_Decrypt(const std::shared_ptr<LWECryptoParams>& params, std::vector<std::vector<LWEPrivateKey>>& sk, ConstLWECiphertext& ct,
    LWEPlaintext* result, int k, LWEPlaintextModulus p) const {

    //cout<<"进入dec"<<endl;
    const NativeInteger& mod = ct->GetModulus();
    if (mod % (p * 2) != 0 && mod.ConvertToInt() & (1 == 0)) {
        std::string errMsg = "ERROR: ciphertext modulus q needs to be divisible by plaintext modulus p*2.";
        OPENFHE_THROW(not_implemented_error, errMsg);
    }

    NativeVector a   = ct->GetA();
    int ki=sk[0].size();
    // std::cout<<"DECRYPT: a = "<<a<<std::endl;
    // std::cout<<"DECRYPT: k = "<<k<<", ki = "<<ki<<std::endl;

    
    std::vector<std::vector<NativeVector>> s(k,std::vector<NativeVector>(ki));
    for(int i=0;i<k;++i){
        for(int j=0;j<ki;++j){
    	    s[i][j]=sk[i][j]->GetElement();
        }
    }
    //std::cout<<"DECRYPT: s[0][0] = "<<s[0][0]<<std::endl;


    NativeInteger mu = mod.ComputeMu();
    for(int i=0;i<k;++i){
    	for(int j=0;j<ki;++j){
    	    s[i][j].SwitchModulus(mod);
        }
    }

    NativeInteger inner(0);
    for(int i = 0; i < k; ++i){
    	for (int j = 0; j < ki; ++j) {
            for(uint32_t l=0;l<params->Getn();++l){
        	    inner += a[params->Getn()*i+l].ModMulFast(s[i][j][l], mod, mu);
    	    }
        }
    }
    //std::cout<<"DECRYPT: inner before ModEq = "<<inner<<std::endl;

    inner.ModEq(mod);                                       //=inner%mod
    //std::cout<<"DECRYPT: inner after ModEq = "<<inner<<std::endl;
    NativeInteger r = ct->GetB();   
    //std::cout<<"DECRYPT: b (ct->GetB()) = "<<r<<std::endl;
    r.ModAddFastEq(inner, mod);                             // r = B + \sum{ai*si} % q
    //std::cout<<"DECRYPT: r after b + inner = "<<r<<std::endl;
    // 使用标准LWE解密的舍入方法，正确处理负数情况
    //std::cout<<"DECRYPT: final r before rounding = "<<r<<std::endl;
    //std::cout<<"DECRYPT: 2*p = "<<(2*p)<<std::endl;
    
    // // 添加 q/(2p) 进行舍入调整
    // r.ModAddFastEq((mod / (p * 2)), mod);
    // std::cout<<"DECRYPT: r after adding q/(2p) = "<<r<<std::endl;
    
    // // 使用整数除法进行舍入
    // *result = ((NativeInteger(p) * r) / mod).ConvertToInt();
    // std::cout<<"DECRYPT: final result = "<<*result<<std::endl;

    if(r>mod/2){
        r=mod-r;
    }
    //std::cout<<"DECRYPT: final result = "<<*result<<std::endl;
    auto res = (r.MultiplyAndRound(NativeInteger(4),mod)).ConvertToInt();

    *result = res;

}





// classical public key LWE encryption
// a = As' + e' of dimension n; with integers mod q
// b = vs' + e" + m floor(q/4) is an integer mod q
LWECiphertext LWEEncryptionScheme::EncryptN(const std::shared_ptr<LWECryptoParams>& params, ConstLWEPublicKey& pk,
                                            LWEPlaintext m, LWEPlaintextModulus p, NativeInteger mod) const {
    if (mod % p != 0 && mod.ConvertToInt() & (1 == 0)) {
        std::string errMsg = "ERROR: ciphertext modulus q needs to be divisible by plaintext modulus p.";
        OPENFHE_THROW(not_implemented_error, errMsg);
    }
    NativeVector bp             = pk->Getv();
    std::vector<NativeVector> A = pk->GetA();

    uint32_t N = bp.GetLength();
    bp.SwitchModulus(mod);  // todo : this is probably not required

    auto dgg        = params->GetDgg();
    NativeInteger b = (m % p) * (mod / p) + dgg.GenerateInteger(mod);

    TernaryUniformGeneratorImpl<NativeVector> tug;
    NativeVector sp = tug.GenerateVector(N, mod);
    NativeVector ep = dgg.GenerateVector(N, mod);

    // compute a in the ciphertext (a, b)
    NativeVector a   = ep;
    NativeInteger mu = mod.ComputeMu();

    for (size_t j = 0; j < N; ++j) {
        // columnwise a = A_1s1 + ... + A_NsN
        a.ModAddEq(A[j].ModMul(sp[j]));
    }

    // compute b in ciphertext (a,b)
    for (size_t i = 0; i < N; ++i) {
        b.ModAddEq(bp[i].ModMulFast(sp[i], mod, mu), mod);
    }

    auto ct = std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(a, b));
    ct->SetptModulus(p);  //wkx: m_p = p
    return ct;
}

// convert ciphertext with modulus Q and dimension N to ciphertext with modulus q and dimension n
LWECiphertext LWEEncryptionScheme::SwitchCTtoqn(const std::shared_ptr<LWECryptoParams>& params,
                                                ConstLWESwitchingKey& ksk, ConstLWECiphertext& ct) const {
    // Modulus switching to a middle step Q'
    auto ctMS = ModSwitch(params->GetqKS(), ct);
    // Key switching
    auto ctKS = KeySwitch(params, ksk, ctMS);
    // Modulus switching
    return ModSwitch(params->Getq(), ctKS);
}

// classical LWE decryption
// m_result = Round(4/q * (b - a*s))
void LWEEncryptionScheme::Decrypt(const std::shared_ptr<LWECryptoParams>& params, ConstLWEPrivateKey& sk,
                                  ConstLWECiphertext& ct, LWEPlaintext* result, LWEPlaintextModulus p) const {
    // TODO in the future we should add a check to make sure sk parameters match
    // the ct parameters

    // Create local variables to speed up the computations
    const NativeInteger& mod = ct->GetModulus();
    if (mod % (p * 2) != 0 && mod.ConvertToInt() & (1 == 0)) {
        std::string errMsg = "ERROR: ciphertext modulus q needs to be divisible by plaintext modulus p*2.";
        OPENFHE_THROW(not_implemented_error, errMsg);
    }

    NativeVector a   = ct->GetA();
    NativeVector s   = sk->GetElement();
    uint32_t n       = s.GetLength();
    NativeInteger mu = mod.ComputeMu();
    s.SwitchModulus(mod);
    NativeInteger inner(0);
    for (size_t i = 0; i < n; ++i) {
        inner += a[i].ModMulFast(s[i], mod, mu);
    }
    inner.ModEq(mod);

    NativeInteger r = ct->GetB();

    r.ModSubFastEq(inner, mod);
    if(r>mod/2){
        r=mod-r;
    }

    // Alternatively, rounding can be done as
    // *result = (r.MultiplyAndRound(NativeInteger(4),q)).ConvertToInt();
    // But the method below is a more efficient way of doing the rounding
    // the idea is that Round(4/q x) = q/8 + Floor(4/q x)
   // r.ModAddFastEq((mod / (p * 2)), mod);
    cout<<"p="<<p<<endl;
    //*result = ((NativeInteger(p) * r) / mod).ConvertToInt();
    *result=(r.MultiplyAndRound(NativeInteger(4),mod)).ConvertToInt();

#if defined(WITH_NOISE_DEBUG)
    double error =
        (static_cast<double>(p) * (r.ConvertToDouble() - mod.ConvertToDouble() / (p * 2))) / mod.ConvertToDouble() -
        static_cast<double>(*result);
    std::cerr << error * mod.ConvertToDouble() / static_cast<double>(p) << std::endl;
#endif
}

void LWEEncryptionScheme::EvalAddEq(LWECiphertext& ct1, ConstLWECiphertext& ct2) const {
    ct1->GetA().ModAddEq(ct2->GetA());
    ct1->GetB().ModAddFastEq(ct2->GetB(), ct1->GetModulus());
}

void LWEEncryptionScheme::EvalAddConstEq(LWECiphertext& ct, NativeInteger cnst) const {
    ct->GetB().ModAddFastEq(cnst, ct->GetModulus());
}

void LWEEncryptionScheme::EvalSubEq(LWECiphertext& ct1, ConstLWECiphertext& ct2) const {
    ct1->GetA().ModSubEq(ct2->GetA());
    ct1->GetB().ModSubFastEq(ct2->GetB(), ct1->GetModulus());
}

void LWEEncryptionScheme::EvalSubEq2(ConstLWECiphertext& ct1, LWECiphertext& ct2) const {
    ct2->GetA() = ct1->GetA().ModSub(ct2->GetA());
    ct2->GetB() = ct1->GetB().ModSubFast(ct2->GetB(), ct1->GetModulus());
}

void LWEEncryptionScheme::EvalSubConstEq(LWECiphertext& ct, NativeInteger cnst) const {
    ct->GetB().ModSubFastEq(cnst, ct->GetModulus());
}

void LWEEncryptionScheme::EvalMultConstEq(LWECiphertext& ct1, NativeInteger cnst) const {
    ct1->GetA().ModMulEq(cnst);
    ct1->GetB().ModMulFastEq(cnst, ct1->GetModulus());
}

// Modulus switching - directly applies the scale-and-round operation RoundQ
LWECiphertext LWEEncryptionScheme::ModSwitch(NativeInteger q, ConstLWECiphertext& ctQ) const {
    auto n = ctQ->GetLength();
    //cout<<"MS n="<<n<<endl;
    auto Q = ctQ->GetModulus();
    NativeVector a(n, q);
    for (size_t i = 0; i < n; ++i)
        a[i] = RoundqQ(ctQ->GetA()[i], q, Q);
    return std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(a), RoundqQ(ctQ->GetB(), q, Q)));
}


// Switching key as described in Section 3 of https://eprint.iacr.org/2014/816
LWESwitchingKey LWEEncryptionScheme::KeySwitchGen(const std::shared_ptr<LWECryptoParams>& params,
                                                  ConstLWEPrivateKey& sk, ConstLWEPrivateKey& skN) const {
    const size_t n(params->Getn());
    const size_t N(params->GetN());
    NativeInteger qKS(params->GetqKS());
    NativeInteger::Integer value{1};
    NativeInteger::Integer baseKS(params->GetBaseKS());
    const auto digitCount =
        static_cast<size_t>(std::ceil(log(qKS.ConvertToDouble()) / log(static_cast<double>(baseKS))));
    std::vector<NativeInteger> digitsKS;
    digitsKS.reserve(digitCount);
    for (size_t i = 0; i < digitCount; ++i) {
        digitsKS.emplace_back(value);
        value *= baseKS;
    }

    // newSK stores negative values using modulus q
    // we need to switch to modulus Q
    NativeVector sv(sk->GetElement());
    sv.SwitchModulus(qKS);

    NativeVector svN(skN->GetElement());
    svN.SwitchModulus(qKS);
    //    NativeVector oldSK(oldSKlargeQ.GetLength(), qKS);
    //    for (size_t i = 0; i < oldSK.GetLength(); i++) {
    //        if ((oldSKlargeQ[i] == 0) || (oldSKlargeQ[i] == 1)) {
    //            oldSK[i] = oldSKlargeQ[i];
    //        }
    //        else {
    //            oldSK[i] = qKS - 1;
    //        }
    //    }

    DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(qKS);

    NativeInteger mu(qKS.ComputeMu());

    std::vector<std::vector<std::vector<NativeVector>>> resultVecA(N);
    std::vector<std::vector<std::vector<NativeInteger>>> resultVecB(N);

// TODO (cpascoe/dsuponit): this pragma needs to be revised as it may have to be removed completely
// #if !defined(__MINGW32__) && !defined(__MINGW64__)
// #pragma omp parallel for num_threads(N)
// #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(N))
// #endif
    for (size_t i = 0; i < N; ++i) {
        std::vector<std::vector<NativeVector>> vector1A;
        vector1A.reserve(baseKS);
        std::vector<std::vector<NativeInteger>> vector1B;
        vector1B.reserve(baseKS);

        for (size_t j = 0; j < baseKS; ++j) {
            std::vector<NativeVector> vector2A;
            vector2A.reserve(digitCount);
            std::vector<NativeInteger> vector2B;
            vector2B.reserve(digitCount);
            for (size_t k = 0; k < digitCount; ++k) {
                vector2A.emplace_back(dug.GenerateVector(n));
                NativeVector& a = vector2A.back();
                NativeInteger b =
                    (params->GetDggKS().GenerateInteger(qKS)).ModAdd(svN[i].ModMul(j * digitsKS[k], qKS), qKS);
#if NATIVEINT == 32
                for (size_t i = 0; i < n; ++i) {
                    b.ModAddFastEq(a[i].ModMulFast(sv[i], qKS, mu), qKS);
                }
#else
                for (size_t i = 0; i < n; ++i) {
                    b += a[i].ModMulFast(sv[i], qKS, mu);
                }
                b.ModEq(qKS);
#endif
                vector2B.emplace_back(b);
            }
            vector1A.push_back(std::move(vector2A));
            vector1B.push_back(std::move(vector2B));
        }
        resultVecA[i] = std::move(vector1A);
        resultVecB[i] = std::move(vector1B);
    }
    return std::make_shared<LWESwitchingKeyImpl>(LWESwitchingKeyImpl(std::move(resultVecA), std::move(resultVecB)));
}

LWESwitchingKey LWEEncryptionScheme::MGKeySwitchGen(const std::shared_ptr<LWECryptoParams>& params,
                                                  ConstLWEPrivateKey& sk, ConstLWEPrivateKey& skN) const {
    const size_t n(params->Getn());
    const size_t N(params->GetN());
    NativeInteger qKS(params->GetqKS());
    NativeInteger::Integer value{1};
    NativeInteger::Integer baseKS(params->GetBaseKS());
    const auto digitCount =
        static_cast<size_t>(std::ceil(log(qKS.ConvertToDouble()) / log(static_cast<double>(baseKS))));
    std::vector<NativeInteger> digitsKS;
    digitsKS.reserve(digitCount);
    for (size_t i = 0; i < digitCount; ++i) {
        digitsKS.emplace_back(value);
        value *= baseKS;
    }

    // newSK stores negative values using modulus q
    // we need to switch to modulus Q
    NativeVector sv(sk->GetElement());
    //取反
    // auto q = sv.GetModulus();
    // for(uint32_t i=0;i<n;++i){
    //     if(sv[i]<(q/2)){
    //         sv[i]=q-sv[i];
    //     }
    // }

    sv.SwitchModulus(qKS);

    NativeVector svN(skN->GetElement());
    svN.SwitchModulus(qKS);
    //    NativeVector oldSK(oldSKlargeQ.GetLength(), qKS);
    //    for (size_t i = 0; i < oldSK.GetLength(); i++) {
    //        if ((oldSKlargeQ[i] == 0) || (oldSKlargeQ[i] == 1)) {
    //            oldSK[i] = oldSKlargeQ[i];
    //        }
    //        else {
    //            oldSK[i] = qKS - 1;
    //        }
    //    }

    DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(qKS);

    NativeInteger mu(qKS.ComputeMu());

    std::vector<std::vector<std::vector<NativeVector>>> resultVecA(N);
    std::vector<std::vector<std::vector<NativeInteger>>> resultVecB(N);

// TODO (cpascoe/dsuponit): this pragma needs to be revised as it may have to be removed completely
// #if !defined(__MINGW32__) && !defined(__MINGW64__)
// #pragma omp parallel for num_threads(N)
// #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(N))
// #endif
    for (size_t i = 0; i < N; ++i) {
        std::vector<std::vector<NativeVector>> vector1A;
        vector1A.reserve(baseKS);
        std::vector<std::vector<NativeInteger>> vector1B;
        vector1B.reserve(baseKS);

        for (size_t j = 0; j < baseKS; ++j) {
            std::vector<NativeVector> vector2A;
            vector2A.reserve(digitCount);
            std::vector<NativeInteger> vector2B;
            vector2B.reserve(digitCount);
            for (size_t k = 0; k < digitCount; ++k) {
                vector2A.emplace_back(dug.GenerateVector(n));
                NativeVector& a = vector2A.back();
                NativeInteger b =
                    (params->GetDggKS().GenerateInteger(qKS)).ModAdd(svN[i].ModMul(j * digitsKS[k], qKS), qKS);
#if NATIVEINT == 32
                for (size_t i = 0; i < n; ++i) {
                    b.ModAddFastEq(a[i].ModMulFast(sv[i], qKS, mu), qKS);
                }
#else
                for (size_t i = 0; i < n; ++i) {
                    b += a[i].ModMulFast(sv[i], qKS, mu);
                }
                b.ModEq(qKS);
#endif
                vector2B.emplace_back(b);
            }
            vector1A.push_back(std::move(vector2A));
            vector1B.push_back(std::move(vector2B));
        }
        resultVecA[i] = std::move(vector1A);
        resultVecB[i] = std::move(vector1B);
    }
    return std::make_shared<LWESwitchingKeyImpl>(LWESwitchingKeyImpl(std::move(resultVecA), std::move(resultVecB)));

    //以上是原版
}

// the key switching operation as described in Section 3 of
// https://eprint.iacr.org/2014/816
LWECiphertext LWEEncryptionScheme::KeySwitch(const std::shared_ptr<LWECryptoParams>& params, ConstLWESwitchingKey& K,
                                             ConstLWECiphertext& ctQN) const {
    const size_t n(params->Getn());
    const size_t N(params->GetN());
    NativeInteger Q(params->GetqKS());
    NativeInteger::Integer baseKS(params->GetBaseKS());
    const auto digitCount = static_cast<size_t>(std::ceil(log(Q.ConvertToDouble()) / log(static_cast<double>(baseKS))));

    NativeVector a(n, Q);
    NativeInteger b(ctQN->GetB());
    for (size_t i = 0; i < N; ++i) {
        auto& refA = K->GetElementsA()[i];
        auto& refB = K->GetElementsB()[i];
        NativeInteger::Integer atmp(ctQN->GetA(i).ConvertToInt());
        for (size_t j = 0; j < digitCount; ++j) {
            const auto a0 = (atmp % baseKS);
            atmp /= baseKS;
            b.ModSubFastEq(refB[a0][j], Q);
            auto& refAj = refA[a0][j];
            for (size_t k = 0; k < n; ++k)
                a[k].ModSubFastEq(refAj[k], Q);
        }
    }
    return std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(a), b));
}


LWECiphertext LWEEncryptionScheme::MGKeySwitch(const std::shared_ptr<LWECryptoParams>& params, ConstLWESwitchingKey& K,
                                             ConstLWECiphertext& ctQN) const {

                                                
    const size_t n(params->Getn());
    const size_t N(params->GetN());
    NativeInteger Q(params->GetqKS());
    NativeInteger::Integer baseKS(params->GetBaseKS());
    const auto digitCount = static_cast<size_t>(std::ceil(log(Q.ConvertToDouble()) / log(static_cast<double>(baseKS))));

    NativeVector a(n, Q);
    NativeInteger b(ctQN->GetB());
    for (size_t i = 0; i < N; ++i) {
        auto& refA = K->GetElementsA()[i];
        auto& refB = K->GetElementsB()[i];
        NativeInteger::Integer atmp(ctQN->GetA(i).ConvertToInt());
        for (size_t j = 0; j < digitCount; ++j) {
            const auto a0 = (atmp % baseKS);
            atmp /= baseKS;
            b.ModSubFastEq(refB[a0][j], Q);
            auto& refAj = refA[a0][j];
            for (size_t k = 0; k < n; ++k)
                a[k].ModAddFastEq(refAj[k], Q);
        }
    }
    return std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(a), b));
}

LWECiphertext LWEEncryptionScheme::MGKeySwitch_new(const std::shared_ptr<LWECryptoParams>& params,
                                             ConstLWECiphertext& ctQN, 
                                             std::vector<std::vector<std::vector<std::vector<NativeInteger>>>>& A_all,
                                             std::vector<std::vector<std::vector<NativeInteger>>>& B_all, int numofgroups) const {
                                           
    const size_t n(numofgroups * params->Getn());
    const size_t N(params->GetN());
    NativeInteger qKS(params->GetqKS());
    NativeInteger::Integer baseKS(params->GetBaseKS());
    const auto digitCount = static_cast<size_t>(std::ceil(log(qKS.ConvertToDouble()) / log(static_cast<double>(baseKS))));


    
    NativeVector a(n, qKS);   //新密文的a向量,长度为numofgroups * n
    NativeInteger b(ctQN->GetB());  //新密文的b元素
    for(int i = 0; i < numofgroups; ++i){
        NativeVector boldai(N,qKS);
        for(size_t x = 0; x < N; ++x){
            boldai[x] = ctQN->GetA()[i*N+x];
        }
        
        // 初始化output，使其有digitCount个向量，每个向量长度为N
        std::vector<NativeVector> output(digitCount, NativeVector(N, qKS));
        CompleteSignedDigitDecomposeNativeVector(params, boldai, output);

    
        // 计算b
        // 对每个j，计算output(d*N)的每个向量与B_all(k*N*d)的内积
        NativeInteger sum(0);
        auto mu = qKS.ComputeMu();
        for(size_t j = 0; j < N; ++j){
            // 计算output每个向量与B_all[i][j]的内积
            for(uint32_t d = 0; d < digitCount; ++d) {
                // 计算output[d][j] * B_all[i][j][d]
                sum.ModAddEq(output[d][j].ModMulFast(B_all[i][j][d],mu),qKS);
                //cout<<sum<<" ";
            }
        }
        
        // for(auto j=0;j<N;++j){
        //     for(auto d=0;d<digitCount;++d){
        //         cout<<output[d][j]<<" ";
        //     }
        // }
        // cout<<endl;
        // cout<<"qKS="<<qKS<<endl;
        // cout<<"sum="<<sum<<endl;
        // cout<<"BEFORE b="<<b<<endl;
        b+=sum;
        //cout<<"AFTER b="<<b<<endl;

        //计算A  
        NativeVector c_i(params->Getn(),qKS);
        for(uint32_t j=0;j<params->Getn();++j) c_i[j]=0;

        for(uint32_t l=0;l<N;++l){
            for(uint32_t j=0;j<params->Getn();++j){
                for(uint32_t d=0;d<digitCount;++d){
                    c_i[j] += output[d][l].ModMulFast(A_all[i][l][d][j],qKS,mu);
                }
            }
        }
        
        //将c_k添加到a向量中
        for(uint32_t j=0;j<params->Getn();++j){
            a[i*params->Getn()+j]=c_i[j];
        }

    }
    return std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(a), b));
}



// noiseless LWE embedding
// a is a zero vector of dimension n; with integers mod q
// b = m floor(q/4) is an integer mod q
LWECiphertext LWEEncryptionScheme::NoiselessEmbedding(const std::shared_ptr<LWECryptoParams>& params,
                                                      LWEPlaintext m) const {
    NativeInteger q(params->Getq());
    NativeInteger b(m * (q >> 2));
    NativeVector a(params->Getn(), q);
    return std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(a), b));
}

void LWEEncryptionScheme::CompleteSignedDigitDecomposeNativeVector(const std::shared_ptr<LWECryptoParams>& params, const NativeVector& input, std::vector<NativeVector>& output) const {
    // 从input中获取模数Q
    NativeInteger Q = input.GetModulus();
    uint32_t N = input.GetLength();
    uint32_t baseKS = params->GetBaseKS();
    //cout<<baseKS<<endl;
    // 计算Q/2的值 - 一致使用input的模数Q
    auto QHalf{Q.ConvertToInt<BasicInteger>() >> 1};
    auto Q_int{Q.ConvertToInt<NativeInteger::SignedNativeInt>()};
    auto gBits{static_cast<NativeInteger::SignedNativeInt>(__builtin_ctz(baseKS))};
    auto gBitsMaxBits{static_cast<NativeInteger::SignedNativeInt>(NativeInteger::MaxBits() - gBits)};
    
    // 计算分解的位数
    const uint32_t digitCount = static_cast<size_t>(std::ceil(log(Q.ConvertToDouble()) / log(static_cast<double>(baseKS))));
    
        for (uint32_t k{0}; k < N; ++k) {
        auto t0{input[k].ConvertToInt<BasicInteger>()};    //拿到多项式的第k项的值,转换成无符号64位整数
        auto d0{static_cast<NativeInteger::SignedNativeInt>(t0 < QHalf ? t0 : t0 - Q_int)}; //将其转变为mod Q_int的值，高位直接截断


        for (uint32_t d{0}; d < digitCount; ++d) {     //d表示分解的第d个多项式，k表示多项式的第k项
            auto r0{(d0 << gBitsMaxBits) >> gBitsMaxBits};      //d0的高位置0，只留下小于gBitMaxBits的位
            d0 = (d0 - r0) >> gBits;
            if (r0 < 0)
                r0 += Q_int;
            output[d][k] += r0;
        } 
    }
}

};  // namespace lbcrypto
