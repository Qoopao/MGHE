//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2023, NJIT, Duality Technologies Inc. and other contributors
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

#include "vntru-acc-xzddf.h"
#include "binfhe-base-scheme.h"
#include "lattice/hal/lat-backend.h"
#include "math/hal/nativeintbackend.h"
#include "secret-sharing.h"
#include "utils/inttypes.h"

#include <cstddef>
#include <string>
#include <vector>
namespace lbcrypto {

//因为LWE解密形式的不同，ek全是反的
VectorNTRUACCKey VectorNTRUAccumulatorXZDDF::KeyGenAcc(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                                       const NativePoly& skNTT, const NativePoly& invskNTT,
                                                       ConstLWEPrivateKey& LWEsk) const {
    auto sv{LWEsk->GetElement()};

    
    auto mod{sv.GetModulus().ConvertToInt<int32_t>()};  //q_ks
    auto modHalf{mod >> 1};
    size_t n{sv.GetLength()};
    auto q{params->Getq().ConvertToInt<size_t>()};
    
    VectorNTRUACCKey ek = std::make_shared<VectorNTRUACCKeyImpl>(1, 2, q - 1 > n + 1 ? q - 1 : n + 1);
    //生成评估秘钥
    auto s{sv[0].ConvertToInt<int32_t>()};
                                             // 有三种取值q-1, 0, 1 
    (*ek)[0][0][0] = KDMKeyGenXZDDF(params, invskNTT, s > modHalf ? mod - s : -s);
    //(*ek)[0][0][0] = KDMKeyGenXZDDF(params, invskNTT, s > modHalf ? -(mod - s) : s);  //改改改
    //第一个evk(KDM-form)，s的存在形式都是正数，-1用mod-1表示，此时s>modHalf,传入1；如果s是0，传入0；如果s是1，传入-1
    //X^{-1}等于X^{N-1}

//#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(n))
    for (size_t i = 1; i < n; ++i) {
        auto s{sv[i].ConvertToInt<int32_t>()};
        (*ek)[0][0][i] = KeyGenXZDDF(params, invskNTT, s > modHalf ?  mod - s : -s);
        //(*ek)[0][0][i] = KeyGenXZDDF(params, invskNTT, s > modHalf ? -(mod - s) : s);   //改改改
        //如果s大于modHalf，则返回s - mod，否则返回s
    }
    auto sums = 0;
    for (size_t i = 0; i < n; ++i) {
        auto s{sv[i].ConvertToInt<int32_t>()};
        sums = sums +s;
    }
    sums %= mod;
    if (sums > modHalf) {
        sums -= mod;
    }
    (*ek)[0][0][n] = KeyGenXZDDF(params, invskNTT, sums);
    //生成自同构秘钥
    int64_t intq = params->Getq().ConvertToInt<int64_t>();  
    int64_t N    = params->GetN();
    for (auto i = 0; i < intq - 1; ++i) {
        (*ek)[0][1][i] = KeyGenAuto(params, skNTT, invskNTT, (2 * N / intq) * (i + 1) + 1);
    }
    return ek;
}

void VectorNTRUAccumulatorXZDDF::EvalAcc(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                         ConstVectorNTRUACCKey& ek, NTRUCiphertext& acc, const NativeVector& a) const {
    size_t n   = a.GetLength();
    uint32_t N = params->GetN();
    int32_t q  = params->Getq().ConvertToInt<int32_t>();
    std::vector<uint32_t> ua(n);
    std::vector<uint32_t> w(n);
    std::vector<uint32_t> invw(n + 1);
    invw[n] = 1;
    std::vector<NativeInteger> NATIVEw(n);  //自同构的次数
    std::vector<uint32_t> invindex(n);      //对应到autk 的index

    for (size_t i = 0; i < n; i++) {
        ua[i]   = a[i].ConvertToInt<int32_t>();       //a
        w[i]    = (2 * N / q) * ua[i] + 1;            //w_i
        invw[i] = ModInverse(w[i], 2 * N) % (2 * N);  //w_inv
    }
    for (size_t i = 0; i < n; i++) {
        NATIVEw[i] = NativeVector::Integer((w[i] * invw[i + 1]) % (2 * N));
        invindex[i] = (NATIVEw[i].ConvertToInt<int32_t>() - 2*N/q -1) / ( 2*N/q);
    }
    for (size_t i = 0; i < n; i++) {
        AddToAccXZDDF(params, (*ek)[0][0][i], acc);  ///evk_{0 ~ n-1}
        if (NATIVEw[i].ConvertToInt<int32_t>() != 1) {
            Automorphism(params, NATIVEw[i], (*ek)[0][1][invindex[i]], acc);
        }
    }
   
    AddToAccXZDDF(params, (*ek)[0][0][n], acc);
}

void VectorNTRUAccumulatorXZDDF::EvalAccMD(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                         ConstVectorNTRUACCKey& ek, NTRUCiphertext& acc, const NativeVector& a, size_t sizeofaperP, int index) const {
    size_t n   =  sizeofaperP; 
    uint32_t N = params->GetN();
    int32_t q  = params->Getq().ConvertToInt<int32_t>();
    std::vector<uint32_t> ua(n);
    std::vector<uint32_t> w(n);
    std::vector<uint32_t> invw(n + 1);
    invw[n] = 1;
    std::vector<NativeInteger> NATIVEw(n);  //自同构的次数
    std::vector<uint32_t> invindex(n);      //对应到autk 的index
    
    //std::cout<<"ua========="<<std::endl;
    for (size_t i = 0; i < n; i++) {
        ua[i]   = a[index*n+i].ConvertToInt<int32_t>();       //a
        //std::cout<<ua[i]<<std::endl;
        w[i]    = (2 * N / q) * ua[i] + 1;            //w_i
        invw[i] = ModInverse(w[i], 2 * N) % (2 * N);  //w_inv
    }
    for (size_t i = 0; i < n; i++) {
        NATIVEw[i] = NativeVector::Integer((w[i] * invw[i + 1]) % (2 * N));     //第9行的变量
        invindex[i] = (NATIVEw[i].ConvertToInt<int32_t>() - 2*N/q -1) / ( 2*N/q);   //ksk的下标
    }

    if(acc==nullptr){
        std::cout<<"acc of index:"<< index <<" is empty"<<std::endl;
    }
    for (size_t i = 0; i < n; i++) {
        AddToAccXZDDF(params, (*ek)[0][0][i], acc);  ///evk_{0 ~ n-1}
        if (NATIVEw[i].ConvertToInt<int32_t>() != 1) {
            Automorphism(params, NATIVEw[i], (*ek)[0][1][invindex[i]], acc);
        }
    }
    AddToAccXZDDF(params, (*ek)[0][0][n], acc);
}

// KDM-form
VectorNTRUEvalKey VectorNTRUAccumulatorXZDDF::KDMKeyGenXZDDF(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                                             const NativePoly& invskNTT, LWEPlaintext m) const {                          
    auto polyParams = params->GetPolyParams();  //(Q,2N)
    auto Gpow       = params->GetGPower();      //
    DiscreteUniformGeneratorImpl<NativeVector> dug;
    NativeInteger Q{params->GetQ()};
    dug.SetModulus(Q);  //确保dug的模数是Q
                        //Reduce mod q (dealing with negative number as well)
    int64_t N  = params->GetN();
    int64_t mm = (((m % N) + N) % N);  // m=-1, mm=N-1;   m=0,mm=0;   m=1,mm=1
    bool isReducedMM{false};
    if (m < 0) {
        isReducedMM = true;
    }
    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG2{(params->GetDigitsG() - 1)};
    std::vector<NativePoly> tempA(digitsG2, NativePoly(dug, polyParams, Format::COEFFICIENT));
    VectorNTRUEvalKeyImpl result(digitsG2);
    for (uint32_t i = 0; i < digitsG2; ++i) {
        result[i] = NativePoly(params->GetDgg(), polyParams, Format::COEFFICIENT);  //采样g
        if (!isReducedMM)//s<0 m>=0 mm=1,0
            result[i][mm].ModAddFastEq(Gpow[i+1],Q);  // g+X^m*G  g是多项式，在第mm位加上B^{i+1}等于g+B^{i+1}*X^{mm}, 第一个ek与其他的不一样，所以第一个ek的生成用KDMKeyGen 
        else    //s>=0 m<0 mm=N-1
            result[i][mm].ModSubFastEq(Gpow[i+1],Q);  // g-X^m*G   result[i][N-m]=(g-B^{i+1})/f
        result[i].SetFormat(Format::EVALUATION);
        result[i] = result[i] * invskNTT;
    }
    return std::make_shared<VectorNTRUEvalKeyImpl>(result);
}
//NO KDM-form
VectorNTRUEvalKey VectorNTRUAccumulatorXZDDF::KeyGenXZDDF(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                                          const NativePoly& invskNTT, LWEPlaintext m) const {
    auto polyParams = params->GetPolyParams();  //(Q,2N)
    auto Gpow       = params->GetGPower();      //
    NativeInteger Q{params->GetQ()};
    int64_t N  = params->GetN();
    int64_t mm = (((m % N) + N) % N);  // s=0 m=0 mm=0; s=q-1 m=1 mm=1; s=1 m=-1 mm=N-1
    //std::cout<<m<<" "<<mm<<std::endl;
    bool isReducedMM{false};
    if (m < 0) {
        isReducedMM = true;
    }
    uint32_t digitsG2{(params->GetDigitsG() - 1)};  //logQ < (B)^k的最小的k
    NativePoly zeroPoly(polyParams, Format::COEFFICIENT);
    zeroPoly.SetValuesToZero();
    std::vector<NativePoly> tempA(digitsG2, zeroPoly);

    VectorNTRUEvalKeyImpl result(digitsG2);
    for (uint32_t i = 0; i < digitsG2; ++i) {
        // result[i][0] = tempA[i];
        tempA[i].SetFormat(Format::COEFFICIENT);
        result[i] = NativePoly(params->GetDgg(), polyParams, Format::COEFFICIENT);  //采样g
        result[i].SetFormat(Format::EVALUATION);
        result[i] = result[i] * invskNTT;  // g/f
        if (!isReducedMM)   //m>=0 mm==m=0,1
            tempA[i][mm].ModAddFastEq(Gpow[i+1], Q);  // X^m*G
        else{                //m<0 mm=N-1
            tempA[i][mm].ModSubFastEq(Gpow[i+1], Q);  // X^m*G      这里要modsub是因为 X^{-1} = - X^{N-1}  mod X^{N} + 1
        } 
        tempA[i].SetFormat(Format::EVALUATION);
        result[i] = result[i] + tempA[i];
    }
    return std::make_shared<VectorNTRUEvalKeyImpl>(result);
}
VectorNTRUEvalKey VectorNTRUAccumulatorXZDDF::KeyGenAuto(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                                         const NativePoly& skNTT, const NativePoly& invskNTT,
                                                         LWEPlaintext k) const {
    //auto polyParams{params->GetPolyParams()};
    // m_polyParams{std::make_shared<ILNativeParams>(2 * N, Q)},
    // auto Gpow{params->GetGPower()};//m_Gpower,是一个3长度vector (0,1024,1048576)，误导！！！！
    auto polyParams = params->GetPolyParams();  //(Q,2N)
    auto Gpow       = params->GetGPower();      //

    //打印Gpow，这个才是对的，如果gadgetbase=1<<16，那么Gpow=1, 65535=2^16, 4294967296=2^32...
    // for(uint32_t x=0;x<params->GetDigitsG();++x){
    //     std::cout<<Gpow[x]<<" ";
    // }
    //std::cout<<std::endl;
    DiscreteUniformGeneratorImpl<NativeVector> dug;
    NativeInteger Q{params->GetQ()};
    dug.SetModulus(Q);
    auto skAuto{skNTT.AutomorphismTransform(k)};  //生成f(X^k)

    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG{params->GetDigitsG() - 1};
    VectorNTRUEvalKeyImpl result(digitsG);

    for (uint32_t i = 0; i < digitsG; ++i) {
        result[i] = NativePoly(params->GetDgg(), polyParams, EVALUATION) + skAuto * Gpow[i+1];//这里本来是i+1
        result[i] = result[i] * invskNTT;
    }
    return std::make_shared<VectorNTRUEvalKeyImpl>(result);
}

void VectorNTRUAccumulatorXZDDF::AddToAccXZDDF(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                               ConstVectorNTRUEvalKey& ek, NTRUCiphertext& acc) const {
    NativePoly ct(acc->GetElements());
    ct.SetFormat(Format::COEFFICIENT);
    
    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG{(params->GetDigitsG() - 1)};
    std::vector<NativePoly> dct(digitsG,NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));  // d-1维N长多项式

    SignedDigitDecompose(params, ct, dct);   //分解acc
    //如果acc只有第一项！=0，那么打印出dct
    // if(ct[0]!=0&&ct[1]==0){
    //     std::cout<<"Q="<<params->GetQ()<<std::endl;
    //     for(uint32_t ss=0;ss<dct.size();++ss){
    //         std::cout<<"第"<<ss<<"个多项式"<<std::endl;
    //         for(uint32_t sss=0;sss<params->GetN();++sss){
    //             std::cout<<dct[ss][sss]<<" ";
    //         }
    //         std::cout<<std::endl;
    //     }
        
    // }
    
                                                         
    // calls digitsG2 NTTs
    NativePoly sum(params->GetPolyParams(), Format::EVALUATION, true);
//#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(digitsG))
    for (uint32_t d = 0; d < digitsG; ++d)
        dct[d].SetFormat(Format::EVALUATION);
    // acc = dct * ek (matrix product);
    const std::vector<NativePoly>& ev = ek->GetElements();
    //std::cout<<"tmp__________________"<<std::endl;
    for (uint32_t d = 0; d < digitsG; ++d){
        auto tmp = dct[d]*ev[d];  //lyh注：两个NTT多项式相乘就是对应项执行模乘法
        // if(ct[0]!=0 && ct[1]==0){
        //     NativePoly eee(params->GetPolyParams(),Format::EVALUATION,true);
        //     eee=tmp;
        //     eee.SetFormat(Format::COEFFICIENT);
        //     std::cout<<"tmp的值"<<std::endl;
        //     for(uint32_t j=0;j<2048;++j){
        //         cout<<eee[j]<<" ";
        //     }
        //     std::cout<<std::endl;


        //     cout<<"tmp解密一下"<<endl;
        //     eee.SetFormat(Format::EVALUATION);
    
        //     std::shared_ptr<NativeVector> fvec=params->m_f;
        //     NativePoly fpoly(params->GetPolyParams());
        //     fpoly.SetValues(*fvec,Format::COEFFICIENT);    //注意不能直接设置EVALUTION
        //     fpoly.SetFormat(Format::EVALUATION);
    
        //     auto rr=fpoly*eee;

        //     rr.SetFormat(Format::COEFFICIENT);
        //     for(uint32_t j=0;j<2048;++j){
        //         NativeInteger x;
        //         x=rr[j];
        //         x=x.MultiplyAndRound(NativeInteger(8),9007199254614017).ConvertToInt();
        //         cout<<x<<" ";
        //     }
        //     std::cout<<std::endl;


        // }
        // for(uint32_t j=0;j<2048;++j){
        //     std::cout<<tmp[j]<<" ";
        // }
        // std::cout<<std::endl;
        // {
        //     if(ct[0]!=0&&ct[1]==0&&d==2){
        //         std::cout<<"dct[2]的值"<<std::endl;
        //         for(uint32_t i=0;i<params->GetN();++i){
        //             std::cout<<dct[2][i]<<" ";
        //         }
        //         std::cout<<std::endl;

        //         std::cout<<"ek[2]的值"<<std::endl;
        //         for(uint32_t i=0;i<digitsG;++i){
        //             std::cout<<"ek["<<i<<"]"<<std::endl;
        //             for(uint32_t j=0;j<2048;++j){
        //                 if(i==2){
        //                     cout<<ev[i][j]<<" ";                    
        //                 }
        //             }
        //             std::cout<<std::endl;
        //         }

        //         std::cout<<"相乘的值"<<endl;
        //         for(uint32_t j=0;j<params->GetN();++j){
        //             std::cout<<tmp[j]<<" ";
        //         }
        //         std::cout<<std::endl;
        //     }
        // }

        sum += tmp;
    }

    // if(ct[0]!=0 && ct[1]==0){       //解密了一下evk的值 确实没错
    //     NativePoly evk0(params->GetPolyParams(),Format::EVALUATION,true);
    //     evk0=ev[2];

    //     std::shared_ptr<NativeVector> fvec=params->m_f;
    //     NativePoly fpoly(params->GetPolyParams());
    //     fpoly.SetValues(*fvec,Format::COEFFICIENT);    //注意不能直接设置EVALUTION
    //     fpoly.SetFormat(Format::EVALUATION);

    //     auto rr=fpoly*evk0;
        

    //     rr.SetFormat(Format::COEFFICIENT);
    //     std::cout<<"fdksajfldkjsafldjsf"<<std::endl;
    //     for(uint32_t j=0;j<2048;++j){
    //         cout<<rr[j]<<" ";
    //     }
    //     std::cout<<std::endl;

    // }

    // if(ct[0]!=0 && ct[1]==0){
    //     //打出acc的值
    //     ct.SetFormat(Format::COEFFICIENT);
    //     std::cout<<"acc的值"<<std::endl;
    //     for(uint32_t j=0;j<2048;++j){
    //         std::cout<<ct[j]<<" ";
    //     }
    //     std::cout<<std::endl;        


    //     // std::cout<<"acc分解后的值"<<std::endl;
    //     // std::cout<<"Q="<<params->GetQ()<<std::endl;
    //     // for(uint32_t ss=0;ss<dct.size();++ss){
    //     //     std::cout<<"第"<<ss<<"个多项式"<<std::endl;
    //     //     dct[ss].SetFormat(Format::COEFFICIENT);
    //     //     for(uint32_t sss=0;sss<params->GetN();++sss){
    //     //         std::cout<<dct[ss][sss]<<" ";
    //     //     }
    //     //     std::cout<<std::endl;
    //     // }

    //     // std::cout<<"ek的值"<<std::endl;
    //     // for(uint32_t i=0;i<digitsG;++i){
    //     //     std::cout<<"ek["<<i<<"]"<<std::endl;
    //     //     for(uint32_t j=0;j<2048;++j){
    //     //         if(i==2){
    //     //             cout<<ek->GetElements()[i][j]<<" ";                    
    //     //         }

    //     //     }
    //     //     std::cout<<std::endl;
    //     // }


    //     NativePoly sumPoly(sum);
    //     //cout<<sumPoly.GetFormat()<<endl;   

    //     std::shared_ptr<NativeVector> fvec=params->m_f;
    //     NativePoly fpoly(params->GetPolyParams());
    //     fpoly.SetValues(*fvec,Format::COEFFICIENT);    //注意不能直接设置EVALUTION

    //     // std::cout<<"------- fpoly --------"<<std::endl;
    //     // for(uint32_t j=0;j<2048;++j){
    //     //     std::cout<<fpoly[j]<<" ";
    //     // }
    //     // std::cout<<std::endl;

    //     fpoly.SetFormat(Format::EVALUATION);

    //     auto result=sumPoly*fpoly;
    //     result.SetFormat(Format::COEFFICIENT);

    //     // std::cout<<"--------"<<std::endl;
    //     // for(uint32_t j=0;j<2048;++j){
    //     //     std::cout<<result[j]<<" ";
    //     // }
    //     // std::cout<<std::endl;

    //     std::cout<<"acc*evk_0"<<std::endl;
    //     for(uint32_t j=0;j<2048;++j){
    //         //std::cout<<result[j]<<" ";
    //         NativeInteger r=result[j];
    //         result[j]=r.MultiplyAndRound(NativeInteger(8),9007199254614017).ConvertToInt();
    //         std::cout<<result[j]<<" ";
    //     }
    //     std::cout<<std::endl;
    // }
    
    sum.SetFormat(Format::EVALUATION);
    acc->GetElements() = sum;
}

void VectorNTRUAccumulatorXZDDF::Automorphism(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                              const NativeInteger& a, ConstVectorNTRUEvalKey& ak,
                                              NTRUCiphertext& acc) const {
    // precompute bit reversal for the automorphism into vec
    uint32_t N{params->GetN()};
    std::vector<usint> vec(N);
    PrecomputeAutoMap(N, a.ConvertToInt<usint>(), &vec);  //
    NativePoly ct(acc->GetElements());
    acc->GetElements().SetValuesToZero();
    ct = ct.AutomorphismTransform(a.ConvertToInt<usint>(), vec);
    ct.SetFormat(COEFFICIENT);
    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG{params->GetDigitsG() - 1};
    std::vector<NativePoly> dct(digitsG, NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));
    SignedDigitDecompose(params, ct, dct);
    NativePoly sum(params->GetPolyParams(), Format::EVALUATION, true);
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(digitsG))
    for (uint32_t d = 0; d < digitsG; ++d)
        dct[d].SetFormat(Format::EVALUATION);
    // acc = dct * input (matrix product);
    const std::vector<NativePoly>& ev = ak->GetElements();
    for (uint32_t d = 0; d < digitsG; ++d)
        sum += (dct[d] * ev[d]);

    acc->GetElements() = sum;
}

//---------------------------MG--------------------------------

VectorNTRUACCKey VectorNTRUAccumulatorXZDDF::MGKeyGenAcc(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                                       const NativePoly& skNTT, const NativePoly& invskNTT,
                                                       ConstLWEPrivateKey& LWEsk) const {
    auto sv{LWEsk->GetElement()};

    
    auto mod{sv.GetModulus().ConvertToInt<int32_t>()};  //q_ks
    auto modHalf{mod >> 1};
    size_t n{sv.GetLength()};
    auto q{params->Getq().ConvertToInt<size_t>()};
    
    VectorNTRUACCKey ek = std::make_shared<VectorNTRUACCKeyImpl>(1, 2, q - 1 > n + 1 ? q - 1 : n + 1);
    //生成评估秘钥
    auto s{sv[0].ConvertToInt<int32_t>()};
                                             // 有三种取值q-1, 0, 1 
    //因为OpenFHE输出的密文（a,b）=(a,△m + e + as)：noise(m)=-(-b+as)=b-as，论文里是（a,b）=(a,as-(△m + e))：noise(m)=-b+as
    //我们需要改成b+as
    (*ek)[0][0][0] = MGKDMKeyGenXZDDF(params, invskNTT, s > modHalf ? s-mod : s);
    //(*ek)[0][0][0] = KDMKeyGenXZDDF(params, invskNTT, s > modHalf ? -(mod - s) : s);  //改改改
    //第一个evk(KDM-form)，s的存在形式都是正数，-1用mod-1表示，此时s>modHalf,传入1；如果s是0，传入0；如果s是1，传入-1
    //X^{-1}等于X^{N-1}

//#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(n))
    for (size_t i = 1; i < n; ++i) {
        auto s{sv[i].ConvertToInt<int32_t>()};
        (*ek)[0][0][i] = MGKeyGenXZDDF(params, invskNTT, s > modHalf ?  s-mod : s);
        //(*ek)[0][0][i] = KeyGenXZDDF(params, invskNTT, s > modHalf ? -(mod - s) : s);   //改改改
        //如果s大于modHalf，则返回s - mod，否则返回s
    }
    auto sums = 0;
    for (size_t i = 0; i < n; ++i) {
        auto s{sv[i].ConvertToInt<int32_t>()};
        sums = sums +s;
    }
    sums %= mod;
    if (sums > modHalf) {
        sums -= mod;
    }
    (*ek)[0][0][n] = MGKeyGenXZDDF(params, invskNTT, sums);
    //生成自同构秘钥
    int64_t intq = params->Getq().ConvertToInt<int64_t>();  
    int64_t N    = params->GetN();
    for (auto i = 0; i < intq - 1; ++i) {
        (*ek)[0][1][i] = MGKeyGenAuto(params, skNTT, invskNTT, (2 * N / intq) * (i + 1) + 1);
    }
    return ek;
}

void VectorNTRUAccumulatorXZDDF::MGEvalAcc(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                         ConstVectorNTRUACCKey& ek, NTRUCiphertext& acc, const NativeVector& a) const {
    size_t n   = a.GetLength();
    uint32_t N = params->GetN();
    int32_t q  = params->Getq().ConvertToInt<int32_t>();
    std::vector<uint32_t> ua(n);
    std::vector<uint32_t> w(n);
    std::vector<uint32_t> invw(n + 1);
    invw[n] = 1;
    std::vector<NativeInteger> NATIVEw(n);  //自同构的次数
    std::vector<uint32_t> invindex(n);      //对应到autk 的index

    for (size_t i = 0; i < n; i++) {
        ua[i]   = a[i].ConvertToInt<int32_t>();       //a
        w[i]    = (2 * N / q) * ua[i] + 1;            //w_i
        invw[i] = ModInverse(w[i], 2 * N) % (2 * N);  //w_inv
    }
    for (size_t i = 0; i < n; i++) {
        NATIVEw[i] = NativeVector::Integer((w[i] * invw[i + 1]) % (2 * N));
        invindex[i] = (NATIVEw[i].ConvertToInt<int32_t>() - 2*N/q -1) / ( 2*N/q);
    }
    
    for (size_t i = 0; i < n; i++) {
        MGAddToAccXZDDF(params, (*ek)[0][0][i], acc);  ///evk_{0 ~ n-1}
        if (NATIVEw[i].ConvertToInt<int32_t>() != 1) {
            MGAutomorphism(params, NATIVEw[i], (*ek)[0][1][invindex[i]], acc);
        }
    }
   
    MGAddToAccXZDDF(params, (*ek)[0][0][n], acc);
}

void VectorNTRUAccumulatorXZDDF::MGEvalAcc_new(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                         ConstVectorNTRUACCKey& ek, NTRUCiphertext& acc, const NativeVector& a, const NativeInteger& lwedimension, int index) const {
    //size_t n   = a.GetLength();
    size_t n = lwedimension.ConvertToInt<size_t>();
    uint32_t N = params->GetN();
    int32_t q  = params->Getq().ConvertToInt<int32_t>();
    std::vector<uint32_t> ua(n);
    std::vector<uint32_t> w(n);
    std::vector<uint32_t> invw(n + 1);
    invw[n] = 1;
    std::vector<NativeInteger> NATIVEw(n);  //自同构的次数
    std::vector<uint32_t> invindex(n);      //对应到autk 的index

    for (size_t i = 0; i < n; i++) {
        ua[i]   = a[index*n+i].ConvertToInt<int32_t>();       //a
        w[i]    = (2 * N / q) * ua[i] + 1;            //w_i
        invw[i] = ModInverse(w[i], 2 * N) % (2 * N);  //w_inv
    }
    for (size_t i = 0; i < n; i++) {
        NATIVEw[i] = NativeVector::Integer((w[i] * invw[i + 1]) % (2 * N));
        invindex[i] = (NATIVEw[i].ConvertToInt<int32_t>() - 2*N/q -1) / ( 2*N/q);
    }
    
    for (size_t i = 0; i < n; i++) {
        MGAddToAccXZDDF(params, (*ek)[0][0][i], acc);  ///evk_{0 ~ n-1}
        if (NATIVEw[i].ConvertToInt<int32_t>() != 1) {
            MGAutomorphism(params, NATIVEw[i], (*ek)[0][1][invindex[i]], acc);
        }
    }
   
    MGAddToAccXZDDF(params, (*ek)[0][0][n], acc);
}

//KDM-form  完全分解
//仅修改了digitG2的值
VectorNTRUEvalKey VectorNTRUAccumulatorXZDDF::MGKDMKeyGenXZDDF(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                                             const NativePoly& invskNTT, LWEPlaintext m) const {                          
    auto polyParams = params->GetPolyParams();  //(Q,2N)
    auto Gpow       = params->GetGPower();      //
    DiscreteUniformGeneratorImpl<NativeVector> dug;
    NativeInteger Q{params->GetQ()};
    dug.SetModulus(Q);  //确保dug的模数是Q
                        //Reduce mod q (dealing with negative number as well)
    int64_t N  = params->GetN();
    int64_t mm = (((m % N) + N) % N);  // m=-1, mm=N-1;   m=0,mm=0;   m=1,mm=1
    bool isReducedMM{false};
    if (m < 0) {
        isReducedMM = true;
    }

    uint32_t digitsG2{(params->GetDigitsG())};
    std::vector<NativePoly> tempA(digitsG2, NativePoly(dug, polyParams, Format::COEFFICIENT));
    VectorNTRUEvalKeyImpl result(digitsG2);
    for (uint32_t i = 0; i < digitsG2; ++i) {
        result[i] = NativePoly(params->GetDgg(), polyParams, Format::COEFFICIENT);  //采样g
        if (!isReducedMM)//s<0 m>=0 mm=1,0
            result[i][mm].ModAddFastEq(Gpow[i],Q);  // g+X^m*G  g是多项式，在第mm位加上B^{i+1}等于g+B^{i+1}*X^{mm}, 第一个ek与其他的不一样，所以第一个ek的生成用KDMKeyGen 
        else    //s>=0 m<0 mm=N-1
            result[i][mm].ModSubFastEq(Gpow[i],Q);  // g-X^m*G   result[i][N-m]=(g-B^{i+1})/f
        result[i].SetFormat(Format::EVALUATION);
        result[i] = result[i] * invskNTT;
    }
    return std::make_shared<VectorNTRUEvalKeyImpl>(result);
}
//NO KDM-form 完全分解
VectorNTRUEvalKey VectorNTRUAccumulatorXZDDF::MGKeyGenXZDDF(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                                          const NativePoly& invskNTT, LWEPlaintext m) const {
    auto polyParams = params->GetPolyParams();  //(Q,2N)
    auto Gpow       = params->GetGPower();      //
    NativeInteger Q{params->GetQ()};
    int64_t N  = params->GetN();
    int64_t mm = (((m % N) + N) % N);  // s=0 m=0 mm=0; s=q-1 m=1 mm=1; s=1 m=-1 mm=N-1
    //std::cout<<m<<" "<<mm<<std::endl;
    bool isReducedMM{false};
    if (m < 0) {
        isReducedMM = true;
    }
    uint32_t digitsG2{(params->GetDigitsG())};  //logQ < (B)^k的最小的k
    NativePoly zeroPoly(polyParams, Format::COEFFICIENT);
    zeroPoly.SetValuesToZero();
    std::vector<NativePoly> tempA(digitsG2, zeroPoly);

    VectorNTRUEvalKeyImpl result(digitsG2);
    for (uint32_t i = 0; i < digitsG2; ++i) {
        // result[i][0] = tempA[i];
        tempA[i].SetFormat(Format::COEFFICIENT);
        result[i] = NativePoly(params->GetDgg(), polyParams, Format::COEFFICIENT);  //采样g
        result[i].SetFormat(Format::EVALUATION);
        result[i] = result[i] * invskNTT;  // g/f
        if (!isReducedMM){   //m>=0 mm==m=0,1
            tempA[i][mm].ModAddFastEq(Gpow[i], Q);  // X^m*G  Gpow[0]=1
        }
        else{                //m<0 mm=N-1
            tempA[i][mm].ModSubFastEq(Gpow[i], Q);  // X^m*G      这里要modsub是因为 X^{-1} = - X^{N-1}  mod X^{N} + 1
        } 
        tempA[i].SetFormat(Format::EVALUATION);
        result[i] = result[i] + tempA[i];
    }
    return std::make_shared<VectorNTRUEvalKeyImpl>(result);
}

VectorNTRUEvalKey VectorNTRUAccumulatorXZDDF::MGKeyGenAuto(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                                         const NativePoly& skNTT, const NativePoly& invskNTT,
                                                         LWEPlaintext k) const {
    //auto polyParams{params->GetPolyParams()};
    // m_polyParams{std::make_shared<ILNativeParams>(2 * N, Q)},
    // auto Gpow{params->GetGPower()};//m_Gpower,是一个3长度vector (0,1024,1048576)，误导！！！！
    auto polyParams = params->GetPolyParams();  //(Q,2N)
    auto Gpow       = params->GetGPower();      //

    //打印Gpow，这个才是对的，如果gadgetbase=1<<16，那么Gpow=1, 65535=2^16, 4294967296=2^32...
    // for(uint32_t x=0;x<params->GetDigitsG();++x){
    //     std::cout<<Gpow[x]<<" ";
    // }
    //std::cout<<std::endl;
    DiscreteUniformGeneratorImpl<NativeVector> dug;
    NativeInteger Q{params->GetQ()};
    dug.SetModulus(Q);
    auto skAuto{skNTT.AutomorphismTransform(k)};  //生成f(X^k)

    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG{params->GetDigitsG()};
    VectorNTRUEvalKeyImpl result(digitsG);

    for (uint32_t i = 0; i < digitsG; ++i) {
        result[i] = NativePoly(params->GetDgg(), polyParams, EVALUATION) + skAuto * Gpow[i];//这里本来是i+1
        result[i] = result[i] * invskNTT;
    }
    return std::make_shared<VectorNTRUEvalKeyImpl>(result);
}

void VectorNTRUAccumulatorXZDDF::MGAddToAccXZDDF(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                               ConstVectorNTRUEvalKey& ek, NTRUCiphertext& acc) const {
    NativePoly ct(acc->GetElements());
    ct.SetFormat(Format::COEFFICIENT);
    
    
    
    uint32_t digitsG{(params->GetDigitsG())};
    std::vector<NativePoly> dct(digitsG,NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));  // d-1维N长多项式

    CompleteSignedDigitDecompose(params, ct, dct);   //分解acc
                                                
    // calls digitsG2 NTTs
    NativePoly sum(params->GetPolyParams(), Format::EVALUATION, true);

//#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(digitsG))
    #pragma omp parallel for schedule(static)   
    for (uint32_t d = 0; d < digitsG; ++d)
        dct[d].SetFormat(Format::EVALUATION);
    // acc = dct * ek (matrix product);
    const std::vector<NativePoly>& ev = ek->GetElements();
    //std::cout<<"tmp__________________"<<std::endl;
    for (uint32_t d = 0; d < digitsG; ++d){
        auto tmp = dct[d]*ev[d];  //lyh注：两个NTT多项式相乘就是对应项执行模乘法
        sum += tmp;
    }

    
    sum.SetFormat(Format::EVALUATION);
    acc->GetElements() = sum;
}

void VectorNTRUAccumulatorXZDDF::MGAutomorphism(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                              const NativeInteger& a, ConstVectorNTRUEvalKey& ak,
                                              NTRUCiphertext& acc) const {
    // precompute bit reversal for the automorphism into vec
    uint32_t N{params->GetN()};
    std::vector<usint> vec(N);
    PrecomputeAutoMap(N, a.ConvertToInt<usint>(), &vec);  //
    NativePoly ct(acc->GetElements());
    acc->GetElements().SetValuesToZero();
    ct = ct.AutomorphismTransform(a.ConvertToInt<usint>(), vec);
    ct.SetFormat(COEFFICIENT);
    
    uint32_t digitsG{params->GetDigitsG()};
    std::vector<NativePoly> dct(digitsG, NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));
    CompleteSignedDigitDecompose(params, ct, dct);
    NativePoly sum(params->GetPolyParams(), Format::EVALUATION, true);
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(digitsG))
    for (uint32_t d = 0; d < digitsG; ++d)
        dct[d].SetFormat(Format::EVALUATION);
    // acc = dct * input (matrix product);
    const std::vector<NativePoly>& ev = ak->GetElements();
    for (uint32_t d = 0; d < digitsG; ++d)
        sum += (dct[d] * ev[d]);

    acc->GetElements() = sum;
}

//---------------------------MG—ek的生成改为秘密共享版本--------------------------------

VectorNTRUACCKey VectorNTRUAccumulatorXZDDF::MGKeyGenAcc_new(BinFHEContext &cc,
    std::vector<NativePoly> F_shares,std::vector<NativePoly> F_inv_shares,
    std::vector<std::vector<NativePoly>>& sk_sum_shares,std::vector<NativePoly> sk_inv_all_sum_shares, int groupidx, int numofparties) const{

    auto params = cc.GetParams()->GetVectorNTRUParams();
    size_t n{sk_sum_shares.size()};
    auto q{params->Getq().ConvertToInt<size_t>()};
    
    VectorNTRUACCKey ek = std::make_shared<VectorNTRUACCKeyImpl>(1, 2, q - 1 > n + 1 ? q - 1 : n + 1);

    AdditiveSecretSharing secretSharing(params, numofparties);
    secretSharing.TriplesGen(cc, groupidx, numofparties);

    (*ek)[0][0][0] = MGKDMKeyGenXZDDF_new(cc, F_inv_shares, sk_sum_shares[0], groupidx, numofparties, secretSharing);
   
//#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(n))
    for (size_t i = 1; i < n; ++i) {
        (*ek)[0][0][i] = MGKeyGenXZDDF_new(cc, F_inv_shares, sk_sum_shares[i], groupidx, numofparties, secretSharing);
    }
    (*ek)[0][0][n] = MGKeyGenXZDDF_new(cc, F_inv_shares, sk_inv_all_sum_shares, groupidx, numofparties, secretSharing);
    
    
    int64_t intq = params->Getq().ConvertToInt<int64_t>();  
    int64_t N    = params->GetN();
    for (auto i = 0; i < intq - 1; ++i) {
        (*ek)[0][1][i] = MGKeyGenAuto_new(cc, F_shares, F_inv_shares, (2 * N / intq) * (i + 1) + 1, groupidx, numofparties, secretSharing);
    }
    return ek;

}

VectorNTRUEvalKey VectorNTRUAccumulatorXZDDF::MGKDMKeyGenXZDDF_new(BinFHEContext &cc,
 std::vector<NativePoly>& invskNTT_share, std::vector<NativePoly>& Xpows_share, int groupidx, int numofparties, AdditiveSecretSharing& secretSharing) const{
    // (g+Gpow*X^m)/f

    auto params = cc.GetParams()->GetVectorNTRUParams();    
    auto polyParams = params->GetPolyParams();  
    auto Gpow       = params->GetGPower();      
    NativeInteger Q{params->GetQ()};
    int64_t N  = params->GetN();
   
    uint32_t digitsG2{(params->GetDigitsG())};
    VectorNTRUEvalKeyImpl result(digitsG2);
    

    for (uint32_t i = 0; i < digitsG2; ++i) {
        //各方选取噪声，并计算其和
        std::vector<NativePoly> noisesum_share(numofparties);
        for(int partyid=0;partyid<numofparties;++partyid){
            NativePoly noise(params->GetDgg(), polyParams, EVALUATION);
            auto noise_share = secretSharing.Split(noise, numofparties);
            if(partyid == 0){
                noisesum_share = noise_share;
            }
            else{
                noisesum_share = secretSharing.Add(noisesum_share, noise_share);
            }
        }
        NativePoly GpowPoly(polyParams, Format::COEFFICIENT, true);
        GpowPoly[0] = Gpow[i];
        GpowPoly.SetFormat(Format::EVALUATION);
        auto GXpows_share = secretSharing.ScalarMult(GpowPoly, Xpows_share);
        auto upper_share = secretSharing.Add(noisesum_share, GXpows_share);
        auto result_share = secretSharing.Mult(upper_share, invskNTT_share, secretSharing.getTriples()[0][0],secretSharing.getTriples()[0][1],secretSharing.getTriples()[0][2]);
        result[i] = secretSharing.Recover(result_share);
    }
    return std::make_shared<VectorNTRUEvalKeyImpl>(result);
}

VectorNTRUEvalKey VectorNTRUAccumulatorXZDDF::MGKeyGenXZDDF_new(BinFHEContext &cc,
 std::vector<NativePoly>& invskNTT_share, std::vector<NativePoly>& Xpows_share,int groupidx, int numofparties, AdditiveSecretSharing& secretSharing) const{
    // (g/f) + Gpow * X^m

    auto params = cc.GetParams()->GetVectorNTRUParams();
    auto polyParams = params->GetPolyParams();  
    auto Gpow       = params->GetGPower();      
    NativeInteger Q{params->GetQ()};
    int64_t N  = params->GetN();

    uint32_t digitsG2{(params->GetDigitsG())};  

    VectorNTRUEvalKeyImpl result(digitsG2);
    for (uint32_t i = 0; i < digitsG2; ++i) {
        //各方选取噪声，并计算其和
        std::vector<NativePoly> noisesum_share(numofparties);
        for(int partyid=0;partyid<numofparties;++partyid){
            NativePoly noise(params->GetDgg(), polyParams, EVALUATION);
            auto noise_share = secretSharing.Split(noise, numofparties);
            if(partyid == 0){
                noisesum_share = noise_share;
            }
            else{
                noisesum_share = secretSharing.Add(noisesum_share, noise_share);
            }
        }

        auto gdivf = secretSharing.Mult(noisesum_share, invskNTT_share,secretSharing.getTriples()[0][0],secretSharing.getTriples()[0][1],secretSharing.getTriples()[0][2]);
        
        NativePoly GpowPoly(polyParams, Format::COEFFICIENT, true);
        GpowPoly[0] = Gpow[i];
        GpowPoly.SetFormat(Format::EVALUATION);
        auto GXpows_share = secretSharing.ScalarMult(GpowPoly, Xpows_share);

        auto result_share = secretSharing.Add(gdivf, GXpows_share);

        result[i] = secretSharing.Recover(result_share);
    }
    return std::make_shared<VectorNTRUEvalKeyImpl>(result);


}

VectorNTRUEvalKey VectorNTRUAccumulatorXZDDF::MGKeyGenAuto_new(BinFHEContext &cc,
    std::vector<NativePoly>& skNTT_share, std::vector<NativePoly>& invskNTT_share, LWEPlaintext k, int groupidx, int numofparties, AdditiveSecretSharing& secretSharing) const{

    auto params = cc.GetParams()->GetVectorNTRUParams();    
    auto polyParams = params->GetPolyParams();  
    auto Gpow       = params->GetGPower();      

    DiscreteUniformGeneratorImpl<NativeVector> dug;
    NativeInteger Q{params->GetQ()};
    dug.SetModulus(Q);
    

    uint32_t kVal = static_cast<uint32_t>(k);
    std::vector<NativePoly> skAuto_share = secretSharing.Auto(skNTT_share, kVal);

    uint32_t digitsG{params->GetDigitsG()};
    VectorNTRUEvalKeyImpl result(digitsG);

    for (uint32_t i = 0; i < digitsG; ++i) {

        //各方选取噪声，并计算其和
        std::vector<NativePoly> noisesum_share(numofparties);
        for(int partyid=0;partyid<numofparties;++partyid){
            NativePoly noise(params->GetDgg(), polyParams, EVALUATION);
            auto noise_share = secretSharing.Split(noise, numofparties);
            if(partyid == 0){
                noisesum_share = noise_share;
            }
            else{
                noisesum_share = secretSharing.Add(noisesum_share, noise_share);
            }
        }
        
        //计算skAuto * Gpow[i]的分片
        NativePoly GpowPoly(polyParams, Format::COEFFICIENT, true);
        GpowPoly[0] = Gpow[i];
        GpowPoly.SetFormat(Format::EVALUATION);
        auto gskAuto_share = secretSharing.ScalarMult(GpowPoly, skAuto_share);

        //计算和
        auto uppersum_share = secretSharing.Add(noisesum_share, gskAuto_share);
        auto result_share = secretSharing.Mult(uppersum_share, invskNTT_share, secretSharing.getTriples()[0][0], secretSharing.getTriples()[0][1], secretSharing.getTriples()[0][2]);
        
        result[i] = secretSharing.Recover(result_share);
    }
    return std::make_shared<VectorNTRUEvalKeyImpl>(result);
}


};  // namespace lbcrypto
